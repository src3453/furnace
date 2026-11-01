/**
 * Furnace Tracker - multi-system chiptune tracker
 * Copyright (C) 2021-2025 tildearrow and contributors
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include "3hs816.h"
#include "../engine.h"
#include <math.h>

#define CHIP_FREQBASE 48000 // DAC output rate; frequency registers are Hz

#define rWrite(a,v) if (!skipRegisterWrites) {chipWrite(a,v); if (dumpWrites) {addWrite(a,v);} }

// Register offsets inside a channel block (see doc/spec_3HS816.md in CPT32).
#define CH_BASE(c) ((c)*hs816::reg::kChannelStride)

// Operators whose volume follows the channel volume, per algorithm. Kept
// identical to the 3HS88PWN4 (S3HS) dispatch so shared instruments match.
static constexpr bool hs816VolumeScaling[13][8]={
  {true,true,true,true,true,true,true,true},
  {true,true,true,true,false,false,false,false},
  {true,true,true,true,false,false,false,false},
  {true,true,false,false,false,false,false,false},
  {true,false,false,false,false,false,false,false},
  {true,true,false,false,true,true,false,false},
  {true,false,true,false,true,false,true,false},
  {true,false,false,false,false,false,false,false},
  {true,true,false,false,false,false,false,false},
  {true,false,false,false,false,false,false,false},
  {true,true,true,true,false,false,false,false},
  {true,true,false,false,false,false,false,false},
  {true,false,false,false,false,false,false,false}
};

void DivPlatform3HS816::chipWrite(unsigned int addr, unsigned char val) {
  chipRegs->write(addr,val);
  if (addr<hs816::reg::kSize) regPool[addr]=chipRegs->read(addr);
}

void DivPlatform3HS816::writeOperatorVolumes(int ch) {
  const int mode=chan[ch].mode;
  for (int op=0; op<8; op++) {
    int volume=chan[ch].opVolumes[op];
    if (mode>=0 && mode<13 && hs816VolumeScaling[mode][op]) {
      volume=(volume*chan[ch].outVol)/255;
    }
    rWrite(CH_BASE(ch)+hs816::reg::kVolume+op,volume);
  }
}

// Instrument parameters, written on note-on before the gate.
void DivPlatform3HS816::writeInstrument(int ch, DivInstrument* ins) {
  const DivInstrumentS3HS& s=ins->s3hs;
  const unsigned char ratioHigh[7]={s.op2fu,s.op3fu,s.op4fu,s.op5fu,s.op6fu,s.op7fu,s.op8fu};
  const unsigned char ratioLow[7]={s.op2fl,s.op3fl,s.op4fl,s.op5fl,s.op6fl,s.op7fl,s.op8fl};
  const unsigned char waves[8]={s.op1w,s.op2w,s.op3w,s.op4w,s.op5w,s.op6w,s.op7w,s.op8w};
  const unsigned char volumes[8]={s.op1v,s.op2v,s.op3v,s.op4v,s.op5v,s.op6v,s.op7v,s.op8v};
  const unsigned char adsr[8][4]={
    {s.op1a,s.op1d,s.op1s,s.op1r},{s.op2a,s.op2d,s.op2s,s.op2r},
    {s.op3a,s.op3d,s.op3s,s.op3r},{s.op4a,s.op4d,s.op4s,s.op4r},
    {s.op5a,s.op5d,s.op5s,s.op5r},{s.op6a,s.op6d,s.op6s,s.op6r},
    {s.op7a,s.op7d,s.op7s,s.op7r},{s.op8a,s.op8d,s.op8s,s.op8r}
  };
  const unsigned int base=CH_BASE(ch);
  for (int op=1; op<8; op++) {
    rWrite(base+hs816::reg::kRatio+(op-1)*2,ratioHigh[op-1]);
    rWrite(base+hs816::reg::kRatio+(op-1)*2+1,ratioLow[op-1]);
  }
  for (int pair=0; pair<4; pair++) {
    rWrite(base+hs816::reg::kWave+pair,((waves[pair*2]&15)<<4)|(waves[pair*2+1]&15));
  }
  for (int op=0; op<8; op++) {
    for (int i=0; i<4; i++) {
      rWrite(base+hs816::reg::kEnvelope+op*4+i,adsr[op][i]);
    }
    chan[ch].opVolumes[op]=volumes[op];
  }
  chan[ch].mode=s.mode;
  rWrite(base+hs816::reg::kAlgorithm,s.mode);
  rWrite(base+hs816::reg::kFeedback,(int)s.fb+128);
  rWrite(hs816::reg::kFeedbackTargetBase+ch,s.fbDest);
  writeOperatorVolumes(ch);
}

// Song wavetables 0..3 are the chip's user waves A..D: 32 unsigned samples.
void DivPlatform3HS816::uploadUserWave(int wave) {
  if (wave<0 || wave>=hs816::kUserWaveCount) return;
  DivWavetable* wt=(wave<parent->song.waveLen)?parent->getWave(wave):NULL;
  for (int i=0; i<hs816::kUserWaveLength; i++) {
    int value=0;
    if (wt!=NULL && wt->len>0 && wt->max>0) {
      value=wt->data[i*wt->len/hs816::kUserWaveLength]*255/wt->max;
      if (value<0) value=0;
      if (value>255) value=255;
    }
    rWrite(hs816::reg::kUserWaveBase+wave*hs816::kUserWaveLength+i,value);
  }
}

void DivPlatform3HS816::acquire(short** buf, size_t len) {
  hs816::ChannelOutputs taps;
  for (int i=0; i<16; i++) {
    oscBuf[i]->begin(len);
    oscTap[i].resize(len*2);
    taps.left[i]=oscTap[i].data();
    taps.right[i]=oscTap[i].data()+len;
  }

  chip->render(buf[0],buf[1],len,&taps);

  for (int i=0; i<16; i++) {
    for (size_t j=0; j<len; j++) {
      // taps are DAC units; both sides summed for a full-scale scope
      int out=oscTap[i][j]+oscTap[i][len+j];
      if (out<-32768) out=-32768;
      if (out>32767) out=32767;
      oscBuf[i]->putSample(j,(short)out);
    }
    oscBuf[i]->end(len);
  }
}

void DivPlatform3HS816::tick(bool sysTick) {
  for (int i=0; i<16; i++) {
    const unsigned int base=CH_BASE(i);
    chan[i].std.next();
    if (chan[i].std.vol.had) {
      chan[i].outVol=(chan[i].vol*MIN(255,chan[i].std.vol.val))/255;
      writeOperatorVolumes(i);
    }
    if (NEW_ARP_STRAT) {
      chan[i].handleArp();
    } else if (chan[i].std.arp.had) {
      if (!chan[i].inPorta) {
        chan[i].baseFreq=chan[i].calcBaseFreq(parent->calcArp(chan[i].note,chan[i].std.arp.val));
      }
      chan[i].freqChanged=true;
    }
    if (chan[i].std.pitch.had) {
      if (chan[i].std.pitch.mode) {
        chan[i].pitch2+=chan[i].std.pitch.val;
        CLAMP_VAR(chan[i].pitch2,-32768,32767);
      } else {
        chan[i].pitch2=chan[i].std.pitch.val;
      }
      chan[i].freqChanged=true;
    }
    if (chan[i].std.panL.had) {
      chan[i].pan=((chan[i].std.panL.val&15)<<4)|(chan[i].pan&15);
      rWrite(base+hs816::reg::kPan,chan[i].pan);
    }
    if (chan[i].std.panR.had) {
      chan[i].pan=(chan[i].pan&0xf0)|(chan[i].std.panR.val&15);
      rWrite(base+hs816::reg::kPan,chan[i].pan);
    }
    if (chan[i].std.phaseReset.had && chan[i].std.phaseReset.val==1 && chan[i].active) {
      // gate=1 restarts envelopes and phases
      rWrite(base+hs816::reg::kGate,1);
    }
    if (chan[i].freqChanged || chan[i].keyOn || chan[i].keyOff) {
      double frequency=chan[i].calcFreq();
      chan[i].freq=(int)CLAMP(round(frequency),0.0,65535.0);
      rWrite(base+hs816::reg::kFrequency,chan[i].freq>>8);
      rWrite(base+hs816::reg::kFrequency+1,chan[i].freq&0xff);
      chan[i].freqChanged=false;
      chan[i].keyOn=false;
      chan[i].keyOff=false;
    }
  }
}

int DivPlatform3HS816::dispatch(DivCommand c) {
  if (c.chan<0 || c.chan>=16) return 0;
  const unsigned int base=CH_BASE(c.chan);
  switch (c.cmd) {
    case DIV_CMD_NOTE_ON: {
      DivInstrument* ins=parent->getIns(chan[c.chan].ins,DIV_INS_S3HS);
      if (c.value!=DIV_NOTE_NULL) {
        chan[c.chan].baseFreq=chan[c.chan].calcBaseFreq(c.value);
        chan[c.chan].freqChanged=true;
        chan[c.chan].note=c.value;
      }
      chan[c.chan].active=true;
      chan[c.chan].keyOn=true;
      chan[c.chan].macroInit(ins);
      if (!parent->song.compatFlags.brokenOutVol && !chan[c.chan].std.vol.will) {
        chan[c.chan].outVol=chan[c.chan].vol;
      }
      chan[c.chan].insChanged=false;
      writeInstrument(c.chan,ins);
      rWrite(base+hs816::reg::kGate,1);
      break;
    }
    case DIV_CMD_NOTE_OFF:
      chan[c.chan].active=false;
      chan[c.chan].keyOff=true;
      chan[c.chan].keyOn=false;
      chan[c.chan].macroInit(NULL);
      rWrite(base+hs816::reg::kGate,0);
      break;
    case DIV_CMD_NOTE_OFF_ENV:
      chan[c.chan].active=false;
      chan[c.chan].keyOff=true;
      chan[c.chan].keyOn=false;
      chan[c.chan].std.release();
      rWrite(base+hs816::reg::kGate,0);
      break;
    case DIV_CMD_ENV_RELEASE:
      chan[c.chan].std.release();
      break;
    case DIV_CMD_INSTRUMENT:
      if (chan[c.chan].ins!=c.value || c.value2==1) {
        chan[c.chan].ins=c.value;
        chan[c.chan].insChanged=true;
      }
      break;
    case DIV_CMD_VOLUME:
      chan[c.chan].vol=MIN(255,c.value);
      if (!chan[c.chan].std.vol.has) {
        chan[c.chan].outVol=chan[c.chan].vol;
      }
      writeOperatorVolumes(c.chan);
      break;
    case DIV_CMD_GET_VOLUME:
      return chan[c.chan].vol;
    case DIV_CMD_GET_VOLMAX:
      return 255;
    case DIV_CMD_PITCH:
      chan[c.chan].pitch=c.value;
      chan[c.chan].freqChanged=true;
      break;
    case DIV_CMD_NOTE_PORTA: {
      int destFreq=chan[c.chan].calcBaseFreq(c.value2);
      bool return2=false;
      if (destFreq>chan[c.chan].baseFreq) {
        chan[c.chan].baseFreq+=c.value;
        if (chan[c.chan].baseFreq>=destFreq) {
          chan[c.chan].baseFreq=destFreq;
          return2=true;
        }
      } else {
        chan[c.chan].baseFreq-=c.value;
        if (chan[c.chan].baseFreq<=destFreq) {
          chan[c.chan].baseFreq=destFreq;
          return2=true;
        }
      }
      chan[c.chan].freqChanged=true;
      if (return2) {
        chan[c.chan].inPorta=false;
        return 2;
      }
      break;
    }
    case DIV_CMD_LEGATO:
      chan[c.chan].baseFreq=chan[c.chan].calcBaseFreq(c.value);
      chan[c.chan].freqChanged=true;
      chan[c.chan].note=c.value;
      break;
    case DIV_CMD_PRE_PORTA:
      if (chan[c.chan].active && c.value2) {
        if (parent->song.compatFlags.resetMacroOnPorta) chan[c.chan].macroInit(parent->getIns(chan[c.chan].ins,DIV_INS_S3HS));
      }
      chan[c.chan].inPorta=c.value;
      break;
    case DIV_CMD_PANNING:
      chan[c.chan].pan=((c.value>>4)<<4)|(c.value2>>4);
      rWrite(base+hs816::reg::kPan,chan[c.chan].pan);
      break;
    case DIV_CMD_S3HS_OP_VOLUME:
      if (c.value>=0 && c.value<8) {
        chan[c.chan].opVolumes[c.value]=c.value2;
        writeOperatorVolumes(c.chan);
      }
      break;
    case DIV_CMD_S3HS_OP_WAVE:
      if (c.value>=0 && c.value<8) {
        const unsigned int addr=base+hs816::reg::kWave+c.value/2;
        int reg=chipRegs->read(addr);
        if (c.value%2==0) {
          reg=((c.value2&15)<<4)|(reg&15);
        } else {
          reg=(reg&0xf0)|(c.value2&15);
        }
        rWrite(addr,reg);
      }
      break;
    case DIV_CMD_S3HS_MODMODE:
      chan[c.chan].mode=c.value;
      rWrite(base+hs816::reg::kAlgorithm,c.value);
      writeOperatorVolumes(c.chan);
      break;
    case DIV_CMD_S3HS_FEEDBACK:
      rWrite(base+hs816::reg::kFeedback,c.value);
      break;
    case DIV_CMD_S3HS_OP_FREQ_FU:
      if (c.value>=0 && c.value<8) rWrite(base+c.value*2,c.value2);
      break;
    case DIV_CMD_S3HS_OP_FREQ_FL:
      if (c.value>=0 && c.value<8) rWrite(base+c.value*2+1,c.value2);
      break;
    default:
      break;
  }
  return 1;
}

void DivPlatform3HS816::muteChannel(int ch, bool mute) {
  isMuted[ch]=mute;
  rWrite(hs816::reg::kMuteBase+ch,mute?1:0);
}

SharedChannel* DivPlatform3HS816::getChanState(int ch) {
  return &chan[ch];
}

DivMacroInt* DivPlatform3HS816::getChanMacroInt(int ch) {
  return &chan[ch].std;
}

DivDispatchOscBuffer* DivPlatform3HS816::getOscBuffer(int ch) {
  return oscBuf[ch];
}

unsigned short DivPlatform3HS816::getPan(int ch) {
  return ((chan[ch].pan>>4)<<8)|(chan[ch].pan&15);
}

unsigned char* DivPlatform3HS816::getRegisterPool() {
  return regPool;
}

int DivPlatform3HS816::getRegisterPoolSize() {
  return hs816::reg::kSize;
}

int DivPlatform3HS816::getOutputCount() {
  return 2;
}

void DivPlatform3HS816::notifyWaveChange(int wave) {
  uploadUserWave(wave);
}

void DivPlatform3HS816::notifyInsDeletion(void* ins) {
  for (int i=0; i<16; i++) {
    chan[i].std.notifyInsDeletion((DivInstrument*)ins);
  }
}

void DivPlatform3HS816::notifyPitchTable(int sample) {
  // Fold the /32 of DivPitchTable's frequency mode into its scale so the
  // table yields Hz (as S3W2 does).
  pitchTable.init(parent->song.tuning,chipClock,CHIP_FREQBASE/32.0,0xffff,false,parent->song.compatFlags.linearPitch);
}

void DivPlatform3HS816::poke(unsigned int addr, unsigned short val) {
  rWrite(addr,val);
}

void DivPlatform3HS816::poke(std::vector<DivRegWrite>& wlist) {
  for (DivRegWrite& i: wlist) rWrite(i.addr,i.val);
}

void DivPlatform3HS816::reset() {
  chipRegs->reset();
  for (unsigned int i=0; i<hs816::reg::kSize; i++) {
    regPool[i]=chipRegs->read(i);
  }
  for (int i=0; i<16; i++) {
    chan[i]=DivPlatform3HS816::Channel(parent->song.compatFlags.linearPitch);
    chan[i].pitchTable=&pitchTable;
    chan[i].std.setEngine(parent);
    // neutral feedback (the chip resets to -1 on OP1)
    rWrite(CH_BASE(i)+hs816::reg::kFeedback,0x80);
    rWrite(hs816::reg::kMuteBase+i,isMuted[i]?1:0);
  }
  for (int w=0; w<hs816::kUserWaveCount; w++) {
    uploadUserWave(w);
  }
}

int DivPlatform3HS816::init(DivEngine* p, int channels, int sugRate, const DivConfig& flags) {
  parent=p;
  dumpWrites=false;
  skipRegisterWrites=false;
  for (int i=0; i<16; i++) {
    isMuted[i]=false;
    oscBuf[i]=new DivDispatchOscBuffer;
    oscBuf[i]->setRate(hs816::kSampleRate);
  }
  rate=hs816::kSampleRate;
  chipClock=hs816::kSampleRate;
  notifyPitchTable();
  chip=new hs816::Engine();
  chipRegs=new hs816::RegisterInterface(*chip);
  reset();
  return 16;
}

void DivPlatform3HS816::quit() {
  for (int i=0; i<16; i++) {
    delete oscBuf[i];
  }
  delete chipRegs;
  delete chip;
}

DivPlatform3HS816::~DivPlatform3HS816() {
}
