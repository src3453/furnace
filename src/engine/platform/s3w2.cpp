
#include "s3w2.h"
#include "../engine.h"
#include <math.h>

#define CHIP_DIVIDER 1 // 本来のクロックは9216000Hzだが、無視する
#define CHIP_FREQBASE 48000 // DAC output rate
#define regBase 0x800

#define rWrite(a,v) {if (!skipRegisterWrites) {chip->writeRegister(a,v);regPool[a]=v;if (dumpWrites) addWrite(a,v); }}
#define rRead(a) (regPool[a])
#define rWrite16(a,v) {rWrite(a+1,(v)&0xff); rWrite((a),((v)>>8)&0xff); } // in big endian
#define rWrite24(a,v) {rWrite(a+2,(v)&0xff); rWrite((a+1),((v)>>8)&0xff); rWrite((a),((v)>>16)&0xff); } // in big endian for pcm addr

void DivPlatformS3W2::acquire(short** buf, size_t len) {
  // オシロスコープバッファを初期化
  for (int i=0; i<8; i++) {
    oscBuf[i]->begin(len);
  }

  // オーディオデータを取得
  std::vector<std::vector<std::vector<int16_t>>> output = chip->clock(len);
  //std::printf("Acquired %zu samples\n", len);
  
  // 各サンプルを処理
  for (size_t i=0; i<len; i++) {
    float outL = 0;
    float outR = 0;

    // 各チャンネルを処理
    for (unsigned char j=0; j<8; j++) {
      // マスターアウトプットを計算
      outL += output[j][0][i]/4;
      outR += output[j][1][i]/4;
      // モノラル出力を計算（オシロスコープ用)
      float chanOut = ((output[j][0][i] + output[j][1][i]) / 2);
      if (chanOut < -32768.0f) chanOut = -32768.0f;
      if (chanOut > 32767.0f) chanOut = 32767.0f;
      // オシロスコープバッファにサンプルを格納
      oscBuf[j]->putSample(i, (short)((float)chanOut / 2));
    }
    
    // 値のクランプ
    if (outL < -32768) outL = -32768;
    if (outL > 32767) outL = 32767;
    if (outR < -32768) outR = -32768;
    if (outR > 32767) outR = 32767;
    
    // 出力バッファに書き込み
    buf[0][i] = outL;
    buf[1][i] = outR;
  }
  
  // オシロスコープバッファの終了処理
  for (int i=0; i<8; i++) {
    oscBuf[i]->end(len);
  }
}

void DivPlatformS3W2::updateWave(int ch) {
  if (chan[ch].pcm) return; // PCMモードでは更新しない
  for (int i=0; i<256; i++) {
    rWrite(ch*256+i,(unsigned char)chan[ch].ws.output[i]);
  }
}

/*
3WS8PN (S3W2) Specification
=========================
Overall: Fantasy wavetable soundchip with 8 of wavetable/PCM/noise channels
Wavetable: Size: 256x256
Wavetable Modulation: with combining two wavetables (ex. CH1+CH2, CH3+CH8)... Phase Modulation, Ring Modulation, Hard Sync, Window, ...
PCM: 8bit, RAM 8Mbit (1MB), 20bit I/0 (also can output 16bit with DMA, no volume control)
Noise: LFSR, 1bit output
Mixer: 16bit Stereo Linear PCM, 48KHz Master Output, each channel has 256 volume steps, 16 panpot (each stereo channel has 16 steps)
(Virtual) Clock: 9.216MHz, Main Sample Clock, divived by 48, 192KHz, all sound frequencies were quantized by this value
Registers (all bytes are in big endian):
0x000-0x7FF: Channel Wave Table Data (8 channels x 256 bytes)
 - +0x00~0xFF: Wave Table Data (256 bytes)
 in PCM Mode:
 - +0x00~0x02: PCM Start Address (MSB 4bit: Reserved, LSB 20bit: Address)
 - +0x03~0x05: PCM End Address (MSB 4bit: Reserved, LSB 20bit: Address)
 - +0x06~0x08: PCM Loop Address (MSB 4bit: Reserved, LSB 20bit: Address)
 - +0x09: PCM Playing Control (bit0: Play/Stop, bit1: Loop On/Off, bit2~7: Reserved)
0x800-0x8FF: Channel Control Registers (8 channels x 32 bytes)
 - +0x00~0x01: Frequency (in Hertz, 16bit, in PCM Mode: Sample Rate (32x Frequency (ex. 64hz -> 0x0002)))
 - +0x02: Waveform Type (0: Wavetable, 1: PCM, 2: Noise, 3: DMA PCM (16bit), 4~: Reserved)
 - +0x03: Volume (0~255)
 - +0x04: Panpot (MSB: Left 4bit, LSB: Right 4bit)
 - +0x05: Waveform Modulation Type (MSB 5bit: Type, LSB 3bit: Target Channel (0~7))
 - +0x06~0x09: Modulation Depth/Parameter (depends on Modulation Type)
 - +0x0A: Any access in this register will reset the channel phase
*/

void DivPlatformS3W2::tick(bool sysTick) {
  for (int i=0; i<8; i++) {
    int addrbase = regBase + 0x20 * i;
    chan[i].std.next();
    if (chan[i].std.vol.had) {
      chan[i].outVol=((chan[i].vol&255)*MIN(255,chan[i].std.vol.val))/255;
      rWrite(addrbase+0x03,chan[i].outVol); // +0x03: Volume (0~255)
    }
    if (NEW_ARP_STRAT) {
      chan[i].handleArp();
    } else if (chan[i].std.arp.had) {
      if (!chan[i].inPorta) {
        chan[i].baseFreq=chan[i].calcBaseFreq(parent->calcArp(chan[i].note,chan[i].std.arp.val));
      }
      chan[i].freqChanged=true;
    }
    if (chan[i].std.wave.had) {
      if (chan[i].wave!=chan[i].std.wave.val || chan[i].ws.activeChanged()) {
        chan[i].wave=chan[i].std.wave.val;
        chan[i].ws.changeWave1(chan[i].wave);
      }
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
    if (chan[i].std.duty.had) 
    {
      if (chan[i].pcm) {
        // PCM mode
        rWrite(addrbase+0x02,1); // +0x02: Waveform Type (PCM)
      } else
      // Wavetable/Noise mode
      rWrite(addrbase+0x02,chan[i].std.duty.val==1?2:0); // +0x02: Waveform Type (Wave/Noise)
      
    }
    if (chan[i].std.phaseReset.had && chan[i].std.phaseReset.val == 1) {
      rWrite(addrbase+0x0A,0); // +0x0A: Phase Reset (any access resets phase)
    }
    if (chan[i].std.ex1.had) {
      // Modulation Type : (MSB 5bit: Type, LSB 3bit: Target Channel (0~7))
      uint8_t currentReg = rRead(addrbase+0x05);
      uint8_t newReg = (chan[i].std.ex1.val << 3) | (currentReg & 0x07);
      rWrite(addrbase+0x05,newReg); // +0x05: Waveform Modulation Type
    }
    if (chan[i].std.ex2.had) {
      // Modulation Target : (MSB 5bit: Type, LSB 3bit: Target Channel (0~7))
      uint8_t currentReg = rRead(addrbase+0x05);
      uint8_t newReg = (currentReg & 0xF8) | (chan[i].std.ex2.val & 0x07);
      rWrite(addrbase+0x05,newReg); // +0x05: Waveform Modulation Type
    }
    if (chan[i].std.ex3.had) {
      rWrite16(addrbase+0x06,chan[i].std.ex3.val); // +0x06~0x07: Modulation Depth/Parameter MSB 16bit
    }
    if (chan[i].std.ex4.had) {
      rWrite16(addrbase+0x08,chan[i].std.ex4.val); // +0x08~0x09: Modulation Depth/Parameter LSB 16bit
    }
    if (chan[i].active) {
      if (chan[i].ws.tick()) {
        updateWave(i);
      }
    }
    if (chan[i].keyOn) {
      if (chan[i].pcm) {
          int sNum=chan[i].sample;
          DivSample* sample=parent->getSample(sNum);
          if (sample!=NULL && sNum>=0 && sNum<parent->song.sampleLen) {
            unsigned int sampleEnd=sampleOffSU[sNum]+(sample->getLoopEndPosition());
            unsigned int off=sampleOffSU[sNum]+chan[i].pcmOffset;
            chan[i].hasOffset=0;
            if (sampleEnd>=getSampleMemCapacity(0)) sampleEnd=getSampleMemCapacity(0)-1;
            rWrite24(i*0x100+0, off); // +0x00~0x02: PCM Start Address
            rWrite24(i*0x100+3, sampleEnd); // +0x03~0x05: PCM End Address
            rWrite24(i*0x100+6, 0x000000); // +0x06~0x08: PCM Loop Address (0x000000 = no loop)
            rWrite(i*0x100+9,0b01); // +0x09: PCM Playing Control (bit0: Play/Stop, bit1: Loop On/Off, bit2~7: Reserved)
            if (sample->isLoopable()) {
              unsigned int sampleLoop=sampleOffSU[sNum]+sample->getLoopStartPosition();
              unsigned int sampleLoopEnd=sampleOffSU[sNum]+sample->getLoopEndPosition();
              if (sampleLoop>=getSampleMemCapacity(0)) sampleLoop=getSampleMemCapacity(0)-1;
              rWrite24(i*0x100+6, sampleLoop); // +0x06~0x08: PCM Loop Address
              rWrite24(i*0x100+3, sampleLoopEnd); // +0x03~0x05: PCM End Address
              rWrite(i*0x100+9,0b11); // +0x09: PCM Playing Control (bit0: Play/Stop, bit1: Loop On/Off, bit2~7: Reserved)
              chan[i].pcmLoop=true;
            } else {
              rWrite(i*0x100+9,0b01); // +0x09: PCM Playing Control (bit0: Play/Stop, bit1: Loop On/Off, bit2~7: Reserved)
              chan[i].pcmLoop=false;
            }
          }
          rWrite(regBase+i*0x20+0x02,1); // +0x02: Waveform Type (PCM)
        }
    }
    if (isMuted[i]) {
      rWrite(addrbase+0x03,0); // +0x03: Volume (0~255)
    }
    if (chan[i].freqChanged || chan[i].keyOn || chan[i].keyOff) {
      //DivInstrument* ins=parent->getIns(chan[i].ins,DIV_INS_SU);
      chan[i].freq=parent->calcFreq(chan[i].baseFreq,chan[i].pitch,chan[i].fixedArp?chan[i].baseNoteOverride:chan[i].arpOff,chan[i].fixedArp,0,0,chan[i].pitch2,chipClock,CHIP_FREQBASE)/32;
      if (chan[i].pcm) {
        DivSample* sample=parent->getSample(chan[i].sample);
        if (sample!=NULL) {
          float off=0.5;
          if (sample->centerRate<1) {
            off=0.5;
          } else {
            off=(float)sample->centerRate/(8363.0*1.0);
          }
          chan[i].freq=(float)chan[i].freq*off;
          rWrite16(addrbase+0x00,(unsigned short)(chan[i].freq)); // +0x00~0x01: Sample Rate (32x Frequency)
        }
        chan[i].freqChanged=false;
      } else {
          rWrite16(addrbase+0x00,(unsigned short)(chan[i].freq)); // +0x00~0x01: Frequency
          chan[i].freqChanged=false;
      }
    }
  }
}

int DivPlatformS3W2::dispatch(DivCommand c) {
  int addrbase = regBase + 0x20 * c.chan;
  switch (c.cmd) {
    case DIV_CMD_NOTE_ON: {
      DivInstrument* ins=parent->getIns(chan[c.chan].ins,DIV_INS_S3W2);
      if (chan[c.chan].pcm && !(ins->type==DIV_INS_AMIGA || ins->amiga.useSample)) {
        chan[c.chan].pcm=(ins->type==DIV_INS_AMIGA || ins->amiga.useSample);
      }
      chan[c.chan].pcm=(ins->type==DIV_INS_AMIGA || ins->amiga.useSample);
      if (chan[c.chan].pcm) {
        if (c.value!=DIV_NOTE_NULL) {
          chan[c.chan].sample=ins->amiga.getSample(c.value);
          chan[c.chan].baseFreq=chan[c.chan].calcBaseFreq(c.value);
          chan[c.chan].freqChanged=true;
        }
        rWrite(addrbase+0x0A,0); // +0x0A: Phase Reset (means trigger PCM play from start)
      } else {
        if (c.value!=DIV_NOTE_NULL) {
          chan[c.chan].baseFreq=chan[c.chan].calcBaseFreq(c.value);
          chan[c.chan].freqChanged=true;
        }
      }
      chan[c.chan].active=true;
      chan[c.chan].macroInit(ins);
      if (!parent->song.compatFlags.brokenOutVol && !chan[c.chan].std.vol.will) {
        chan[c.chan].outVol=chan[c.chan].vol;
      }
      if (chan[c.chan].wave<0) {
        chan[c.chan].wave=0;
        chan[c.chan].ws.changeWave1(chan[c.chan].wave);
      }
      rWrite(addrbase+0x03,chan[c.chan].outVol); // +0x03: Volume (0~255)
      chan[c.chan].ws.init(ins,256,255,chan[c.chan].insChanged);
      chan[c.chan].insChanged=false;
      chan[c.chan].keyOn=true;
      break;
      
    }
    case DIV_CMD_NOTE_OFF:
      chan[c.chan].active=false;
      chan[c.chan].keyOff=true;
      chan[c.chan].keyOn=false;
      rWrite(addrbase+0x03,0); // +0x03: Volume (0~255)
      chan[c.chan].macroInit(NULL);
      break;
    case DIV_CMD_NOTE_OFF_ENV:
    case DIV_CMD_ENV_RELEASE:
      chan[c.chan].std.release();
      break;
    case DIV_CMD_INSTRUMENT:
      if (chan[c.chan].ins!=c.value || c.value2==1) {
        chan[c.chan].ins=c.value;
      }
      break;
    case DIV_CMD_VOLUME:
      if (chan[c.chan].vol!=c.value) {
        chan[c.chan].vol=c.value;
        if (!chan[c.chan].std.vol.has) {
          chan[c.chan].outVol=c.value;
          rWrite(addrbase+0x03,c.value); // +0x03: Volume (0~255)

        }
      }
      break;
    case DIV_CMD_GET_VOLUME:
      if (chan[c.chan].std.vol.has) {
        return chan[c.chan].vol;
      }
      return chan[c.chan].outVol;
      break;
    case DIV_CMD_PITCH:
      chan[c.chan].pitch=c.value;
      chan[c.chan].freqChanged=true;
      break;
    case DIV_CMD_WAVE:
      chan[c.chan].wave=c.value;
      chan[c.chan].ws.changeWave1(chan[c.chan].wave);
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
      chan[c.chan].baseFreq=chan[c.chan].calcBaseFreq(c.value+((HACKY_LEGATO_MESS)?(chan[c.chan].std.arp.val):(0)));
      chan[c.chan].freqChanged=true;
      chan[c.chan].note=c.value;
      break;
    case DIV_CMD_PRE_PORTA:
      if (chan[c.chan].active && c.value2) {
        if (parent->song.compatFlags.resetMacroOnPorta) chan[c.chan].macroInit(parent->getIns(chan[c.chan].ins,DIV_INS_SCC));
      }
      if (!chan[c.chan].inPorta && c.value && !parent->song.compatFlags.brokenPortaArp && chan[c.chan].std.arp.will && !NEW_ARP_STRAT) chan[c.chan].baseFreq=chan[c.chan].calcBaseFreq(chan[c.chan].note);
      chan[c.chan].inPorta=c.value;
      break;
    case DIV_CMD_GET_VOLMAX:
      return 255;
      break;
    case DIV_CMD_MACRO_OFF:
      chan[c.chan].std.mask(c.value,true);
      break;
    case DIV_CMD_MACRO_ON:
      chan[c.chan].std.mask(c.value,false);
      break;
    case DIV_CMD_MACRO_RESTART:
      chan[c.chan].std.restart(c.value);
      break;
    default:
      break;
  }
  return 1;
}

void DivPlatformS3W2::muteChannel(int ch, bool mute) {
  isMuted[ch]=mute;
}

/*void DivPlatformS3W2::forceIns() {
  for (int i=0; i<5; i++) {
    chan[i].insChanged=true;
    chan[i].freqChanged=true;
    chan[i].freqInit=false;
    if (isPlus || i<3) {
      updateWave(i);
    }
    rWrite(regBase+10+i,chan[i].outVol);
  }
  if (!isPlus) {
    if (lastUpdated34>=3) {
      updateWave(lastUpdated34);
    }
  }
}*/

SharedChannel* DivPlatformS3W2::getChanState(int ch) {
  return &chan[ch];
}

DivMacroInt* DivPlatformS3W2::getChanMacroInt(int ch) {
  return &chan[ch].std;
}

DivDispatchOscBuffer* DivPlatformS3W2::getOscBuffer(int ch) {
  return oscBuf[ch];
}

unsigned char* DivPlatformS3W2::getRegisterPool() {
  return (unsigned char*)regPool;
}

int DivPlatformS3W2::getRegisterPoolSize() {
  return 0x900;
}

void DivPlatformS3W2::reset() {
  memset(regPool,0,0x900);
  chip->reset();
  for (int i=0; i<8; i++) {
    chan[i]=DivPlatformS3W2::Channel();
    chan[i].pitchTable=&pitchTable;
    chan[i].std.setEngine(parent);
    chan[i].ws.setEngine(parent,128);
    chan[i].ws.init(NULL,256,255,false);
    chan[i].vol=255;
    chan[i].outVol=255;
  }
  for (int i=0; i<sampleMemSize; i++) {
    chip->writePCMRAM(0x0+i,(uint8_t)((sampleMem[i]+128)&0xff));
  }
}

int DivPlatformS3W2::getOutputCount() {
  return 2;
}

bool DivPlatformS3W2::hasAcquireDirect() {
  return false;
}

const void* DivPlatformS3W2::getSampleMem(int index) {
  return (index==0)?sampleMem:NULL;
}

size_t DivPlatformS3W2::getSampleMemCapacity(int index) {
  return (index==0)?((sampleMemSize)):0;
}

size_t DivPlatformS3W2::getSampleMemUsage(int index) {
  return (index==0)?sampleMemLen:0;
}

const DivMemoryComposition* DivPlatformS3W2::getMemCompo(int index) {
  if (index!=0) return NULL;
  return &memCompo;
}

bool DivPlatformS3W2::isSampleLoaded(int index, int sample) {
  if (index!=0) return false;
  if (sample<0 || sample>255) return false;
  return sampleLoaded[sample];
}

void DivPlatformS3W2::renderSamples(int sysID) {
  memset(sampleMem,0,sampleMemSize);
  memset(sampleOffSU,0,256*sizeof(unsigned int));
  memset(sampleLoaded,0,256*sizeof(bool));

  memCompo=DivMemoryComposition();
  memCompo.name="Sample RAM";

  size_t memPos=0;
  for (int i=0; i<parent->song.sampleLen; i++) {
    DivSample* s=parent->song.sample[i];
    if (s->data8==NULL) {
      continue;
    };
    if (!s->renderOn[0][sysID]) {
      sampleOffSU[i]=0;
      continue;
    }
    
    int len=s->getLoopEndPosition(DIV_SAMPLE_DEPTH_8BIT);
    if (len == -1) {
      len=s->getEndPosition(DIV_SAMPLE_DEPTH_8BIT);
    }
    int paddedLen=MIN((int)(getSampleMemCapacity(0)-memPos),len);
    if (memPos>=getSampleMemCapacity(0)) {
      logW("out of PCM memory for sample %d!",i);
      break;
    }
    if (memPos+paddedLen>=getSampleMemCapacity(0)) {
      memcpy(sampleMem+memPos,s->data8,getSampleMemCapacity(0)-memPos-1);
      logW("out of PCM memory for sample %d!",i);
    } else {
      memcpy(sampleMem+memPos,s->data8,paddedLen);
      sampleLoaded[i]=true;
    }
    sampleOffSU[i]=memPos;
    memCompo.entries.push_back(DivMemoryEntry(DIV_MEMORY_SAMPLE,"Sample",i,memPos,memPos+paddedLen));
    memPos+=paddedLen;
  }
  sampleMemLen=memPos;

  for (int i=0; i<sampleMemSize; i++) {
    chip->writePCMRAM(0x0+i,(uint8_t)((sampleMem[i]+128)&0xff));
  }

  memCompo.used=sampleMemLen;
  memCompo.capacity=sampleMemSize;
}


void DivPlatformS3W2::notifyWaveChange(int wave) {
  for (int i=0; i<8; i++) {
    if (chan[i].wave==wave) {
      chan[i].ws.changeWave1(chan[i].wave);
      if (chan[i].active) {
        updateWave(i);
      }
    }
  }
}

void DivPlatformS3W2::notifyInsDeletion(void* ins) {
  for (int i=0; i<8; i++) {
    chan[i].std.notifyInsDeletion((DivInstrument*)ins);
  }
}

void DivPlatformS3W2::poke(unsigned int addr, unsigned short val) {
  rWrite(addr,val);
}

void DivPlatformS3W2::poke(std::vector<DivRegWrite>& wlist) {
  for (DivRegWrite& i: wlist) rWrite(i.addr,i.val);
}

int DivPlatformS3W2::init(DivEngine* p, int channels, int sugRate, const DivConfig& flags) {
  parent=p;
  dumpWrites=false;
  skipRegisterWrites=false;
  writeOscBuf=0;
  for (int i=0; i<8; i++) {
    isMuted[i]=false;
    oscBuf[i]=new DivDispatchOscBuffer;
    oscBuf[i]->rate=CHIP_FREQBASE;
  }
  rate=48000;
  chipClock=48000;
  pitchTable.init(parent->song.tuning,chipClock,CHIP_FREQBASE,0xffff,false,parent->song.compatFlags.linearPitch);
  sampleMemSize=0x100000;// 1024KB (1048576 bytes)
  sampleMem=new unsigned char[sampleMemSize];
  memset(sampleMem,0,sampleMemSize); 
  sampleMemLen=0;
  chip = new S3W2_Sound;
  reset();
  return 8;
}

void DivPlatformS3W2::quit() {
  for (int i=0; i<8; i++) {
    delete oscBuf[i];
  }
  if (chip!=NULL) {
    delete chip;
  }
}

DivPlatformS3W2::~DivPlatformS3W2() {
}
