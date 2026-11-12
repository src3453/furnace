#include "34si22.h"
#include "../engine.h"
#include <cmath>
#include <cstring>

void DivPlatform34SI22::write(unsigned a,unsigned v) {
  if (skipRegisterWrites || a>=unsigned(getRegisterPoolSize())) return;
  chip.write(a,v); registerPool[a]=chip.read(a);
  if (dumpWrites) addWrite(a,v);
}
void DivPlatform34SI22::controlWrite(unsigned ch,unsigned cmd) {
  unsigned commands=0x0f; commands=(commands&~(3u<<(ch*2)))|(cmd<<(ch*2));
  write(0x0c,control|commands);
}
void DivPlatform34SI22::address(unsigned a,unsigned v) { write(a,v>>16); write(a+1,v>>8); write(a+2,v); }
void DivPlatform34SI22::updateVolume() {
  // The FIFO DAC ignores PLAY_VOL (fixed 0 dB); its channels scale the PCM they stream.
  int v[2];
  for (int i=0;i<2;i++) v[i]=muted[i]?0:(chan[i].fifoMode?15:chan[i].outVol);
  write(0x0d,v[0]|(v[1]<<4));
}

// Next 16-bit FIFO sample of an AMIGA instrument: no interpolation, loops as in
// Generic PCM DAC, linear volume 0..15. 0 after a one-shot sample ends.
int16_t DivPlatform34SI22::fifoSample(Channel& c) {
  DivSample* s=(c.sample>=0 && c.sample<parent->song.sampleLen)?parent->getSample(c.sample):NULL;
  if (s==NULL || s->samples<=0 || s->data16==NULL) { c.sample=-1; return 0; }
  int out=(c.audPos>=0 && c.audPos<(int)s->samples)?s->data16[c.audPos]:0;
  c.audSub+=c.step;
  while (c.audSub>=0x10000) {
    c.audSub-=0x10000;
    c.audPos+=c.audDir?-1:1;
    if (s->isLoopable()) {
      if (c.audDir) {
        if (c.audPos<s->loopStart) {
          if (s->loopMode==DIV_SAMPLE_LOOP_BACKWARD) c.audPos=s->loopEnd-1-(s->loopStart-c.audPos);
          else { c.audPos=s->loopStart+(s->loopStart-c.audPos); c.audDir=false; }
        }
      } else if (c.audPos>=s->loopEnd) {
        if (s->loopMode==DIV_SAMPLE_LOOP_FORWARD) c.audPos=(c.audPos+s->loopStart)-s->loopEnd;
        else { c.audPos=s->loopEnd-1-(c.audPos-s->loopEnd+1); c.audDir=true; }
      }
    } else if (c.audPos>=(int)s->samples) {
      c.sample=-1;
      break;
    }
  }
  return int16_t((out*CLAMP(c.outVol,0,15))/15);
}
void DivPlatform34SI22::updatePan() {
  control=0;
  for (int i=0;i<2;i++) { if (!(chan[i].pan&1)) control|=0x10<<(i*2); if (!(chan[i].pan&2)) control|=0x20<<(i*2); }
  write(0x0c,control|0x0f);
}
void DivPlatform34SI22::acquire(short** buf,size_t len) {
  for (int c=0;c<2;c++) osc[c]->begin(len);
  for (size_t i=0;i<len;i++) {
    // The driver keeps the FIFO fed one sample ahead (DIV 0: one sample per tick).
    for (int c=0;c<2;c++) {
      if (!chan[c].fifoPlaying || chip.fifoCount(c)) continue;
      int16_t pcm=fifoSample(chan[c]);
      uint16_t v=muted[c]?0:uint16_t(pcm); // PLAY_VOL cannot mute the FIFO DAC
      chip.write(0x25+c*2,v>>8); chip.write(0x26+c*2,v&0xff);
    }
    auto out=chip.tick(); buf[0][i]=out[0]; buf[1][i]=out[1];
    for (int c=0;c<2;c++) osc[c]->putSample(i,chip.channelOutput(c));
  }
  for (int c=0;c<2;c++) osc[c]->end(len);
}
void DivPlatform34SI22::tick(bool) {
  for (int i=0;i<2;i++) {
    auto& c=chan[i]; c.std.next();
    if (c.fifoMode) {
      if (NEW_ARP_STRAT) {
        c.handleArp();
      } else if (c.std.arp.had && !c.rawFreq) {
        if (!c.inPorta) c.baseFreq=c.calcBaseFreq(parent->calcArp(c.note,c.std.arp.val));
        c.freqChanged=true;
      }
      if (c.std.pitch.had) {
        if (c.std.pitch.mode) {
          c.pitch2+=c.std.pitch.val;
          CLAMP_VAR(c.pitch2,-32768,32767);
        } else {
          c.pitch2=c.std.pitch.val;
        }
        c.freqChanged=true;
      }
    }
    if (c.std.vol.had) { c.outVol=(c.vol*CLAMP(c.std.vol.val,0,15))/15; updateVolume(); }
    if (c.std.duty.had && !c.fifoMode) { c.directDivider=true; write(0x0e + i,CLAMP(c.std.duty.val,0,15)); }
    if (c.std.panL.had) { c.pan=((c.std.panL.val&1)<<1)|((c.std.panL.val&2)>>1); updatePan(); }
    if (c.std.phaseReset.had && c.std.phaseReset.val && c.active) c.keyOn=true;
    if (c.keyOn && c.fifoMode) {
      // FIFO DAC: stop (empties the FIFO), PLAY_FORMAT 5, DIV 0, restart.
      controlWrite(i,0);
      unsigned f=chip.read(0x21), shift=i*4; write(0x21,(f&~(7u<<shift))|(5u<<shift));
      write(0x0e + i,0);
      c.audPos=0; c.audSub=0; c.audDir=false; c.freqChanged=true;
      c.fifoPlaying=c.sample>=0;
      controlWrite(i,1);
      updateVolume();
      c.keyOn=false;
    }
    if (c.keyOn) {
      controlWrite(i,0);
      if (c.sample>=0 && c.sample<32768 && loaded[c.sample]) {
        address(i*6,offsets[c.sample]); address(i*6+3,offsets[c.sample]+lengths[c.sample]);
        unsigned f=chip.read(0x21), shift=i*4; write(0x21,(f&~(7u<<shift))|(formats[c.sample]<<shift));
        controlWrite(i,1);
      }
      c.keyOn=false;
    }
    if (c.keyOff) { controlWrite(i,0); c.fifoPlaying=false; c.keyOff=false; }
    if (c.freqChanged && c.fifoMode) {
      // Use Generic PCM DAC's 8.16 pitch units directly, including raw notes.
      c.freq=c.calcFreq();
      if (!c.rawFreq && c.freq>16777215) c.freq=16777215;
      c.step=c.freq;
    } else if (c.freqChanged && !c.directDivider && !c.dividerPending) {
      double freq=c.calcFreq(); auto* s=parent->getSample(c.sample);
      double center=(s && s->centerRate>0)?s->centerRate:parent->getCenterRate();
      if (!c.rawFreq) freq*=center/(440.0*std::pow(2.0,-9.0/12.0));
      int div=(freq>0)?int(std::round(48000.0/freq))-1:15;
      write(0x0e + i,CLAMP(div,0,15));
    }
    c.freqChanged=false;
    if (c.dividerPending) {
      // 20xx sets the hardware divider of memory playback; FIFO DAC channels stay at DIV 0.
      if (!c.fifoMode) { c.directDivider=true; write(0x0e + i,c.pendingDivider); }
      c.dividerPending=false;
    }
  }
}
int DivPlatform34SI22::dispatch(DivCommand cmd) {
  Channel& c=chan[cmd.chan];
  switch (cmd.cmd) {
    case DIV_CMD_NOTE_ON: {
      DivInstrument* ins=parent->getIns(c.ins,DIV_INS_34SI22);
      bool fifo=variantA && ins->type==DIV_INS_AMIGA;
      if (fifo!=c.fifoMode) { c.fifoMode=fifo; updateVolume(); }
      int note=cmd.value;
      if (note!=DIV_NOTE_NULL) {
        c.directDivider=false;
        c.sample=ins->amiga.getSample(note); c.sampleNote=note;
        note=ins->amiga.getFreq(note); c.sampleNoteDelta=note-c.sampleNote;
      } else if (fifo && c.sampleNote!=DIV_NOTE_NULL) {
        c.sample=ins->amiga.getSample(c.sampleNote);
        note=ins->amiga.getFreq(c.sampleNote);
      }
      c.pitchTable=fifo?fifoPitchTable.get(c.sample):&pitchTable;
      if (note!=DIV_NOTE_NULL) { c.note=note; c.baseFreq=c.calcBaseFreq(note); c.freqChanged=true; }
      c.macroInit(ins); if (!c.std.vol.will) c.outVol=c.vol;
      c.active=true; c.keyOn=true; c.keyOff=false; updateVolume(); break;
    }
    case DIV_CMD_NOTE_OFF: c.active=false; c.keyOff=true; c.keyOn=false; c.macroInit(NULL); break;
    case DIV_CMD_NOTE_OFF_ENV: case DIV_CMD_ENV_RELEASE: c.std.release(); break;
    case DIV_CMD_INSTRUMENT: c.ins=cmd.value; c.insChanged=true; break;
    case DIV_CMD_VOLUME: c.vol=CLAMP(cmd.value,0,15); if (!c.std.vol.has) { c.outVol=c.vol; updateVolume(); } break;
    case DIV_CMD_GET_VOLUME: return c.vol;
    case DIV_CMD_GET_VOLMAX: return 15;
    case DIV_CMD_PANNING: c.pan=(cmd.value?1:0)|(cmd.value2?2:0); updatePan(); break;
    case DIV_CMD_SAMPLE_FREQ: c.pendingDivider=cmd.value&15; c.dividerPending=true; break;
    case DIV_CMD_PITCH: c.pitch=cmd.value; c.freqChanged=true; break;
    case DIV_CMD_LEGATO: {
      // HACKY_LEGATO_MESS uses a command named c; spell it out for this channel reference.
      bool relativeArp=!c.rawFreq && c.std.arp.will &&
        ((c.std.arp.val&0xc0000000)==0 || (c.std.arp.val&0xc0000000)==0xc0000000) && !NEW_ARP_STRAT;
      c.baseFreq=c.calcBaseFreq(cmd.value+c.sampleNoteDelta+((c.fifoMode && relativeArp)?c.std.arp.val:0));
      c.note=c.fifoMode?cmd.value:cmd.value+c.sampleNoteDelta;
      c.freqChanged=true; break;
    }
    case DIV_CMD_NOTE_PORTA: { int dest=c.calcBaseFreq(cmd.value2+c.sampleNoteDelta); if (c.baseFreq<dest) c.baseFreq=std::min(c.baseFreq+cmd.value,dest); else c.baseFreq=std::max(c.baseFreq-cmd.value,dest); c.freqChanged=true; if (c.baseFreq==dest) { c.inPorta=false; return 2; } break; }
    case DIV_CMD_PRE_PORTA:
      if (c.fifoMode && c.active && cmd.value2 && parent->song.compatFlags.resetMacroOnPorta) c.macroInit(parent->getIns(c.ins,DIV_INS_AMIGA));
      c.inPorta=cmd.value; break;
    case DIV_CMD_MACRO_OFF: c.std.mask(cmd.value,true); break;
    case DIV_CMD_MACRO_ON: c.std.mask(cmd.value,false); break;
    case DIV_CMD_MACRO_RESTART: c.std.restart(cmd.value); break;
    default: break;
  }
  return 1;
}
void DivPlatform34SI22::reset() {
  chip.reset(); memset(registerPool,0,sizeof(registerPool)); control=0;
  for (int i=0;i<2;i++) { chan[i]=Channel(); chan[i].pitchTable=&pitchTable; chan[i].std.setEngine(parent); }
  updateVolume();
}
void DivPlatform34SI22::notifyPitchTable(int sample) {
  pitchTable.init(parent->song.tuning,48000,48000.0,0xffff,false,parent->song.compatFlags.linearPitch);
  fifoPitchTable.update<Channel>(chan,2,parent->song.tuning,48000,65536,0xffffff,false,parent->song.compatFlags.linearPitch,sample);
}
void DivPlatform34SI22::setFlags(const DivConfig& flags) {
  // chipType 1: 34SI22A (FIFO DAC, SGU v3.2); 0: 34SI22.
  variantA=flags.getInt("chipType",0)==1;
  chip.setVariantA(variantA);
  if (!variantA) for (auto& c:chan) { c.fifoMode=false; c.fifoPlaying=false; c.pitchTable=&pitchTable; }
  size_t capacity=flags.getBool("largeSampleMemory",false)?0x1000000:0x100000;
  if (memory.size()!=capacity) {
    memory.resize(capacity);
    used=0;
    memset(loaded,0,sizeof(loaded));
    composition=DivMemoryComposition();
    composition.name="Sample RAM";
    composition.capacity=capacity;
  }
}
void DivPlatform34SI22::notifyInsDeletion(void* ins) { for (auto& c:chan) c.std.notifyInsDeletion((DivInstrument*)ins); }
void DivPlatform34SI22::muteChannel(int ch,bool mute) { muted[ch]=mute; updateVolume(); }
unsigned char* DivPlatform34SI22::getRegisterPool() { for (int i=0;i<getRegisterPoolSize();i++) registerPool[i]=chip.read(i); return registerPool; }
void DivPlatform34SI22::poke(std::vector<DivRegWrite>& list) { for (auto& w:list) write(w.addr,w.val); }
void DivPlatform34SI22::renderSamples(int sysID) {
  std::fill(memory.begin(),memory.end(),0); memset(loaded,0,sizeof(loaded)); used=0; composition=DivMemoryComposition(); composition.name="Sample RAM";
  for (int i=0;i<parent->song.sampleLen && i<32768;i++) {
    DivSample* s=parent->song.sample[i]; if (!s->renderOn[0][sysID]) continue;
    std::vector<uint8_t> data; uint8_t format=0;
    if (divSampleIsSQV4(s->depth)) { data=s->getSQV4(s->samples); format=2; }
    else if (s->depth==DIV_SAMPLE_DEPTH_16BIT) { format=1; data.reserve(size_t(s->samples)*2); for (unsigned j=0;j<s->samples;j++) { uint16_t x=uint16_t(s->data16[j]); data.push_back(x>>8); data.push_back(x); } }
    else { data.reserve(s->samples); for (unsigned j=0;j<s->samples;j++) data.push_back(uint8_t(int(s->data8[j])+128)); }
    if (data.empty() || data.size()>memory.size()-used || used+data.size()>0xffffff) continue;
    offsets[i]=unsigned(used); lengths[i]=unsigned(data.size()); formats[i]=format; loaded[i]=true;
    std::copy(data.begin(),data.end(),memory.begin()+used); composition.entries.push_back(DivMemoryEntry(DIV_MEMORY_SAMPLE,"Sample",i,used,used+data.size())); used+=data.size();
  }
  composition.capacity=memory.size(); composition.used=used;
}
int DivPlatform34SI22::init(DivEngine* p,int,int,const DivConfig& flags) {
  parent=p; fifoPitchTable.init(parent); dumpWrites=false; skipRegisterWrites=false; rate=48000; chipClock=48000;
  for (int i=0;i<2;i++) { osc[i]=new DivDispatchOscBuffer; osc[i]->setRate(rate); muted[i]=false; }
  notifyPitchTable(); setFlags(flags); reset(); return 2;
}
void DivPlatform34SI22::quit() { fifoPitchTable.destroy<Channel>(chan,2); for (auto* o:osc) delete o; }
