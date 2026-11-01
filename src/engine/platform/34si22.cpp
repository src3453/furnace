#include "34si22.h"
#include "../engine.h"
#include <cmath>
#include <cstring>

void DivPlatform34SI22::write(unsigned a,unsigned v) {
  if (skipRegisterWrites || a>=0x29) return;
  chip.write(a,v); registerPool[a]=chip.read(a);
  if (dumpWrites) addWrite(a,v);
}
void DivPlatform34SI22::controlWrite(unsigned ch,unsigned cmd) {
  unsigned commands=0x0f; commands=(commands&~(3u<<(ch*2)))|(cmd<<(ch*2));
  write(0x0c,control|commands);
}
void DivPlatform34SI22::address(unsigned a,unsigned v) { write(a,v>>16); write(a+1,v>>8); write(a+2,v); }
void DivPlatform34SI22::updateVolume() { write(0x0d,(muted[0]?0:chan[0].outVol)|(muted[1]?0:chan[1].outVol)<<4); }
void DivPlatform34SI22::updatePan() {
  control=0;
  for (int i=0;i<2;i++) { if (!(chan[i].pan&1)) control|=0x10<<(i*2); if (!(chan[i].pan&2)) control|=0x20<<(i*2); }
  write(0x0c,control|0x0f);
}
void DivPlatform34SI22::acquire(short** buf,size_t len) {
  for (int c=0;c<2;c++) osc[c]->begin(len);
  for (size_t i=0;i<len;i++) {
    auto out=chip.tick(); buf[0][i]=out[0]; buf[1][i]=out[1];
    for (int c=0;c<2;c++) osc[c]->putSample(i,chip.channelOutput(c));
  }
  for (int c=0;c<2;c++) osc[c]->end(len);
}
void DivPlatform34SI22::tick(bool) {
  for (int i=0;i<2;i++) {
    auto& c=chan[i]; c.std.next();
    if (c.std.vol.had) { c.outVol=(c.vol*CLAMP(c.std.vol.val,0,15))/15; updateVolume(); }
    if (c.std.duty.had) { c.directDivider=true; write(0x0e+i,CLAMP(c.std.duty.val,0,15)); }
    if (c.std.panL.had) { c.pan=((c.std.panL.val&1)<<1)|((c.std.panL.val&2)>>1); updatePan(); }
    if (c.std.phaseReset.had && c.std.phaseReset.val && c.active) c.keyOn=true;
    if (c.keyOn) {
      controlWrite(i,0);
      if (c.sample>=0 && c.sample<32768 && loaded[c.sample]) {
        address(i*6,offsets[c.sample]); address(i*6+3,offsets[c.sample]+lengths[c.sample]);
        unsigned f=chip.read(0x21), shift=i*4; write(0x21,(f&~(7u<<shift))|(formats[c.sample]<<shift));
        controlWrite(i,1);
      }
      c.keyOn=false;
    }
    if (c.keyOff) { controlWrite(i,0); c.keyOff=false; }
    if (c.freqChanged && !c.directDivider && !c.dividerPending) {
      double freq=c.calcFreq(); auto* s=parent->getSample(c.sample);
      double center=(s && s->centerRate>0)?s->centerRate:parent->getCenterRate();
      if (!c.rawFreq) freq*=center/(440.0*std::pow(2.0,-9.0/12.0));
      int div=(freq>0)?int(std::round(48000.0/freq))-1:15;
      write(0x0e+i,CLAMP(div,0,15));
    }
    c.freqChanged=false;
    if (c.dividerPending) { c.directDivider=true; write(0x0e+i,c.pendingDivider); c.dividerPending=false; }
  }
}
int DivPlatform34SI22::dispatch(DivCommand cmd) {
  Channel& c=chan[cmd.chan];
  switch (cmd.cmd) {
    case DIV_CMD_NOTE_ON: {
      DivInstrument* ins=parent->getIns(c.ins,DIV_INS_34SI22);
      if (cmd.value!=DIV_NOTE_NULL) { c.directDivider=false; c.sample=ins->amiga.getSample(cmd.value); c.sampleNote=cmd.value; c.note=ins->amiga.getFreq(cmd.value); c.sampleNoteDelta=c.note-c.sampleNote; c.baseFreq=c.calcBaseFreq(c.note); c.freqChanged=true; }
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
    case DIV_CMD_LEGATO: c.note=cmd.value+c.sampleNoteDelta; c.baseFreq=c.calcBaseFreq(c.note); c.freqChanged=true; break;
    case DIV_CMD_NOTE_PORTA: { int dest=c.calcBaseFreq(cmd.value2+c.sampleNoteDelta); if (c.baseFreq<dest) c.baseFreq=std::min(c.baseFreq+cmd.value,dest); else c.baseFreq=std::max(c.baseFreq-cmd.value,dest); c.freqChanged=true; if (c.baseFreq==dest) { c.inPorta=false; return 2; } break; }
    case DIV_CMD_PRE_PORTA: c.inPorta=cmd.value; break;
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
void DivPlatform34SI22::notifyPitchTable(int) { pitchTable.init(parent->song.tuning,48000,48000.0,0xffff,false,parent->song.compatFlags.linearPitch); }
void DivPlatform34SI22::setFlags(const DivConfig& flags) {
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
unsigned char* DivPlatform34SI22::getRegisterPool() { for (int i=0;i<0x29;i++) registerPool[i]=chip.read(i); return registerPool; }
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
  parent=p; dumpWrites=false; skipRegisterWrites=false; rate=48000; chipClock=48000;
  for (int i=0;i<2;i++) { osc[i]=new DivDispatchOscBuffer; osc[i]->setRate(rate); muted[i]=false; }
  notifyPitchTable(); setFlags(flags); reset(); return 2;
}
void DivPlatform34SI22::quit() { for (auto* o:osc) delete o; }
