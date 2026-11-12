#include "../../src/engine/engine.h"
#include "../../src/engine/platform/34si22.h"
#include "../../src/engine/platform/pcmdac.h"
#include <cassert>
#include <cstdio>
#ifdef main
#undef main
#endif
void reportError(String what) { fprintf(stderr,"%s\n",what.c_str()); }
int main() {
  static DivEngine engine;
  assert(!engine.load(new unsigned char[32](),32));
  engine.song.systemLen=1; engine.song.system[0]=DIV_SYSTEM_34SI22;
  engine.song.systemChans[0]=engine.song.chans=2;
  auto* sample=new DivSample; assert(sample->init(4));
  sample->data16[0]=0x1234; sample->data16[1]=-1200; sample->data16[2]=2000; sample->data16[3]=0;
  sample->depth=DIV_SAMPLE_DEPTH_16BIT; sample->centerRate=48000;
  engine.song.sample.push_back(sample); engine.song.sampleLen=1;
  engine.renderSamples();
  assert(sample->dataSQV4H==nullptr && sample->dataSQV4L==nullptr && sample->dataSQV4LPlus==nullptr);
  auto* instrument=new DivInstrument; instrument->type=DIV_INS_34SI22; instrument->amiga.initSample=0;
  engine.song.ins.push_back(instrument); engine.song.insLen=1;
  assert(instrument->save("instrument.fui",nullptr,true));
  auto imported=engine.instrumentFromFile("instrument.fui",false,true);
  assert(imported.size()==1 && imported[0]->type==DIV_INS_34SI22 && imported[0]->amiga.initSample==0);
  delete imported[0]; remove("instrument.fui");
  static DivPlatform34SI22 dispatch; dispatch.init(&engine,2,48000,DivConfig()); dispatch.renderSamples(0);
  assert(dispatch.isSampleLoaded(0,0)); assert(dispatch.getSampleMemUsage(0)==8);
  assert(dispatch.getRegisterPoolSize()==0x29);
  for (unsigned a=0x25;a<=0x28;a++) { dispatch.poke(a,0xff); assert(dispatch.getRegisterPool()[a]==0); }
  dispatch.dispatch(DivCommand(DIV_CMD_INSTRUMENT,0,0)); dispatch.dispatch(DivCommand(DIV_CMD_NOTE_ON,0,48)); dispatch.tick();
  auto* regs=dispatch.getRegisterPool(); assert(regs[0x21]==1); assert(regs[0x0e]==0); assert((regs[0x0c]&3)==1);
  short left[4],right[4]; short* output[]={left,right}; dispatch.acquire(output,4);
  assert(left[0]==0x1234 && left[1]==-1200 && left[2]==2000 && left[3]==0);
  for (int i=0;i<4;i++) assert(left[i]==right[i]);
  // Macro bit 0 is right and bit 1 is left, matching Furnace panBits.
  instrument->std.panLMacro.len=1;
  for (int pan=0; pan<4; ++pan) {
    instrument->std.panLMacro.val[0]=pan;
    dispatch.dispatch(DivCommand(DIV_CMD_NOTE_ON,0,48)); dispatch.tick();
    dispatch.acquire(output,1);
    assert(left[0]==((pan&2)?0x1234:0));
    assert(right[0]==((pan&1)?0x1234:0));
  }
  // Exercise real 34SI22 sample upload and SQV4 playback across 63-sample boundaries.
  instrument->std.panLMacro.len=0;
  auto* compressed=new DivSample; compressed->centerRate=48000; assert(compressed->init(130));
  for (int i=0;i<130;i++) compressed->data16[i]=i*100-6500;
  compressed->convert(DIV_SAMPLE_DEPTH_SQV4L_PLUS,1u<<DIV_SAMPLE_DEPTH_8BIT);
  auto image=compressed->getSQV4(); sqv4::Header sqHeader; std::vector<int16_t> expected;
  assert(sqv4::decode(image.data(),image.size(),expected,sqHeader));
  engine.song.sample.push_back(compressed); engine.song.sampleLen=2;
  instrument->amiga.initSample=1; dispatch.renderSamples(0); dispatch.reset();
  dispatch.dispatch(DivCommand(DIV_CMD_INSTRUMENT,0,0)); dispatch.dispatch(DivCommand(DIV_CMD_NOTE_ON,0,48)); dispatch.tick();
  short sqLeft[130],sqRight[130]; short* sqOut[]={sqLeft,sqRight}; dispatch.acquire(sqOut,130);
  for (int i=0;i<130;i++) { assert(sqLeft[i]==expected[i]); assert(sqRight[i]==sqLeft[i]); }
  engine.song.sample.pop_back(); engine.song.sampleLen=1; delete compressed; instrument->amiga.initSample=0;
  dispatch.renderSamples(0); dispatch.reset(); dispatch.dispatch(DivCommand(DIV_CMD_INSTRUMENT,0,0));
  instrument->std.panLMacro.len=0;
  dispatch.dispatch(DivCommand(DIV_CMD_NOTE_ON,0,48)); dispatch.tick();
  assert(dispatch.getRegisterPool()[0x0e]==0);
  sample->centerRate=1000;
  dispatch.dispatch(DivCommand(DIV_CMD_NOTE_ON,0,48)); dispatch.tick();
  assert(dispatch.getRegisterPool()[0x0e]==15); // automatic pitch clamps instead of wrapping
  sample->centerRate=48000;
  dispatch.dispatch(DivCommand(DIV_CMD_SAMPLE_FREQ,0,0x10)); dispatch.tick();
  assert(dispatch.getRegisterPool()[0x0e]==0); // reserved high bits ignored
  dispatch.dispatch(DivCommand(DIV_CMD_SAMPLE_FREQ,0,255)); dispatch.tick();
  assert(dispatch.getRegisterPool()[0x0e]==15);
  dispatch.dispatch(DivCommand(DIV_CMD_PITCH,0,100)); dispatch.tick();
  assert(dispatch.getRegisterPool()[0x0e]==15);
  dispatch.dispatch(DivCommand(DIV_CMD_NOTE_ON,0,48)); dispatch.tick();
  assert(dispatch.getRegisterPool()[0x0e]!=15);
  instrument->std.dutyMacro.len=1; instrument->std.dutyMacro.val[0]=2;
  dispatch.dispatch(DivCommand(DIV_CMD_NOTE_ON,0,48));
  dispatch.dispatch(DivCommand(DIV_CMD_SAMPLE_FREQ,0,1)); dispatch.tick();
  assert(dispatch.getRegisterPool()[0x0e]==1); // effect wins over macro on the same tick
  instrument->std.dutyMacro.len=2; instrument->std.dutyMacro.val[1]=4;
  dispatch.dispatch(DivCommand(DIV_CMD_INSTRUMENT,1,0));
  dispatch.dispatch(DivCommand(DIV_CMD_NOTE_ON,1,48)); dispatch.tick();
  assert(dispatch.getRegisterPool()[0x0f]==2);
  dispatch.tick(); assert(dispatch.getRegisterPool()[0x0f]==4);
  assert(dispatch.getRegisterPool()[0x0e]==1);
  auto* def=engine.getSystemDef(DIV_SYSTEM_34SI22);
  assert(def->postEffectHandlers.count(0x20) && def->postEffectHandlers.at(0x20).dispatchCmd==DIV_CMD_SAMPLE_FREQ);
  auto* large=new DivSample; assert(large->init(524288)); large->depth=DIV_SAMPLE_DEPTH_16BIT;
  engine.song.sample.push_back(large); engine.song.sampleLen=2;
  dispatch.renderSamples(0); assert(dispatch.getSampleMemCapacity(0)==0x100000); assert(!dispatch.isSampleLoaded(0,1));
  DivConfig flags; flags.set("largeSampleMemory",1); dispatch.setFlags(flags); dispatch.renderSamples(0);
  assert(dispatch.getSampleMemCapacity(0)==0x1000000); assert(dispatch.isSampleLoaded(0,1));
  assert(dispatch.getSampleMemUsage(0)==0x100008);
  dispatch.setFlags(DivConfig()); dispatch.renderSamples(0);
  assert(dispatch.getSampleMemCapacity(0)==0x100000 && !dispatch.isSampleLoaded(0,1));
  engine.song.sample.pop_back(); engine.song.sampleLen=1; delete large;
  engine.song.systemFlags[0].set("largeSampleMemory",1);
  dispatch.quit();
  {
    // 34SI22A: a Generic Sample (AMIGA) instrument streams 16-bit PCM through the FIFO DAC,
    // resampled and scaled by the platform; 34SI22 instruments still play from memory.
    DivConfig a; a.set("chipType",1);
    static DivPlatform34SI22 fifo; fifo.init(&engine,2,48000,a); fifo.renderSamples(0);
    assert(fifo.getRegisterPoolSize()==0x2d);
    auto* generic=new DivInstrument; generic->type=DIV_INS_AMIGA; generic->amiga.initSample=0;
    engine.song.ins.push_back(generic); engine.song.insLen=2;
    sample->centerRate=48000;
    fifo.dispatch(DivCommand(DIV_CMD_INSTRUMENT,0,1)); fifo.dispatch(DivCommand(DIV_CMD_NOTE_ON,0,DIV_NOTE_RAW_FLAG|65536)); fifo.tick();
    auto* r=fifo.getRegisterPool();
    assert((r[0x21]&7)==5 && r[0x0e]==0 && (r[0x0d]&15)==15 && (r[0x0c]&3)==1);
    short fl[6],fr[6]; short* fo[]={fl,fr}; fifo.acquire(fo,6);
    const short plain[6]={0x1234,-1200,2000,0,0,0};
    for (int i=0;i<6;i++) assert(fl[i]==plain[i] && fr[i]==fl[i]);
    fifo.dispatch(DivCommand(DIV_CMD_VOLUME,0,5)); fifo.dispatch(DivCommand(DIV_CMD_NOTE_ON,0,DIV_NOTE_RAW_FLAG|65536)); fifo.tick();
    assert((fifo.getRegisterPool()[0x0d]&15)==15); // the chip stays at 0 dB
    fifo.acquire(fo,2); assert(fl[0]==0x1234*5/15 && fl[1]==-1200*5/15);
    fifo.dispatch(DivCommand(DIV_CMD_VOLUME,0,15));
    sample->centerRate=24000; fifo.notifyPitchTable(); // raw notes ignore center rate
    fifo.dispatch(DivCommand(DIV_CMD_NOTE_ON,0,DIV_NOTE_RAW_FLAG|32768)); fifo.tick(); fifo.acquire(fo,6);
    const short half[6]={0x1234,0x1234,-1200,-1200,2000,2000};
    for (int i=0;i<6;i++) assert(fl[i]==half[i]);
    sample->centerRate=48000; sample->loop=true; sample->loopMode=DIV_SAMPLE_LOOP_FORWARD;
    sample->loopStart=1; sample->loopEnd=3;
    fifo.dispatch(DivCommand(DIV_CMD_NOTE_ON,0,DIV_NOTE_RAW_FLAG|65536)); fifo.tick(); fifo.acquire(fo,6);
    const short looped[6]={0x1234,-1200,2000,-1200,2000,-1200};
    for (int i=0;i<6;i++) assert(fl[i]==looped[i]);
    fifo.muteChannel(0,true); fifo.acquire(fo,1); assert(fl[0]==0); fifo.muteChannel(0,false);
    fifo.dispatch(DivCommand(DIV_CMD_NOTE_OFF,0,0)); fifo.tick();
    assert((fifo.getRegisterPool()[0x0c]&3)==0 && fifo.getRegisterPool()[0x29]==0 && fifo.getRegisterPool()[0x2a]==0);
    fifo.acquire(fo,1); assert(fl[0]==0);
    sample->loop=false;
    fifo.dispatch(DivCommand(DIV_CMD_INSTRUMENT,0,0)); fifo.dispatch(DivCommand(DIV_CMD_NOTE_ON,0,48)); fifo.tick();
    assert((fifo.getRegisterPool()[0x21]&7)==1); // DIV_INS_34SI22: memory playback (PCM16)
    // Compare actual dispatches at the same output rate, with every pitch mode.
    static DivPlatformPCMDAC pcm;
    DivConfig pcmFlags; pcmFlags.set("rate",48000);
    pcm.init(&engine,1,48000,pcmFlags);
    auto both=[&](DivCommand command) { fifo.dispatch(command); pcm.dispatch(command); };
    auto tickBoth=[&]() {
      fifo.tick(); pcm.tick();
      auto* f=fifo.getChanState(0); auto* p=pcm.getChanState(0);
      assert(f->baseFreq==p->baseFreq && f->freq==p->freq);
      assert(f->pitch2==p->pitch2 && f->arpOff==p->arpOff && f->fixedArp==p->fixedArp);
    };
    both(DivCommand(DIV_CMD_INSTRUMENT,0,1));
    for (float tuning: {440.0f,480.0f}) for (int linear: {0,1,2}) for (bool oldArp: {false,true}) for (bool oldCenter: {false,true}) {
      engine.song.tuning=tuning;
      engine.song.compatFlags.linearPitch=linear;
      engine.song.compatFlags.oldArpStrategy=oldArp;
      engine.song.compatFlags.oldCenterRate=oldCenter;
      for (int center: {0,8363,22050,48000}) {
        sample->centerRate=center;
        fifo.notifyPitchTable(); pcm.notifyPitchTable();
        for (int note: {12,48,69,108,int(DIV_NOTE_RAW_FLAG|32768),int(DIV_NOTE_RAW_FLAG|65536)}) {
          generic->std.arpMacro.len=0; generic->std.pitchMacro.len=0;
          both(DivCommand(DIV_CMD_PITCH,0,0));
          both(DivCommand(DIV_CMD_NOTE_ON,0,note)); tickBoth();
          for (int pitch: {-128,0,128}) { both(DivCommand(DIV_CMD_PITCH,0,pitch)); tickBoth(); }
        }
      }
      for (bool relative: {false,true}) {
        generic->std.arpMacro.len=3;
        generic->std.arpMacro.val[0]=0; generic->std.arpMacro.val[1]=7; generic->std.arpMacro.val[2]=0x40000000|12;
        generic->std.pitchMacro.len=3; generic->std.pitchMacro.mode=relative;
        generic->std.pitchMacro.val[0]=64; generic->std.pitchMacro.val[1]=-32; generic->std.pitchMacro.val[2]=128;
        generic->amiga.useNoteMap=true;
        generic->amiga.noteMap[48].map=0; generic->amiga.noteMap[48].freq=55;
        both(DivCommand(DIV_CMD_NOTE_ON,0,48)); tickBoth();
        both(DivCommand(DIV_CMD_LEGATO,0,50)); tickBoth();
        tickBoth();
        engine.song.compatFlags.resetMacroOnPorta=true;
        both(DivCommand(DIV_CMD_PRE_PORTA,0,1,1));
        for (int i=0;i<8;i++) { both(DivCommand(DIV_CMD_NOTE_PORTA,0,64,60)); tickBoth(); }
        both(DivCommand(DIV_CMD_NOTE_ON,0,DIV_NOTE_NULL)); tickBoth();
        generic->amiga.useNoteMap=false;
      }
    }
    generic->std.arpMacro.len=0; generic->std.pitchMacro.len=0;
    // Check the calculated step reaches the FIFO resampler in rendered audio.
    sample->centerRate=22050; sample->loop=true; sample->loopStart=0; sample->loopEnd=4;
    sample->loopMode=DIV_SAMPLE_LOOP_FORWARD;
    fifo.notifyPitchTable(); pcm.notifyPitchTable();
    both(DivCommand(DIV_CMD_PITCH,0,0));
    both(DivCommand(DIV_CMD_NOTE_ON,0,48)); tickBoth();
    short renderedL[4096],renderedR[4096]; short* rendered[]={renderedL,renderedR};
    fifo.acquire(rendered,4096);
    unsigned phase=0;
    for (int i=0;i<4096;i++) {
      assert(renderedL[i]==sample->data16[(phase>>16)%4]);
      assert(renderedR[i]==renderedL[i]);
      phase+=pcm.getChanState(0)->freq;
    }
    sample->loop=false;
    pcm.quit();
    fifo.quit();
    // On the 34SI22 a Generic Sample instrument also plays from memory.
    static DivPlatform34SI22 base; base.init(&engine,2,48000,DivConfig()); base.renderSamples(0);
    base.dispatch(DivCommand(DIV_CMD_INSTRUMENT,0,1)); base.dispatch(DivCommand(DIV_CMD_NOTE_ON,0,48)); base.tick();
    assert((base.getRegisterPool()[0x21]&7)==1);
    base.quit();
    engine.song.ins.pop_back(); engine.song.insLen=1; delete generic;
    engine.song.systemFlags[0].set("chipType",1);
  }
  auto* saved=engine.saveFur(); assert(saved); size_t size=saved->size(); auto* data=new unsigned char[size]; memcpy(data,saved->getFinalBuf(),size); saved->finish(); delete saved;
  static DivEngine reload; assert(reload.load(data,size)); assert(reload.song.system[0]==DIV_SYSTEM_34SI22); assert(reload.song.systemFlags[0].getBool("largeSampleMemory",false)); assert(reload.song.systemFlags[0].getInt("chipType",0)==1); assert(reload.song.ins[0]->type==DIV_INS_34SI22); assert(reload.song.ins[0]->amiga.initSample==0); assert(reload.song.ins[0]->std.dutyMacro.len==2 && reload.song.ins[0]->std.dutyMacro.val[1]==4);
  puts("34SI22/34SI22A engine playback and .fur round trip passed");
}

