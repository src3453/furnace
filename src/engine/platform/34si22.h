#ifndef FURNACE_PLATFORM_34SI22_H
#define FURNACE_PLATFORM_34SI22_H
#include "../dispatch.h"
#include "sound/34si22/playback.h"
#include <vector>

class DivPlatform34SI22: public DivDispatch, private Si22Playback::Memory {
  struct Channel: SharedChannel {
    int sample=-1; unsigned char pan=3;
    bool directDivider=false, dividerPending=false;
    unsigned char pendingDivider=0;
    // 34SI22A FIFO DAC: an AMIGA (Generic Sample) instrument is resampled and
    // scaled here, like Generic PCM DAC, and streamed as 16-bit PCM at 48 kHz.
    bool fifoMode=false, fifoPlaying=false, audDir=false;
    int audPos=0;
    unsigned audSub=0, step=0;
    Channel(): SharedChannel(15,true) {}
  } chan[2];
  bool variantA=false;
  DivDispatchOscBuffer* osc[2]={};
  bool muted[2]={};
  std::vector<uint8_t> memory;
  unsigned offsets[32768]={}, lengths[32768]={};
  uint8_t formats[32768]={};
  bool loaded[32768]={};
  size_t used=0;
  DivMemoryComposition composition;
  DivPitchTable pitchTable;
  Si22Playback chip;
  uint8_t control=0;
  uint8_t read(uint32_t a) const override { return a<memory.size()?memory[a]:0; }
  void write(unsigned a, unsigned v);
  void controlWrite(unsigned ch, unsigned cmd);
  void updateVolume();
  void updatePan();
  void address(unsigned reg, unsigned value);
  int16_t fifoSample(Channel& c);
public:
  DivPlatform34SI22(): chip(*this) {}
  void acquire(short** buf, size_t len) override;
  int dispatch(DivCommand c) override;
  void tick(bool sysTick=true) override;
  void reset() override;
  int init(DivEngine* p,int channels,int sugRate,const DivConfig& flags) override;
  void quit() override;
  void setFlags(const DivConfig&) override;
  void notifyPitchTable(int sample=-1) override;
  void notifyInsDeletion(void* ins) override;
  void muteChannel(int ch,bool mute) override;
  int getOutputCount() override { return 2; }
  SharedChannel* getChanState(int ch) override { return &chan[ch]; }
  DivMacroInt* getChanMacroInt(int ch) override { return &chan[ch].std; }
  DivDispatchOscBuffer* getOscBuffer(int ch) override { return osc[ch]; }
  unsigned char* getRegisterPool() override;
  int getRegisterPoolSize() override { return variantA?0x2d:0x29; }
  void poke(unsigned a,unsigned short v) override { write(a,v); }
  void poke(std::vector<DivRegWrite>& writes) override;
  const void* getSampleMem(int index) override { return index==0?memory.data():nullptr; }
  size_t getSampleMemCapacity(int index) override { return index==0?memory.size():0; }
  size_t getSampleMemUsage(int index) override { return index==0?used:0; }
  bool isSampleLoaded(int index,int sample) override { return index==0 && sample>=0 && sample<32768 && loaded[sample]; }
  const DivMemoryComposition* getMemCompo(int index) override { return index==0?&composition:nullptr; }
  void renderSamples(int sysID) override;
private:
  unsigned char registerPool[0x2d]={};
};
#endif
