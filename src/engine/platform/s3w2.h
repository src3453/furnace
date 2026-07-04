#include "sound/3ws8pn/sound.hpp"
#include "../waveSynth.h"
#include "../dispatch.h"
  
class DivPlatformS3W2: public DivDispatch {
  struct Channel: public SharedChannel {
    bool freqInit;
    bool pcm;
    signed short wave;
    signed int sample;
    int hasOffset, pcmOffset;
    signed char waveROM[256] = {0}; // 8 bit signed waveform
    bool pcmLoop;
    DivWaveSynth ws;
    Channel():
      SharedChannel(15,true),
      freqInit(false),
      pcm(false),
      sample(-1),
      wave(-1),
      hasOffset(0),
      pcmOffset(0),
      pcmLoop(false) {}
  };
  Channel chan[8];
  DivDispatchOscBuffer* oscBuf[8];
  bool isMuted[8];
  DivMemoryComposition memCompo;
  unsigned char writeOscBuf;
  unsigned int sampleMemSize;
  unsigned char ilCtrl, ilSize, fil1;
  unsigned char initIlCtrl, initIlSize, initFil1;
  bool sampleLoaded[256];
  unsigned char* sampleMem;
  size_t sampleMemLen;
  unsigned int sampleOffSU[256];
  S3W2_Sound* chip;
  DivPitchTable pitchTable;
  unsigned char regPool[0x900]; // 0x000-0x8FF
  void updateWave(int ch);
  friend void putDispatchChip(void*,int);
  friend void putDispatchChan(void*,int,int);
  public:
    void acquire(short** buf, size_t len);
    int dispatch(DivCommand c);
    SharedChannel* getChanState(int chan);
    DivMacroInt* getChanMacroInt(int ch);
    DivDispatchOscBuffer* getOscBuffer(int chan);
    unsigned char* getRegisterPool();
    int getRegisterPoolSize();
    void reset();
    void tick(bool sysTick=true);
    void muteChannel(int ch, bool mute);
    int getOutputCount();
    bool hasAcquireDirect();
    void notifyWaveChange(int wave);
    void notifyInsDeletion(void* ins);
        const void* getSampleMem(int index);
    size_t getSampleMemCapacity(int index);
    size_t getSampleMemUsage(int index);
    const DivMemoryComposition* getMemCompo(int index);
    bool isSampleLoaded(int index, int sample);
    void renderSamples(int chipID);
    void poke(unsigned int addr, unsigned short val);
    void poke(std::vector<DivRegWrite>& wlist);
    int init(DivEngine* parent, int channels, int sugRate, const DivConfig& flags);
    void quit();
    ~DivPlatformS3W2();
};