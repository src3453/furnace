#include "sound/3hs/3hs88pwn4/include/hs88/registers.hpp"
#include <vector>
#include "../waveSynth.h"
#include "../dispatch.h"

class DivPlatformS3HS: public DivDispatch {
  struct Channel : public SharedChannel {
    int freq, baseFreq, pitch;
    int wave, sample;
    unsigned short pos;
    bool active, freqChanged;
    bool volumeChanged;
    bool waveChanged, waveUpdated;
    bool pcm;
    unsigned char vol;
    unsigned char actualVol;
    unsigned char outVol, resVol;
    signed char amp;
    bool pcmLoop;
    int hasOffset;
    int pan;
    int modmode;
    int opvols[8];
    int opwaves[8];
    int opfreqs[8];
    int feedback;
    int modmode_o;
    int pcmOffset;
    int directWriteAddr;
    bool opvoloverride[8];
    bool opwaveoverride[8];
    bool opfreqoverride[8];
    bool modmodeoverride;
    bool feedbackoverride;
    DivWaveSynth ws;
    Channel(): 
      SharedChannel(0,true),
      freq(0), 
      baseFreq(0), 
      pitch(0), 
      wave(-1),
      sample(0),
      pos(0), 
      active(false), 
      freqChanged(false),  
      volumeChanged(false),
      waveChanged(false),
      waveUpdated(false),
      pcm(false),
      vol(0),
      outVol(0),
      resVol(0),
      amp(0),
      pcmLoop(false),
      hasOffset(0),
      pan(0),
      pcmOffset(0){}
  };
  bool enableComp;
  bool enableEQ;
  int compThresh;
  int compRatio;
  int compVolume;
  int EQlo;
  int EQmid;
  int EQhi;
  hs88::Engine* chip;
  hs88::RegisterInterface* chipRegs;
  std::vector<int> oscTap[12];
  DivPitchTable pitchTable;
  unsigned char regPool[1024];
  DivMemoryComposition memCompo;
  Channel chan[12];
  DivDispatchOscBuffer* oscBuf[12];
  bool isMuted[12];
  int muteGain[12];
  unsigned char chans;  
  unsigned int sampleMemSize;
  unsigned char ilCtrl, ilSize, fil1;
  unsigned char initIlCtrl, initIlSize, initFil1;
  bool sampleLoaded[256];
  unsigned char* sampleMem;
  size_t sampleMemLen;
  unsigned int sampleOffSU[256];
  friend void putDispatchChip(void*,int);
  friend void putDispatchChan(void*,int,int);
  int volmax=0;
  // Shadow copies of core state that the vendored engine keeps private
  // (envelope stage/level, PCM cursor), replayed from the same register
  // writes and sample counts so the visualizer can show them.
  hs816::Envelope visEnv[8][8];
  hs816::EnvelopeRates visRates[8][8];
  unsigned int visIdle[8][8];
  bool visGate[8];
  unsigned long long visPcmPos[4];
  unsigned int visPcmFrac[4];
  unsigned int visTrigger[12];
  void visReset();
  void visWrite(unsigned int offset, unsigned char data);
  void visAdvance(size_t len);
  public:
    // Read by the GUI without locking; values are for display only.
    const hs88::Engine* getChip() const { return chip; }
    // Envelope stage (hs816::Envelope::Stage) and level (Q16) of an operator.
    int getVisStage(int ch, int op) const { return (int)visEnv[ch][op].stage(); }
    int getVisLevel(int ch, int op) const { return visEnv[ch][op].level(); }
    // PCM/WT cursor in table/PCM samples since the last sync.
    unsigned long long getVisPcmPos(int ch) const { return visPcmPos[ch]; }
    // Increments on every key-on/sync of the channel.
    unsigned int getVisTrigger(int ch) const { return visTrigger[ch]; }
    struct QueuedWrite {
    unsigned char addr;
    unsigned char val;
    QueuedWrite(): addr(0), val(0) {}
    QueuedWrite(unsigned char a, unsigned char v): addr(a), val(v) {}
   };
    unsigned char* getRegisterPool();
    void updateWave(int ch);
    DivChannelModeHints getModeHints(int ch);
    int getRegisterPoolSize();
    void acquire(short** buf, size_t len);
    void muteChannel(int ch, bool mute);
    int dispatch(DivCommand c);
    int NewFunction(DivCommand &c);
    DivMacroInt *getChanMacroInt(int ch);
      SharedChannel* getChanState(int chan);
    DivDispatchOscBuffer* getOscBuffer(int chan);
    void notifyWaveChange(int wave);
    void notifyInsDeletion(void* ins);
    int getOutputCount();
    unsigned short getPan(int ch);
    void reset();
    void tick(bool sysTick=true);
    int init(DivEngine* parent, int channels, int sugRate, const DivConfig& flags);
    void quit();
    void poke(unsigned int addr, unsigned short val);
    void poke(std::vector<DivRegWrite>& wlist);
    const void* getSampleMem(int index);
    size_t getSampleMemCapacity(int index);
    size_t getSampleMemUsage(int index);
    const DivMemoryComposition* getMemCompo(int index);
    bool isSampleLoaded(int index, int sample);
    void renderSamples(int chipID);
    void doWrite(unsigned int addr, unsigned char data);
    void writeOpVolumes(int channel, int channelVolume, int mode, int* opVolumes);
    void loadSampleMemory();
    void setFlags(const DivConfig& flags);
    ~DivPlatformS3HS();
};
