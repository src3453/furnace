#include "s3hs.h"
#include "../engine.h"
#include <stdio.h>
#include <math.h>
#include "furIcons.h"
#include "IconsFontAwesome4.h"
#include <chrono>

#define CHIP_FREQBASE 48000
#define rWrite(a,v) if (!skipRegisterWrites) {if(a>=0x400000 && a<0x400400) {doWrite(a,v); regPool[((a-0x400000)%0x400)]=v;} if (dumpWrites && false) {addWrite(a,v);} }
#define rRead(a) cpt->ram_peek(cpt->ram,(int)a);

unsigned int s3hs_chanaddrs_freq[12] = {0x400000,0x400040,0x400080,0x4000c0,0x400100,0x400140,0x400180,0x4001c0,0x400200,0x400230,0x400260,0x400290};
unsigned int s3hs_chanaddrs_volume[12] = {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x400202,0x400232,0x400262,0x400292};

void DivPlatformS3HS::doWrite(unsigned int addr, unsigned char data) {
  //printf("S3W: %lld %08x : %02x\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), addr, (int)data);
  cpt->ram_poke(cpt->ram,(int)addr,(Byte)data);
}

void DivPlatformS3HS::updateWave(int ch) {
  if (ch>=8 && ch<=11)
  {
    for (int j=0;j<32;j++) {
        rWrite(0x400210+0x30*(ch-8)+j,chan[ch].ws.output[j]);
    }
  }
}

void DivPlatformS3HS::acquire(short** buf, size_t len) {
  // オシロスコープバッファを初期化
  for (int i=0; i<chans; i++) {
    oscBuf[i]->begin(len);
  }

  // S3HSエミュレータからオーディオデータを取得
  std::vector<std::vector<std::vector<float_t>>> output = cpt->AudioCallBack(len);
  
  // 各サンプルを処理
  for (size_t i=0; i<len; i++) {
    float outL = 0;
    float outR = 0;

    // 各チャンネルを処理
    for (unsigned char j=0; j<chans; j++) {
      // モノラル出力を計算（オシロスコープ用）
      float chanOut = ((output[0][j][i] + output[1][j][i]) / 2);
      if (chanOut < -32768.0f) chanOut = -32768.0f;
      if (chanOut > 32767.0f) chanOut = 32767.0f;

      // オシロスコープバッファにサンプルを格納
      oscBuf[j]->putSample(i, (short)((float)chanOut / 2));
      
      // チャンネル位置を更新
      chan[j].pos += chan[j].freq;
    }
    
    // マスター出力（チャンネル12）を取得
    outL = (int)output[0][12][i];
    outR = (int)output[1][12][i];
    
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
  for (int i=0; i<chans; i++) {
    oscBuf[i]->end(len);
  }
}

unsigned short DivPlatformS3HS::getPan(int ch) {
  return (chan[ch].pan>>4)<<8|(chan[ch].pan&0xf);
}

int DivPlatformS3HS::getOutputCount() {
  return 2;
}

void DivPlatformS3HS::muteChannel(int ch, bool mute) {
  isMuted[ch]=mute;
  rWrite(0x4002D0+ch,mute?1:0);
}

void DivPlatformS3HS::setFlags(const DivConfig& flags) {
  enableComp = flags.getBool("enableComp",false);
  enableEQ = flags.getBool("enableEQ",false);
  compThresh = flags.getInt("compThresh",0);
  compRatio = flags.getInt("compRatio",0);
  compVolume = flags.getInt("compVolume",32);
  EQlo = flags.getInt("EQlo",0);
  EQmid = flags.getInt("EQmid",0);
  EQhi = flags.getInt("EQhi",0);
  rWrite(0x4002C0,enableComp?1:0);
  rWrite(0x4002C1,enableEQ?1:0);
  rWrite(0x4002C2,compThresh);
  rWrite(0x4002C3,compRatio);
  rWrite(0x4002C4,compVolume);
  rWrite(0x4002C5,EQlo);
  rWrite(0x4002C6,EQmid);
  rWrite(0x4002C7,EQhi);
}

void DivPlatformS3HS::tick(bool sysTick) {
  rWrite(0x4002C0,enableComp?1:0);
  rWrite(0x4002C1,enableEQ?1:0);
  rWrite(0x4002C2,compThresh);
  rWrite(0x4002C3,compRatio);
  rWrite(0x4002C4,compVolume);
  rWrite(0x4002C5,EQlo);
  rWrite(0x4002C6,EQmid);
  rWrite(0x4002C7,EQhi);
  for (unsigned char i=0; i<chans; i++) {
    chan[i].std.next();
    
    if (chan[i].std.vol.had) {
      int volume = chan[i].std.vol.val*chan[i].vol/256;
      if(i>=8 && i<=11) {
        chan[i].outVol=(MIN(255,volume)*(chan[i].vol))/(256);
        if (chan[i].outVol<0) chan[i].outVol=0;
        if (chan[i].outVol>255) chan[i].outVol=255;
      } else {
        chan[i].outVol=(MIN(255,volume)*parent->getIns(chan[i].ins,DIV_INS_S3HS)->cpt.op1v)/(256);
        if (chan[i].outVol<0) chan[i].outVol=0;
        if (chan[i].outVol>255) chan[i].outVol=255;
      }
      if (chan[i].resVol!=chan[i].outVol) {
        chan[i].resVol=chan[i].outVol;
        if (!isMuted[i]) {
          chan[i].volumeChanged=true;
        }
      }
      if (i < 8) {
        rWrite(0x400010 + 64*i,(unsigned char)((float)(chan[i].opvols[0]*volume)/256));
        rWrite(0x400011 + 64*i,(unsigned char)((float)(chan[i].opvols[1]*(chan[i].modmode>3?255:volume))/256));
        rWrite(0x400012 + 64*i,(unsigned char)((float)(chan[i].opvols[2]*(chan[i].modmode>2?255:volume))/256));
        rWrite(0x400013 + 64*i,(unsigned char)((float)(chan[i].opvols[3]*(chan[i].modmode>2?255:volume))/256));
        rWrite(0x400014 + 64*i,(unsigned char)((float)(chan[i].opvols[4]*(chan[i].modmode>3||chan[i].modmode==1?255:volume))/256));
        rWrite(0x400015 + 64*i,(unsigned char)((float)(chan[i].opvols[5]*(chan[i].modmode>3||chan[i].modmode==1?255:volume))/256));
        rWrite(0x400016 + 64*i,(unsigned char)((float)(chan[i].opvols[6]*(chan[i].modmode>3||chan[i].modmode==1?255:volume))/256));
        rWrite(0x400017 + 64*i,(unsigned char)((float)(chan[i].opvols[7]*(chan[i].modmode>3||chan[i].modmode==1?255:volume))/256));
      }
    }
    if (NEW_ARP_STRAT) {
      chan[i].handleArp();
    } else if (chan[i].std.arp.had) {
      if (!chan[i].inPorta) {
        chan[i].baseFreq=NOTE_FREQUENCY(parent->calcArp(chan[i].note,chan[i].std.arp.val));
      }
      chan[i].freqChanged=true;
    }
    if (chan[i].std.wave.had) {
      if (chan[i].wave!=chan[i].std.wave.val || chan[i].ws.activeChanged()) {
        chan[i].wave=chan[i].std.wave.val;
        chan[i].ws.changeWave1(chan[i].wave);
        chan[i].waveUpdated=true;
        
      }
    }
    if (chan[i].std.pitch.had) {
      if (chan[i].std.pitch.mode) {
        chan[i].pitch+=chan[i].std.pitch.val;
        CLAMP_VAR(chan[i].pitch,-32768,32767);
      } else {
        chan[i].pitch=chan[i].std.pitch.val;
      }
      chan[i].freqChanged=true;
    }
    if (chan[i].std.duty.had) {
      if (i>=8 && i<=11) {
        rWrite(0x400203+48*(i-8),chan[i].std.duty.val+1);
      }
    }
    if (chan[i].std.ex1.had) {
      if (i>=8 && i<=11) {
        rWrite(0x400204+(i-8)*0x30,chan[i].std.ex1.val);
      }
    }
    if (chan[i].std.ex2.had) {
      if (i>=8 && i<=11) {
        rWrite(0x400205+(i-8)*0x30,chan[i].std.ex2.val);
      }
    }
    if (chan[i].std.ex3.had) {
      if (i>=8 && i<=11) {
        rWrite(0x400206+(i-8)*0x30,chan[i].std.ex3.val);
      }
    }
    if (chan[i].std.panL.had) {
      if (i>=0 && i<=7) {
        int panRegV = rRead(0x40001d+i*0x40);
        //std::cout << std::hex << panRegV << std::endl;
        panRegV = chan[i].std.panL.val*16+panRegV%16;
        //std::cout << std::hex << panRegV << std::endl;
        rWrite(0x40001d+i*0x40,panRegV);
      } else {
        int panRegV = rRead(0x400207+(i-8)*0x30);
        panRegV = chan[i].std.panL.val*16+panRegV%16;
        rWrite(0x400207+(i-8)*0x30,panRegV);
      }
    }
    if (chan[i].std.panR.had) {
      if (i>=0 && i<=7) {
        int panRegV = rRead(0x40001d+i*0x40);
        //std::cout << std::hex << panRegV << std::endl;
        panRegV = chan[i].std.panR.val+panRegV/16*16;
        //std::cout << std::hex << panRegV << std::endl;
        rWrite(0x40001d+i*0x40,panRegV);
      } else {
        int panRegV = rRead(0x400207+(i-8)*0x30);
        panRegV = chan[i].std.panR.val+panRegV/16*16;
        rWrite(0x400207+(i-8)*0x30,panRegV);
      }
    }
    if (chan[i].std.phaseReset.had) {
      if (i>=0 && i<=7) {
        if (chan[i].std.phaseReset.val == 1){
          cpt->resetGate(i);
        }
      } else {
        if (chan[i].std.phaseReset.val == 1){
          cpt->wtSync(i-8);
        }
      }
    }
    if (chan[i].waveUpdated) {
      if (i>=8 && i<=11 && !chan[i].pcm)
      {
        for (int j=0;j<32;j++) {
        rWrite(0x400210+0x30*(i-8)+j,chan[i].ws.output[j]);
        }
      }
      if (chan[i].active) {
        if (chan[i].ws.tick()) {
          updateWave(i);
        }
        //chan[i].freqChanged=true;
      }
      chan[i].waveChanged=false;
    }
    if (chan[i].freqChanged || chan[i].keyOn || chan[i].keyOff) {
      //DivInstrument* ins=parent->getIns(chan[i].ins,DIV_INS_SU);
      chan[i].freq=parent->calcFreq(chan[i].baseFreq,chan[i].pitch,chan[i].fixedArp?chan[i].baseNoteOverride:chan[i].arpOff,chan[i].fixedArp,0,2,chan[i].pitch2,chipClock,CHIP_FREQBASE)/16;
      if (chan[i].pcm && i>=8 && i<=11) {
        DivSample* sample=parent->getSample(chan[i].sample);
        if (sample!=NULL) {
          float off=0.5;
          if (sample->centerRate<1) {
            off=0.5;
          } else {
            off=(float)sample->centerRate/(8363.0*2);
          }
          chan[i].freq=(float)chan[i].freq*off;
        }
        chan[i].freqChanged=false;
      } else {
          chan[i].freqChanged=false;
      }
      if (chan[i].freq<0) chan[i].freq=0;
      if (chan[i].freq>65535) chan[i].freq=65535;
      rWrite(s3hs_chanaddrs_freq[i],chan[i].freq>>8);
      rWrite(s3hs_chanaddrs_freq[i]+1,chan[i].freq&0xff);
    }
    if (chan[i].keyOn) {
      if (i>=8) {
        if (chan[i].pcm) {
          int sNum=chan[i].sample;
          DivSample* sample=parent->getSample(sNum);
          if (sample!=NULL && sNum>=0 && sNum<parent->song.sampleLen) {
            unsigned int sampleEnd=sampleOffSU[sNum]+(sample->getLoopEndPosition());
            unsigned int off=sampleOffSU[sNum]+chan[i].pcmOffset;
            chan[i].hasOffset=0;
            if (sampleEnd>=getSampleMemCapacity(0)) sampleEnd=getSampleMemCapacity(0)-1;
            rWrite(0x400210+(48*(i-8)),(off>>16));
            rWrite(0x400211+(48*(i-8)),off>>8&0xff);
            rWrite(0x400212+(48*(i-8)),off&0xff);
            rWrite(0x400213+(48*(i-8)),(sampleEnd>>16));
            rWrite(0x400214+(48*(i-8)),sampleEnd>>8&0xff);
            rWrite(0x400215+(48*(i-8)),sampleEnd&0xff);
            rWrite(0x400216+(48*(i-8)),255);
            rWrite(0x400217+(48*(i-8)),255);
            rWrite(0x400218+(48*(i-8)),255);
            if (sample->isLoopable()) {
              unsigned int sampleLoop=sampleOffSU[sNum]+sample->getLoopStartPosition();
              unsigned int sampleLoopEnd=sampleOffSU[sNum]+sample->getLoopEndPosition();
              if (sampleLoop>=getSampleMemCapacity(0)) sampleLoop=getSampleMemCapacity(0)-1;
              rWrite(0x400216+(48*(i-8)),(sampleLoop>>16));
              rWrite(0x400217+(48*(i-8)),sampleLoop>>8&0xff);
              rWrite(0x400218+(48*(i-8)),sampleLoop&0xff);
              rWrite(0x400213+(48*(i-8)),(sampleLoopEnd>>16));
              rWrite(0x400214+(48*(i-8)),sampleLoopEnd>>8&0xff);
              rWrite(0x400215+(48*(i-8)),sampleLoopEnd&0xff);
              chan[i].pcmLoop=true;
            } else {
              chan[i].pcmLoop=false;
            }
          }
        }
        rWrite(0x400204+(i-8)*0x30,chan[i].std.ex1.val);
      }
    }
    if (chan[i].pcm && i>=8 && i<=11) {
      rWrite(0x400203+48*(i-8),0);
    } else {
      //rWrite(0x400203+48*(i-8),chan[i].std.duty.val+1);
    }
    if (chan[i].keyOff) {
      rWrite(0x40001e + 64*i,0);
      chan[i].keyOff=false;
    }
    if (chan[i].volumeChanged){
      if (i<8) {
        rWrite(0x400010 + 64*i,(unsigned char)((float)(chan[i].opvols[0]*chan[i].outVol)/256));
        rWrite(0x400011 + 64*i,(unsigned char)((float)(chan[i].opvols[1]*(chan[i].modmode>3?255:chan[i].outVol))/256));
        rWrite(0x400012 + 64*i,(unsigned char)((float)(chan[i].opvols[2]*(chan[i].modmode>2?255:chan[i].outVol))/256));
        rWrite(0x400013 + 64*i,(unsigned char)((float)(chan[i].opvols[3]*(chan[i].modmode>2?255:chan[i].outVol))/256));
        rWrite(0x400014 + 64*i,(unsigned char)((float)(chan[i].opvols[4]*(chan[i].modmode>3||chan[i].modmode==1?255:chan[i].outVol))/256));
        rWrite(0x400015 + 64*i,(unsigned char)((float)(chan[i].opvols[5]*(chan[i].modmode>3||chan[i].modmode==1?255:chan[i].outVol))/256));
        rWrite(0x400016 + 64*i,(unsigned char)((float)(chan[i].opvols[6]*(chan[i].modmode>3||chan[i].modmode==1?255:chan[i].outVol))/256));
        rWrite(0x400017 + 64*i,(unsigned char)((float)(chan[i].opvols[7]*(chan[i].modmode>3||chan[i].modmode==1?255:chan[i].outVol))/256));
      } else {
        rWrite(s3hs_chanaddrs_volume[i],chan[i].outVol);
      }
    }
    for (int j=0;j<8;j++) {
      if (chan[i].opvoloverride[j]) {
        rWrite(0x400010 + 64*i+j,chan[i].opvols[j]);
        chan[i].opvoloverride[j]=false;
      }
    }
  }
}

void* DivPlatformS3HS::getChanState(int ch) {
  return &chan[ch];
}

DivDispatchOscBuffer* DivPlatformS3HS::getOscBuffer(int ch) {
  return oscBuf[ch];
}

unsigned char* DivPlatformS3HS::getRegisterPool() {
  return regPool;
}

int DivPlatformS3HS::getRegisterPoolSize() {
  return 1024;
}

int DivPlatformS3HS::dispatch(DivCommand c) {
  switch (c.cmd) {
    case DIV_CMD_NOTE_ON:
      {
      DivInstrument* ins = parent->getIns(chan[c.chan].ins,DIV_INS_S3HS);
      if (chan[c.chan].pcm && !(ins->type==DIV_INS_AMIGA || ins->amiga.useSample)) {
        chan[c.chan].pcm=(ins->type==DIV_INS_AMIGA || ins->amiga.useSample);
      }
      chan[c.chan].pcm=(ins->type==DIV_INS_AMIGA || ins->amiga.useSample);
      if (chan[c.chan].pcm) {
        if (c.value!=DIV_NOTE_NULL) {
          chan[c.chan].sample=ins->amiga.getSample(c.value);
          chan[c.chan].baseFreq=NOTE_FREQUENCY(c.value);
          chan[c.chan].freqChanged=true;
        }
      } else {
        if (c.value!=DIV_NOTE_NULL) {
          chan[c.chan].baseFreq=NOTE_FREQUENCY(c.value);
          chan[c.chan].freqChanged=true;
        }
      }
      chan[c.chan].note = c.value;
      chan[c.chan].macroInit(ins);
      if (c.chan<8)
      {
        cpt->resetGate(c.chan);
        rWrite(0x400002 + 64*c.chan,ins->s3hs.op2fu);
        rWrite(0x400003 + 64*c.chan,ins->s3hs.op2fl);
        rWrite(0x400004 + 64*c.chan,ins->s3hs.op3fu);
        rWrite(0x400005 + 64*c.chan,ins->s3hs.op3fl);
        rWrite(0x400006 + 64*c.chan,ins->s3hs.op4fu);
        rWrite(0x400007 + 64*c.chan,ins->s3hs.op4fl);
        rWrite(0x400008 + 64*c.chan,ins->s3hs.op5fu);
        rWrite(0x400009 + 64*c.chan,ins->s3hs.op5fl);
        rWrite(0x40000a + 64*c.chan,ins->s3hs.op6fu);
        rWrite(0x40000b + 64*c.chan,ins->s3hs.op6fl);
        rWrite(0x40000c + 64*c.chan,ins->s3hs.op7fu);
        rWrite(0x40000d + 64*c.chan,ins->s3hs.op7fl);
        rWrite(0x40000e + 64*c.chan,ins->s3hs.op8fu);
        rWrite(0x40000f + 64*c.chan,ins->s3hs.op8fl);
        rWrite(0x400010 + 64*c.chan,(unsigned char)((float)(ins->s3hs.op1v*chan[c.chan].vol)/256));
        rWrite(0x400011 + 64*c.chan,(unsigned char)((float)(ins->s3hs.op2v*(ins->s3hs.mode>3?255:chan[c.chan].vol))/256));
        rWrite(0x400012 + 64*c.chan,(unsigned char)((float)(ins->s3hs.op3v*(ins->s3hs.mode>2?255:chan[c.chan].vol))/256));
        rWrite(0x400013 + 64*c.chan,(unsigned char)((float)(ins->s3hs.op4v*(ins->s3hs.mode>2?255:chan[c.chan].vol))/256));
        rWrite(0x400014 + 64*c.chan,(unsigned char)((float)(ins->s3hs.op5v*(ins->s3hs.mode>3||ins->s3hs.mode==1?255:chan[c.chan].vol))/256));
        rWrite(0x400015 + 64*c.chan,(unsigned char)((float)(ins->s3hs.op6v*(ins->s3hs.mode>3||ins->s3hs.mode==1?255:chan[c.chan].vol))/256));
        rWrite(0x400016 + 64*c.chan,(unsigned char)((float)(ins->s3hs.op7v*(ins->s3hs.mode>3||ins->s3hs.mode==1?255:chan[c.chan].vol))/256));
        rWrite(0x400017 + 64*c.chan,(unsigned char)((float)(ins->s3hs.op8v*(ins->s3hs.mode>3||ins->s3hs.mode==1?255:chan[c.chan].vol))/256));
        rWrite(0x400018 + 64*c.chan,ins->s3hs.op1w*16+ins->s3hs.op2w);
        rWrite(0x400019 + 64*c.chan,ins->s3hs.op3w*16+ins->s3hs.op4w);
        rWrite(0x40001a + 64*c.chan,ins->s3hs.op5w*16+ins->s3hs.op6w);
        rWrite(0x40001b + 64*c.chan,ins->s3hs.op7w*16+ins->s3hs.op8w);
        rWrite(0x400020 + 64*c.chan,ins->s3hs.op1a);
        rWrite(0x400021 + 64*c.chan,ins->s3hs.op1d);
        rWrite(0x400022 + 64*c.chan,ins->s3hs.op1s);
        rWrite(0x400023 + 64*c.chan,ins->s3hs.op1r);
        rWrite(0x400024 + 64*c.chan,ins->s3hs.op2a);
        rWrite(0x400025 + 64*c.chan,ins->s3hs.op2d);
        rWrite(0x400026 + 64*c.chan,ins->s3hs.op2s);
        rWrite(0x400027 + 64*c.chan,ins->s3hs.op2r);
        rWrite(0x400028 + 64*c.chan,ins->s3hs.op3a);
        rWrite(0x400029 + 64*c.chan,ins->s3hs.op3d);
        rWrite(0x40002a + 64*c.chan,ins->s3hs.op3s);
        rWrite(0x40002b + 64*c.chan,ins->s3hs.op3r);
        rWrite(0x40002c + 64*c.chan,ins->s3hs.op4a);
        rWrite(0x40002d + 64*c.chan,ins->s3hs.op4d);
        rWrite(0x40002e + 64*c.chan,ins->s3hs.op4s);
        rWrite(0x40002f + 64*c.chan,ins->s3hs.op4r);
        rWrite(0x400030 + 64*c.chan,ins->s3hs.op5a);
        rWrite(0x400031 + 64*c.chan,ins->s3hs.op5d);
        rWrite(0x400032 + 64*c.chan,ins->s3hs.op5s);
        rWrite(0x400033 + 64*c.chan,ins->s3hs.op5r);
        rWrite(0x400034 + 64*c.chan,ins->s3hs.op6a);
        rWrite(0x400035 + 64*c.chan,ins->s3hs.op6d);
        rWrite(0x400036 + 64*c.chan,ins->s3hs.op6s);
        rWrite(0x400037 + 64*c.chan,ins->s3hs.op6r);
        rWrite(0x400038 + 64*c.chan,ins->s3hs.op7a);
        rWrite(0x400039 + 64*c.chan,ins->s3hs.op7d);
        rWrite(0x40003a + 64*c.chan,ins->s3hs.op7s);
        rWrite(0x40003b + 64*c.chan,ins->s3hs.op7r);
        rWrite(0x40003c + 64*c.chan,ins->s3hs.op8a);
        rWrite(0x40003d + 64*c.chan,ins->s3hs.op8d);
        rWrite(0x40003e + 64*c.chan,ins->s3hs.op8s);
        rWrite(0x40003f + 64*c.chan,ins->s3hs.op8r);
        rWrite(0x40001c + 64*c.chan,ins->s3hs.mode);
        rWrite(0x40001f + 64*c.chan,(int)(ins->s3hs.fb)+128);
        chan[c.chan].modmode = ins->s3hs.mode;
        chan[c.chan].opvols[0] = ins->s3hs.op1v;
        chan[c.chan].opvols[1] = ins->s3hs.op2v;
        chan[c.chan].opvols[2] = ins->s3hs.op3v;
        chan[c.chan].opvols[3] = ins->s3hs.op4v;
        chan[c.chan].opvols[4] = ins->s3hs.op5v;
        chan[c.chan].opvols[5] = ins->s3hs.op6v;
        chan[c.chan].opvols[6] = ins->s3hs.op7v;
        chan[c.chan].opvols[7] = ins->s3hs.op8v;
      } else {
      
        cpt->wtSync(c.chan-8);
      }
      if (chan[c.chan].insChanged) {
        if (!parent->song.brokenOutVol && !chan[c.chan].std.vol.will) {
          chan[c.chan].outVol=chan[c.chan].vol;
        }
        if (chan[c.chan].wave<0) {
          chan[c.chan].wave=0;
          chan[c.chan].ws.changeWave1(chan[c.chan].wave);
          chan[c.chan].waveUpdated=true;
        }
        
        chan[c.chan].insChanged=false;
      }
      chan[c.chan].ws.init(ins,32,255,chan[c.chan].insChanged);
      chan[c.chan].active=true;
      chan[c.chan].keyOn=true;
      chan[c.chan].amp=(signed char)255;
      if (c.chan<8)
      {
        rWrite(0x40001e + 64*c.chan,1);
      }  
      break;
      }
    case DIV_CMD_NOTE_OFF:
      chan[c.chan].active=false;
      chan[c.chan].keyOff=true;
      if (c.chan<8)
      {
        rWrite(0x40001e + 64*c.chan,0);
      }  
      break;
    case DIV_CMD_VOLUME:
      chan[c.chan].vol=c.value;
      //chan[c.chan].vol=255;
      if (chan[c.chan].vol>255) chan[c.chan].vol=255;
      if (c.chan < 8 ) {
        rWrite(0x400010 + 64*c.chan,(unsigned char)((float)(chan[c.chan].opvols[0]*chan[c.chan].vol)/256));
        rWrite(0x400011 + 64*c.chan,(unsigned char)((float)(chan[c.chan].opvols[1]*(chan[c.chan].modmode>3?255:chan[c.chan].vol))/256));
        rWrite(0x400012 + 64*c.chan,(unsigned char)((float)(chan[c.chan].opvols[2]*(chan[c.chan].modmode>2?255:chan[c.chan].vol))/256));
        rWrite(0x400013 + 64*c.chan,(unsigned char)((float)(chan[c.chan].opvols[3]*(chan[c.chan].modmode>2?255:chan[c.chan].vol))/256));
        rWrite(0x400014 + 64*c.chan,(unsigned char)((float)(chan[c.chan].opvols[4]*(chan[c.chan].modmode>3||chan[c.chan].modmode==1?255:chan[c.chan].vol))/256));
        rWrite(0x400015 + 64*c.chan,(unsigned char)((float)(chan[c.chan].opvols[5]*(chan[c.chan].modmode>3||chan[c.chan].modmode==1?255:chan[c.chan].vol))/256));
        rWrite(0x400016 + 64*c.chan,(unsigned char)((float)(chan[c.chan].opvols[6]*(chan[c.chan].modmode>3||chan[c.chan].modmode==1?255:chan[c.chan].vol))/256));
        rWrite(0x400017 + 64*c.chan,(unsigned char)((float)(chan[c.chan].opvols[7]*(chan[c.chan].modmode>3||chan[c.chan].modmode==1?255:chan[c.chan].vol))/256));
      } else {
        rWrite(0x400202 + 48*(c.chan-8),chan[c.chan].vol); 
      }
      //printf("vol %d %d\n",c.chan,c.value);
      break;
    case DIV_CMD_GET_VOLUME:
      return chan[c.chan].vol;
      break;
    case DIV_CMD_INSTRUMENT:
      if (chan[c.chan].ins!=c.value || c.value2==1) {
        chan[c.chan].ins=c.value;
        chan[c.chan].insChanged=true;
      }
      break;
    case DIV_CMD_PITCH:
      chan[c.chan].pitch=c.value;
      chan[c.chan].freqChanged=true;
      break;
    case DIV_CMD_NOTE_PORTA: {
      int destFreq=NOTE_FREQUENCY(c.value2);
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
      if (return2) return 2;
      break;
    }
    case DIV_CMD_LEGATO:
      chan[c.chan].baseFreq=NOTE_FREQUENCY(c.value);
      chan[c.chan].freqChanged=true;
      break;
    case DIV_CMD_GET_VOLMAX:
      return 255;
      break;
    case DIV_CMD_WAVE:
      if (chan[c.chan].wave!=c.value) {
        chan[c.chan].wave=c.value;
        chan[c.chan].ws.changeWave1(chan[c.chan].wave);
        chan[c.chan].waveUpdated=true;
      }
      break;
    case DIV_CMD_PANNING: {
      chan[c.chan].pan=(c.value>>4)*16+(c.value2>>4);
      rWrite((c.chan<8)?0x40001d+64*c.chan:0x400207+48*(c.chan-8),(c.value>>4)*16+(c.value2>>4));
      break;
    }
    case DIV_CMD_SAMPLE_POS: {
      chan[c.chan].pcmOffset=c.value;
      break;
    }
    case DIV_CMD_ENV_RELEASE:
      chan[c.chan].std.release();
      break;
    case DIV_CMD_NOTE_OFF_ENV:
      chan[c.chan].keyOff=true;
      chan[c.chan].keyOn=false;
      chan[c.chan].active=false;
      chan[c.chan].std.release();
      break;
    case DIV_CMD_S3HS_OP_VOLUME: {
      chan[c.chan].opvols[c.value]=c.value2;
      rWrite(0x400010 + 64*c.chan+c.value,c.value2);
      chan[c.chan].opvoloverride[c.value]=true;
      //printf("opvol %d %d %d\n",c.chan,c.value,c.value2);
      break;
    }
    case DIV_CMD_S3HS_OP_WAVE: {
      if (c.value % 2 == 0) {
        int waveRegV = rRead(0x400018+c.chan*0x40+c.value/2);
        //std::cout << std::hex << panRegV << std::endl;
        waveRegV = c.value*16+waveRegV%16;
        //std::cout << std::hex << panRegV << std::endl;
        rWrite(0x400018+c.chan*0x40+c.value/2,waveRegV);
      } else {
        int waveRegV = rRead(0x400018+c.chan*0x40+c.value/2);
        //std::cout << std::hex << panRegV << std::endl;
        waveRegV = c.value+waveRegV/16*16;
        //std::cout << std::hex << panRegV << std::endl;
        rWrite(0x400018+c.chan*0x40+c.value/2,waveRegV);
      }
      chan[c.chan].opwaveoverride[c.value]=true;
      break;
    }
    case DIV_CMD_S3HS_MODMODE: {
      rWrite(0x40001c+c.chan*0x40,c.value);
      break;
    }
    case DIV_CMD_S3HS_FEEDBACK: {
      rWrite(0x40001f+c.chan*0x40,c.value);
      break;
    }
    case DIV_CMD_S3HS_OP_FREQ_FU: {
      rWrite(0x400000 + 64*c.chan+c.value*2,c.value2);
      chan[c.chan].opfreqoverride[c.value]=true;
      break;
    }
    case DIV_CMD_S3HS_OP_FREQ_FL: {
      rWrite(0x400001 + 64*c.chan+c.value*2,c.value2);
      chan[c.chan].opfreqoverride[c.value]=true;
      break;
    }
    case DIV_CMD_S3HS_FILTER: {
      if (c.chan >= 8) {
        rWrite(0x400204+48*(c.chan-8)+c.value,c.value2);
      }
      break;
    }
    default:
      break;
  }
  return 1;
}

DivChannelModeHints DivPlatformS3HS::getModeHints(int ch) {
  DivChannelModeHints ret;
  if (ch<8) return ret;
  ret.count=1;
  ret.hint[0]=ICON_FA_VOLUME_UP;
  ret.type[0]=0;

  if (chan[ch].pcm) ret.type[0]=4;
  
  return ret;
}

void DivPlatformS3HS::notifyWaveChange(int wave) {
  for (int i=8; i<12; i++) {
    if (chan[i].wave==wave) {
      chan[i].ws.changeWave1(wave);
    }
    chan[i].waveUpdated=true;
  }
}

void DivPlatformS3HS::notifyInsDeletion(void* ins) {
  for (int i=0; i<12; i++) {
    chan[i].std.notifyInsDeletion((DivInstrument*)ins);
  }
}

DivMacroInt* DivPlatformS3HS::getChanMacroInt(int ch) {
  return &chan[ch].std;
}

void DivPlatformS3HS::poke(unsigned int addr, unsigned short val) {
  rWrite(addr,val);
}

void DivPlatformS3HS::poke(std::vector<DivRegWrite>& wlist) {
  for (DivRegWrite& i: wlist) rWrite(i.addr,i.val);
}

const void* DivPlatformS3HS::getSampleMem(int index) {
  return (index==0)?sampleMem:NULL;
}

size_t DivPlatformS3HS::getSampleMemCapacity(int index) {
  return (index==0)?((sampleMemSize)):0;
}

size_t DivPlatformS3HS::getSampleMemUsage(int index) {
  return (index==0)?sampleMemLen:0;
}

const DivMemoryComposition* DivPlatformS3HS::getMemCompo(int index) {
  if (index!=0) return NULL;
  return &memCompo;
}

bool DivPlatformS3HS::isSampleLoaded(int index, int sample) {
  if (index!=0) return false;
  if (sample<0 || sample>255) return false;
  return sampleLoaded[sample];
}

void DivPlatformS3HS::renderSamples(int sysID) {
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
    cpt->ram_poke(cpt->ram,0x0+i,(Byte)((sampleMem[i]+128)&0xff));
  }

  memCompo.used=sampleMemLen;
  memCompo.capacity=sampleMemSize;
}

void DivPlatformS3HS::reset() {
  cpt->ram_boot(cpt->ram);
  cpt->initSound();
  for (int i=0; i<chans; i++) {
    chan[i]=DivPlatformS3HS::Channel();
    chan[i].vol=0xff;
    chan[i].pan=0xff;
    chan[i].std.setEngine(parent);
    chan[i].ws.setEngine(parent);
    chan[i].ws.init(NULL,32,255,false);
    chan[i].feedbackoverride=false;
    chan[i].modmodeoverride=false;
    for (int j=0; j<8; j++) {
      chan[i].opvoloverride[j]=false;
      chan[i].opfreqoverride[j]=false;
      chan[i].opwaveoverride[j]=false;
    }
  }
  for (int i=0; i<sampleMemSize; i++) {
    cpt->ram_poke(cpt->ram,0x0+i,(Byte)((sampleMem[i]+128)&0xff));
  }
  for (int i=0;i<8;i++) {
    rWrite(0x40001f+64*i,0x80);
  }
  for (int i=8;i<12;i++) {
    rWrite(0x400203+48*(i-8),1);
  }
  for (int i=0;i<12;i++) {
    rWrite(0x4002D0+i,isMuted[i]?1:0);
  }
}

int DivPlatformS3HS::init(DivEngine* p, int channels, int sugRate, const DivConfig& flags) {
  parent=p;
  dumpWrites=false;
  skipRegisterWrites=false;
  for (int i=0; i<channels; i++) {
    isMuted[i]=false;
    if (i<channels) {
      oscBuf[i]=new DivDispatchOscBuffer;
      oscBuf[i]->rate=48000;
    }
  }
  enableComp = false;
  enableEQ = false;
  compThresh = 0;
  compRatio = 0;
  compVolume = 128;
  EQlo = 0;
  EQmid = 0;
  EQhi = 0;
  sampleMemSize=0x400000;// 4096KB (4194304 bytes)
  sampleMem=new unsigned char[sampleMemSize];
  memset(sampleMem,0,sampleMemSize); 
  sampleMemLen=0;
  rate=48000;
  chipClock=48000;
  chans=channels;
  cpt = new S3HS_sound();
  memset(regPool,0,1024);
  reset();
  setFlags(flags);
  return channels;
}

void DivPlatformS3HS::quit() {
  delete[] sampleMem;
  for (int i=0; i<chans; i++) {
    delete oscBuf[i];
  }
  delete cpt;
}

DivPlatformS3HS::~DivPlatformS3HS() {
}
