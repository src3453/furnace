// 34SI22 playback engine. Memory addresses are 24-bit byte offsets.
// The host owns memory and decides how out-of-range addresses are decoded.
// With setVariantA(true) it is the 34SI22A (SGU v3.2): writes to REC_x_ADC push
// 16-bit samples into a 2048-sample FIFO per playback channel, FIFO_x_COUNT
// (0x29..0x2C) reads the fill level, and PLAY_FORMAT 5 plays the FIFO at
// 48000 / (divider + 1) Hz (CPT32 doc/spec_34SI22.md section 11).
#ifndef FURNACE_34SI22_H
#define FURNACE_34SI22_H
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include "../../../sqv4.h"

class Si22Playback {
public:
  struct Memory {
    virtual ~Memory() = default;
    virtual uint8_t read(uint32_t address) const = 0;
  };
  static constexpr unsigned FIFO_CAPACITY=2048;
  explicit Si22Playback(const Memory& m): memory(m) { reset(); }
  void reset() { reg.fill(0); voice = {}; channelOut = {}; adc = {}; fifo = {}; }
  void setVariantA(bool a) { variantA=a; if (!a) fifo={}; }
  bool isVariantA() const { return variantA; }
  unsigned fifoCount(unsigned ch) const { return fifo[ch].count; }
  uint8_t read(unsigned a) const {
    if (variantA && a >= 0x29 && a < 0x2d) {
      unsigned count=fifo[(a-0x29)/2].count;
      return ((a-0x29)&1)?uint8_t(count):uint8_t(count>>8);
    }
    if (a >= 0x29) return 0;
    if (a >= 0x25) {
      uint16_t sample=uint16_t(adc[(a-0x25)/2]);
      return ((a-0x25)&1)?uint8_t(sample):uint8_t(sample>>8);
    }
    if (a == 0x0c) return (reg[a]&0xf0) | voice[0].status | (voice[1].status<<2);
    return reg[a];
  }
  void write(unsigned a, uint8_t value) {
    if (a >= 0x25) {
      if (!variantA || a > 0x28) return;
      Fifo& f=fifo[(a-0x25)/2];
      if (!((a-0x25)&1)) { f.high=value; return; }
      if (f.count<FIFO_CAPACITY) { f.data[(f.head+f.count)%FIFO_CAPACITY]=int16_t((uint16_t(f.high)<<8)|value); f.count++; }
      return;
    }
    if (a==0x0e || a==0x0f || a==0x1e || a==0x1f) value &= 0x0f;
    if (a == 0x20) { reg[a]=(value&3)|(reg[a]&~(value&0xf4)&0xf4); return; }
    if (a != 0x0c) { reg[a]=value; return; }
    reg[a]=value&0xf0;
    for (unsigned ch=0; ch<2; ch++) {
      unsigned cmd=(value>>(ch*2))&3;
      Voice& v=voice[ch];
      if (!cmd) {
        if (v.status==1 && isFifo(v)) { fifo[ch].count=0; fifo[ch].head=0; }
        v.status=0; v.held=0;
      }
      else if (cmd==1 && (v.status==0 || v.status==3)) {
        v=Voice(); v.status=1;
        v.start=address(ch*6); v.stop=address(ch*6+3); v.pos=v.start;
        v.format=(reg[0x21]>>(ch*4))&7;
      }
    }
  }
  int16_t channelOutput(unsigned ch) const { return channelOut[ch]; }
  std::array<int16_t,2> tick(int16_t adcL=0, int16_t adcR=0) {
    adc={adcL,adcR};
    static const int vol[16]={0,583,777,1036,1382,1843,2457,3277,4370,5827,7771,10362,13818,18427,24573,32768};
    int out[2]={0,0};
    for (unsigned ch=0; ch<2; ch++) {
      Voice& v=voice[ch];
      if (v.status==1) {
        if (v.counter) --v.counter;
        else if (isFifo(v)) {
          // one sample per divider period; when empty, hold and retry next tick
          Fifo& f=fifo[ch];
          if (f.count) { v.held=f.data[f.head]; f.head=(f.head+1)%FIFO_CAPACITY; f.count--; v.counter=reg[0x0e + ch]; }
        }
        else { v.counter=reg[0x0e + ch]; if (!sample(ch)) done(ch,false); }
      }
      // The FIFO DAC has a fixed 0 dB gain (PLAY_VOL ignored); panning still applies.
      int gain=isFifo(v)?32768:vol[(reg[0x0d]>>(ch*4))&15];
      int scaled=(int(v.held)*gain)>>15;
      channelOut[ch]=int16_t(scaled);
      if (!(reg[0x0c]&(0x10<<(ch*2)))) out[0]+=scaled;
      if (!(reg[0x0c]&(0x20<<(ch*2)))) out[1]+=scaled;
    }
    return {int16_t(std::max(-32768,std::min(32767,out[0]))),int16_t(std::max(-32768,std::min(32767,out[1])))};
  }
private:
  struct Voice {
    uint8_t status=0, counter=0, format=0;
    uint32_t start=0, stop=0, pos=0, sampleIndex=0;
    int16_t held=0;
    bool headerRead=false;
    sqv4::Header header;
    sqv4::State state;
  };
  struct Fifo {
    std::array<int16_t,FIFO_CAPACITY> data{};
    unsigned head=0, count=0;
    uint8_t high=0;
  };
  const Memory& memory;
  bool variantA=false;
  std::array<Fifo,2> fifo{};
  bool isFifo(const Voice& v) const { return variantA && v.format==5; }
  std::array<uint8_t,0x25> reg{};
  std::array<Voice,2> voice{};
  std::array<int16_t,2> channelOut{};
  std::array<int16_t,2> adc{};
  uint32_t address(unsigned a) const { return (uint32_t(reg[a])<<16)|(uint32_t(reg[a+1])<<8)|reg[a+2]; }
  void done(unsigned ch, bool error) { voice[ch].status=3; voice[ch].held=0; reg[0x20]|=(0x10<<ch)|(error?4:0); }
  bool sample(unsigned ch) {
    Voice& v=voice[ch];
    if (v.format==1) {
      if (uint64_t(v.pos)+2>v.stop) return false;
      v.held=int16_t((uint16_t(memory.read(v.pos))<<8)|memory.read(v.pos+1)); v.pos+=2; return true;
    }
    if (v.format<2 || v.format>4) {
      if (v.pos>=v.stop) return false;
      v.held=int16_t((int(memory.read(v.pos++))-128)*256); return true;
    }
    if (!v.headerRead) {
      uint8_t h[16];
      if (uint64_t(v.start)+16>v.stop) { done(ch,true); return true; }
      for (unsigned i=0;i<16;i++) h[i]=memory.read(v.start+i);
      if (std::memcmp(h,"SQV4",4) || h[4]!=3 || h[5]>2 || h[6]!=1 || h[7]) { done(ch,true); return true; }
      v.header.version=h[4]; v.header.variant=h[5]; v.header.channels=1;
      v.header.rate=sqv4::read32(h+8); v.header.samples=sqv4::read32(h+12); v.headerRead=true;
    }
    if (v.sampleIndex>=v.header.samples) return false;
    uint64_t frame=uint64_t(v.start)+16+uint64_t(v.sampleIndex/63)*v.header.stride();
    if (frame+v.header.stride()>v.stop) return false;
    unsigned local=v.sampleIndex%63, bits=v.header.bits(), code=0;
    if (local==0) { v.state=sqv4::State(); v.state.index=std::min(88,int(memory.read(uint32_t(frame+1)))); }
    for (unsigned b=0;b<bits;b++) code|=((memory.read(uint32_t(frame+2+(local*bits+b)/8))>>((local*bits+b)%8))&1)<<b;
    v.held=int16_t(sqv4::applyVolume(sqv4::decodeCode(v.header,v.state,code),int(memory.read(uint32_t(frame))),3));
    ++v.sampleIndex; return true;
  }
};
#endif
