// SqueezVOX 4 v1/v2 container and ADPCM adapter.
// Contract: CPT32 doc/squeezvox/squeezvox4_spec_stub.md.
#ifndef FURNACE_SQV4_H
#define FURNACE_SQV4_H
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace sqv4 {
struct Header {
  unsigned version=2, variant=0, channels=1;
  uint32_t rate=0, samples=0;
  size_t stride() const { return variant==0?33:25; }
  unsigned bits() const { return variant==0?4:3; }
};
struct State { int predictor=0, index=0; };
inline uint32_t read32(const uint8_t* p) {
  return uint32_t(p[0])|(uint32_t(p[1])<<8)|(uint32_t(p[2])<<16)|(uint32_t(p[3])<<24);
}
inline void write32(uint8_t* p, uint32_t n) {
  for (unsigned i=0;i<4;i++) p[i]=n>>(8*i);
}
inline size_t fileSize(uint32_t samples, unsigned channels, unsigned variant) {
  uint64_t size=16+((uint64_t(samples)+63)/64)*channels*(variant==0?33:25);
  return size>std::numeric_limits<size_t>::max()?0:static_cast<size_t>(size);
}
inline unsigned codeAt(const uint8_t* payload, unsigned sample, unsigned bits) {
  unsigned code=0;
  for (unsigned b=0;b<bits;b++) code|=((payload[(sample*bits+b)/8]>>((sample*bits+b)%8))&1)<<b;
  return code;
}
inline bool inspect(const uint8_t* data, size_t size, Header& h) {
  if (size<16 || !data || std::memcmp(data,"SQV4",4) || (data[4]!=1 && data[4]!=2) ||
      data[5]>2 || (data[6]!=1 && data[6]!=2) || data[7] || !read32(data+8)) return false;
  h.version=data[4]; h.variant=data[5]; h.channels=data[6];
  h.rate=read32(data+8); h.samples=read32(data+12);
  size_t expected=fileSize(h.samples,h.channels,h.variant);
  if (!expected || size!=expected) return false;
  if (h.samples%64) {
    size_t offset=size-h.channels*h.stride();
    for (unsigned ch=0;ch<h.channels;ch++)
      for (unsigned i=h.samples%64;i<64;i++)
        if (codeAt(data+offset+ch*h.stride()+1,i,h.bits())) return false;
  }
  return true;
}
inline int step(int index) {
  static const int table[89]={7,8,9,10,11,12,13,14,16,17,19,21,23,25,28,31,34,37,41,45,50,55,60,66,73,80,88,97,107,118,130,143,157,173,190,209,230,253,279,307,337,371,408,449,494,544,598,658,724,796,876,963,1060,1166,1282,1411,1552,1707,1878,2066,2272,2499,2749,3024,3327,3660,4026,4428,4871,5358,5894,6484,7132,7845,8630,9493,10442,11487,12635,13899,15289,16818,18500,20350,22385,24623,27086,29794,32767};
  return table[index];
}
inline int delta(const Header& h, const State& s, unsigned code) {
  int v=step(s.index), diff;
  if (!h.variant) {
    diff=(v>>3)+((code&1)?v>>2:0)+((code&2)?v>>1:0)+((code&4)?v:0);
    return (code&8)?-diff:diff;
  }
  static const int levels[3][4]={{1,3,5,7},{10,11,12,14},{1,2,3,14}};
  diff=(levels[h.variant==1?0:(h.version==1?1:2)][code&3]*v+4)/8;
  return (code&4)?-diff:diff;
}
inline int decodeCode(const Header& h, State& s, unsigned code) {
  static const int hi[8]={-1,-1,-1,-1,2,4,6,8}, lo[4]={-1,2,4,6};
  s.predictor=std::max(-32768,std::min(32767,s.predictor+delta(h,s,code)));
  s.index=std::max(0,std::min(88,s.index+(!h.variant?hi[code&7]:lo[code&3])));
  return s.predictor;
}
inline unsigned encodeCode(const Header& h, State& s, int input) {
  unsigned code=0;
  if (!h.variant) { // IMA ADPCM threshold quantizer; no WAV block header in SQV4.
    int diff=input-s.predictor, v=step(s.index);
    if (diff<0) { code=8; diff=-diff; }
    for (int bit=2;bit>=0;bit--,v>>=1) if (diff>=v) { code|=1<<bit; diff-=v; }
  } else {
    int best=std::numeric_limits<int>::max();
    for (unsigned c=0;c<8;c++) {
      int reconstructed=std::max(-32768,std::min(32767,s.predictor+delta(h,s,c)));
      int error=input-reconstructed; if (error<0) error=-error;
      if (error<best) { best=error; code=c; }
    }
  }
  decodeCode(h,s,code);
  return code;
}
inline std::vector<uint8_t> encode(const int16_t* pcm, uint32_t samples, uint32_t rate, unsigned variant) {
  if (!rate || variant>2 || (samples && !pcm)) return {};
  Header h; h.variant=variant; h.rate=rate; h.samples=samples;
  size_t size=fileSize(samples,1,variant); if (!size) return {};
  std::vector<uint8_t> out(size,0);
  std::memcpy(out.data(),"SQV4",4); out[4]=2; out[5]=variant; out[6]=1;
  write32(out.data()+8,rate); write32(out.data()+12,samples);
  State state;
  for (uint64_t first=0;first<samples;first+=64) {
    uint8_t* fragment=out.data()+16+(first/64)*h.stride(); fragment[0]=255;
    for (unsigned i=0;i<64 && first+i<samples;i++) {
      unsigned code=encodeCode(h,state,pcm[first+i]);
      for (unsigned b=0;b<h.bits();b++) fragment[1+(i*h.bits()+b)/8]|=((code>>b)&1)<<((i*h.bits()+b)%8);
    }
  }
  return out;
}
inline bool decode(const uint8_t* data, size_t size, std::vector<int16_t>& out, Header& h) {
  if (!inspect(data,size,h) || uint64_t(h.samples)*h.channels>16777215) return false;
  out.resize(size_t(h.samples)*h.channels);
  State state[2];
  for (uint32_t i=0;i<h.samples;i++) for (unsigned ch=0;ch<h.channels;ch++) {
    const uint8_t* f=data+16+(size_t(i/64)*h.channels+ch)*h.stride();
    out[size_t(i)*h.channels+ch]=decodeCode(h,state[ch],codeAt(f+1,i%64,h.bits()))*f[0]/255;
  }
  return true;
}
}
#endif
