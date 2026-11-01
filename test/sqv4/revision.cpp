#include "../../src/engine/platform/sound/3ws8pn/sound.hpp"
#include <cassert>
#include <cstdio>

static void push(S3W2_Sound& s, int16_t sample, int ch=0) {
  s.writeRegister(0x80C+ch*32,static_cast<uint16_t>(sample)>>8);
  s.writeRegister(0x80D+ch*32,sample&255);
}

static void address(S3W2_Sound& s, uint32_t reg, uint32_t value) {
  s.writeRegister(reg,value>>16);
  s.writeRegister(reg+1,value>>8);
  s.writeRegister(reg+2,value);
}

static void sqv4Interpolation() {
  std::vector<int16_t> pcm(130);
  for (int i=0;i<130;i++) pcm[i]=i*40-3000;
  for (unsigned variant=0;variant<3;variant++) {
    auto bytes=sqv4::encode(pcm.data(),pcm.size(),8000,variant);
    sqv4::Header header;
    assert(sqv4::inspect(bytes.data(),bytes.size(),header));
    // Distinct fragment gains verify that interpolation follows fragment gain.
    bytes[16+header.stride()]=128;
    bytes[16+2*header.stride()]=64;
    for (unsigned version=3;version<=3;version++) {
      bytes[4]=version;
      std::vector<int16_t> decoded;
      assert(sqv4::decode(bytes.data(),bytes.size(),decoded,header));
      for (int revision=0;revision<2;revision++) {
        for (int loop=0;loop<2;loop++) {
          S3W2_Sound s;
          s.setRevision(revision);
          for (size_t i=0;i<bytes.size();i++) s.writePCMRAM(256+i,bytes[i]);
          address(s,0x810,256);
          address(s,0x813,256+bytes.size());
          address(s,0x816,256+16+header.stride()); // loop to decoded sample 63
          s.writeRegister(0x800,1500>>8);
          s.writeRegister(0x801,1500&255); // exactly 0.25 sample per call
          s.writeRegister(0x802,4);
          s.writeRegister(0x803,173); // retain Q8 fraction through channel volume
          s.writeRegister(0x819,loop?3:1);
          auto expected=[&](unsigned index,unsigned frac) {
            int32_t q8=int32_t(decoded[index])*256;
            if (revision && index+1<decoded.size())
              q8=int32_t(decoded[index])*(256-frac)+int32_t(decoded[index+1])*frac;
            return q8*173/(255*256*4);
          };
          // Includes both 62->63 and 125->126 frame crossings, and a partial
          // final fragment. Repeated fractional reads must keep both frames.
          for (unsigned i=0;i<130;i++) for (unsigned frac=0;frac<256;frac+=64)
            assert(s.generateSample(0)==expected(i,frac));
          if (loop) {
            for (unsigned i=63;i<130;i++) for (unsigned frac=0;frac<256;frac+=64)
              assert(s.generateSample(0)==expected(i,frac));
          } else {
            assert(s.generateSample(0)==0);
          }
          s.writeRegister(0x80A,0);
          s.writeRegister(0x819,1);
          assert(s.generateSample(0)==expected(0,0));
          // Changing revision after lookahead invalidates the decoder cache.
          s.setRevision(1-revision);
          int32_t q8=int32_t(decoded[0])*256;
          if (!revision) q8=int32_t(decoded[0])*192+int32_t(decoded[1])*64;
          assert(s.generateSample(0)==q8*173/(255*256*4));
          // RAM writes must invalidate cached lookahead as well.
          s.writePCMRAM(256,0);
          assert(s.generateSample(0)==0);
        }
      }
    }
  }
}

int main() {
  sqv4Interpolation();
  for (int revision=0;revision<2;revision++) {
    S3W2_Sound s;
    s.setRevision(revision);
    s.writeRegister(0,0);
    s.writeRegister(1,128);
    s.writeRegister(0x801,24);
    s.writeRegister(0x803,4);
    int16_t sample=0;
    for (int i=0;i<16;i++) sample=s.generateSample(0);
    assert(revision?(sample>=-64 && sample<=-60):sample==-128);
    // Modulation parameters must round-trip independently.
    for (int i=0;i<4;i++) s.writeRegister(0x806+i,0x12+i*0x11);
    for (int i=0;i<4;i++) assert(s.readRegister(0x806+i)==0x12+i*0x11);
    s.writeRegister(0x805,0);

    s.reset(); // revision survives reset
    s.writeRegister(0x802,1);
    s.writeRegister(0x803,4);
    s.writeRegister(0x800,1500>>8);
    s.writeRegister(0x801,1500&255); // 0.25 PCM byte/internal tick
    s.writeRegister(0x815,2);
    s.writeRegister(0x819,1);
    s.writePCMRAM(0,0);
    s.writePCMRAM(1,100);
    s.writePCMRAM(2,200);
    s.generateSample(0); s.generateSample(0);
    assert(s.generateSample(0)==(revision?-78:-128)); // byte 0.5
    for (int i=0;i<3;i++) s.generateSample(0);
    assert(s.generateSample(0)==-28); // byte 1.5: one-shot holds last tap
    s.writeRegister(0x819,3);
    s.writeRegister(0x80A,0);
    for (int i=0;i<6;i++) s.generateSample(0);
    assert(s.generateSample(0)==(revision?-78:-28)); // loop tap at byte 1.5

    s.reset();
    s.writeRegister(0x802,5);
    push(s,1000);
    assert(s.readRegister(0x802)==(revision?5:0));
    assert(s.readRegister(0x80F)==revision);
    if (!revision) continue;
    push(s,-2000); push(s,3000);
    s.writeRegister(0x801,1); // decimation: two output steps/sample
    const int expected[]={1000,1000,-2000,-2000,3000,3000,3000};
    for (int value: expected) {
      for (int j=0;j<4;j++) assert(s.generateSample(0)==value);
    }
    assert(s.readRegister(0x80F)==0);
    push(s,32767);
    s.writeRegister(0x80A,0);
    assert(s.readRegister(0x80F)==0 && s.generateSample(0)==0);
    for (int i=0;i<2060;i++) push(s,32767);
    assert(s.readRegister(0x80E)==8 && s.readRegister(0x80F)==0);
    s.writeRegister(0x80A,0);
    push(s,32767);
    s.writeRegister(0x801,0);
    s.writeRegister(0x804,0xF0);
    auto out=s.clock(2);
    assert(out[0][0][0]==32767 && out[0][1][0]==0);
    s.setRevision(0);
    assert(s.readRegister(0x802)==0 && s.readRegister(0x80F)==0);
    s.setRevision(1);
    assert(s.readRegister(0x80F)==0);
  }
  // Wavetable interpolation must wrap from entry 255 to entry 0.
  S3W2_Sound s;
  s.setRevision(1);
  s.writeRegister(255,0); s.writeRegister(0,200);
  s.writeRegister(0x803,4);
  s.writeRegister(0x801,125); // exactly 1/6 entry per tick
  for (int i=0;i<1532;i++) s.generateSample(0);
  int sample=s.generateSample(0); // approximately entry 255.5
  assert(sample>=-30 && sample<=-28);
  puts("S3W2 wave/PCM/SQV4 revision interpolation and FIFO checks passed");
}
