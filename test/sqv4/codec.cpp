#include "../../src/engine/sqv4.h"
#include "../../src/engine/sample.h"
#include "../../src/engine/platform/sound/3ws8pn/sound.hpp"
#include "../../src/ta-log.h"
#include <cassert>
#include <cstdio>
#include "v3_golden.h"

// Standalone sample tests use the real codecs/serialization, with quiet logging.
int writeLog(int, const char*, fmt::printf_args) { return 0; }
FILE* ps_fopen(const char* name, const char* mode) { return fopen(name,mode); }

void address(S3W2_Sound& sound, unsigned reg, uint32_t value) {
  sound.writeRegister(reg,value>>16); sound.writeRegister(reg+1,value>>8); sound.writeRegister(reg+2,value);
}
void configure(S3W2_Sound& sound, const std::vector<uint8_t>& bytes, unsigned loop=128, unsigned control=1) {
  for (size_t i=0;i<bytes.size();i++) sound.writePCMRAM(128+i,bytes[i]);
  sound.writeRegister(0x800,0x17); sound.writeRegister(0x801,0x70); // one sample per 192 kHz call
  sound.writeRegister(0x802,4); sound.writeRegister(0x803,255);
  address(sound,0x810,128); address(sound,0x813,128+bytes.size()); address(sound,0x816,loop);
  sound.writeRegister(0x80a,0); sound.writeRegister(0x819,control);
}

int main() {
  // Compare encoder bytes against the upstream Rust v3 encoder at three fragment boundaries.
  std::vector<int16_t> fixture(130);
  for (int i=0;i<130;i++) fixture[i]=int(i*1777%65536)-32768;
  for (unsigned variant=0;variant<3;variant++) {
    auto bytes=sqv4::encode(fixture.data(),fixture.size(),8000,variant);
    assert(bytes==upstreamGolden[variant]);
    sqv4::Header h; std::vector<int16_t> decoded;
    assert(sqv4::decode(bytes.data(),bytes.size(),decoded,h));
    assert(decoded==upstreamDecoded[variant]);
  }
  std::vector<int16_t> ramp(65);
  for (unsigned i=0;i<65;i++) ramp[i]=i+1;
  auto h=sqv4::encode(ramp.data(),65,8000,0);
  const int16_t lplus[]={12,35,75,147,275};
  auto plus=sqv4::encode(lplus,5,8000,2);
  sqv4::Header header;
  std::vector<int16_t> decoded;
  assert(h.size()==16+2*34 && h[4]==3);
  assert(sqv4::decode(h.data(),h.size(),decoded,header));
  std::vector<int16_t> hDecoded=decoded;
  assert(sqv4::decode(plus.data(),plus.size(),decoded,header));
  std::vector<int16_t> plusDecoded=decoded;
  auto legacy=h;
  for (unsigned version: {1u,2u}) { legacy[4]=version; assert(!sqv4::inspect(legacy.data(),legacy.size(),header)); }
  auto badIndex=h; badIndex[17]=89; assert(!sqv4::inspect(badIndex.data(),badIndex.size(),header));
  auto badPadding=h; badPadding[16+33]|=0x80; assert(!sqv4::inspect(badPadding.data(),badPadding.size(),header));
  // Empty, truncated, trailing bytes, header fields and final padding.
  auto empty=sqv4::encode(nullptr,0,8000,1);
  assert(empty.size()==16 && sqv4::inspect(empty.data(),empty.size(),header));
  for (unsigned field: {0u,4u,5u,6u,7u,8u}) {
    auto bad=h;
    if (field==8) std::fill(bad.begin()+8,bad.begin()+12,0);
    else bad[field]=255;
    assert(!sqv4::inspect(bad.data(),bad.size(),header));
  }
  assert(!sqv4::inspect(h.data(),15,header));
  assert(!sqv4::inspect(h.data(),h.size()-1,header));
  auto bad=h; bad.push_back(0); assert(!sqv4::inspect(bad.data(),bad.size(),header));
  bad=h; bad.back()=1; assert(!sqv4::inspect(bad.data(),bad.size(),header));

  // Exercise all 3-bit packing positions, state continuity, partial fragments,
  // editor conversion, PCM preview, FUR persistence and undo/redo.
  for (unsigned count: {0u,1u,5u,63u,64u,65u,129u,1024u}) {
    for (auto depth: {DIV_SAMPLE_DEPTH_SQV4H,DIV_SAMPLE_DEPTH_SQV4L,DIV_SAMPLE_DEPTH_SQV4L_PLUS}) {
      DivSample sample;
      sample.centerRate=8000; assert(sample.init(count));
      for (unsigned i=0;i<count;i++) sample.data16[i]=int(i*1777%65536)-32768;
      sample.prepareUndo(true);
      sample.convert(depth,(1u<<DIV_SAMPLE_DEPTH_8BIT));
      assert(sample.samples==count && sample.depth==depth);
      auto bytes=sample.getSQV4();
      assert(sqv4::decode(bytes.data(),bytes.size(),decoded,header));
      assert(header.variant==divSampleSQV4Variant(depth));
      assert(std::equal(decoded.begin(),decoded.end(),sample.data16));
      SafeWriter writer; writer.init(); sample.putSampleData(&writer);
      SafeReader reader(writer.getFinalBuf(),writer.size());
      DivSample loaded; assert(loaded.readSampleData(reader,248)==DIV_DATA_SUCCESS);
      loaded.render(1u<<DIV_SAMPLE_DEPTH_8BIT);
      assert(loaded.depth==depth && loaded.getCurBufLen()==sample.getCurBufLen());
      assert(!memcmp(loaded.getCurBuf(),sample.getCurBuf(),sample.getCurBufLen()));
      assert(std::equal(decoded.begin(),decoded.end(),loaded.data16));
      assert(sample.undo()==2 && sample.depth==DIV_SAMPLE_DEPTH_16BIT);
      assert(sample.redo()==2 && sample.depth==depth);
      if (count) {
        auto trimmed=sample.getSQV4(count-1);
        assert(sqv4::inspect(trimmed.data(),trimmed.size(),header) && header.samples==count-1);
      }
      writer.finish();
    }
  }
  // IMA remains a distinct WAV-style block; explicit .sqv export adapts it to H.
  DivSample ima; assert(ima.init(65)); std::copy(ramp.begin(),ramp.end(),ima.data16);
  ima.convert(DIV_SAMPLE_DEPTH_IMA_ADPCM,1u<<DIV_SAMPLE_DEPTH_8BIT);
  assert(memcmp(ima.getCurBuf(),"SQV4",4)!=0);
  auto imaImage=ima.getSQV4(); assert(sqv4::inspect(imaImage.data(),imaImage.size(),header) && header.variant==0);

  for (int revision=0;revision<2;revision++) {
  S3W2_Sound sound;
  sound.setRevision(revision);
  configure(sound,h,128+16+34,3);
  for (int i=0;i<65;i++) assert(sound.generateSample(0)==hDecoded[i]/4);
  assert(sound.generateSample(0)==hDecoded[63]/4);
  sound.writeRegister(0x80a,0); assert(sound.generateSample(0)==0);
  configure(sound,plus);
  for (int expected: plusDecoded) assert(sound.generateSample(0)==expected/4);
  assert(sound.generateSample(0)==0 && sound.generateSample(0)==0);
  configure(sound,h,128+51); assert(sound.generateSample(0)==0); // invalid loop
  bad=h; bad[0]=0; configure(sound,bad); assert(sound.generateSample(0)==0);
  configure(sound,h); address(sound,0x813,PCM_RAM_SIZE+8); assert(sound.generateSample(0)==0);
  // PCM RAM mutation invalidates decoder state.
  configure(sound,h); sound.writePCMRAM(128,0); assert(sound.generateSample(0)==0);
  }
  std::puts("SQV4 codec, sample conversion/persistence/undo and S3W2 playback passed");
}
