// Real-engine integration harness; link in place of src/main.cpp.
#include "../../src/engine/engine.h"
#include "../../src/engine/sqv4.h"
#include "../../src/engine/platform/s3w2.h"
#include <cassert>
#include <cstdio>
#ifdef main
#undef main
#endif

void reportError(String what) { fprintf(stderr,"%s\n",what.c_str()); }

void write(const char* path, const std::vector<uint8_t>& data) {
  FILE* f=fopen(path,"wb"); assert(f);
  assert(fwrite(data.data(),1,data.size(),f)==data.size()); assert(!fclose(f));
}
int main() {
  static DivEngine engine;
  engine.song.systemLen=1; engine.song.system[0]=DIV_SYSTEM_S3W2;
  const char* path="import.sqv4";
  const int16_t pcm[]={12,35,75,147,275};
  for (unsigned variant=0;variant<3;variant++) {
    auto bytes=sqv4::encode(pcm,5,16000,variant);
    bytes[16]=128; // gains and exact H codes must survive without re-encoding
    write(path,bytes);
    auto imported=engine.sampleFromFile(path);
    assert(imported.size()==1);
    DivSample* sample=imported[0];
    assert(sample->samples==5 && sample->centerRate==16000);
    assert(sample->depth==(variant==0?DIV_SAMPLE_DEPTH_SQV4H:
      variant==1?DIV_SAMPLE_DEPTH_SQV4L:DIV_SAMPLE_DEPTH_SQV4L_PLUS));
    sample->render(1u<<DIV_SAMPLE_DEPTH_8BIT);
    assert(sample->getCurBufLen()==bytes.size());
    assert(!memcmp(sample->getCurBuf(),bytes.data(),bytes.size()));
    assert(sample->saveRaw("export.sqv4"));
    auto reimported=engine.sampleFromFile("export.sqv4"); assert(reimported.size()==1);
    assert(reimported[0]->samples==5 && reimported[0]->depth==sample->depth);
    assert(reimported[0]->getCurBufLen()==bytes.size());
    assert(!memcmp(reimported[0]->getCurBuf(),bytes.data(),bytes.size()));
    delete reimported[0];
    auto legacy=bytes; legacy[4]=1; write(path,legacy);
    auto old=engine.sampleFromFile(path); assert(old.size()==1 && old[0]->depth==sample->depth);
    assert(old[0]->getCurBufLen()==legacy.size());
    assert(!memcmp(old[0]->getCurBuf(),legacy.data(),legacy.size())); delete old[0];
    engine.song.sample.push_back(sample); engine.song.sampleLen=engine.song.sample.size();
  }
  // Ordinary IMA remains its own format and uses the S3W2 PCM fallback.
  auto* ima=new DivSample; assert(ima->init(5)); std::copy(pcm,pcm+5,ima->data16);
  ima->convert(DIV_SAMPLE_DEPTH_IMA_ADPCM,1u<<DIV_SAMPLE_DEPTH_8BIT);
  engine.song.sample.push_back(ima); engine.song.sampleLen=engine.song.sample.size();
  // Stereo is rejected, including an empty, otherwise valid stereo file.
  auto stereo=sqv4::encode(nullptr,0,8000,0); stereo[6]=2; write(path,stereo);
  assert(engine.sampleFromFile(path).empty());
  assert(engine.getLastError().find("Stereo")!=String::npos);
  auto bad=sqv4::encode(pcm,5,8000,1); bad.push_back(0); write(path,bad);
  assert(engine.sampleFromFile(path).empty());
  assert(engine.sampleFromFileRaw(path,DIV_SAMPLE_DEPTH_SQV4L,1,false,false,false,8000)==nullptr);
  // The signed RAM image must round-trip compressed bytes through the legacy +128 export path.
  DivPlatformS3W2 dispatch; dispatch.init(&engine,8,48000,DivConfig()); dispatch.renderSamples(0);
  const uint8_t* memory=(const uint8_t*)dispatch.getSampleMem(0);
  size_t offset=0;
  for (int i=0;i<3;i++) {
    auto bytes=engine.song.sample[i]->getSQV4();
    assert(dispatch.isSampleLoaded(0,i));
    for (size_t j=0;j<bytes.size();j++) assert(uint8_t(memory[offset+j]+128)==bytes[j]);
    offset+=bytes.size();
  }
  assert(dispatch.isSampleLoaded(0,3));
  for (unsigned i=0;i<ima->samples;i++) assert(memory[offset+i]==uint8_t(ima->data8[i]));
  assert(dispatch.getSampleMemUsage(0)==offset+ima->samples);
  DivInstrument* ins=new DivInstrument; ins->type=DIV_INS_AMIGA; ins->amiga.initSample=2;
  engine.song.ins.push_back(ins); engine.song.insLen=1;
  dispatch.dispatch(DivCommand(DIV_CMD_INSTRUMENT,0,0));
  dispatch.dispatch(DivCommand(DIV_CMD_NOTE_ON,0,48)); dispatch.tick();
  auto regs=dispatch.getRegisterPool(); assert(regs[0x802]==4 && regs[0x819]==1);
  assert(regs[0x816]==regs[0x810] && regs[0x817]==regs[0x811] && regs[0x818]==regs[0x812]);
  engine.song.sample[2]->loop=true; engine.song.sample[2]->loopStart=0; engine.song.sample[2]->loopEnd=5;
  dispatch.dispatch(DivCommand(DIV_CMD_NOTE_ON,0,48)); dispatch.tick();
  assert(regs[0x802]==4 && regs[0x819]==3);
  dispatch.quit();
  remove(path); remove("export.sqv4");
  puts("SQV4 mono import, stereo/error rejection, save/reimport and S3W2 RAM upload passed");
}
