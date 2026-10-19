#include "../../src/engine/export/sgubWriter.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
static int checkIndex=0;
static void check(bool b) { ++checkIndex; if (!b) { std::fprintf(stderr,"check %d failed\n",checkIndex); std::abort(); } }
static unsigned be16(const std::vector<uint8_t>& v, int n) { return (unsigned(v[n])<<8)|v[n+1]; }
static unsigned be32(const std::vector<uint8_t>& v, int n) {
  return (unsigned(v[n])<<24)|(unsigned(v[n+1])<<16)|(unsigned(v[n+2])<<8)|v[n+3];
}
template<class F> static void rejects(F f) {
  bool caught=false; try { f(); } catch (const std::runtime_error&) { caught=true; } check(caught);
}
int main() {
  check(divSGUBModulation((3<<3)|7)==((5<<3)|7));
  check(divSGUBModulation((4<<3)|2)==((6<<3)|2));
  check(divSGUBModulation((31<<3)|3)==3);
  check(divSGUBModulation((2<<3)|1)==((2<<3)|1));
  DivSGUBWriter w(true);
  std::vector<uint8_t> wave(256,128), pcm={128,129,127};
  auto waveId=w.addPool(1,wave.data(),wave.size());
  auto pcmId=w.addPool(2,pcm.data(),pcm.size());
  w.initialFill(2,0,0x100000,128);
  w.initialCopy(2,0,pcmId,0,pcm.size());
  w.copy(0,0,waveId,0,256); w.write(0,0x80a,0); w.wait(400);
  w.write(0,0x80a,0); w.wait(400); // repeated strobes retain their order
  w.wait(800); // silent frame
  w.copy(1,0x700,waveId,0,256); w.wait(1); // final partial frame
  const auto& b=w.finish();
  check(b.size()==404 && b[4]==2 && b[5]==16 && b[6]==1 && b[7]==0);
  check(be32(b,8)==3 && be32(b,12)==2 && be32(b,16)==2 && be32(b,32)==b.size());
  check(be32(b,20)==40 && be32(b,24)==331 && be32(b,28)==356 && be16(b,36)==40);
  check(be16(b,40)==0 && b[42]==1 && be32(b,44)==256 && be32(b,52)==72);
  check(be16(b,56)==1 && b[58]==2 && be32(b,60)==3 && be32(b,64)==3 && be32(b,68)==328);
  check(b[331]==1 && be16(b,356)==3); // initial FILL followed by three ordered frame operations
  check(be16(b,356)==3); check(b[358]==2); check(b[373]==0); check(b[379]==0);
  check(be16(b,385)==0 && be16(b,387)==1 && b[389]==2);
  rejects([&]() { w.initial(0,0,0); });
  DivSGUBWriter invalid;
  rejects([&]() { invalid.write(0,0x900,0); });
  rejects([&]() { invalid.write(2,0x100000,0); });
  rejects([&]() { invalid.write(3,0,0); });
  rejects([&]() { invalid.initialFill(0,0x7ff,2,0); });
  rejects([&]() { invalid.addPool(1,pcm.data(),pcm.size()); });
  rejects([&]() { invalid.initialCopy(2,0,0,0,1); });
  DivSGUBWriter overflow;
  for (int i=0;i<65535;i++) overflow.write(0,0,0);
  rejects([&]() { overflow.write(0,0,0); });
  DivSGUBWriter size;
  rejects([&]() { for (;;) size.wait(800); });
  DivSGUBWriter drift;
  for (int i=0;i<125;i++) drift.wait(384);
  check(be32(drift.finish(),8)==60);
  DivSGUBWriter empty;
  check(be32(empty.finish(),8)==1);
  std::cout << "SGUB v2 writer checks passed\n";
}
