#include "../../src/engine/platform/sound/34si22/playback.h"
#include <cassert>
#include <vector>
#include "../sqv4/v3_golden.h"
struct Ram: Si22Playback::Memory { std::vector<uint8_t> b=std::vector<uint8_t>(1024); uint8_t read(uint32_t a) const override { return b.at(a); } };
static void address(Si22Playback& p, unsigned a, unsigned n) { p.write(a,n>>16); p.write(a+1,n>>8); p.write(a+2,n); }
static void push(Si22Playback& p, unsigned ch, int16_t v) { p.write(0x25+ch*2,uint16_t(v)>>8); p.write(0x26+ch*2,uint16_t(v)&0xff); }
static unsigned count(const Si22Playback& p, unsigned ch) { return (p.read(0x29+ch*2)<<8)|p.read(0x2a+ch*2); }
// 34SI22A FIFO DAC (CPT32 doc/spec_34SI22.md section 11).
static void variantA(Ram& r) {
 Si22Playback p(r); p.setVariantA(true);
 push(p,0,0x1234); push(p,0,-2); push(p,0,0x0100); push(p,1,0x7fff);
 assert(count(p,0)==3 && count(p,1)==1 && p.read(0x2d)==0);
 assert(p.read(0x25)==0); // reads still show the ADC
 p.write(0x21,0x55); p.write(0x0d,0xff); p.write(0x0e,1); p.write(0x0c,0x05|0x40|0x20); // A left, B right
 const int16_t expect[8][2]={{0x1234,0x7fff},{0x1234,0x7fff},{-2,0x7fff},{-2,0x7fff},{0x100,0x7fff},{0x100,0x7fff},{0x100,0x7fff},{0x100,0x7fff}};
 for (auto& e: expect) { auto o=p.tick(); assert(o[0]==e[0] && o[1]==e[1]); }
 push(p,0,0x4000); assert(p.tick()[0]==0x4000); // underrun retries every tick
 p.write(0x0d,0x00); { auto o=p.tick(); assert(o[0]==0x4000 && o[1]==0x7fff); } // PLAY_VOL ignored
 p.write(0x0c,0x0f|0x10|0x20|0x40|0x80); { auto o=p.tick(); assert(o[0]==0 && o[1]==0); } // panning applies
 p.write(0x0c,0x0f|0x40|0x20); p.write(0x0d,0xff);
 assert((p.read(0x0c)&15)==5 && p.read(0x20)==0); // never finishes
 for (unsigned i=0;i<Si22Playback::FIFO_CAPACITY+5;i++) push(p,1,1);
 assert(count(p,1)==2048);
 p.write(0x0c,0x0f); assert(count(p,1)==2048 && (p.read(0x0c)&15)==5); // 3 = keep
 p.write(0x0c,0x03); assert(count(p,1)==0 && (p.read(0x0c)&15)==1); // stopping B empties it
 push(p,1,9); p.write(0x0c,0x00); assert(count(p,1)==1); // stopping a stopped channel keeps it
 // The 34SI22 ignores the pushes and plays format 5 as PCM8.
 Si22Playback b(r); push(b,0,0x1234); assert(count(b,0)==0);
 r.b[0]=0xc0; address(b,3,1); b.write(0x21,5); b.write(0x0d,15); b.write(0x0c,0x0d); assert(b.tick()[0]==0x4000);
}
int main() {
 Ram r; Si22Playback p(r);
 variantA(r);
 // Raw ADC readings are read-only, update on every tick, and reset to zero.
 for (unsigned a=0x25;a<=0x28;a++) assert(p.read(a)==0);
 p.tick(0x1234,-2); assert(p.read(0x25)==0x12 && p.read(0x26)==0x34 && p.read(0x27)==0xff && p.read(0x28)==0xfe);
 for (unsigned a=0x25;a<=0x28;a++) p.write(a,0xaa);
 assert(p.read(0x25)==0x12 && p.read(0x26)==0x34 && p.read(0x27)==0xff && p.read(0x28)==0xfe);
 p.tick(-32768,32767); assert(p.read(0x25)==0x80 && p.read(0x26)==0 && p.read(0x27)==0x7f && p.read(0x28)==0xff);
 p.write(0x29,0xff); assert(p.read(0x29)==0 && p.read(0x3f)==0);
 p.reset(); for (unsigned a=0x25;a<=0x28;a++) assert(p.read(a)==0);
 // Divider registers expose only their low four bits, including record registers.
 for (unsigned a: {0x0e,0x0f,0x1e,0x1f}) { p.write(a,0xf0); assert(p.read(a)==0); p.write(a,0xff); assert(p.read(a)==15); }
 p.reset(); r.b[0]=0xff; r.b[1]=0; address(p,0,0); address(p,3,2); p.write(0x0d,15); p.write(0x0e,1); p.write(0x0c,0x0d);
 assert(p.tick()[0]==32512); assert(p.tick()[0]==32512); assert(p.tick()[0]==-32768); assert(p.tick()[0]==-32768); assert(p.tick()[0]==0); assert((p.read(0x0c)&3)==3); assert(p.read(0x20)&0x10); p.write(0x20,0x10); assert(!(p.read(0x20)&0x10));
 p.reset(); r.b[1]=0x12; r.b[2]=0x34; address(p,0,1); address(p,3,3); p.write(0x21,1); p.write(0x0d,15); p.write(0x0c,0x0d); assert(p.tick()[0]==0x1234); assert(p.tick()[0]==0);
 p.reset(); address(p,0,0); address(p,3,15); p.write(0x21,2); p.write(0x0c,0x0d); p.tick(); assert((p.read(0x20)&0x14)==0x14);
 // DIV 15 holds each sample for 16 ticks, including a write with reserved bits set.
 p.reset(); r.b[0]=0xff; r.b[1]=0; address(p,0,0); address(p,3,2);
 p.write(0x0d,15); p.write(0x0e,0xff); p.write(0x0c,0x0d);
 for (int i=0;i<16;i++) assert(p.tick()[0]==32512);
 assert(p.tick()[0]==-32768);
 // Every specified volume coefficient, including signed rounding and mute.
 const int gain[16]={0,583,777,1036,1382,1843,2457,3277,4370,5827,7771,10362,13818,18427,24573,32768};
 for (int volume=0;volume<16;volume++) {
   p.reset(); r.b[0]=0x7f; r.b[1]=0xff; r.b[2]=0x80; r.b[3]=0;
   address(p,0,0); address(p,3,4); p.write(0x21,1); p.write(0x0d,volume); p.write(0x0c,0x0d);
   auto positive=p.tick(); auto negative=p.tick();
   assert(positive[0]==(gain[volume]?gain[volume]-1:0)); assert(positive[1]==positive[0]);
   assert(negative[0]==-gain[volume]); assert(negative[1]==negative[0]);
 }
 // Two simultaneous channels, latched start/stop, output gating and saturation.
 p.reset(); r.b[0]=0xff; r.b[1]=0xff; address(p,0,0); address(p,3,2); address(p,6,0); address(p,9,2);
 p.write(0x0d,0xff); p.write(0x0c,0x05); address(p,3,0); assert(p.tick()[0]==32767);
 p.write(0x0c,0x2f); auto stereo=p.tick(); assert(stereo[0]==32767 && stereo[1]==32512);
 p.tick(); assert((p.read(0x0c)&15)==15); assert((p.read(0x20)&0x30)==0x30);
 // Upstream-generated v3 H/L/L+ containers, across two fragment boundaries.
 for (size_t variant=0;variant<upstreamGolden.size();variant++) {
   const auto& fixture=upstreamGolden[variant];
   p.reset(); for (size_t i=0;i<fixture.size();i++) r.b[i]=fixture[i];
   address(p,0,0); address(p,3,unsigned(fixture.size())); p.write(0x21,2); p.write(0x0d,15); p.write(0x0c,0x0d);
   std::vector<int16_t> expected; sqv4::Header header; assert(sqv4::decode(fixture.data(),fixture.size(),expected,header));
   assert(expected==upstreamDecoded[variant]);
   for (auto sample:upstreamDecoded[variant]) assert(p.tick()[0]==sample);
   assert(p.tick()[0]==0 && (p.read(0x20)&0x10));
 }
 // Hardware clamps an out-of-range step index; the strict file parser rejects it.
 auto invalid=upstreamGolden[0]; invalid[17]=255; sqv4::Header invalidHeader;
 assert(!sqv4::inspect(invalid.data(),invalid.size(),invalidHeader));
 for (size_t i=0;i<invalid.size();i++) r.b[i]=invalid[i];
 p.reset(); address(p,0,0); address(p,3,unsigned(invalid.size())); p.write(0x21,2); p.write(0x0d,15); p.write(0x0c,0x0d);
 auto clamped=p.tick(); invalid[17]=88; for (size_t i=0;i<invalid.size();i++) r.b[i]=invalid[i];
 p.reset(); address(p,0,0); address(p,3,unsigned(invalid.size())); p.write(0x21,2); p.write(0x0d,15); p.write(0x0c,0x0d); assert(p.tick()==clamped);
 std::vector<int16_t> in={1000,-2000,3000}; auto data=sqv4::encode(in.data(),3,48000,1); for (size_t i=0;i<data.size();i++) r.b[i]=data[i]; p.reset(); address(p,0,0); address(p,3,unsigned(data.size())); p.write(0x21,2); p.write(0x0d,15); p.write(0x0c,0x0d); std::vector<int16_t> expected; sqv4::Header h; assert(sqv4::decode(data.data(),data.size(),expected,h)); for (auto v:expected) assert(p.tick()[0]==v); assert(p.tick()[0]==0); return 0;
}
