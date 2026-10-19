/**
 * Furnace Tracker - multi-system chiptune tracker
 * Copyright (C) 2021-2026 tildearrow and contributors
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 */
#ifndef _DIV_SGUB_WRITER_H
#define _DIV_SGUB_WRITER_H
#include <cstdint>
#include <stdexcept>
#include <vector>

inline uint8_t divSGUBModulation(uint8_t value) {
  unsigned mode=value>>3;
  if (mode==3 || mode==4) mode+=2;
  else if (mode>4) mode=0;
  return (mode<<3)|(value&7);
}

class DivSGUBWriter {
  static constexpr size_t limit=0x00e00000;
  struct Pool { uint8_t type; std::vector<uint8_t> data; };
  std::vector<uint8_t> bytes, payloads, initialOps, frame, frameData;
  std::vector<Pool> pools;
  uint32_t frames=0, initialCount=0;
  unsigned samples=0, frameOps=0;
  bool started=false, repeat;
  static void be16(std::vector<uint8_t>& v, uint16_t n) { v.push_back(n>>8); v.push_back(n); }
  static void be24(std::vector<uint8_t>& v, unsigned n) { v.push_back(n>>16); v.push_back(n>>8); v.push_back(n); }
  static void be32(std::vector<uint8_t>& v, uint32_t n) { v.push_back(n>>24); v.push_back(n>>16); v.push_back(n>>8); v.push_back(n); }
  void room(size_t extra) const {
    if (40+pools.size()*16+payloads.size()+initialOps.size()+frameData.size()+frame.size()+extra+2>limit)
      throw std::runtime_error("SGUB exceeds the 14 MiB limit");
  }
  static void targetAddress(unsigned target, unsigned addr, size_t length=1) {
    if (target>2 || !length || addr>(target==2?0xfffff:0x8ff) || length-1>(target==2?0xfffff:0x8ff)-addr)
      throw std::runtime_error("Invalid SGUB register address");
  }
  static void appendWrite(std::vector<uint8_t>& out, unsigned target, unsigned addr, unsigned value) {
    targetAddress(target,addr); out.push_back(0); out.push_back(target); be24(out,addr); out.push_back(value);
  }
  void appendInitial(const std::vector<uint8_t>& op) {
    if (started) throw std::runtime_error("SGUB initialization after playback began");
    room(op.size()); initialOps.insert(initialOps.end(),op.begin(),op.end()); initialCount++;
  }
  void appendFrame(const std::vector<uint8_t>& op) {
    started=true; if (frameOps>=65535) throw std::runtime_error("Too many SGUB operations in one 60 Hz frame");
    room(op.size()); frame.insert(frame.end(),op.begin(),op.end()); frameOps++;
  }
  void flush() {
    room(2); be16(frameData,(uint16_t)frameOps); frameData.insert(frameData.end(),frame.begin(),frame.end());
    frame.clear(); frameOps=0; frames++;
  }
public:
  explicit DivSGUBWriter(bool doRepeat=false): repeat(doRepeat) {}
  void reserve(size_t n) { if (n>limit) throw std::runtime_error("SGUB exceeds the 14 MiB limit"); bytes.reserve(n); }
  uint16_t addPool(uint8_t type, const uint8_t* data, size_t size) {
    if ((type==1 && size!=256) || (type==2 && (!size || size>0x100000)) || type<1 || type>2)
      throw std::runtime_error("Invalid SGUB pool payload");
    if (pools.size()>=65536) throw std::runtime_error("Too many SGUB pool entries");
    room(size); uint16_t id=(uint16_t)pools.size();
    pools.push_back({type,std::vector<uint8_t>(data,data+size)});
    payloads.insert(payloads.end(),data,data+size); return id;
  }
  void initial(unsigned target, unsigned addr, unsigned value) {
    std::vector<uint8_t> op; appendWrite(op,target,addr,value); appendInitial(op);
  }
  void initialFill(unsigned target, unsigned addr, uint32_t length, unsigned value) {
    targetAddress(target,addr,length);
    if (target<2 && (addr>=0x800 || length>0x800-addr)) throw std::runtime_error("Invalid SGUB wave fill");
    std::vector<uint8_t> op{1,(uint8_t)target}; be24(op,addr); be32(op,length); op.push_back(value); appendInitial(op);
  }
  void initialCopy(unsigned target, unsigned addr, uint16_t pool, uint32_t source, uint32_t length) {
    if (pool>=pools.size() || !length || source>pools[pool].data.size() || length>pools[pool].data.size()-source)
      throw std::runtime_error("Invalid SGUB pool copy");
    targetAddress(target,addr,length);
    if (pools[pool].type==1 && (target==2 || addr>=0x800 || length>0x800-addr || (addr&255) || (source&255) || (length&255))) throw std::runtime_error("Invalid SGUB wave copy");
    if (pools[pool].type==2 && target!=2) throw std::runtime_error("Invalid SGUB PCM copy");
    std::vector<uint8_t> op{2,(uint8_t)target}; be24(op,addr); be16(op,pool); be32(op,source); be32(op,length); appendInitial(op);
  }
  void write(unsigned target, unsigned addr, unsigned value) {
    std::vector<uint8_t> op; appendWrite(op,target,addr,value); appendFrame(op);
  }
  void copy(uint8_t target, unsigned addr, uint16_t pool, uint32_t source, uint32_t length) {
    if (pool>=pools.size() || !length || source>pools[pool].data.size() || length>pools[pool].data.size()-source)
      throw std::runtime_error("Invalid SGUB pool copy");
    targetAddress(target,addr,length);
    if (pools[pool].type==1 && (target==2 || addr>=0x800 || length>0x800-addr || (addr&255) || (source&255) || (length&255))) throw std::runtime_error("Invalid SGUB wave copy");
    if (pools[pool].type==2 && target!=2) throw std::runtime_error("Invalid SGUB PCM copy");
    std::vector<uint8_t> op{2,target}; be24(op,addr); be16(op,pool); be32(op,source); be32(op,length); appendFrame(op);
  }
  void wait(unsigned n) {
    started=true;
    while (n>=800-samples) { n-=800-samples; samples=0; flush(); }
    samples+=n;
  }
  const std::vector<uint8_t>& finish() {
    if (samples || frameOps || !frames) flush();
    bytes.clear(); bytes.reserve(40+pools.size()*16+payloads.size()+initialOps.size()+frameData.size());
    bytes.insert(bytes.end(),{'S','G','U','B',2,16,(uint8_t)repeat,0});
    be32(bytes,frames); be32(bytes,(uint32_t)pools.size()); be32(bytes,initialCount);
    uint32_t directory=40, payloadStart=directory+(uint32_t)pools.size()*16;
    uint32_t initStart=payloadStart+(uint32_t)payloads.size();
    uint32_t frameStart=initStart+(uint32_t)initialOps.size();
    be32(bytes,directory); be32(bytes,initStart); be32(bytes,frameStart); be32(bytes,0); be16(bytes,40); be16(bytes,0);
    uint32_t offset=payloadStart;
    for (size_t i=0;i<pools.size();i++) {
      const auto& pool=pools[i];
      be16(bytes,(uint16_t)i); bytes.push_back(pool.type); bytes.push_back(0);
      be32(bytes,(uint32_t)pool.data.size()); be32(bytes,(uint32_t)pool.data.size()); be32(bytes,offset); offset+=(uint32_t)pool.data.size();
    }
    bytes.insert(bytes.end(),payloads.begin(),payloads.end());
    bytes.insert(bytes.end(),initialOps.begin(),initialOps.end());
    bytes.insert(bytes.end(),frameData.begin(),frameData.end());
    uint32_t total=(uint32_t)bytes.size();
    for (int i=0;i<4;i++) bytes[32+i]=(total>>(24-8*i))&255;
    if (bytes.size()>limit) throw std::runtime_error("SGUB exceeds the 14 MiB limit");
    return bytes;
  }
};
#endif
