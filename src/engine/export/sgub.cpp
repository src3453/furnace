/**
 * Furnace Tracker - multi-system chiptune tracker
 * Copyright (C) 2021-2026 tildearrow and contributors
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include "sgub.h"
#include "sgubWriter.h"
#include "../engine.h"
#include <algorithm>
#include <map>
#include <memory>
#include <stdexcept>

void DivExportSGUB::run() {
  try {
    if (!e->isROMExportViable(DIV_ROM_SGUB)) throw std::runtime_error("SGUB requires one S3HS chip or one/two S3W2 chips and no other systems");
    const bool s3hs=e->song.systemLen==1 && e->song.system[0]==DIV_SYSTEM_S3HS;
    DivSGUBWriter stream(conf.getBool("loop",false));
    if (s3hs) stream.setVersion(3);
    stream.reserve(s3hs?0x500000:0x100000*5+0x10000);
    e->stop();
    e->setOrder(0);
    e->synchronizedSoft([&]() {
      const double origRate=e->got.rate;
      const bool origRepeat=e->repeatPattern;
      const bool origExporting=e->exporting;
      // Always restore engine/dump state, including on cancellation or size errors.
      auto cleanup=[&]() {
        e->got.rate=origRate;
        e->repeatPattern=origRepeat;
        e->exporting=origExporting;
        for (int i=0; i<e->song.systemLen; i++) {
          e->disCont[i].dispatch->toggleRegisterDump(false);
          e->disCont[i].dispatch->getRegisterWrites().clear();
        }
        e->remainingLoops=-1;
        e->playing=false;
        e->freelance=false;
        e->extValuePresent=false;
      };
      try {
        e->got.rate=48000;
        e->repeatPattern=false;
        e->exporting=true;
        e->curOrder=0;
        e->freelance=false;
        e->playing=false;
        e->remainingLoops=-1;
        e->playSub(false);
        if (s3hs) {
          auto* dispatch=e->disCont[0].dispatch;
          const auto* signedPCM=static_cast<const int8_t*>(dispatch->getSampleMem(0));
          const size_t pcmSize=dispatch->getSampleMemUsage(0);
          if (!signedPCM || pcmSize>0x400000) throw std::runtime_error("Invalid S3HS PCM RAM image");
          if (pcmSize) {
            std::vector<uint8_t> pcm(pcmSize);
            for (size_t i=0;i<pcmSize;i++) pcm[i]=(uint8_t)(signedPCM[i]+128);
            const uint16_t pool=stream.addPool(3,pcm.data(),pcm.size());
            stream.initialCopy(3,0,pool,0,(uint32_t)pcmSize);
          }
          dispatch->getRegisterWrites().clear();
          dispatch->toggleRegisterDump(true);
          logAppend("Recording S3HS writes at 48 kHz; quantizing to 60 Hz SGUB v3 frames.");
          e->calcSongTimestamps();
          const uint64_t totalTicks=std::max<uint64_t>(1,e->curSubSong->ts.totalTicks);
          uint64_t recordedTicks=0;
          while (true) {
            if (mustAbort) throw std::runtime_error("SGUB export cancelled");
            if (e->nextTick(false,true) || !e->playing) break;
            recordedTicks++;
            auto& writes=dispatch->getRegisterWrites();
            for (const auto& write:writes) {
              if (write.addr>=0x400) throw std::runtime_error("Invalid S3HS register write address");
              stream.write(3,write.addr,(unsigned)write.val&255);
            }
            writes.clear();
            if (e->cycles>0) stream.wait(e->cycles);
            progress=0.2f+0.75f*std::min(1.0f,(float)recordedTicks/(float)totalTicks);
          }
          cleanup();
          return;
        }
        const uint8_t* pcm[2]={nullptr,nullptr};
        size_t used[2]={0,0};
        for (int i=0; i<e->song.systemLen; i++) {
          auto* dispatch=e->disCont[i].dispatch;
          pcm[i]=static_cast<const uint8_t*>(dispatch->getSampleMem(0));
          used[i]=dispatch->getSampleMemUsage(0);
          if (!pcm[i] || dispatch->getSampleMemCapacity(0)!=0x100000 || used[i]>0x100000) {
            throw std::runtime_error("Invalid S3W2 PCM RAM");
          }
        }
        unsigned pcmBase[2]={0,0};
        std::vector<uint8_t> sharedPCM(0x100000,0);
        bool conflicts=pcm[1] && !std::equal(pcm[0],pcm[0]+std::min(used[0],used[1]),pcm[1]);
        if (conflicts) {
          // Keep a silent guard after chip 0: the stub reads the end address
          // before stopping, so adjacent banks must not share that byte.
          pcmBase[1]=(used[0]+1+15)&~size_t(15);
          if (pcmBase[1]+used[1]>=sharedPCM.size()) {
            throw std::runtime_error("Combined S3W2 PCM banks exceed the shared 1 MiB RAM (including end guards)");
          }
          std::copy(pcm[0],pcm[0]+used[0],sharedPCM.begin());
          std::copy(pcm[1],pcm[1]+used[1],sharedPCM.begin()+pcmBase[1]);
          logAppend("Relocating the second S3W2 PCM bank into shared CPT32 RAM.");
        } else {
          int source=used[1]>used[0]?1:0;
          std::copy(pcm[source],pcm[source]+0x100000,sharedPCM.begin());
        }
        unsigned pcmAddress[2][8][3]={};
        // Keep silence implicit and store only used PCM bytes in one raw pool.
        // Furnace's signed samples map to CPT32's unsigned representation.
        unsigned pcmSize=conflicts?pcmBase[1]+used[1]:std::max(used[0],used[1]);
        std::vector<uint8_t> packedPCM(pcmSize);
        for (unsigned i=0;i<pcmSize;i++) packedPCM[i]=(sharedPCM[i]+128)&255;
        uint16_t pcmPool=0;
        if (pcmSize) pcmPool=stream.addPool(2,packedPCM.data(),packedPCM.size());
        stream.initialFill(2,0,0x100000,0x80);
        if (pcmSize) stream.initialCopy(2,0,pcmPool,0,pcmSize);
        progress=0.2f;
        const std::vector<uint8_t> resetWave(256,0x80);
        const uint16_t resetWavePool=stream.addPool(1,resetWave.data(),resetWave.size());
        // Reset both SGUs in frame zero so whole-song repetition also resets them.
        // regPool is NOT a reset snapshot (pan and wavetable defaults differ).
        for (unsigned chip=0; chip<2; chip++) {
          for (unsigned ch=0; ch<8; ch++) {
            unsigned base=0x800+ch*32;
            stream.write(chip,base+3,0);
            stream.write(chip,base+0x19,0);
            for (unsigned reg=0; reg<32; reg++) stream.write(chip,base+reg,reg==4?0xff:0);
            stream.copy(chip,ch*256,resetWavePool,0,256);
            // The play bit gates wave/noise too. Match chip reset's active=true
            // after clearing PCM state; volume remains zero until a note-on.
            stream.write(chip,base+0x19,1);
          }
        }
        for (int i=0; i<e->song.systemLen; i++) {
          e->disCont[i].dispatch->getRegisterWrites().clear();
          e->disCont[i].dispatch->toggleRegisterDump(true);
        }
        logAppend("Recording S3W2 writes at 48 kHz; quantizing to 60 Hz SGUB frames.");
        std::map<std::vector<uint8_t>,uint16_t> wavePools;
        wavePools.emplace(resetWave,resetWavePool);
        e->calcSongTimestamps();
        uint64_t totalTicks=std::max<uint64_t>(1,e->curSubSong->ts.totalTicks);
        uint64_t recordedTicks=0;
        progress=0.2f;
        bool warnedReservedMode=false;
        while (true) {
          if (mustAbort) throw std::runtime_error("SGUB export cancelled");
          // The terminal tick is loop lookahead; do not export its note-ons.
          if (e->nextTick(false,true) || !e->playing) break;
          recordedTicks++;
          for (int i=0; i<e->song.systemLen; i++) {
            auto& writes=e->disCont[i].dispatch->getRegisterWrites();
            for (size_t writeIndex=0; writeIndex<writes.size();) {
              const auto& write=writes[writeIndex];
              if (write.addr<0x800 && (write.addr&255)==0 && writeIndex+256<=writes.size()) {
                bool complete=true;
                std::vector<uint8_t> wave(256);
                for (unsigned pos=0;pos<256;pos++) {
                  const auto& sampleWrite=writes[writeIndex+pos];
                  if (sampleWrite.addr!=write.addr+pos) { complete=false; break; }
                  wave[pos]=sampleWrite.val&255;
                }
                if (complete) {
                  auto found=wavePools.find(wave);
                  uint16_t poolId;
                  if (found==wavePools.end()) {
                    poolId=stream.addPool(1,wave.data(),wave.size());
                    wavePools.emplace(std::move(wave),poolId);
                  } else poolId=found->second;
                  stream.copy(i,write.addr,poolId,0,256);
                  writeIndex+=256;
                  continue;
                }
              }
              unsigned value=write.val&255;
              if (write.addr>=0x800 && (write.addr&31)==5) {
                if ((value>>3)>4 && !warnedReservedMode) {
                  logAppend("WARNING: reserved Furnace S3W2 stub modulation mode exported as None.");
                  warnedReservedMode=true;
                }
                value=divSGUBModulation(value);
              }
              unsigned reg=write.addr&31;
              if (pcmBase[i] && write.addr>=0x800 && reg>=0x10 && reg<=0x18) {
                unsigned ch=(write.addr-0x800)/32;
                unsigned field=(reg-0x10)/3;
                unsigned byte=(reg-0x10)%3;
                unsigned shift=(2-byte)*8;
                unsigned& address=pcmAddress[i][ch][field];
                address=(address&~(255u<<shift))|(value<<shift);
                unsigned relocated=address+pcmBase[i];
                if (relocated>=0x100000) throw std::runtime_error("Relocated S3W2 PCM pointer exceeds shared RAM");
                unsigned base=0x800+ch*32+0x10+field*3;
                // Re-emit the complete address to propagate carries between
                // bytes, regardless of the source's byte-write order.
                stream.write(i,base+2,relocated&255);
                stream.write(i,base+1,(relocated>>8)&255);
                stream.write(i,base,(relocated>>16)&255);
              } else {
                stream.write(i,write.addr,value);
              }
              writeIndex++;
            }
            writes.clear();
          }
          if (e->cycles>0) stream.wait(e->cycles);
          progress=0.2f+0.75f*std::min(1.0f,(float)recordedTicks/(float)totalTicks);
        }
      } catch (...) {
        cleanup();
        throw;
      }
      cleanup();
    });
    const auto& bytes=stream.finish();
    progress=0.95f;
    std::unique_ptr<SafeWriter> writer(new SafeWriter);
    writer->init();
    writer->write(bytes.data(),bytes.size());
    output.push_back(DivROMExportOutput("song.sgub",writer.release()));
    progress=1.0f;
    logAppend("SGUB export finished.");
  } catch (const std::exception& error) {
    e->lastError=error.what();
    logAppend(String("ERROR: ")+error.what());
    failed=true;
  }
  running=false;
}

bool DivExportSGUB::go(DivEngine* eng) {
  if (running) return false;
  wait();
  e=eng;
  failed=false;
  mustAbort=false;
  progress=0.0f;
  running=true;
  exportThread=std::thread(&DivExportSGUB::run,this);
  return true;
}
void DivExportSGUB::wait() {
  if (exportThread.joinable()) exportThread.join();
}
void DivExportSGUB::abort() {
  mustAbort=true;
  wait();
}
DivROMExportProgress DivExportSGUB::getProgress(int index) {
  return {index==0?String("Initialize and record S3W2"):String(""),index==0?progress.load():0.0f};
}
