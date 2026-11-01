/**
 * Furnace Tracker - multi-system chiptune tracker
 * Copyright (C) 2021-2025 tildearrow and contributors
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

#ifndef _3HS816_H
#define _3HS816_H

#include <vector>
#include "../dispatch.h"
#include "sound/3hs/3hs816/include/hs816/registers.hpp"

// 3HS816: 16 harmonic channels of the 3HS88PWN4, no PCM, no effects. Uses
// the 3HS88PWN4 (S3HS) instrument type. Song wavetables 0..3 become the
// chip's user waves A..D.
class DivPlatform3HS816: public DivDispatch {
  struct Channel: public SharedChannel {
    int pan;
    int mode;
    int opVolumes[8];
    Channel(bool linear=true):
      SharedChannel(255,linear),
      pan(0xff),
      mode(0),
      opVolumes{0,0,0,0,0,0,0,0} {}
  };
  Channel chan[16];
  DivDispatchOscBuffer* oscBuf[16];
  bool isMuted[16];
  hs816::Engine* chip;
  hs816::RegisterInterface* chipRegs;
  std::vector<int> oscTap[16];
  DivPitchTable pitchTable;
  unsigned char regPool[hs816::reg::kSize];

  void chipWrite(unsigned int addr, unsigned char val);
  void writeOperatorVolumes(int ch);
  void writeInstrument(int ch, DivInstrument* ins);
  void uploadUserWave(int wave);

  friend void putDispatchChip(void*,int);
  friend void putDispatchChan(void*,int,int);
  public:
    void acquire(short** buf, size_t len);
    int dispatch(DivCommand c);
    SharedChannel* getChanState(int chan);
    DivMacroInt* getChanMacroInt(int ch);
    DivDispatchOscBuffer* getOscBuffer(int chan);
    unsigned short getPan(int chan);
    unsigned char* getRegisterPool();
    int getRegisterPoolSize();
    int getOutputCount();
    void reset();
    void tick(bool sysTick=true);
    void muteChannel(int ch, bool mute);
    void notifyWaveChange(int wave);
    void notifyInsDeletion(void* ins);
    void notifyPitchTable(int sample=-1);
    void poke(unsigned int addr, unsigned short val);
    void poke(std::vector<DivRegWrite>& wlist);
    int init(DivEngine* parent, int channels, int sugRate, const DivConfig& flags);
    void quit();
    ~DivPlatform3HS816();
};

#endif
