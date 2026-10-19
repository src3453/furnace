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

#ifndef _DIV_EXPORT_SGUB_H
#define _DIV_EXPORT_SGUB_H
#include "../export.h"
#include <atomic>
#include <thread>

class DivExportSGUB: public DivROMExport {
  DivEngine* e=nullptr;
  std::thread exportThread;
  std::atomic<bool> running{false}, failed{false}, mustAbort{false};
  std::atomic<float> progress{0.0f};
  void run();
public:
  bool go(DivEngine* eng);
  bool isRunning() { return running.load(); }
  bool hasFailed() { return failed.load(); }
  void abort();
  void wait();
  DivROMExportProgress getProgress(int index=0);
  ~DivExportSGUB() { abort(); }
};
#endif
