# Development notes

## S3W2 SGUB export (2026-10-18)

The S3W2 core in this fork is a **stub**, not the normative CPT32 hardware
implementation. Consult the CPT32 checkout's `doc/spec_SGU.md` and
`doc/spec_SGUB.md`; never assume matching enum numbers imply matching sound.

- ROM export accepts only one or two S3W2 chips. Explicit `.sgub` CLI filenames
  cannot fall back to ZSM or another target when the chip set is incompatible.
- Accumulate engine ticks at 48 kHz and quantize to 800-sample / 60 Hz frames.
  Preserve write order and duplicate phase strobes. SGUB repeats from frame zero;
  its initialization section is applied once. PCM_CONTROL gates wave/noise too,
  so reset must re-enable the play bit at zero volume or wave output is silent.
- Stub modulation IDs 3/4 are Slew; map to CPT32 5/6. CPT32 keeps Hard Sync/Window
  at 3/4. Reserved stub values are no-ops and map to None with a warning.
- Convert signed PCM with +128. Compatible banks share offsets; conflicting
  bank 1 is relocated with an end guard and all start/end/loop pointer byte
  carries are rewritten. Banks plus guards must fit 1 MiB. No per-sample dedup
  is attempted across conflicting banks.
- Pitch uses SharedChannel::calcFreq and Hz-scaled tables (fold /32 into table
  initialization). Applying /32 after non-linear macro offsets changes their
  units. PCM alone needs centerRate/(standard C4 Hz * 32); retain song tuning.
  Clear key transition flags after ticking to avoid repeated PCM triggers.
- Tests: `test/sgub/CMakeLists.txt` / CTest for the standalone encoder and
  `python test/sgub/smoke.py build/Release/furnace.exe` for real-engine export.
  The smoke creates temporary documented FUR fixtures, validates full SGUB
  structures, checks pitch values, Slew translation and PCM relocation, and
  rejects combined banks exceeding capacity and incompatible chip sets.
- Windows GUI-subsystem executables may reopen CONOUT$ when attaching to the
  parent; subprocess capture can then be empty. Failure smoke also requests
  text output to prove the source fixture loaded before checking no ROM exists.
- The initial PCM RAM image contains 1,048,576 records. Avoid calling the
  per-record vector insert path for this bulk section: resize once per chunk
  and fill records directly. Keep export progress moving during initialization
  and song capture; a constant zero progress bar made a healthy worker look hung.

Existing uncommitted S3W2/GUI/core changes were retained. Release build and
focused checks passed; the sibling CPT32 tests also validate real exported
wave audio, repeat resets, core Slew semantics and existing SGUB playback.

## SGUB v2 pooled ROM export (2026-10-18)
- The SGUB writer emits v2: 40-byte header, pool directory and payloads, initial FILL/COPY_POOL operations, and variable-length frame operations. The 14 MiB bound includes all sections; per-frame operation count remains u16.
- Store used PCM bytes after signed-to-unsigned conversion, initialize the shared 1 MiB PCM RAM with one FILL(0x80), then COPY_POOL the used range. Pool repeated 256-byte waveforms by content; COPY_POOL may target only aligned 256-byte wave SRAM ranges.
- PCM bank relocation remains required before packing when two chips conflict. Coalesce captured wave uploads only when all 256 writes are contiguous and aligned; partial or interleaved writes must preserve their original WRITE order.
- `test/sgub` checks the v2 wire layout, pooling, timing, repeated writes, silence frames, size and address rejection. Rebuild Furnace and rerun it whenever operation sizes or section offsets change.

- CPT32 guest v2 playback was repaired: all directory, init, frame and payload offsets must be rebased from file offsets to STREAM_BASE at 0x00200000. A complete user Furnace export (4,995 bytes, two wave pools, 1 MiB PCM silence FILL, 131 Hz note) now passes Rust audio and compiled PipelineCpu guest tests. Build the updated guest with -O speed; stale source/binary copies beside existing SGUB files must be replaced, but the SGUB itself is valid and needs no re-export.

## SQV4 mono samples and S3W2 (2026-10-19)
- The normative SQV4 v1/v2 contract is in the CPT32 sibling checkout's
  doc/squeezvox/squeezvox4_spec_stub.md. Predictor/index starts at zero once
  per file, continues across 64-sample fragments, and gains do not feed back.
  L+ must select [10,11,12,14] for v1 and [1,2,3,14] for v2. New encoders emit v2.
- Furnace imports mono only; reject stereo rather than silently dropping a channel.
  IMA ADPCM is the existing editor type for H. Its WAV-style initial-state block
  is not SQV4 wire data: H import/upload adapts via decoded PCM with re-encoding.
  L/L+ use fork depths 17/18, storing the complete container in SMP2 and undo data.
- S3W2 waveform 4 consumes a complete mono file with an exclusive byte end.
  Validate exact length and padding before playback. A loop points to a fragment's
  volume byte; replay predictor history, never initialize a decoder at the loop.
  Loop starts round down to 64 samples; truncate/pad loop ends canonically.
- The existing sample RAM API/export adds 128 to signed bytes. Store compressed
  containers XOR 0x80 in that image, so both core upload and SGUB export recover
  the original bytes. Do not add 128 directly to an unadapted SQV4 header/payload.
  Use start as the compressed no-loop pointer, so two-chip relocation stays valid.
- Focused tests: test/sqv4 CMake/CTest covers golden H/L+ vectors, v1 codebook,
  malformed inputs, conversions, FUR persistence, undo and incremental playback.
  test/sqv4/msvc_engine.py reuses a Release MSVC build for actual import/export,
  stereo rejection and S3W2 uploads. Use the matching VS developer toolset.
  test/sqv4/sgub_smoke.py validates compressed banks/loop relocation and runs
  the existing SGUB matrix through a v2 parser. Release build and these checks passed.
