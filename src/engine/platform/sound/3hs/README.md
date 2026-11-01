# 3HS816 / 3HS88PWN4 integer models (vendored)

These are copies of the CPT32 repository's `chips/3hs816` and
`chips/3hs88pwn4` (`include/` and `src/` only). The CPT32 checkout is the
source of truth, together with its `doc/spec_3HS816.md`; change the model
there, run its tests (including the comparison with the old float core), and
copy the two directories here again. Do not edit these files in place.

- `3hs816/`: the harmonic core (`hs816::Channel`, envelopes, waves) and the
  16-channel 3HS816 (`hs816::Engine`, `hs816::RegisterInterface`).
- `3hs88pwn4/`: the 3HS88PWN4 (`hs88::Engine`, `hs88::RegisterInterface`),
  built on the 3HS816 core plus PCM/WT channels, filters, EQ and compressor.

Both are integer-only. `src/generated_tables.hpp` in each comes from that
model's `tools/gen_tables.py` in the CPT32 repository.
