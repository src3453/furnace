"""Real-engine CLI regression: python test/sgub/smoke.py build/Release/furnace.exe.

Fixtures use the documented INF2/SNG2/INS2/SMP2/PATN formats and are made in a
TemporaryDirectory. No user song or config is modified.
"""
import math
import struct as st
import subprocess
import sys
import tempfile
from pathlib import Path

VERSION = 248

def pack(fmt, *values):
    return st.pack("<" + fmt, *values)

def block(tag, data):
    return tag + pack("I", len(data)) + data

def song(chips=1, hz=60, macro=True, pcm=False, relative=False, conflict=False, tuning=440, linear=True, center=8363, slew=0, mixed=False, pcm_size=1024):
    channels = chips * 8
    info = b"SGUB regression\0" + b"\0" * 7 + pack("fBfHH", tuning, 1, 1., channels, chips)
    for i in range(chips):
        info += pack("HHfff", 0x82 if mixed and i else 0xfa, 8, 1., 0., 0.)
    info += pack("IB", 0, 1)
    speed = 1
    sub = pack("fBBHHBBHHB", hz, 1, 1, 8, 1, 4, 16, 1, 1, 1)
    sub += pack("16H", speed, *([0] * 15)) + b"\0\0"
    sub += bytes(channels) + bytes([1] * channels) + bytes(channels * 8)
    ins = pack("HH", VERSION, 4 if pcm else 69)
    if pcm:
        ins += b"SM" + pack("HhBB", 4, 0, 2, 255)
    if macro:
        values = [0, 128, -128, 0] if not relative else [0, 64, 64, -128]
        ma = pack("H8B", 8, 4, 4, 255, 255, int(relative), 0x80, 0, 1) + pack("4h", *values) + b"\xff"
        ins += b"MA" + pack("H", len(ma)) + ma
    if slew:
        ma = pack("H9B", 8, 5, 1, 255, 255, 0, 0, 0, 1, slew) + b"\xff"
        ins += b"MA" + pack("H", len(ma)) + ma
    ins += b"EN" + pack("H", 0)
    wave = b"square\0" + pack("III", 256, 0, 255) + pack("256I", *([0] * 128 + [255] * 128))
    elements = [(1, block(b"SNG2", sub)), (4, block(b"INS2", ins)), (5, block(b"WAVE", wave))]
    if not linear:
        elements.append((8, block(b"CFLG", b"linearPitch=0\n\0")))
    if pcm:
        data = (bytes([128, 192, 0, 64]) * ((pcm_size + 3) // 4))[:pcm_size]
        mask = 1 if conflict else (1 << chips) - 1
        sample = b"PCM\0" + pack("III4Bii4I", len(data), 8363, center, 8, 0, 0, 0, 0, len(data), mask, 0, 0, 0) + data
        elements.append((6, block(b"SMP2", sample)))
        if conflict:
            other = b"other\0" + pack("III4Bii4I", len(data), 8363, center, 8, 0, 0, 0, 0, len(data), 2, 0, 0, 0) + bytes([1] * len(data))
            elements.append((6, block(b"SMP2", other)))
            second_ins = pack("HH", VERSION, 4) + b"SM" + pack("HhBB", 4, 1, 2, 255) + b"EN" + pack("H", 0)
            elements.append((4, block(b"INS2", second_ins)))
    for chip in range(chips):
        # A4 is note 117 in PATN, internal note 57. Hold for 8 ticks.
        pat = pack("BHH", 0, chip * 8, 0) + b"\0" + bytes([7, 117, 1 if conflict and chip else 0, 255, 255])
        elements.append((7, block(b"PATN", pat)))
    # One entry per element is supported by the reader only for unique type groups.
    groups = {}
    for typ, data in elements:
        groups.setdefault(typ, []).append(data)
    table_size = sum(5 + 4 * len(group) for group in groups.values()) + 1
    offset = 32 + 8 + len(info) + table_size
    table = b""
    payload = b""
    for typ, group in groups.items():
        table += pack("BI", typ, len(group))
        for data in group:
            table += pack("I", offset)
            offset += len(data)
            payload += data
    header = b"-Furnace module-" + pack("HHI", VERSION, 0, 32) + bytes(8)
    return header + block(b"INF2", info + table + b"\0") + payload

def parse(path):
    data = path.read_bytes()
    assert data[:8] in (b"SGUB\x01\x10\0\0", b"SGUB\x01\x10\x01\0")
    frames, initial, size = st.unpack_from(">III", data, 8)
    assert size == len(data) <= 0xe00000 and frames > 0
    offset = 20
    def records(count):
        nonlocal offset
        out = []
        for _ in range(count):
            target, a, b, c, value = data[offset:offset + 5]
            address = a << 16 | b << 8 | c
            assert target <= 2 and address <= (0xfffff if target == 2 else 0x8ff)
            out.append((target, address, value))
            offset += 5
        return out
    init = records(initial)
    result = []
    for _ in range(frames):
        count, = st.unpack_from(">H", data, offset)
        offset += 2
        result.append(records(count))
    assert offset == size
    return data, init, result

def main():
    exe = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="sgub-smoke-") as tmp:
        root = Path(tmp)
        for name, kwargs in [
            ("wave", {}), ("plain", {"macro": False}),
            ("relative", {"relative": True}),
            ("nonlinear", {"linear": False}),
            ("highpcm", {"pcm": True, "center": 16726}),
            ("dual50", {"chips": 2, "hz": 50}),
            ("fast125", {"hz": 125}),
            ("pcm", {"pcm": True}),
            ("dualpcm", {"chips": 2, "pcm": True}),
            ("tuning", {"tuning": 432}),
            ("pcmtuning", {"tuning": 432, "pcm": True}),
            ("slew", {"slew": 3, "macro": False}),
            ("slewstep", {"slew": 4, "macro": False}),
            ("reserved", {"slew": 31, "macro": False}),
            ("mixed", {"chips": 2, "mixed": True}),
            ("three", {"chips": 3}),
            ("overcapacity", {"chips": 2, "pcm": True, "conflict": True, "pcm_size": 600000}),
            ("conflict", {"chips": 2, "pcm": True, "conflict": True}),
        ]:
            source = root / (name + ".fur")
            dest = root / (name + ".sgub")
            source.write_bytes(song(**kwargs))
            text = root / (name + ".txt")
            run = subprocess.run([str(exe), "-txtout", str(text), "-loglevel", "debug", "-romout", str(dest), str(source)], cwd=root, capture_output=True, timeout=60)
            log = (run.stdout + run.stderr).decode(errors="replace")
            if name in ("overcapacity", "mixed", "three"):
                assert not dest.exists(), log
                assert text.exists(), "Conflict fixture failed to load"
                if log and name == "overcapacity":
                    assert "shared 1 MiB RAM" in log, log
                print(name, "rejected")
                continue
            assert dest.exists(), log
            raw, initial, frames = parse(dest)
            assert len(initial) == 0x100000
            if kwargs.get("pcm"):
                assert [r[2] for r in initial[:4]] == [0, 64, 128, 192]
            registers = [{}, {}]
            pitches = []
            for frame in frames:
                for target, address, value in frame:
                    if target < 2:
                        registers[target][address] = value
                pitches.append((registers[0].get(0x800, 0) << 8) | registers[0].get(0x801, 0))
            active = [pitch for pitch in pitches if pitch]
            base = (kwargs.get("center", 8363) / (440 * 2 ** (-9 / 12) * 32)) * kwargs.get("tuning", 440) if kwargs.get("pcm") else kwargs.get("tuning", 440)
            if kwargs.get("hz", 60) <= 60:
                assert abs(active[0] - base) <= 2, (name, pitches, base)
            assert len(frames) == math.ceil(8 * 60 / kwargs.get("hz", 60)), (name, len(frames))
            if name == "wave":
                assert abs(active[1] - 440 * 2 ** (1/12)) <= 2, pitches
                assert abs(active[2] - 440 * 2 ** (-1/12)) <= 2, pitches
            if name == "nonlinear":
                assert active[:4] == [440, 568, 312, 440], pitches
            if name == "relative":
                assert active[0] < active[1] < active[2] and abs(active[3] - active[0]) <= 1, pitches
            if kwargs.get("chips") == 2:
                assert registers[1].get(0x803) == 255, registers[1]
            if name == "conflict":
                relocated = sum(registers[1].get(0x810 + j, 0) << (16 - j*8) for j in range(3))
                assert relocated == 1040, relocated
                assert initial[relocated][2] == 129
                assert initial[1024][2] == 128  # silent end guard
            if kwargs.get("slew"):
                expected = {3: 5, 4: 6, 31: 0}[kwargs["slew"]]
                assert registers[0].get(0x805, 0) >> 3 == expected
            print(name, "frames", len(frames), "pitches", pitches)
            if len(sys.argv) > 2 and name == "wave":
                Path(sys.argv[2]).write_bytes(raw)
    print("SGUB engine smoke checks passed")

if __name__ == "__main__":
    main()
