"""Real-engine render check for S3HS (3HS88PWN4) and 3HS816.

    python test/3hs/render_check.py build/Release/furnace.exe [--reference OLD_furnace.exe]

Builds FUR fixtures in a TemporaryDirectory (INF2/SNG2/INS2/WAVE/PATN, as
test/sgub/smoke.py does), renders them with `-output`, and checks that:

- both systems sound, at the expected pitch, on every channel,
- an S3HS instrument sounds the same on 3HS816 (shared instrument type),
- 3HS816 panning and its user wave A (from song wavetable 0) work,
- with --reference, the S3HS render matches an older Furnace build (the
  float 3HS88PWN4 core) closely. The integer model is not bit-exact; the
  CPT32 repository's hs88_legacy_compare holds the detailed tolerances.
"""
import argparse
import math
import struct as st
import subprocess
import sys
import tempfile
import wave
from pathlib import Path

VERSION = 248
S3HS = (0xFC, 12)
HS816 = (0xFB, 16)
INS_S3HS = 68
NOTE_A4 = 117  # PATN note number (internal note 57)


def pack(fmt, *values):
    return st.pack("<" + fmt, *values)


def block(tag, data):
    return tag + pack("I", len(data)) + data


def feature_3h(algorithm=0, wave=0, volumes=(255, 0, 0, 0, 0, 0, 0, 0), ratios=(0x1000,) * 7):
    """The INS2 "3H" feature: ratios, waves, op volume/ADSR, mode, fb, fbDest."""
    data = bytes(r >> 8 for r in ratios) + bytes(r & 255 for r in ratios)
    data += bytes([wave] * 8)
    for volume in volumes:
        data += bytes([volume, 0, 0, 255, 0])  # instant attack, full sustain, instant release
    data += bytes([algorithm]) + pack("b", 0) + bytes([0])
    return b"3H" + pack("H", len(data)) + data


def pcm_instrument():
    """An Amiga-type (sample) instrument using sample 0, as test/sgub/smoke.py."""
    return b"SM" + pack("HhBB", 4, 0, 2, 255)


def song(system, notes, instrument, wavetable=None, pan=None, sample=None, tuning=440):
    """One chip, one pattern per channel; notes = {channel: note}."""
    chip_id, channels = system
    info = b"render check\0" + b"\0" * 7 + pack("fBfHH", tuning, 1, 1., channels, 1)
    info += pack("HHfff", chip_id, channels, 1., 0., 0.)
    info += pack("IB", 0, 1)
    sub = pack("fBBHHBBHHB", 60, 1, 1, 16, 1, 4, 16, 1, 1, 1)
    sub += pack("16H", 6, *([0] * 15)) + b"\0\0"
    sub += bytes(channels) + bytes([1] * channels) + bytes(channels * 8)
    ins_type = 4 if instrument.startswith(b"SM") else INS_S3HS
    ins = pack("HH", VERSION, ins_type) + instrument + b"EN" + pack("H", 0)
    elements = [(1, block(b"SNG2", sub)), (4, block(b"INS2", ins))]
    if wavetable is not None:
        data = b"user A\0" + pack("III", len(wavetable), 0, 255) + pack(f"{len(wavetable)}I", *wavetable)
        elements.append((5, block(b"WAVE", data)))
    if sample is not None:
        data = b"PCM\0" + pack("III4Bii4I", len(sample), 8363, 8363, 8, 0, 0, 0, 0, len(sample), 1, 0, 0, 0) + sample
        elements.append((6, block(b"SMP2", data)))
    for channel in range(channels):
        row = b""
        if channel in notes:
            # note, instrument 0, volume FF (+ effect 08xx panning)
            row = bytes([7, notes[channel], 0, 255])
            if pan is not None:
                row = bytes([0x1F, notes[channel], 0, 255, 0x08, pan])
        pat = pack("BHH", 0, channel, 0) + b"\0" + row + b"\xff"
        elements.append((7, block(b"PATN", pat)))
    groups = {}
    for typ, data in elements:
        groups.setdefault(typ, []).append(data)
    table_size = sum(5 + 4 * len(group) for group in groups.values()) + 1
    offset = 32 + 8 + len(info) + table_size
    table = payload = b""
    for typ, group in groups.items():
        table += pack("BI", typ, len(group))
        for data in group:
            table += pack("I", offset)
            offset += len(data)
            payload += data
    header = b"-Furnace module-" + pack("HHI", VERSION, 0, 32) + bytes(8)
    return header + block(b"INF2", info + table + b"\0") + payload


def render(exe, fur, out):
    out.unlink(missing_ok=True)
    run = subprocess.run([str(exe), "-loops", "0", "-output", str(out), str(fur)],
                         capture_output=True, timeout=120)
    if run.returncode != 0 or not out.exists():
        raise SystemExit(f"render failed ({run.returncode}): {fur.name}\n{run.stdout.decode(errors='replace')}"
                         f"{run.stderr.decode(errors='replace')}")
    with wave.open(str(out)) as w:
        assert w.getnchannels() == 2 and w.getsampwidth() == 2, "expected 16-bit stereo"
        rate = w.getframerate()
        raw = w.readframes(w.getnframes())
    samples = st.unpack(f"<{len(raw) // 2}h", raw)
    return rate, samples[0::2], samples[1::2]


def window(samples, rate, start=0.05, length=0.2):
    return samples[int(start * rate):int((start + length) * rate)]


def peak(samples):
    return max((abs(s) for s in samples), default=0)


def frequency(samples, rate):
    """Rising zero crossings, interpolated, from the first to the last."""
    xs = [i + (-a) / (b - a) for i, (a, b) in enumerate(zip(samples, samples[1:])) if a < 0 <= b]
    return (len(xs) - 1) * rate / (xs[-1] - xs[0]) if len(xs) > 1 else 0.0


def snr(reference, test):
    signal = sum(r * r for r in reference)
    error = sum((r - t) ** 2 for r, t in zip(reference, test))
    return math.inf if error == 0 else 10 * math.log10(signal / error)


def shape_snr(reference, test):
    """SNR after fitting the test's gain: compares waveform, not level."""
    energy = sum(t * t for t in test)
    gain = sum(r * t for r, t in zip(reference, test)) / energy if energy else 0.0
    return snr(reference, [gain * t for t in test])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("furnace")
    parser.add_argument("--reference", help="older furnace.exe with the float S3HS core")
    args = parser.parse_args()
    exe = Path(args.furnace).resolve()
    failures = []

    def check(condition, message):
        print(("ok   " if condition else "FAIL ") + message)
        if not condition:
            failures.append(message)

    sine = feature_3h()
    fm = feature_3h(algorithm=1, volumes=(255, 0, 0, 0, 200, 0, 0, 0),
                    ratios=(0x1000, 0x1000, 0x1000, 0x2000, 0x1000, 0x1000, 0x1000))
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)

        def run(name, data, binary=exe):
            fur = tmp / f"{name}.fur"
            fur.write_bytes(data)
            return render(binary, fur, tmp / f"{name}-{'new' if binary == exe else 'ref'}.wav")

        rate, left, _ = run("s3hs", song(S3HS, {0: NOTE_A4}, sine))
        w = window(left, rate)
        check(peak(w) > 2000, f"S3HS sine sounds (peak {peak(w)})")
        check(abs(frequency(w, rate) / 440 - 1) < 0.003, f"S3HS A4 at 440 Hz ({frequency(w, rate):.2f})")

        rate, left, _ = run("hs816", song(HS816, {0: NOTE_A4}, sine))
        w816 = window(left, rate)
        check(peak(w816) > 2000, f"3HS816 sine sounds (peak {peak(w816)})")
        check(abs(frequency(w816, rate) / 440 - 1) < 0.003, f"3HS816 A4 at 440 Hz ({frequency(w816, rate):.2f})")

        # The frequency register is exact Hz. The retired 192 kHz period
        # quantization played A7 (3520 Hz) at 192000 / 54 = 3555.6 Hz.
        for system, name in ((S3HS, "S3HS"), (HS816, "3HS816")):
            rate, left, _ = run(f"{name}-a7", song(system, {0: NOTE_A4 + 36}, sine))
            hz = frequency(window(left, rate), rate)
            check(abs(hz / 3520 - 1) < 0.003, f"{name} A7 at 3520 Hz, not 3555.6 ({hz:.1f})")
        # S3HS's note-on scales operator volume by vol/256, 3HS816 by vol/255.
        check(shape_snr(w, w816) > 60, f"shared instrument: S3HS and 3HS816 match ({shape_snr(w, w816):.1f} dB)")

        for channel in (7, 15):
            rate, left, _ = run(f"hs816-ch{channel}", song(HS816, {channel: NOTE_A4}, sine))
            cw = window(left, rate)
            check(snr(w816, cw) > 60, f"3HS816 channel {channel + 1} matches channel 1 ({snr(w816, cw):.1f} dB)")

        rate, left, right = run("hs816-pan", song(HS816, {0: NOTE_A4}, sine, pan=0xF0))
        check(peak(window(right, rate)) == 0 and peak(window(left, rate)) > 2000,
              f"3HS816 hard-left panning (L {peak(window(left, rate))}, R {peak(window(right, rate))})")

        rate, left, _ = run("hs816-fm", song(HS816, {0: NOTE_A4}, fm))
        check(snr(w816, window(left, rate)) < 20, "3HS816 FM algorithm changes the timbre")

        # User wave A (wave 10) from song wavetable 0. A half-high, half-low
        # table plays (nearly) the built-in square, wave 3. Furnace's output
        # removes DC, so a constant table cannot be checked directly.
        rate, left, _ = run("hs816-square", song(HS816, {0: NOTE_A4}, feature_3h(wave=3)))
        square = window(left, rate)
        rate, left, _ = run("hs816-userwave", song(HS816, {0: NOTE_A4}, feature_3h(wave=10),
                                                    wavetable=[255] * 16 + [0] * 16))
        uw = window(left, rate)
        check(shape_snr(square, uw) > 10 and snr(w816, uw) < 10,
              f"3HS816 user wave A from wavetable 0 ({shape_snr(square, uw):.1f} dB vs square)")

        if args.reference:
            reference = Path(args.reference).resolve()
            # 8-bit signed PCM: a decaying 2-partial tone
            pcm = bytes(int(90 * math.exp(-i / 6000) * (math.sin(i * 0.11) + 0.4 * math.sin(i * 0.37))) & 255
                        for i in range(16384))
            for name, notes, instrument, sample, limit in (
                    ("ab-sine", {0: NOTE_A4, 3: NOTE_A4 - 12}, sine, None, 40),
                    ("ab-fm", {1: NOTE_A4}, fm, None, 20),
                    ("ab-pcm", {8: NOTE_A4 - 24}, pcm_instrument(), pcm, 30)):
                # The float core quantized pitch to 192000 / floor(192000 / F);
                # with A4 = 480 Hz the notes used (480, 240, 120 Hz) divide
                # 192000, so both cores play the same pitch.
                data = song(S3HS, notes, instrument, sample=sample, tuning=480)
                _, new_left, _ = run(name, data)
                _, old_left, _ = run(name, data, reference)
                # Identical output is possible: after resampling to the
                # output rate and 16-bit rounding the cores rarely differ.
                value = snr(old_left, new_left)
                check(peak(old_left) > 2000 and value > limit,
                      f"S3HS {name} vs float core: {value:.1f} dB over {len(new_left)} frames (> {limit})")

    if failures:
        print(f"{len(failures)} check(s) failed")
        return 1
    print("all render checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
