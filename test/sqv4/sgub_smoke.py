"""SQV4 real-engine SGUB v2 checks; also runs the existing SGUB smoke matrix.

python test/sqv4/sgub_smoke.py build/Release/furnace.exe
"""
from pathlib import Path
import importlib.util
import struct
import subprocess
import sys
import tempfile

spec = importlib.util.spec_from_file_location("sgub_smoke", Path(__file__).parents[1] / "sgub/smoke.py")
sgub = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sgub)

def parse(path):
    data = path.read_bytes()
    assert data[:6] == b"SGUB\x02\x10"
    frames, pools, initial, directory, init_offset, frame_offset, size = struct.unpack_from(">7I", data, 8)
    assert size == len(data) and directory == 40
    payloads = []
    for i in range(pools):
        ident, typ, flags, raw, stored, offset = struct.unpack_from(">HBBIII", data, directory + i * 16)
        assert ident == i and typ in (1, 2) and flags == 0 and raw == stored
        payloads.append(data[offset:offset + stored])
    def operations(offset, count):
        out = []
        for _ in range(count):
            op, target = data[offset:offset + 2]
            address = int.from_bytes(data[offset + 2:offset + 5], "big")
            if op == 0:
                values = data[offset + 5:offset + 6]; offset += 6
            elif op == 1:
                length = int.from_bytes(data[offset + 5:offset + 9], "big")
                values = data[offset + 9:offset + 10] * length; offset += 10
            else:
                assert op == 2
                pool, source, length = struct.unpack_from(">HII", data, offset + 5)
                values = payloads[pool][source:source + length]; offset += 15
                assert len(values) == length
            out.extend((target, address + i, value) for i, value in enumerate(values))
        return out, offset
    init, offset = operations(init_offset, initial)
    assert offset == frame_offset
    memory = bytearray(0x100000)
    for target, address, value in init:
        assert target == 2
        memory[address] = value
    result = []
    for _ in range(frames):
        count, = struct.unpack_from(">H", data, offset); offset += 2
        records, offset = operations(offset, count); result.append(records)
    assert offset == len(data)
    return data, [(2, i, value) for i, value in enumerate(memory)], result

def replace_samples(data, depth, container):
    # Read and rebuild the small FUR element table while replacing SMP2 blocks.
    end = data.index(b"\0", 40) + 1 + 7
    channels, chips = struct.unpack_from("<HH", data, end + 9)
    table = end + 13 + chips * 16 + 5
    groups = []; pos = table
    while data[pos]:
        typ = data[pos]; count, = struct.unpack_from("<I", data, pos + 1); pos += 5
        blocks = []
        for _ in range(count):
            ptr, = struct.unpack_from("<I", data, pos); pos += 4
            length, = struct.unpack_from("<I", data, ptr + 4)
            block = data[ptr:ptr + length + 8]
            if typ == 6:
                name = block.index(b"\0", 8) + 1
                header = bytearray(block[8:name + 40])
                struct.pack_into("<I", header, name - 8, 65)
                header[name - 8 + 12] = depth
                struct.pack_into("<ii", header, name - 8 + 16, 64, 65)
                image=bytearray(container)
                if blocks: image[16]=129  # force incompatible bank layouts
                block = sgub.block(b"SMP2", header + image)
            blocks.append(block)
        groups.append((typ, blocks))
    payload = bytearray(); entries = bytearray(); ptr = pos + 1
    for typ, blocks in groups:
        entries += struct.pack("<BI", typ, len(blocks))
        for block in blocks:
            entries += struct.pack("<I", ptr); payload += block; ptr += len(block)
    return data[:32] + sgub.block(b"INF2", data[40:table] + entries + b"\0") + payload

def main():
    exe = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="sqv4-sgub-") as temp:
        root = Path(temp)
        for variant, depth in [(0, 19), (1, 17), (2, 18)]:
            # Two fragments; non-unity fragment gain must survive export exactly.
            container = b"SQV4" + bytes([2, variant, 1, 0]) + struct.pack("<II", 8363, 65)
            payload=32 if variant==0 else 24
            container += bytes([128]) + bytes(payload) + bytes([255, 1]) + bytes(payload-1)
            for chips in (1, 2):
                fixture = replace_samples(sgub.song(chips=chips, pcm=True, conflict=chips==2), depth, container)
                source = root / "source.fur"; dest = root / "output.sgub"
                source.write_bytes(fixture)
                if dest.exists(): dest.unlink()
                subprocess.run([str(exe), "-romout", str(dest), str(source)], cwd=root, capture_output=True, timeout=60)
                assert dest.exists()
                _, init, frames = parse(dest)
                memory = bytes(value for _, _, value in init)
                assert memory[:len(container)] == container
                registers = [{}, {}]
                for frame in frames:
                    for target, addr, value in frame:
                        if target < 2: registers[target][addr] = value
                for chip in range(chips):
                    regs = registers[chip]
                    def addr(reg): return sum(regs[reg+i] << (16-i*8) for i in range(3))
                    assert regs[0x802] == 4 and regs[0x819] == 3
                    start, end, loop = addr(0x810), addr(0x813), addr(0x816)
                    assert end-start == len(container) and loop-start == 17+payload
                    expected=bytearray(container)
                    if chip: expected[16]=129
                    assert memory[start:end] == expected
                    if chip: assert start>0
                print("SQV4", variant, "chips", chips, "SGUB bytes and loop pointers passed")
    sgub.parse = parse
    sgub.main()

if __name__ == "__main__":
    main()
