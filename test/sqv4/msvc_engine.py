"""Run real-engine integration using an existing MSVC Release Furnace build.

Run from a VS developer shell, after cmake --build build --config Release:
    python test/sqv4/msvc_engine.py build
Reuses the build's recorded compiler/linker settings and objects, replacing
only main.cpp. The application target and user configuration are untouched.
"""
from pathlib import Path
import re
import subprocess
import sys

root = Path(__file__).resolve().parents[2]
build = Path(sys.argv[1]).resolve()
tlog = build / "furnace.dir/Release/furnace.tlog"
out = build / "sqv4-engine"
out.mkdir(exist_ok=True)
source = root / "test/sqv4/engine.cpp"
obj = out / "engine_check.obj"
exe = out / "engine_check.exe"
lines = (tlog / "CL.command.1.tlog").read_text(encoding="utf-16").splitlines()
for i, line in enumerate(lines):
    if line.startswith("^") and line.endswith("\\SRC\\MAIN.CPP"):
        command = lines[i + 1]
        original = line[1:]
        break
else:
    raise RuntimeError("main.cpp compiler command not found")
command = command.replace(original, f'"{source}"')
command = re.sub(r'/Fo"[^"]*"', lambda _: f'/Fo"{obj}"', command)
command = command.replace('/D NDEBUG', '/UNDEBUG')
# The harness has assertions enabled and must not reuse the app PCH's macros.
command = re.sub(r'/Y[uc]"[^"]*"|/Fp"[^"]*"|/FI(?:"[^"]*"|\S+)', '', command)
subprocess.run("cl.exe " + command, cwd=build, check=True)
link = (tlog / "link.command.1.tlog").read_text(encoding="utf-16").splitlines()[1:]
command = "\n".join(link)
command = re.sub(r'/OUT:"[^"]*"', lambda _: f'/OUT:"{exe}"', command)
command = command.replace('/SUBSYSTEM:WINDOWS', '/SUBSYSTEM:CONSOLE')
command = re.sub(r'/PDB:"[^"]*"', lambda _: f'/PDB:"{out / "engine_check.pdb"}"', command)
command = re.sub(r'/IMPLIB:"[^"]*"', lambda _: f'/IMPLIB:"{out / "engine_check.lib"}"', command)
command = command.replace('FURNACE.DIR\\RELEASE\\MAIN.OBJ', f'"{obj}"')
rsp = out / "link.rsp"
rsp.write_text(command, encoding="utf-8")
subprocess.run(["link.exe", "@" + str(rsp)], cwd=build, check=True)
subprocess.run([str(exe)], cwd=out, check=True)
