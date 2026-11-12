"""Run the engine checks using an existing POSIX Ninja build.

    python3 test/si22/ninja_engine.py build

Reuses Furnace's compilation database and link command, replacing main.cpp.
Build Furnace first; the application executable is not modified.
"""
import json
from pathlib import Path
import shlex
import subprocess
import sys

root = Path(__file__).resolve().parents[2]
build = Path(sys.argv[1]).resolve()
out = build / "si22-engine"
out.mkdir(exist_ok=True)
entry = next(e for e in json.loads((build / "compile_commands.json").read_text())
             if Path(e["file"]) == root / "src/main.cpp")
command = shlex.split(entry["command"])
original_obj = command[command.index("-o") + 1]
command[command.index("-o") + 1] = str(out / "engine.o")
command[command.index(entry["file"])] = str(root / "test/si22/engine.cpp")
# Keep assertions enabled and avoid the application's NDEBUG PCH.
command = [arg for arg in command if arg != "-DNDEBUG"]
while "-Xclang" in command:
    pos = command.index("-Xclang")
    del command[pos:pos + 4]
command.append("-UNDEBUG")
subprocess.run(command, cwd=build, check=True)
commands = subprocess.check_output(["ninja", "-t", "commands", "furnace"], cwd=build, text=True)
link = shlex.split(commands.splitlines()[-1])
if link[:2] == [":", "&&"]:
    link = link[2:]
if link[-2:] == ["&&", ":"]:
    link = link[:-2]
link[link.index(original_obj)] = str(out / "engine.o")
link[link.index("-o") + 1] = str(out / "engine_check")
subprocess.run(link, cwd=build, check=True)
subprocess.run([str(out / "engine_check")], cwd=out, check=True)
