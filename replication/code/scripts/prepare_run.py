"""Prepare a fresh, isolated study directory; never launch a numerical run."""
from pathlib import Path
import re
import shutil
import sys

root = Path(__file__).resolve().parent.parent
name = sys.argv[1] if len(sys.argv) == 2 else ""
if not re.fullmatch(r"[A-Za-z0-9_-]+", name):
    raise SystemExit("Supply an alphanumeric run name, optionally containing - or _.")
dest = root / "work" / name
if dest.exists():
    raise SystemExit("Refusing to overwrite " + str(dest))
study = root / "study"
for src in study.rglob("*"):
    if not src.is_file():
        continue
    rel = src.relative_to(study)
    s = rel.as_posix()
    frozen = "/production" in s or s.startswith("manuscript/generated/")
    if frozen or s.startswith("output/") or s.endswith((".pdf", ".eps", ".png")):
        continue
    target = dest / rel
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, target)
print("Prepared fresh source and input tree:", dest.relative_to(root))
print("Next: bash run.sh acquire " + name + "; see README.md for analysis inputs.")
