"""Copy six checksum-verified inputs from an existing local study directory."""
from pathlib import Path
import csv
import hashlib
import os
import shutil
import sys
import tempfile

root = Path(__file__).resolve().parent.parent
if len(sys.argv) != 2:
    raise SystemExit("Supply the existing study directory.")
source = Path(sys.argv[1]).expanduser().resolve(strict=True)
dest = Path.cwd().resolve()
if source == dest or not (dest / "real_data/R/arfc_real_data_v1.R").is_file():
    raise SystemExit("Run from a newly prepared study directory.")
with (root / "inputs.csv").open(newline="") as stream:
    rows = list(csv.DictReader(stream))
if len(rows) != 6 or len({row["file"] for row in rows}) != 6:
    raise SystemExit("Expected six unique analysis inputs.")

copies = []
for row in rows:
    rel = Path(row["file"])
    if rel.is_absolute() or ".." in rel.parts:
        raise SystemExit("Invalid input path.")
    src = (source / rel).resolve(strict=True)
    target = dest / rel
    if not src.is_relative_to(source) or not target.resolve().is_relative_to(dest):
        raise SystemExit("Input path leaves its study directory.")
    payload = src.read_bytes()
    if len(payload) != int(row["bytes"]) or hashlib.md5(payload).hexdigest() != row["md5"]:
        raise SystemExit("Input checksum mismatch: " + str(rel))
    if target.exists():
        raise SystemExit("Refusing to overwrite: " + str(rel))
    copies.append((src, target, row))

for src, target, row in copies:
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=target.name + ".", dir=target.parent)
    try:
        with os.fdopen(fd, "wb") as stream, src.open("rb") as original:
            shutil.copyfileobj(original, stream)
        payload = Path(temporary).read_bytes()
        if len(payload) != int(row["bytes"]) or hashlib.md5(payload).hexdigest() != row["md5"]:
            raise SystemExit("Copy verification failed: " + row["file"])
        # A hard link installs the completed file without replacing a concurrent file.
        os.link(temporary, target)
    finally:
        Path(temporary).unlink(missing_ok=True)
print("PASS: six existing inputs copied with matching checksums; source files unchanged.")
