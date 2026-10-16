"""Verify all delivered payloads before running saved-output checks."""
from pathlib import Path
import csv
import hashlib

root = Path(__file__).resolve().parent.parent
with (root / "MANIFEST.csv").open(newline="") as stream:
    rows = list(csv.DictReader(stream))
seen = set()
for row in rows:
    rel = Path(row["file"])
    if rel.is_absolute() or ".." in rel.parts or str(rel) in seen:
        raise SystemExit("Invalid manifest path: " + str(rel))
    seen.add(str(rel))
    path = root / rel
    data = path.read_bytes()
    if len(data) != int(row["bytes"]) or hashlib.sha256(data).hexdigest() != row["sha256"]:
        raise SystemExit("Integrity failure: " + str(rel))
print("PASS:", len(rows), "SHA-256 payload entries")
