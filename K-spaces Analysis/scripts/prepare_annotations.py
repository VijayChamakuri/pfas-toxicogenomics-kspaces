"""Validate and decompress the vendored WormBase GAF snapshot."""

from __future__ import annotations

import gzip
import hashlib
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data_external"
COMPRESSED = DATA / "wb.gaf.gz"
DECOMPRESSED = DATA / "wb.gaf"
EXPECTED_SHA256 = "bda92f6642cc68a16e485ba32eb5890320855895185b9cbcd448f52cf698fed2"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    observed = sha256(COMPRESSED)
    if observed != EXPECTED_SHA256:
        raise RuntimeError(f"WormBase GAF checksum mismatch: {observed}")
    with gzip.open(COMPRESSED, "rb") as source, DECOMPRESSED.open("wb") as target:
        shutil.copyfileobj(source, target)
    print(f"Prepared {DECOMPRESSED.relative_to(ROOT)} ({DECOMPRESSED.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()

