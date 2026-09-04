from __future__ import annotations

import hashlib
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# Only locally reviewed project artifacts belong in the grounded corpus. Generated model
# output, the legacy notebook, and external web text are intentionally excluded.
VERIFIED_SOURCES = {
    "findings": ROOT / "K-spaces Analysis" / "FINDINGS.md",
    "scientific_audit": ROOT / "integrated_system" / "docs" / "SCIENTIFIC_AUDIT.md",
    "reproducibility": ROOT / "integrated_system" / "docs" / "BASELINE_AND_REPRODUCIBILITY.md",
    "dataset_schema": ROOT / "integrated_system" / "docs" / "DATASET_SCHEMA.md",
    "interpretation_rules": ROOT / "integrated_system" / "docs" / "INTERPRETATION_RULES.md",
    "scientific_sources": ROOT / "integrated_system" / "docs" / "SCIENTIFIC_SOURCES.md",
    "methodology": ROOT / "integrated_system" / "docs" / "METHODOLOGY.md",
}


@dataclass(frozen=True)
class Chunk:
    chunk_id: str
    source_id: str
    relative_path: str
    heading: str
    text: str
    source_sha256: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_corpus(max_chars: int = 1800) -> list[Chunk]:
    chunks: list[Chunk] = []
    for source_id, path in VERIFIED_SOURCES.items():
        if not path.is_file():
            raise FileNotFoundError(f"Verified corpus source is missing: {path}")
        source_hash = sha256(path)
        text = path.read_text(encoding="utf-8")
        sections = re.split(r"(?m)(?=^#{1,6}\s+)", text)
        for section_index, section in enumerate(sections):
            section = section.strip()
            if not section:
                continue
            first = section.splitlines()[0]
            heading = re.sub(r"^#{1,6}\s+", "", first).strip() or "Introduction"
            for part_index, start in enumerate(range(0, len(section), max_chars)):
                body = section[start : start + max_chars].strip()
                if body:
                    chunks.append(Chunk(
                        chunk_id=f"{source_id}:{section_index}:{part_index}",
                        source_id=source_id,
                        relative_path=str(path.relative_to(ROOT)),
                        heading=heading,
                        text=body,
                        source_sha256=source_hash,
                    ))
    return chunks


def write_manifest(output: Path) -> Path:
    chunks = build_corpus()
    payload = {
        "schema_version": "1.0",
        "policy": "Local reviewed documents only. Hash changes require corpus revalidation.",
        "chunks": [asdict(chunk) for chunk in chunks],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return output
