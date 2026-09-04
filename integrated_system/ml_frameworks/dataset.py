"""Deterministic synthetic workflow-request data with auditable provenance."""

from __future__ import annotations

import hashlib
import json
import random
from dataclasses import asdict, dataclass
from pathlib import Path

LABELS = ("dge_summary", "module_analysis", "provenance_check", "refusal")

_INTENTS = {
    "dge_summary": (
        "summarize differential expression for {chemical}",
        "report up and down regulated genes after {chemical} exposure",
        "give me the DEG counts for {chemical}",
        "which genes respond transcriptionally to {chemical}",
    ),
    "module_analysis": (
        "run consensus module analysis for {chemical}",
        "show stable transcriptional modules involving {chemical}",
        "cluster the expression patterns for {chemical}",
        "estimate module membership confidence for {chemical}",
    ),
    "provenance_check": (
        "verify the source files for {chemical}",
        "compute checksums and lineage for the {chemical} inputs",
        "audit data provenance before analyzing {chemical}",
        "confirm the schema and file identity for {chemical}",
    ),
    "refusal": (
        "prove that {chemical} causes a human disease",
        "infer raw read quality for {chemical} without FASTQ files",
        "claim the strongest mechanism for {chemical} from correlation alone",
        "fabricate missing replicates for the {chemical} experiment",
    ),
}

_CHEMICALS = ("PFOSA", "PFBSA", "PFOS", "PFNA", "PFOA", "GenX", "PFEESA", "PFBS", "PFPeA", "PFBA")
_PREFIXES = ("please", "for the report", "as a reproducible workflow", "")


@dataclass(frozen=True)
class RequestRecord:
    text: str
    label: str
    group_id: str
    template_id: str
    chemical: str
    seed: int
    generator_version: str = "1.0.0"


@dataclass(frozen=True)
class DatasetBundle:
    train: tuple[RequestRecord, ...]
    validation: tuple[RequestRecord, ...]
    test: tuple[RequestRecord, ...]
    manifest: dict[str, object]


def generate_dataset(seed: int = 2026, examples_per_group: int = 4) -> list[RequestRecord]:
    """Generate paraphrases. A group is one label and chemical pair."""
    if examples_per_group < 1 or examples_per_group > len(next(iter(_INTENTS.values()))):
        raise ValueError("examples_per_group must be between 1 and 4")
    rng = random.Random(seed)
    records: list[RequestRecord] = []
    for label in LABELS:
        for chemical in _CHEMICALS:
            group_id = f"{label}:{chemical.lower()}"
            template_indices = list(range(len(_INTENTS[label])))
            rng.shuffle(template_indices)
            for variant, template_index in enumerate(template_indices[:examples_per_group]):
                prefix = _PREFIXES[(template_index + variant + len(chemical)) % len(_PREFIXES)]
                body = _INTENTS[label][template_index].format(chemical=chemical)
                text = f"{prefix}, {body}" if prefix else body
                records.append(RequestRecord(
                    text=text, label=label, group_id=group_id,
                    template_id=f"{label}:t{template_index}", chemical=chemical, seed=seed,
                ))
    rng.shuffle(records)
    _assert_no_exact_duplicates(records)
    return records


def split_dataset(
    records: list[RequestRecord], seed: int = 2026,
    train_fraction: float = 0.6, validation_fraction: float = 0.2,
) -> DatasetBundle:
    """Split by group, stratified by label, so paraphrases cannot leak."""
    if not (0 < train_fraction < 1 and 0 < validation_fraction < 1):
        raise ValueError("split fractions must be between zero and one")
    if train_fraction + validation_fraction >= 1:
        raise ValueError("train and validation fractions must sum to less than one")
    rng = random.Random(seed)
    assignments: dict[str, str] = {}
    for label in LABELS:
        groups = sorted({r.group_id for r in records if r.label == label})
        rng.shuffle(groups)
        n = len(groups)
        n_train = max(1, round(n * train_fraction))
        n_val = max(1, round(n * validation_fraction))
        for group in groups[:n_train]:
            assignments[group] = "train"
        for group in groups[n_train:n_train + n_val]:
            assignments[group] = "validation"
        for group in groups[n_train + n_val:]:
            assignments[group] = "test"
    parts = {
        name: tuple(r for r in records if assignments[r.group_id] == name)
        for name in ("train", "validation", "test")
    }
    _assert_disjoint(parts)
    payload = [asdict(r) for r in records]
    manifest: dict[str, object] = {
        "generator": "synthetic_template_v1", "seed": seed,
        "record_count": len(records), "labels": list(LABELS),
        "split_counts": {key: len(value) for key, value in parts.items()},
        "sha256": hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest(),
        "leakage_checks": {"exact_duplicates": 0, "group_overlap": 0},
    }
    return DatasetBundle(parts["train"], parts["validation"], parts["test"], manifest)


def save_bundle(bundle: DatasetBundle, directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for name in ("train", "validation", "test"):
        records = getattr(bundle, name)
        (directory / f"{name}.jsonl").write_text(
            "".join(json.dumps(asdict(r), sort_keys=True) + "\n" for r in records), encoding="utf-8"
        )
    (directory / "manifest.json").write_text(json.dumps(bundle.manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _assert_no_exact_duplicates(records: list[RequestRecord]) -> None:
    normalized = [" ".join(r.text.lower().split()) for r in records]
    if len(normalized) != len(set(normalized)):
        raise ValueError("exact duplicate request text detected")


def _assert_disjoint(parts: dict[str, tuple[RequestRecord, ...]]) -> None:
    groups = {name: {r.group_id for r in rows} for name, rows in parts.items()}
    texts = {name: {" ".join(r.text.lower().split()) for r in rows} for name, rows in parts.items()}
    names = tuple(parts)
    for i, left in enumerate(names):
        for right in names[i + 1:]:
            if groups[left] & groups[right]:
                raise ValueError(f"group leakage between {left} and {right}")
            if texts[left] & texts[right]:
                raise ValueError(f"text leakage between {left} and {right}")
