from __future__ import annotations

import hashlib
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

SCHEMA_VERSION = "1.0"
LICENSE = "Project-generated synthetic data; no third-party text or biological labels."


def _record(group: str, request: str, status: str, intent: str = "unknown") -> dict:
    target = {"status": status, "intent": intent, "requires_expert_review": True}
    fingerprint = hashlib.sha256(request.strip().lower().encode()).hexdigest()
    return {
        "schema_version": SCHEMA_VERSION,
        "group_id": group,
        "request": request,
        "target": target,
        "provenance": {"kind": "synthetic_template", "generator": "training_data.py", "license": LICENSE},
        "fingerprint": fingerprint,
    }


def generate(seed: int = 2026) -> list[dict]:
    rng = random.Random(seed)
    chemicals = ["PFOSA", "PFBSA", "PFOS", "PFNA", "PFOA", "GenX", "PFEESA", "PFBS", "PFPeA", "PFBA"]
    templates = {
        "deg": ["Summarize differential expression for {chemical}.", "Report DEGs and effect direction for {chemical}.", "What is the transcriptomic response to {chemical}?"],
        "pair": ["Compare corrected response patterns for {chemical} and {other}.", "What expression patterns are shared by {chemical} and {other}?"],
        "go": ["Review GO enrichment for consensus module {module}.", "Test ontology enrichment for module {module} using the measured-gene background."],
        "gene": ["Find consensus membership for WBGene{gene:08d}.", "Report module and confidence for WBGene{gene:08d}."],
        "raw": ["Run FastQ alignment for {chemical}.", "Perform raw-read batch correction for {chemical}."],
        "causal": ["Prove {chemical} causes toxicity through module {module}.", "Confirm the definitive mechanism of {chemical}."],
        "column": ["Stratify {chemical} response by patient_id.", "Compare {chemical} by sex and age."],
        "ambiguous": ["Analyze the biology.", "Find the best result.", "Do a complete analysis."],
    }
    qualifiers = [
        " Include assumptions.",
        " Return a structured workflow.",
        " State the required inputs.",
        " Flag important limitations.",
        " Include an expert-review warning.",
    ]
    audiences = [
        " This is for a methods review.",
        " This is for a reproducibility check.",
        " This is for a workflow planning exercise.",
        " This is for an analyst handoff.",
    ]
    records = []
    for index in range(800):
        kind = list(templates)[index % len(templates)]
        template = templates[kind][(index // len(templates)) % len(templates[kind])]
        chemical = chemicals[index % len(chemicals)]
        other = chemicals[(index * 3 + 1) % len(chemicals)]
        audience_index = (index // 7) % len(audiences)
        request = template.format(
            chemical=chemical, other=other, module=index % 5 + 1, gene=10000 + index
        ) + qualifiers[(index // len(templates)) % len(qualifiers)] + audiences[audience_index]
        group = f"template-{kind}-{templates[kind].index(template)}-audience-{audience_index}"
        if kind == "deg":
            records.append(_record(group, request, "accepted", "differential_expression_summary"))
        elif kind == "pair":
            records.append(_record(group, request, "accepted", "cross_chemical_comparison"))
        elif kind == "go":
            records.append(_record(group, request, "accepted", "go_enrichment_review"))
        elif kind == "gene":
            records.append(_record(group, request, "accepted", "consensus_module_lookup"))
        elif kind in {"raw", "column"}:
            records.append(_record(group, request, "incompatible"))
        elif kind == "causal":
            records.append(_record(group, request, "abstain"))
        else:
            records.append(_record(group, request, "needs_clarification"))
    # Deduplicate before splitting. Repeated template families remain grouped to prevent paraphrase leakage.
    unique = {record["fingerprint"]: record for record in records}
    result = list(unique.values())
    rng.shuffle(result)
    return result


def validate(records: list[dict]) -> dict:
    required = {"schema_version", "group_id", "request", "target", "provenance", "fingerprint"}
    errors = []
    fingerprints: Counter[str] = Counter()
    for index, row in enumerate(records):
        missing = required - set(row)
        if missing:
            errors.append(f"row {index}: missing {sorted(missing)}")
        expected = hashlib.sha256(row.get("request", "").strip().lower().encode()).hexdigest()
        if row.get("fingerprint") != expected:
            errors.append(f"row {index}: invalid fingerprint")
        fingerprints[str(row.get("fingerprint", ""))] += 1
    duplicates = sorted(key for key, count in fingerprints.items() if count > 1)
    if duplicates:
        errors.append(f"duplicate fingerprints: {len(duplicates)}")
    return {"valid": not errors, "errors": errors, "record_count": len(records), "groups": len({r.get('group_id') for r in records})}


def group_split(records: list[dict], seed: int = 2026) -> dict[str, list[dict]]:
    groups: dict[str, list[dict]] = defaultdict(list)
    for row in records:
        groups[row["group_id"]].append(row)
    strata: dict[tuple[str, str], list[str]] = defaultdict(list)
    for key, rows in groups.items():
        target = rows[0]["target"]
        strata[(target["status"], target["intent"])].append(key)
    assignments: dict[str, list[str]] = {"train": [], "validation": [], "test": []}
    rng = random.Random(seed)
    for keys in strata.values():
        keys.sort()
        rng.shuffle(keys)
        n = len(keys)
        validation_n = max(1, round(n * 0.15))
        test_n = max(1, round(n * 0.15))
        if validation_n + test_n >= n:
            raise ValueError(f"not enough groups for three-way stratification: {n}")
        assignments["validation"].extend(keys[:validation_n])
        assignments["test"].extend(keys[validation_n:validation_n + test_n])
        assignments["train"].extend(keys[validation_n + test_n:])
    splits = {name: [row for key in selected for row in groups[key]] for name, selected in assignments.items()}
    split_groups = [{row["group_id"] for row in split} for split in splits.values()]
    if any(split_groups[i] & split_groups[j] for i in range(3) for j in range(i + 1, 3)):
        raise AssertionError("group leakage across splits")
    fingerprints = [{row["fingerprint"] for row in split} for split in splits.values()]
    if any(fingerprints[i] & fingerprints[j] for i in range(3) for j in range(i + 1, 3)):
        raise AssertionError("exact duplicate leakage across splits")
    return splits


def write_dataset(output_dir: Path, seed: int = 2026) -> dict:
    records = generate(seed)
    quality = validate(records)
    if not quality["valid"]:
        raise ValueError(quality["errors"])
    splits = group_split(records, seed)
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, rows in splits.items():
        path = output_dir / f"{name}.jsonl"
        path.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in rows), encoding="utf-8")
    manifest = {"schema_version": SCHEMA_VERSION, "seed": seed, "quality": quality,
        "split_counts": {key: len(value) for key, value in splits.items()},
        "limitations": ["Synthetic templates teach request routing, not biological truth.", "Template-group splitting reduces but cannot eliminate distribution similarity."]}
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest
