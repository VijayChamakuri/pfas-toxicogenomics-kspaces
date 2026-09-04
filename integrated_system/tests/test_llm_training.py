import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1]))

from llm_system.finetune import dependency_report, model_card
from llm_system.training_data import generate, group_split, validate, write_dataset


def test_llm_synthetic_data_schema_quality_and_provenance():
    records = generate()
    result = validate(records)
    assert result["valid"], result["errors"]
    assert len(records) >= 250
    assert all(row["provenance"]["kind"] == "synthetic_template" for row in records)
    assert all(row["target"]["requires_expert_review"] for row in records)


def test_llm_group_split_has_no_group_or_exact_leakage():
    splits = group_split(generate())
    group_sets = [{row["group_id"] for row in rows} for rows in splits.values()]
    fingerprints = [{row["fingerprint"] for row in rows} for rows in splits.values()]
    assert all(group_sets[i].isdisjoint(group_sets[j]) for i in range(3) for j in range(i + 1, 3))
    assert all(fingerprints[i].isdisjoint(fingerprints[j]) for i in range(3) for j in range(i + 1, 3))
    assert all(splits.values())


def test_llm_dataset_writer_and_manifest(tmp_path):
    manifest = write_dataset(tmp_path)
    assert manifest["quality"]["valid"]
    assert sum(manifest["split_counts"].values()) == manifest["quality"]["record_count"]
    for split in ("train", "validation", "test"):
        rows = [json.loads(line) for line in (tmp_path / f"{split}.jsonl").read_text().splitlines()]
        assert len(rows) == manifest["split_counts"][split]


def test_llm_dependency_report_and_model_card_are_truthful():
    report = dependency_report()
    assert report["ready"] == (not report["missing"])
    assert "pip install" in report["install_command"]
    card = model_card("test/model")
    assert "not_run" in card
    assert "not a biological prediction model" in card
