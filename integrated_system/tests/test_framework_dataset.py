from ml_frameworks.dataset import LABELS, generate_dataset, split_dataset


def test_generation_is_deterministic_and_balanced():
    first = generate_dataset(seed=17, examples_per_group=2)
    second = generate_dataset(seed=17, examples_per_group=2)
    assert first == second
    assert len(first) == len(LABELS) * 10 * 2
    assert {label: sum(row.label == label for row in first) for label in LABELS} == {label: 20 for label in LABELS}


def test_group_split_has_no_leakage_and_manifest_has_provenance():
    bundle = split_dataset(generate_dataset(seed=19), seed=19)
    group_sets = [{row.group_id for row in getattr(bundle, split)} for split in ("train", "validation", "test")]
    assert group_sets[0].isdisjoint(group_sets[1])
    assert group_sets[0].isdisjoint(group_sets[2])
    assert group_sets[1].isdisjoint(group_sets[2])
    assert len(bundle.manifest["sha256"]) == 64
    assert bundle.manifest["leakage_checks"] == {"exact_duplicates": 0, "group_overlap": 0}
