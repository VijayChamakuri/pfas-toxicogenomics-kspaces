"""Framework-neutral features, metrics, and benchmark utilities."""

from __future__ import annotations

import json
import re
import time
import tracemalloc
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Self

import numpy as np

from .dataset import LABELS, RequestRecord

TOKEN_RE = re.compile(r"[a-z0-9]+")


@dataclass
class BenchmarkTimer:
    started: float = 0.0

    def __enter__(self) -> Self:
        tracemalloc.start()
        self.started = time.perf_counter()
        return self

    def __exit__(self, *_: object) -> None:
        self.seconds = time.perf_counter() - self.started
        _, self.peak_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()


def build_vocabulary(records: tuple[RequestRecord, ...], min_count: int = 1) -> dict[str, int]:
    counts = Counter(token for r in records for token in TOKEN_RE.findall(r.text.lower()))
    return {token: i for i, (token, count) in enumerate(sorted(counts.items())) if count >= min_count}


def vectorize(records: tuple[RequestRecord, ...], vocabulary: dict[str, int]) -> tuple[np.ndarray, np.ndarray]:
    features = np.zeros((len(records), len(vocabulary)), dtype=np.float32)
    labels = np.empty(len(records), dtype=np.int64)
    for row, record in enumerate(records):
        for token in TOKEN_RE.findall(record.text.lower()):
            if token in vocabulary:
                features[row, vocabulary[token]] += 1.0
        norm = np.linalg.norm(features[row])
        if norm:
            features[row] /= norm
        labels[row] = LABELS.index(record.label)
    return features, labels


def classification_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, object]:
    confusion = np.zeros((len(LABELS), len(LABELS)), dtype=int)
    for actual, predicted in zip(y_true, y_pred, strict=True):
        confusion[int(actual), int(predicted)] += 1
    per_class = {}
    f1_values = []
    for i, label in enumerate(LABELS):
        tp = confusion[i, i]
        fp = confusion[:, i].sum() - tp
        fn = confusion[i, :].sum() - tp
        precision = float(tp / (tp + fp)) if tp + fp else 0.0
        recall = float(tp / (tp + fn)) if tp + fn else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        f1_values.append(f1)
        per_class[label] = {"precision": precision, "recall": recall, "f1": f1, "support": int(confusion[i].sum())}
    return {
        "accuracy": float((y_true == y_pred).mean()), "macro_f1": float(np.mean(f1_values)),
        "confusion_matrix": confusion.tolist(), "per_class": per_class,
    }


def error_analysis(records: tuple[RequestRecord, ...], truth: np.ndarray, predicted: np.ndarray) -> list[dict[str, str]]:
    return [
        {"text": record.text, "expected": LABELS[int(actual)], "predicted": LABELS[int(pred)]}
        for record, actual, pred in zip(records, truth, predicted, strict=True) if actual != pred
    ]


def write_report(report: dict[str, object], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
