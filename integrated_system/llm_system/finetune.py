from __future__ import annotations

import argparse
import importlib.util
import json
import platform
import random
import time
from pathlib import Path

from .training_data import write_dataset

OPTIONAL = ("torch", "transformers", "datasets", "peft", "accelerate")


def dependency_report() -> dict:
    missing = [name for name in OPTIONAL if importlib.util.find_spec(name) is None]
    return {"ready": not missing, "missing": missing, "python": platform.python_version(),
        "install_command": "python -m pip install 'torch>=2.4,<3' 'transformers>=4.45,<5' 'datasets>=3,<5' 'peft>=0.13,<1' 'accelerate>=1,<2'"}


def train(data_dir: Path, output_dir: Path, model_name: str, max_steps: int, seed: int) -> dict:
    report = dependency_report()
    if not report["ready"]:
        raise RuntimeError("Optional fine-tuning dependencies missing: " + ", ".join(report["missing"]) + ". Install with: " + report["install_command"])
    import numpy as np
    import torch
    from datasets import load_dataset
    from peft import LoraConfig, TaskType, get_peft_model
    from transformers import AutoModelForCausalLM, AutoTokenizer, Trainer, TrainingArguments

    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    dataset = load_dataset("json", data_files={split: str(data_dir / f"{split}.jsonl") for split in ("train", "validation", "test")})
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token

    def encode(row):
        target = json.dumps(row["target"], sort_keys=True)
        text = f"Request: {row['request']}\nWorkflow JSON: {target}"
        prompt = f"Request: {row['request']}\nWorkflow JSON:"
        encoded = tokenizer(text, truncation=True, max_length=256, padding="max_length")
        prompt_length = len(tokenizer(prompt, truncation=True, max_length=256)["input_ids"])
        encoded["labels"] = [
            token if index >= prompt_length and mask else -100
            for index, (token, mask) in enumerate(zip(encoded["input_ids"], encoded["attention_mask"]))
        ]
        return encoded

    tokenized = dataset.map(encode, remove_columns=dataset["train"].column_names)
    base = AutoModelForCausalLM.from_pretrained(model_name)
    model = get_peft_model(base, LoraConfig(task_type=TaskType.CAUSAL_LM, r=8, lora_alpha=16, lora_dropout=0.05, target_modules="all-linear"))
    args = TrainingArguments(output_dir=str(output_dir), max_steps=max_steps, per_device_train_batch_size=2,
        per_device_eval_batch_size=2, gradient_accumulation_steps=4, learning_rate=2e-4, warmup_ratio=0.05,
        eval_strategy="steps", eval_steps=max(1, max_steps // 4), save_strategy="steps", save_steps=max(1, max_steps // 4),
        load_best_model_at_end=True, metric_for_best_model="eval_loss", greater_is_better=False, save_total_limit=2,
        seed=seed, data_seed=seed, report_to="none", use_cpu=not torch.cuda.is_available())
    started = time.time()
    trainer = Trainer(model=model, args=args, train_dataset=tokenized["train"], eval_dataset=tokenized["validation"])
    train_result = trainer.train()
    test_metrics = trainer.evaluate(tokenized["test"], metric_key_prefix="test")
    trainer.save_model(output_dir / "adapter")
    metrics = {**train_result.metrics, **test_metrics, "runtime_total_seconds": time.time() - started,
        "model_name": model_name, "seed": seed, "max_steps": max_steps}
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "metrics.json").write_text(json.dumps(metrics, indent=2, default=float) + "\n", encoding="utf-8")
    (output_dir / "MODEL_CARD.md").write_text(model_card(model_name, metrics), encoding="utf-8")
    return metrics


def model_card(model_name: str, metrics: dict | None = None) -> str:
    return f"""# Synthetic workflow router LoRA adapter

Base model: `{model_name}`

Purpose: educational parameter-efficient fine-tuning for structured request routing. It is not a biological prediction model and must not be used to make scientific conclusions.

Training data: project-generated synthetic requests with group-aware train, validation, and held-out test splits. No project results are treated as model labels.

Metrics: `{json.dumps(metrics or {'status': 'not_run'}, sort_keys=True)}`

Limitations: synthetic templates have limited linguistic and domain coverage. Exact JSON validity, routing accuracy, abstention quality, and expert review remain required. Synthetic examples do not replace expert biological review.
"""


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate data or run a bounded LoRA smoke experiment.")
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).with_name("generated_data"))
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).with_name("training_output"))
    parser.add_argument("--model", default="HuggingFaceTB/SmolLM2-135M-Instruct")
    parser.add_argument("--max-steps", type=int, default=20)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    if args.max_steps < 1 or args.max_steps > 500:
        parser.error("--max-steps must be between 1 and 500")
    manifest = write_dataset(args.data_dir, args.seed)
    if args.prepare_only:
        print(json.dumps({"dataset": manifest, "dependencies": dependency_report()}, indent=2))
        return
    print(json.dumps(train(args.data_dir, args.output_dir, args.model, args.max_steps, args.seed), indent=2, default=float))


if __name__ == "__main__":
    main()
