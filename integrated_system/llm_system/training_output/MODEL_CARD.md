# Synthetic workflow router LoRA adapter

Base model: `HuggingFaceTB/SmolLM2-135M-Instruct`

Purpose: educational parameter-efficient fine-tuning for structured request routing. It is not a biological prediction model and must not be used to make scientific conclusions.

Training data: project-generated synthetic requests with group-aware train, validation, and held-out test splits. No project results are treated as model labels.

Metrics: `{"epoch": 0.2119205298013245, "max_steps": 8, "model_name": "HuggingFaceTB/SmolLM2-135M-Instruct", "runtime_total_seconds": 50.8690972328186, "seed": 2026, "test_loss": 2.7038605213165283, "test_runtime": 6.1027, "test_samples_per_second": 9.504, "test_steps_per_second": 4.752, "total_flos": 10680306499584.0, "train_loss": 2.8772830963134766, "train_runtime": 44.4824, "train_samples_per_second": 1.439, "train_steps_per_second": 0.18}`

Limitations: synthetic templates have limited linguistic and domain coverage. Exact JSON validity, routing accuracy, abstention quality, and expert review remain required. Synthetic examples do not replace expert biological review.
