# Grounded workflow and fine-tuning extension

This directory implements the project workflow-planning and validation layer. It does not alter the PFAS scientific results.

## Grounded retrieval

`corpus.py` chunks four reviewed local documents and records each source path and SHA-256 hash. `rag.py` uses transparent lexical retrieval because the corpus is small and auditability matters more than embedding complexity. Accepted proposals contain structured steps, inputs, assumptions, warnings, and source citations. Missing evidence, unavailable data, causal claims, unknown chemicals, fake columns, and the retracted module-correlation method trigger clarification, incompatibility, or abstention.

Run the deterministic comparison:

```bash
cd pfas-toxicogenomics-kspaces/integrated_system
PYTHONPATH=. python -m llm_system.evaluate --output llm_system/evaluation_report.json
```

The evaluator compares a fixed prompt-only rules baseline with retrieval grounding. It measures exact status and intent routing, citation presence, workflow structure, and absence of a small set of forbidden unsupported-claim phrases. These are deterministic engineering metrics, not a substitute for expert answer review. Use `HUMAN_REVIEW_RUBRIC.md` for biological and methodological quality.

## Fine-tuning

`training_data.py` creates 300 candidate synthetic workflow requests, removes exact duplicates, and splits by template family. Group-aware separation prevents paraphrases from one template family appearing across train, validation, and test. Provenance, schema version, and content fingerprints are stored per row.

Prepare and validate data without model dependencies:

```bash
PYTHONPATH=. python -m llm_system.finetune --prepare-only
```

Install the optional stack and run the bounded default LoRA experiment:

```bash
python -m pip install 'torch>=2.4,<3' 'transformers>=4.45,<5' 'datasets>=3,<5' 'peft>=0.13,<1' 'accelerate>=1,<2'
PYTHONPATH=. python -m llm_system.finetune --max-steps 20
```

The default base is `HuggingFaceTB/SmolLM2-135M-Instruct`. Model download requires network access and is subject to the model repository terms. The output contains adapter checkpoints, Trainer logs, held-out loss, and a model card. Loss alone does not establish workflow correctness. A production experiment must also parse generated JSON and run the same held-out routing, abstention, adversarial, and human-review evaluations.

No fine-tuned model is used to produce scientific conclusions. RAG is used only to constrain workflow proposals to reviewed project context. Fine-tuning is used only to learn the synthetic request-to-routing schema.
