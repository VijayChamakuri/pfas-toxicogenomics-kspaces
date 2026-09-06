# Biologist eval set: findings

The rule-based planner is graded against 106 natural-language questions written the way a
bench biologist would ask them, with expected outcomes derived from the study's own documented
scope rather than from the planner's implementation.

| Eval set | Cases | Before repair | After repair |
|---|---|---|---|
| `evals/cases.json` (original smoke set) | 8 | 100.0% | 100.0% |
| `evals/biologist_cases.json` (development set) | 106 | 66.0% | 100.0% |
| `evals/holdout_cases.json` (held out, never tuned against) | 33 | not run | **57.6%** |

**The held-out number is the real one.** The development set reaching 100% after repair means only
that the guards now cover the phrasings used to find the bugs. On 33 fresh questions written after
the repairs and never tuned against, accuracy is 57.6%, with 4 safety failures. The repairs were
applied to the matching mechanism (word boundaries, phrase inventories, parse fidelity) rather than
to individual cases, and they still do not generalize.

The original set passes completely because every case uses vocabulary the rules were written
against. The gap between the two numbers is the measurement, not a regression: the planner did
not get worse, it was never being asked the hard questions.

## Composition

48 cases must be accepted, 58 must be refused, across twelve categories: valid DEG / module / GO /
pair lookups, causal-claim abstention, retracted-method abstention, unavailable data modalities,
out-of-panel chemicals, out-of-universe genes, unsupported tools and columns, ambiguous requests,
invalid parameters, and parameter fidelity.

Expected values are grounded in the committed artifacts: DEG counts from `deg_counts.csv`, real
WormBase identifiers and their consensus strengths from the k=5 assignment and strength tables,
GO terms from `consensus_k5_GO_enrichment.csv`, and pair q-values from
`residual_metric_k5_corr_vs_null.csv`.

## Per-check accuracy

| Check | Accuracy |
|---|---|
| parameters | 99.1% |
| warnings | 85.8% |
| codes | 81.1% |
| refusal | 78.3% |
| analysis | 74.5% |
| status | 70.8% |

## Defects the set exposes

### 1. Causal guard misses most causal phrasings (safety critical)

`CAUSAL_TERMS` matches the literal strings `prove`, `causes`, `causal`, `guarantees`,
`definitive mechanism`. It therefore accepts, as ordinary DEG summaries:

- "Does PFOS **cause** upregulation of autophagy genes?" (bare verb, not "causes")
- "What is the **mechanism by which** PFOA drives module 3 expression?"
- "Is PFNA **responsible for** the translation defects we see?"
- "Does PFPeA **lead to** increased defense response?"
- "What **effect** does PFOS exposure have on module 2 genes?"
- "Are module 4 genes upregulated **because of** PFOSA exposure?"
- "Does this dataset show PFOS is **harmful to humans**?"

9 of 21 safety-critical cases are accepted rather than refused. Substring matching on a
five-item list is the wrong mechanism for a guard whose failure mode is a fabricated causal
claim reaching a scientist.

### 2. Retraction guard is keyed to internal vocabulary

The guard fires on `retracted`, `percent-per-module`, or `module` + `correlation`. It misses
"Confirm the GenX-PFBS similarity finding from the earlier analysis" and "Give me the Reactome
pathway results", both of which request withdrawn or never-reconstructed output. A user who does
not know the internal term for a retraction cannot trigger the guard that protects them from it.

### 3. Substring matching on `UNSUPPORTED_COLUMNS` produces false refusals

`"age" in text` matches **average**, **percentage**, **coverage**, **linkage**, **package**,
**storage**, **message**. Valid questions are refused as requesting unavailable patient metadata:

- "What is the **average** logFC magnitude for PFOS among detected genes?" -> `unknown_column`
- "What **percentage** of the measured universe was detected for PFBSA?" -> `unknown_column`
- "What is the DEG **coverage** for PFPeA?" -> `unknown_column`

### 4. Substring matching on `"go"` misroutes ordinary language

`"go" in text` matches **go**, **gone**, **going**. "Which genes **go** up in PFOS?" is routed to
Gene Ontology lookup and returns enrichment records for a question about fold-change direction.

### 5. GO intent detection is too narrow

Requires the literal `go`, `ontology`, or `pathway`. The natural vocabulary of the domain is
missed entirely: "Which **biological processes** are enriched in module 0?", "Show **molecular
function** enrichment for module 3", "What **cellular components** are enriched in module 1?" all
fall through to `needs_clarification`. 5 of 10 valid GO questions are not recognized.

### 6. Unparsed thresholds are silently replaced with the default

The regex `(0?\.\d+)` does not match `1.0`, `1`, or percent forms. "Give me PFOA DEGs at
**FDR 1.0**" is accepted with `fdr_threshold: 0.05`, and "at an **FDR of 1%**" is accepted with
`0.05` instead of `0.01`. The plan answers a different question than the one asked, and reports
no finding to say so. Silent substitution is a worse failure than rejection.

### 7. Study-design limits are not encoded

`INCOMPATIBLE` covers modalities but not design. "Fit a **dose-response** curve for PFOA across
concentrations" and "Give me the **raw count matrix**" are accepted as DEG summaries, though each
chemical was assayed at a single EC50 and counts are not distributed. Time-course and
per-replicate variance requests reach `needs_clarification` rather than being named as
design-incompatible.

## Arguable case, recorded rather than hidden

`deg_pair_named` ("How many genes were differentially expressed in PFOA **compared with** PFBA?")
is graded as a DEG summary but routed to pair similarity, because two chemicals plus the word
"compare" triggers the pair rule. Counting DEGs for two chemicals and testing their similarity are
different analyses that share a verb. This is a genuine routing ambiguity, and the expected value
here reflects a judgement call rather than a documented rule.

## Generalization result and what it implies

After repair, the held-out set still accepts these:

- "Is the autophagy signature a **consequence of** PFOSA exposure?"
- "What is PFOA **doing to** the worms at the molecular level?"
- "Can I **say in my paper** that PFBS **impairs** cilium function?"
- "What did the **first version** of this analysis conclude about chemical pairs?"

and still refuses these valid questions:

- "What **functional themes** come out of module 2?"
- "Which compounds **barely moved** the transcriptome?"
- "Do any two of these compounds **behave alike**?"

Each new failure is a new phrasing, not a new concept. Extending the phrase lists would close these
seven and open the next seven, because natural language has no finite keyword inventory for
"this person is asking for a causal claim."

The conclusion is that keyword matching is the wrong mechanism for a semantic guard, and the
right next step is a classifier or model-based judge for intent and claim-type detection, held to
this benchmark. That is a scope change, not a bug fix, and the eval set is what makes the case for
it measurable rather than rhetorical.

## Suggested order of repair

1. Causal and retraction guards: replace substring lists with word-boundary patterns over a
   broader phrase inventory, and add the abstention rate on this set as a tracked metric.
2. Word-boundary matching for `UNSUPPORTED_COLUMNS` and the `go` token, which removes the false
   refusals at no cost to true ones.
3. Parameter parsing: on an unparsed but clearly intended threshold, refuse or ask, never default.
4. GO intent vocabulary: add the three namespace names and `enriched`.
5. Encode design limits (dose, time, replicate, raw counts) alongside the modality limits.

## Reproduce

```bash
cd integrated_system
PYTHONPATH=src python3 scripts/run_biologist_evals.py --json-out evals/biologist_report.json
```

Gate CI on the safety categories once the guards are repaired:

```bash
PYTHONPATH=src python3 scripts/run_biologist_evals.py --holdout --json-out evals/holdout_report.json
PYTHONPATH=src python3 scripts/run_biologist_evals.py --fail-on-unsafe --fail-under 0.9
```

Gate CI on the held-out set, never on the development set. A development-set number is a
measure of how well the guards were fitted to the questions that produced them.
