# ICHNOS — Ablation Study

## Repository layout

```
ichnos-ablation/
├── run_all.py                  one command: snapshots -> checks -> arms
├── protocol.py                 single source of truth for the protocol
├── ichnos_ablation.py          shared core: simulation, metrics, gate, lesions
├── requirements.txt
├── models/
│   ├── build_snapshots.py      clones upstream repos, merges, records hashes
│   ├── MANIFEST.json           commit hashes + SHA-256 of every snapshot
│   └── merged_{ox,er}_{FINAL,PULSE}.sbml
├── checks/
│   ├── test_downstream_identical.py   is B2 still a controlled comparison?
│   └── test_protocol_invariance.py    do conclusions survive protocol choices?
├── arms/
│   ├── run_A_lesions.py        Step 2  - gated lesion ladder
│   ├── run_A3_sensitivity.py   Step 3  - % change per ablated part
│   ├── run_A4_kd_spec.py       extension - functional Kd window
│   ├── run_A5_decoding.py      Cramer-Rao bound on (dose, t) + t* sweep
│   ├── run_B2_structure.py     static vs adaptive sensor, same protocol
│   └── run_B1_aicc.py          model selection on data (awaiting data)
├── data/README.md              provenance spec for B1 input
├── docs/wiki_ablation_study.md full wiki text
└── results/                    CSV + PNG, named with upstream commit hash
```

Reproducible ablation study for the ICHNOS biosensor.
This repo does **not** contain models: it generates them as snapshots from the
upstream repos and records commit hash + SHA-256, so every number is traceable
to a specific version.

## Quick start

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python run_all.py            # builds snapshots, runs checks, runs every arm
```

Total runtime is a few minutes. Individual steps:

```bash
python models/build_snapshots.py        # once, and after any upstream change
python checks/test_downstream_identical.py
python checks/test_protocol_invariance.py
python arms/run_A_lesions.py PULSE ox        # A step 2: lesion ladder
python arms/run_A_lesions.py PULSE er
python arms/run_A3_sensitivity.py PULSE ox   # A step 3: sensitivity comparison
python arms/run_A3_sensitivity.py PULSE er
python arms/run_A4_kd_spec.py PULSE ox       # A step 4: design specification
python arms/run_B2_structure.py ox
python arms/run_B2_structure.py er
```

## Mapping to the iGEM Engineering Success framework

| Framework step | Script | Output |
|---|---|---|
| 1. Baseline | any arm (the `full` row) | primary predictive curve |
| 2. Parameter ablation | `arms/run_A_lesions.py` | gated lesion ladder |
| 3. Sensitivity comparison | `arms/run_A3_sensitivity.py` | % change per ablated part |
| (extension) design spec | `arms/run_A4_kd_spec.py` | functional Kd window |
| (extension) decoding bound | `arms/run_A5_decoding.py` | sigma_t, sigma_dose per lesion |

## The three arms — do NOT merge them into one table

| | Question | Method | Does AICc apply? |
|---|---|---|---|
| **A** | What does the sensor lose without circuit X? | lesion + performance metrics | **No** |
| **B2** | What does sensor structure change in the readout? | same protocol, two structures | **No** |
| **B1** | Which structure describes the data? | fit to Delaunay 2000 | **Yes** |

A and B2 do not fit anything — they compare performance, so there is no
likelihood and AIC/AICc is undefined. `aicc()` lives exclusively inside
`arms/run_B1_aicc.py` so that it cannot be called from elsewhere.

## Metrics

**Cooperativity is measured by fitting, not by a slope trick.** `hill_fit()`
fits `base + (top-base)*d^n/(d^n+EC50^n)` on the wide grid `P.FIT_DOSES` and
reports `n_H` with its standard error. The previous `n_eff` (logit slope on a
min-max normalised 6-point grid) returned **2.03 for a true Hill n=1 and 2.03
for a purely linear response** -- it was the estimator's fixed point, not a
property of the circuit. `linearity_r2()` complements it: plain straight-line
R^2 inside `OPERATING_WINDOW` only.

**The decoding bound is the primary figure of merit.** `decoding_crlb()`
returns the Cramer-Rao lower bound on a joint estimate of (time since onset,
log dose) from both observables under multiplicative noise of size
`NOISE_CV`. No decoder can beat it, so "this lesion costs X hours of timing
resolution" is a hard statement. It also cannot be gamed by a dead circuit:
the Fisher matrix goes singular and the bound diverges.

## The viability gate

`ratio_spread` is **minimised by a dead circuit**, and `n_eff` stays well
defined on top of negligible signal (measured: at `Kd = 5000 nM`,
spread = 0.0004 while fold = 1.003). Shape metrics without a gate reward the
dead circuit. The gate is therefore **binary and comes first**, and shape
metrics are recorded only for what passes it.

Two folds, because they diagnose different things:

| gate | `fold_range` | `fold_vs_zero` | meaning |
|---|---|---|---|
| `viable` | ≥ θ | ≥ θ | functional |
| `saturated` | < θ | ≥ θ | responds strongly, **cannot discriminate doses** |
| `dead` | < θ | < θ | does not respond |

`saturated` is a **result**, not a rejection.

The gate now **classifies** (dead vs saturated) rather than **scores**.
Fitted `n_H` and `lin_r2` are reported for every lesion regardless of gate --
a saturated circuit has a perfectly meaningful Hill coefficient, and that is
exactly the headline finding for `no_TetR_feedback`. Only the shape metrics
a dead circuit can win (`ratio_spread`) stay gated.

**θ must be tied to the image-pipeline noise floor.** On the current grid,
θ = 3.0 rejects even the intact circuit (ox) — always run the sweep.

## Protocol

Defined **only** in `protocol.py`; captions are generated from it.
Critical: pre-equilibration ≥ 20 h at `S = 0`. Without it the dose information
is underestimated ~5×, because filling the reporter pools from zero initial
conditions has the same timescale as the stress response itself.

Known pitfall: `Observed_Green` and `Measured_Ratio_RG` are assignment-rule
**parameters** — roadrunner does not return them in the default selections
(you would read back the initial value 1.0). The core recomputes them from the
species.

## Framing

Arms A and B2 run on in silico data: they are **methodology validation and
power analysis, not biological proof**. B1 is not a substitute for wet-lab
validation.
