# data/

Digitised primary data for arm B1, with the code that produced it.

| File | What it is | Re-executable? |
|---|---|---|
| `delaunay_fig2b_timecourse.csv` | Fig. 2B — Yap1 oxidation time course, 400 µM | **Yes** |
| `gel_fig2b_nonreducing.png` | the crop the 2B numbers came from | — |
| `delaunay_fig2b_densitometry.py` | the densitometry that produced them | **Yes** |
| `delaunay_fig2c_oxfraction.csv` | Fig. 2C — dose response at 5 min | No, values frozen |
| `fig2c_crop_as_archived.png` | the 2C crop as it survives on disk | shape only  |
| `delaunay_fig2c_densitometry.py` | the 2C densitometry | shape only |
| `fit_Kact_from_fig2c.m` | how `K_act_ox` and `n_ox` were fitted to 2C | — |

**Verified 2026-09-25:** re-running `delaunay_fig2b_densitometry.py` on
`gel_fig2b_nonreducing.png` reproduces `ox_fraction_window` and
`lane_total_intensity` exactly, to three decimals. The `ox_fraction_unmix`
column does **not** reproduce, because SciPy changed its `nnls`
implementation — which is why the window method is the primary column and the
unmix column is only a sensitivity check (`python arms/run_B1_aicc.py unmix`).

The 2C source crop was archived at a smaller scale than the one actually used,
so those values re-derive to within ~0.03 but not exactly. They are therefore
frozen. This costs nothing here: 2C calibrated `K_act_ox` and `n_ox`, which
arm B1 holds **fixed**, and B1 tests on 2B only — see below.

## Mandatory provenance header

Every CSV starts with `#` comments recording source, figure and panel, units,
dose, who digitised it, when, with what tool, and whether the input image is
archived. `arms/run_B1_aicc.py` refuses to run without `source:`, `figure:`,
`digitised_by:` and `tool:`.

**Why:** four separate source-mixing incidents were found in the POC —
parameters traced to secondary descriptions rather than to primary
measurements, including a beta ratio off by 4.7×. The header exists so a
fifth cannot happen silently.

## Training data vs test data

`K_act_ox` and `n_ox` were fitted to **Fig. 2C** (`fit_Kact_from_fig2c.m`).
Refitting the sensor on 2C would be circular, so 2C is training data.
**Fig. 2B was never used for calibration**, which makes it a genuine held-out
test — and it is the time axis, the only axis on which a static and an
adaptive sensor differ at all.

Note the loop recorded in `fit_Kact_from_fig2c.m`: `k_off = 150` and
`d_x = 0.6` were themselves calibrated under `K_act = 271`, so refitting
`K_act` invalidates them. B1 avoids this entirely by refitting the kinetics on
2B and holding only `K_act`/`n` fixed.

## Source hierarchy for the oxidative module

| Source | Use | Must NOT be used for |
|---|---|---|
| Delaunay 2000 | direct Yap1 oxidation state — calibration of `A_ox` | — |
| Dacquay & McMillen 2021 | accumulated reporter at 3 h — EC50 target | timing |
| Gasch 2000 | TRX2 mRNA kinetics | parameter fitting — prediction checking only |
| Kuge & Jones 1994 | single endpoint (1 h, 1 mM) — consistency check | time-course claims |