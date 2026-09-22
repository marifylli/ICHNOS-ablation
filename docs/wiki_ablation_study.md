# Ablation Study

*Generated from commit `a4af617f` (PULSE) and `caf14266` (FINAL). Every number
below is reproduced by `python run_all.py`; nothing is quoted from memory.*

---

## 1. What this study asks, and what it deliberately does not

An ablation study removes one part of a system at a time and measures what is
lost. In machine learning this is routine [1]; in synthetic biology it is rarer,
because "removing a part" usually means a new strain and three weeks of wet
lab. A validated model lets us do it in seconds — provided we are honest that
the result is a statement about **the model**, and therefore about **our
design reasoning**, not a biological measurement.

We therefore separate three questions that are often blurred together, and we
keep them in separate scripts so that a metric belonging to one cannot leak
into another.

| Arm | Question | Method | Is AICc meaningful? |
|---|---|---|---|
| **A** | What does the biosensor lose if circuit *X* is removed? | Lesion + performance metrics | **No** |
| **B2** | What does the sensor's internal structure change in the readout? | One protocol, two structures | **No** |
| **B1** | Which structure actually describes the data? | Fit to Delaunay 2000 | **Yes** |

The distinction matters because information criteria [2–4] require a likelihood, and
a likelihood requires fitting something to data. Arms A and B2 fit nothing —
they compare the *performance* of circuits under an identical protocol. Quoting
an AIC there would be a category error. Accordingly, `aicc()` is defined inside
`arms/run_B1_aicc.py` and is deliberately not importable from the shared core,
so it cannot be called from the arms where it is meaningless.

**Framing, stated once and meant throughout:** arms A and B2 run on *in
silico* data. They are **methodology validation and power analysis, not
biological proof**. They tell us which parts our design depends on and which
experiments are worth doing. They do not tell us what the cell does.

---

## 2. The system under ablation

ICHNOS is a two-layer device. A **sensing layer** converts a stress dose into
production of TIP, a novel inducer molecule. A **reporting layer** uses TIP to
relieve TetR repression of a tandem fluorescent timer (sfGFP + mCherry). The
division of labour is explicit: **the sensor measures HOW MUCH (dose); the
tandem timer measures WHEN (time since onset)**. Timing information comes from
differential chromophore maturation, not from the sensor's own dynamics.

### 2.1 Reporting layer (shared by every model in this study)

Free TetR represses the reporter promoter with a Hill function, and TIP
sequesters TetR into an inactive complex:

$$
P = P_{\min} + (1-P_{\min})\,\frac{1}{1+\left(\dfrac{[\mathrm{TetR_{active}}]}{K_R}\right)^{n}}
$$

$$
\frac{d[\mathrm{TetR_{active}}]}{dt} = a_{\mathrm{TetR}}\,P - (\mu + k_{\deg,\mathrm{TetR}})[\mathrm{TetR_{active}}] - b\,[\mathrm{TIP}][\mathrm{TetR_{active}}] + u_w[\mathrm{TetR{:}TIP}]
$$

The term $a_{\mathrm{TetR}}P$ is the **autoregulatory negative feedback**:
TetR represses its own promoter. That negative autoregulation linearises a
dose response was shown experimentally in exactly this yeast system by
Nevozhay et al. [5], whose parameter set we adopt. This is the part lesioned as
`no_TetR_feedback`, and it is the mechanism to which we attribute
linearisation.

The timer [6] is a fusion protein with **two parallel, independent maturation
chains** — not a sequential dark→green→red pathway. Green matures in one step
at rate $m$; red matures in two steps at rates $m_1, m_2$, which is what makes
the red/green ratio a clock:

$$
\text{green: } \mathrm{dark} \xrightarrow{m} \mathrm{green}
\qquad
\text{red: } \mathrm{dark} \xrightarrow{m_1} \mathrm{intermediate} \xrightarrow{m_2} \mathrm{red}
$$

Two observables are computed by assignment rules, with $E$ the FRET
bleed-through coefficient:

$$
\mathrm{Observed\_Green} = [\mathrm{green}]\left(1 - E\,\frac{[\mathrm{red}]}{[\mathrm{red}]+[\mathrm{int}]+[\mathrm{dark\,red}]}\right),
\qquad
\mathrm{Ratio} = \frac{[\mathrm{red}]}{\mathrm{Observed\_Green}+\varepsilon}
$$

### 2.2 Sensing layer — this is where FINAL and PULSE differ

**FINAL (static Hill).** TIP production is an instantaneous, memoryless
function of the dose:

$$
\frac{d[\mathrm{TIP}]}{dt} = \beta_{\mathrm{basal}} + (\beta_{\max}-\beta_{\mathrm{basal}})\frac{(S/EC_{50})^{n}}{1+(S/EC_{50})^{n}} - (k_{\deg,\mathrm{TIP}}+\mu)[\mathrm{TIP}]
$$

with $EC_{50,\mathrm{ox}} = 271\ \mu\mathrm{M}$, $n_{\mathrm{ox}} = 1.7$,
$\beta_{\max} = 25$. Under a constant dose this reaches a **sustained
plateau** and stays there.

**PULSE (adaptive sensor).** The static Hill is replaced by an explicit
negative-feedback motif — Yap1 → Trx2 → Yap1 for the oxidative module — with
$A$ the active transcription-factor fraction and $X$ the adaptation mediator:

$$
\frac{dA}{dt} = k_{\mathrm{on}}\,\frac{(S/K_{\mathrm{act}})^{n}}{1+(S/K_{\mathrm{act}})^{n}}\,(1-A) \;-\; k_{\mathrm{off}}\,A\,X
$$

$$
\frac{dX}{dt} = k_x A - d_x X,
\qquad
\frac{d[\mathrm{TIP}]}{dt} = \beta_{\mathrm{basal}} + (\beta_{\max}-\beta_{\mathrm{basal}})A - (k_{\deg,\mathrm{TIP}}+\mu)[\mathrm{TIP}]
$$

The factor $(1-A)$ encodes a fixed total Yap1 pool (Delaunay Fig. 2B lower
panel shows total Yap1 is constant [7]; Kuge & Jones show activation is
post-translational [8]). The term $-k_{\mathrm{off}}AX$ is the thioredoxin-mediated
reduction documented *in vitro* by Delaunay Fig. 7D [7]. This structure was chosen
from **biological argument**, not from curve fitting: Ma et al. 2009 showed by
exhaustive three-node topology search that only two families achieve
adaptation [9], and negative feedback with a buffer node is one of them — ours.

The parameter $d_x$ sets the **adaptation depth**: $d_x = 0$ gives a perfect
adaptor (the signal returns fully to baseline); $d_x > 0$ gives partial
adaptation. This single parameter is the reason the two modules behave
differently.

### 2.3 Why the ox and er modules are not the same experiment

| | **ox (oxidative)** | **er (ER stress)** |
|---|---|---|
| Pathway | Yap1 → Trx2 → Yap1 [7, 8] | Ire1 → Hac1 → Kar2/BiP → Ire1 [12, 13] |
| Promoter | pTRX2 | pKAR2 |
| Stressor | H₂O₂ | DTT / tunicamycin |
| $K_{\mathrm{act}}$ (PULSE) | 208 µM | 2345 (arb.) |
| $n$ (PULSE) | 2.75 | 3.97 |
| $d_x$ | **0.5 — partial adaptor** | **0.0 — perfect adaptor** |
| $\beta_{\max}$ FINAL → PULSE | 25 → 118 | 25 → 25 (**unchanged**) |
| Calibration status | Calibrated against Delaunay 2000 | Structure transferred, parameters provisional |

Two consequences follow directly, and both show up in the results:

1. **The ox sensor adapts only partially, so under a constant dose
   `Observed_Green` rises to a plateau — a step. The er sensor is a perfect
   adaptor, so its readout genuinely pulses.** The same reporter therefore
   behaves qualitatively differently on the two modules. This is a property of
   $d_x$, not an artefact.
2. **`perfect_adaptation` is not an applicable lesion for er**, because
   $d_{x,\mathrm{er}}$ is already zero. The code now detects this and skips the
   row. Previously it printed an exact `+0.0%` line that read like a null
   result when it was in fact a non-applicable one.

A caveat we state rather than bury: **the er parameters are not calibrated to
the same standard as ox.** The ox module has a documented source hierarchy
(§7). The er module inherited the structure and has provisional values, so er
numbers are reported as *descriptive*, and no design decision rests on them
alone.

---

## 3. Protocol

Defined once in `protocol.py`; figure captions are generated from it, so code
and captions cannot drift apart.

**Pre-equilibration: 50 h at $S = 0$ before the dose is applied.** All species
start at zero in the merged SBML. That is not a cell, it is an empty vessel.
Without pre-equilibration, the first hours of the ratio are dominated by the
*filling* of the reporter pools rather than by the stress response — and the
two have the same timescale, so they are not separable. Measured: without
pre-equilibration the dose information is **underestimated roughly fivefold**
(spread 1.8 % instead of 8.0 %). One would conclude that the ratio is nearly
dose-invariant, which is the wrong conclusion and leads to a wrong decoder
design. `checks/test_protocol_invariance.py` verifies convergence for
20/50/200 h and asserts that the cold start does underestimate.

**Readout $t^{*} = 3$ h**, matching Dacquay & McMillen's accumulated-reporter
measurement [10], so that model output and calibration target are read at the same
time. See §6.4 — this choice is now tested rather than assumed.

**Dose grids.** Two grids, for two different purposes:

- `DOSES` — the **operating grid** (ox: 50–600 µM; er: 500–6000). Linearity is
  judged only here. The upper bound is biological: Dacquay Fig. 3A shows signal
  collapse above ~600 µM from toxicity [10].
- `FIT_DOSES` — the **characterisation grid** (14 points, geometric, 5–5000 µM
  for ox). Fitting a four-parameter curve on six points that sample neither
  plateau is not identifiable. The extra doses exist only to pin the baseline
  and the top asymptote; they are **not claims about the cell**.

**Reproducibility.** `models/build_snapshots.py` clones both upstream repos at
**pinned commits**, builds the merged SBML, and records commit hash, SHA-256 of
every snapshot, and the solver versions in `MANIFEST.json`. Every runner calls
`verify()` at start-up and halts if any file has changed. Result filenames
carry the upstream commit hash.

---

## 4. Metrics, and one metric we had to throw away

### 4.1 What went wrong with the first version

The original study reported an "effective Hill coefficient" $n_{\mathrm{eff}}$,
computed as the slope of $\mathrm{logit}$ of the min–max-normalised response
against $\log_{10}$ dose on the six-point operating grid. Tested against
synthetic curves with known answers, this estimator returns:

| True curve | Returned $n_{\mathrm{eff}}$ |
|---|---|
| Hill, $n = 1$ | 2.03 |
| **Purely linear, $y = d$** | **2.03** |
| Hill, $n = 4$, $EC_{50} = 208$ | 4.29 |
| Hill, $n = 4$, $EC_{50} = 1000$ | 8.85 |

It returns the same number for a true Hill $n=1$ and for a straight line, and
it changes by a factor of two for the *same* cooperativity at a different
$EC_{50}$. The headline claim "$n_{\mathrm{eff}} \approx 2.0$, therefore
linearised" was reporting the estimator's fixed point for any gently curving
response, not a property of the circuit.

We found this while preparing the study for review, and we report it because
the correction is more informative than the original claim: **the conclusion
survived, but the evidence for it did not.**

### 4.2 Fitted Hill coefficient — replacement

`hill_fit()` performs a four-parameter least-squares fit of the Hill equation
[14, 15] on the wide grid, using SciPy's bounded trust-region-reflective
solver [16, 17]:

$$
y(d) = y_{\mathrm{base}} + (y_{\mathrm{top}} - y_{\mathrm{base}})\,\frac{d^{\,n_H}}{d^{\,n_H} + EC_{50}^{\,n_H}}
$$

and reports $n_H$ **with its standard error** from the covariance matrix of the
fit (the asymptotic Gauss–Newton approximation $\sigma^2(J^{\mathsf T}J)^{-1}$), plus $EC_{50}$ and $R^2$. On a grid that samples neither plateau the
parameters are not jointly identifiable and the standard errors blow up — which
is the point of reporting them.

### 4.3 Linearity $R^2$ — the direct test of the design claim

$n_H$ answers "how cooperative", which is related to but not identical with
"how usable". `linearity_r2()` fits a plain straight line inside the operating
window and reports $R^2$. No normalisation and no transform, so it cannot
manufacture a number on a dead signal. For a design claim phrased as
"the feedback linearises the response", this is the most direct evidence we
can offer.

### 4.4 Viability gate

`ratio_spread` — the dose-induced spread of the red/green ratio — is
**minimised by a dead circuit**, and shape metrics remain well defined on top
of negligible signal (measured: at $K_d = 5000$ nM, spread = 0.0004 while fold
change = 1.003). Optimising a shape metric would therefore select a
non-functional device. The gate is binary and comes first, and it distinguishes
two failure modes that must not be conflated:

| Gate | `fold_range` | `fold_vs_zero` | Meaning |
|---|---|---|---|
| `viable` | ≥ θ | ≥ θ | Functional |
| `saturated` | < θ | ≥ θ | Responds strongly but **cannot discriminate doses** |
| `dead` | < θ | < θ | Does not respond |

`saturated` is a **result**, not a rejection: it is exactly what removing the
linearising feedback is predicted to cause.

The gate now **classifies** rather than **scores**. Fitted $n_H$ and $R^2$ are
reported for every lesion regardless of gate — a saturated circuit has a
perfectly meaningful Hill coefficient, and that is the headline finding for
`no_TetR_feedback`. Only `ratio_spread`, the metric a dead circuit can win,
stays gated.

**θ is provisional and must be tied to the image-pipeline noise floor.** We
report the whole sweep θ ∈ {1.5, 2.0, 3.0} rather than a single value, and note
that θ = 3.0 rejects even the intact ox circuit.

### 4.5 Decoding bound — the metric the device is actually judged on

Every metric above is a proxy. The deliverable is a decoder returning
(dose, time-since-onset) per cell, so the question a reviewer will ask is *how
much worse does the estimate get without part X* — in hours, not in units of
ratio spread.

`decoding_crlb()` computes the **Cramér–Rao lower bound** [18–20] on a joint
estimate of $\theta = (t, \ln d)$ from both observables. Fisher-information
analysis of this kind is established for identifiability in biochemical
kinetics [21] and connects directly to the information-theoretic view of
noisy signalling [22]. With multiplicative lognormal
measurement noise of coefficient of variation $c$, the log-observables carry
additive noise of constant standard deviation, so with
$\mathbf{y} = (\ln \mathrm{Green}, \ln \mathrm{Red})$:

$$
J = \frac{\partial \mathbf{y}}{\partial (t,\ \ln d)},
\qquad
\Sigma = \mathrm{diag}(c_g^2,\ c_r^2),
\qquad
\mathcal{F} = J^{\mathsf T}\Sigma^{-1}J
$$

$$
\mathrm{Cov}(\hat\theta) \succeq \mathcal{F}^{-1}
\quad\Longrightarrow\quad
\sigma_t \ge \sqrt{[\mathcal{F}^{-1}]_{11}},
\qquad
\sigma_{\ln d} \ge \sqrt{[\mathcal{F}^{-1}]_{22}}
$$

Three properties make this the right primary metric:

- **It is a bound, not an estimate.** No decoder, however clever, beats it. So
  "this lesion costs 6 hours of timing resolution" is a hard statement about
  the device, not about our software.
- **It cannot be gamed by a dead circuit.** As the response flattens,
  $J \to 0$, $\mathcal{F}$ becomes singular, and the bound diverges. The time derivative
  in $J$ is taken on the simulated trajectory; the dose derivative is a central
  finite difference with step 0.05 in $\ln d$. The
  pathology that forced the viability gate does not arise here.
- **Because the parameter is $\ln d$, the dose error comes out relative** —
  the natural unit for a dose-response device.

Two assumptions we state rather than hide. First, $\Sigma$ is diagonal: the two
channels are treated as independently noisy. Shared extrinsic noise (cell size,
expression capacity) is correlated and would make the true bound slightly
worse. Second, the CRLB is a **local** bound, valid when the estimator is
approximately unbiased; when $\sigma_{\ln d}$ comes out above roughly 0.3 the
correct reading is **"dose is not identifiable at this operating point"**, not
"the error is N %".

---

## 5. Lesions

| Lesion | Implementation | What it tests |
|---|---|---|
| `full` | — | Reference |
| `no_TetR_feedback` | Freeze $P$ in the TetR production reaction only | The linearisation mechanism |
| `no_adaptation` | $k_{\mathrm{off}} = 0$ | The sensor's own negative feedback |
| `perfect_adaptation` | $d_x = 0$ | Adaptation depth (ox only; non-applicable for er) |
| `no_timer_split` | $m_1 = m_2 = 2m$ | **Positive control** — red matched to green |
| `no_leaky_expression` | $\beta_{\mathrm{basal}} = 0$ | Basal TIP production |
| `no_FRET_correction` | $E = 0$ | *Observation model*, not a circuit part |

Three design decisions worth defending explicitly:

**The TetR feedback lesion is surgical.** $P$ is frozen *only* in the TetR
production reaction, leaving the reporter's $P$ untouched. Otherwise we would
be ablating the signal path, not the autoregulation. The reaction is located
**topologically** — product is TetR_active *and* the rate law contains $P$ —
not by name, because names change on merge while topology does not. If two such
reactions are found the code raises rather than guesses.

**The frozen value is $P$ at zero stress**, so both circuits start from the
same baseline TetR level. Otherwise part of the measured difference would
simply be a different amount of repressor rather than the feedback itself.

**`no_FRET_correction` is separated into its own category.** $E$ is a parameter
of *how we observe*, not of *what the device does*. Listing it beside
`no_leaky_expression` invites a reader to treat an optical correction as a
biological part.

---

## 6. Results

### 6.1 Arm A — lesion ladder, oxidative module

```
              lesion        gate    fold      vs0        n_H    EC50   linR2   spr3h
                full      viable   2.687    4.737  1.97±0.07     112   0.907   0.080
    no_TetR_feedback   saturated   1.221   12.155  3.01±0.06      32   0.414     -
  perfect_adaptation      viable   2.591    4.157  1.98±0.07     117   0.918   0.068
       no_adaptation   saturated   1.249    8.536  3.19±0.05      34   0.420     -
      no_timer_split      viable   2.666    4.669  1.97±0.07     111   0.907   0.001
```

**The linearisation claim holds, and now has an error bar.** Removing the TetR
autoregulation raises $n_H$ from 1.97 ± 0.07 to 3.01 ± 0.06 — a separation of
roughly fifteen combined standard errors — and drops linearity $R^2$ from 0.91
to 0.41. The gate independently classifies the lesioned circuit as
`saturated`: it responds twelvefold versus zero dose but only 1.2-fold across
the dose range. That is precisely the predicted failure — the device still
knows that *something happened* but no longer knows *how much*.

**`no_adaptation` behaves almost identically to `no_TetR_feedback`** (3.19
versus 3.01, both saturated). With $k_{\mathrm{off}} = 0$ the sensor integrates
without restraint, $A \to 1$, TIP saturates, and the downstream consequence is
the same. Two different mechanisms, one shared failure mode.

**`perfect_adaptation` changes almost nothing** (1.98 versus 1.97). Adaptation
*depth* is not what the dose-response steepness depends on; the *presence* of
adaptation is. This is a useful negative result: it tells the wet lab that
$d_x$ is a tolerant parameter.

**Threshold sensitivity is a warning, not a footnote.** At θ = 3.0 *no* lesion
passes, including the intact circuit. Until the image pipeline delivers a
noise-floor estimate, θ = 2.0 is a choice and we present it as one.

### 6.2 Arm A — lesion ladder, ER module

```
              lesion        gate    fold      vs0        n_H    EC50   linR2   spr3h
                full      viable   3.511    3.993  2.79±0.08    1497   0.936   0.058
    no_TetR_feedback      viable  20.290   31.135  5.56±0.09     988   0.602   0.245
       no_adaptation      viable   6.072    7.054  4.57±0.08    1149   0.737   0.282
      no_timer_split      viable   3.487    3.943  2.81±0.08    1498   0.935   0.002
```

**The same qualitative result, with a quantitatively different signature.**
$n_H$ falls from 5.56 to 2.79 and $R^2$ rises from 0.60 to 0.94 — a larger
absolute reduction than in ox, consistent with the steeper sensor
($n_{\mathrm{er}} = 3.97$ versus $n_{\mathrm{ox}} = 2.75$). The feedback has
more non-linearity to remove, and removes it.

**Nothing is gated here, and that itself is the difference between the
modules.** Because the er sensor is a perfect adaptor, the lesioned circuit
retains a 20-fold dose range instead of collapsing to 1.2-fold. The failure
mode of the ox module (saturation) does not reproduce in er. **Conclusions
about gate behaviour are therefore module-specific and must not be stated as
properties of "the ICHNOS circuit".**

### 6.3 Arm A3 — sensitivity, and Arm A4 — design specification

Percent change in the primary output per ablated part, at a reference dose near
$K_{\mathrm{act}}$. This is a **one-at-a-time, knockout-style** analysis, not a
global sensitivity analysis [23]: it tells us what each part contributes, not
how parameter uncertainties interact. It is at $t^{*}$ rather than at true steady state —
because the decoder operates at 3 h and the circuit is still transient there.
Steady-state numbers would describe a regime the device is never read in.

For ox: `no_TetR_feedback` +207.5 %, `no_adaptation` +116.0 %,
`perfect_adaptation` −13.0 %, `no_leaky_expression` −11.7 %, `no_timer_split`
−6.4 %, and (observation model) `no_FRET_correction` +12.0 %. The ordering is
the actionable output: **the two feedback loops dominate; basal expression and
the FRET correction are tolerant at the ~10 % level.** For the wet lab that
means construct identity of the feedback elements must be verified, while
promoter leakiness need not be tuned precisely.

Arm A4 sweeps $K_d^{\mathrm{TIP\text{-}TetR}}$, which has **no literature
baseline** — TIP is a novel molecule, so this parameter is exploratory by
necessity. The sweep does two jobs:

```
  Kd (nM)      gate     fold      vs0   spread@3h      n_H  sigma_t(h)
     0.05    viable    2.711    4.847     0.09698    1.974       2.611
     0.25    viable    2.687    4.737     0.07975    1.966       2.808
        1    viable    2.553    4.214     0.05730    1.962       3.352
        5    viable    2.093    2.847     0.03895    1.979       4.513
       25      dead    1.454    1.618     0.02546    2.012       7.368
      100      dead    1.143    1.179     0.01294    2.025      16.838
     5000      dead    1.003    1.004     0.00039    2.029     635.446
```

First, it is the cleanest possible demonstration of §4.4: **as the circuit
dies, `ratio_spread` falls monotonically toward zero while $\sigma_t$ diverges
by three orders of magnitude.** The two curves move in opposite directions.
Anyone optimising the shape metric would drive the design straight into a
non-functional device; the decoding bound is immune.

Second, it yields an actual **specification for the wet lab**: the device is
functional for $K_d \lesssim 5$ nM at θ = 2.0, and timing resolution is already
degraded by 60 % at 5 nM relative to sub-nanomolar binding. The recommendation
we pass on is $K_d \le 1$ nM, which is a tighter and more useful statement than
the gate alone gives.

### 6.4 Arm A5 — decoding bound

```
                lesion  sigma_t (h)  sigma_dose    corr  vs intact
                  full        3.040       0.915  -0.859      1.00x
      no_TetR_feedback        9.652     186.133  -0.997      3.17x
    perfect_adaptation        2.886       0.530  -0.502      0.95x
         no_adaptation      145.833    2421.717  -1.000     47.97x
        no_timer_split      270.293      40.223  -1.000     88.90x
```

**The positive control works, and that validates the metric.** Matching red
maturation to green destroys timing resolution by a factor of 89 (96 for er)
while leaving every dose-response metric essentially untouched ($n_H$ 1.97
versus 1.97). This is exactly what a control should do: hit one axis and leave
the other alone. Under the old `ratio_spread` metric the same lesion produced a
"+63 %" that was a units artefact. **The tandem-timer design is now supported by
a quantitative statement rather than by construction.**

**The two channels are jointly, not separately, informative.** The correlation
between $\hat t$ and $\widehat{\ln d}$ is −0.86 in the intact circuit and
approaches −1.000 in every ablated one. A correlation at −1 means the two
observables have become effectively one: dose and time trade off along a
degenerate direction and cannot be separated. This is the quantitative form of
a claim the project previously made qualitatively — *both observables are
structurally necessary for the decoder*.

**The decoder reads cultures, not single cells, and that is what makes the
bound usable.** The observable is a mean over $N$ cells, whose noise has two
parts that behave completely differently:

$$
\sigma^2_{\text{measured}}(N) = \frac{c^2_{\text{cell}}}{N} + c^2_{\text{shared}}
$$

$c_{\text{cell}}$ (expression noise, per-cell photon statistics, segmentation
error) averages away as $1/\sqrt{N}$. $c_{\text{shared}}$ — illumination drift,
exposure, flat-field, background subtraction, the crosstalk coefficient — is
identical for every cell in the field, so averaging does nothing to it. The
consequence is a **floor**:

| Shared error | N=1 | N=10 | N=100 | N=1000 | Floor ($N\to\infty$) | N for $\sigma_t \le 0.5$ h |
|---|---|---|---|---|---|---|
| 0 % | 3.04 h | 0.96 h | 0.30 h | 0.10 h | 0 | 37 |
| 2 % | 3.05 h | 0.99 h | 0.39 h | 0.26 h | 0.24 h | 49 |
| 5 % | 3.10 h | 1.14 h | 0.68 h | 0.62 h | **0.61 h** | **never** |
| 10 % | 3.28 h | 1.55 h | 1.25 h | 1.22 h | **1.22 h** | **never** |

Two actionable numbers come out of this. The crossover
$N^{*} = (c_{\text{cell}}/c_{\text{shared}})^2$ — 156 cells at 2 %, 25 at 5 %,
6 at 10 % — is the point beyond which imaging more cells buys almost nothing.
And at 5 % shared error the target is **unreachable at any N**. The
requirement this places on the image pipeline is therefore not "segment more
cells" but **"keep systematic error below ~2 %"**, which is a statement about
calibration, not about throughput.

**Timing resolution depends on the age of the event, and that is not a choice
of readout time.** An earlier version of this study reported that $t^{*}=3$ h
"is not the optimum". That framing was wrong: it treats the readout time as
something we select, when in deployment the time since onset is precisely the
**unknown being estimated**. The correct object is $\sigma_t(\text{age})$ —
how precisely an event can be dated, given how long ago it happened
(N = 100 cells):

| Age of event | Timing error | Dose error | Timing at 5 % shared |
|---|---|---|---|
| 0.5 h | 9 min | **132 %** | 20 min |
| 1 h | 4 min | 20 % | 8 min |
| 2 h | 7 min | 10 % | 15 min |
| 3 h | 18 min | 9 % | 41 min |
| 4 h | 49 min | 10 % | 110 min |
| 6 h | 5.0 h | 14 % | 11 h |

**The two axes are complementary, which is the quantitative form of the design
claim.** Timing degrades with age, because the ratio flattens as the
population equilibrates. Dose *improves* with age, because signal accumulates:
below ~0.75 h the dose is not identifiable at all. The sensor answering HOW
MUCH and the timer answering WHEN are not just conceptually separate — they
have measurably different windows of competence.

**This narrows a claim the project has been making.** The memory window is
described elsewhere as 0–4 h. The defensible version, at N = 100 with shared
error under 2 %, is: **events 0.5–3 h old are dated to better than 30 minutes,
and dose is recovered to ~10 % for events older than 1 h**; beyond 4 h the
device reports "some hours ago" rather than a time. That is a smaller claim
than 0–4 h, and it is one we can support.

$t^{*} = 3$ h is retained throughout this study as the **calibration** readout,
matching Dacquay's measurement, so that every comparison is internally
consistent. It is not a claim about when the device is read. One consequence
for the wet lab: a validation experiment must span the age grid above rather
than sit at a single time point, and 3 h is the part of that range where
timing is *weakest*.

### 6.5 Arm B2 — sensor structure, FINAL versus PULSE

This arm is controlled by construction: `reporter_module_v2.sbml` and
`TIP_TetR_binding.sbml` are **byte-identical across the two upstream repos**,
verified by SHA-256 in `MANIFEST.json` and asserted by
`checks/test_downstream_identical.py`. Only the sensing module differs. If that
ever stops being true the check fails and the arm does not run.

```
=== ox
               arm    fold     vs0     n_H    spr3h    spr6h   spr12h
             FINAL   2.228   2.466   1.694   0.1306   0.0259   0.0021
             PULSE   2.687   4.737   1.966   0.0798   0.0091   0.0007
     FINAL_beta118   4.689   7.441   1.703   0.1387   0.0216   0.0016
```

**The confound must be disclosed, not omitted.** $\beta_{\max,\mathrm{ox}}$ is
25 in FINAL and 118 in PULSE, following the correction of the
$\beta_{\max}/\beta_{\mathrm{basal}}$ ratio from 5 to 23.6 against Dacquay
Fig. 3A. So "structural difference" is also a gain difference, and the third
row re-runs FINAL at PULSE's gain to separate them.

**The honest reading is uncomfortable and we report it as such: on
dose-response performance, the adaptive sensor does not win.** At matched gain,
FINAL delivers a *larger* dose range (4.69 versus 2.69) and a *lower* $n_H$
(1.70 versus 1.97). The static sensor is the better dose meter.

There is a clean explanation. FINAL's sensing layer already has
$n_{\mathrm{ox}} = 1.7$, and its fitted output $n_H = 1.694$ — the promoter
output reproduces the sensor's own cooperativity almost exactly. PULSE's sensor
has $n_{\mathrm{ox}} = 2.75$, and the circuit brings it down to 1.97. **The
feedback is doing more work in PULSE precisely because it has more
non-linearity to remove.** The non-linearity already decreases along the chain
before the TetR loop is involved at all — a point that must be stated, or the
circuit will be credited with linearisation that pre-existed it.

**So why keep PULSE?** Because the two models are not competing on the same
criterion. A static Hill is **structurally incapable** of producing the
transient activation that Delaunay et al. measured directly: under constant
$S$ it gives a sustained plateau, full stop. PULSE reproduces the transient and
is calibrated against primary data ($K_{\mathrm{act}} = 208$ µM, $n = 2.75$
from densitometry of Delaunay Figs. 2B/2C). **The adaptive sensor is chosen for
biological fidelity, and that choice costs measured dose-response performance.
Naming the price is the point of running this arm.**

```
=== er
  ! gain-confound control SKIPPED: override beta_max_er=25 is a NO-OP
    (current value is already 25).

               arm    fold     vs0     n_H    spr3h    spr6h   spr12h
             FINAL   2.380   2.614   1.994   0.1393   0.0273   0.0022
             PULSE   3.511   3.993   2.787   0.0582   0.0413   0.0242
```

**The er comparison is descriptive only, and the gain control does not exist
for it.** Two independent reasons:

1. $\beta_{\max,\mathrm{er}}$ is 25 in *both* repos, so the confound-control
   arm would be an exact copy of the FINAL row. An earlier version of this
   study printed it anyway, as three rows of which two were identical,
   presented as a controlled comparison. A control that controls nothing is
   worse than no control, because it looks like evidence.
   `assert_override_bites()` now raises on any no-op override and the arm is
   skipped with the reason printed.
2. **The er sensing module is not a single-variable change.** Structure,
   $n_{\mathrm{er}}$ (2.0 → 3.97) and $\beta_{\mathrm{basal}}$ (5 → 0.77) all
   change at once between repos. Even with a gain control, attribution to
   "structure" would be unsupported.

One er observation is worth flagging because it distinguishes the modules
sharply: the ratio spread in PULSE er *stays high at 12 h* (0.024, versus
0.0007 for ox). Perfect adaptation keeps the sensor pulsing, so the ratio
retains dose information far longer. **The usable memory window is a property
of $d_x$, and therefore differs between the two stress modules.** The decoder
cannot assume one window for both.

### 6.6 Arm B1 — model selection on data

**Not yet run.** It requires `data/delaunay2000_fig2b.csv`, a digitisation of
Delaunay et al. 2000 Fig. 2B with a mandatory provenance header recording
figure, panel, axis units, dose, who digitised it, when, and with which tool.
The script halts cleanly rather than producing a placeholder result.

That requirement is not bureaucracy. Four separate source-mixing incidents were
found during model construction — parameters traced to secondary descriptions
rather than to primary measurements — including one where a
$\beta_{\max}/\beta_{\mathrm{basal}}$ ratio was off by a factor of 4.7. The
provenance header exists so that a fifth cannot happen silently.

Rules the implementation will follow, recorded now so the analysis cannot drift:

- **AICc, not AIC** [3, 4], since $n/k < 40$ with few timepoints:
  $\mathrm{AICc} = n\ln(\mathrm{RSS}/n) + 2k + \dfrac{2k(k+1)}{n-k-1}$
- **Akaike weights**, not a declared winner: $\Delta < 2$ means
  indistinguishable, $\Delta > 10$ strong preference.
- **Same data and same error model for both structures**, or AICc compares
  error models rather than biology.
- $k$ counts only the free parameters of *this* fit, plus one for $\sigma$.
- **Most important:** a static Hill *cannot* structurally produce a transient.
  If it loses, that is **structural inadequacy, not a statistical victory** —
  the outcome is close to predetermined, and presenting a large $\Delta$AICc as
  evidence would be overclaiming. We will report it as such.

A stronger use of the same data, which we recommend over AICc alone, is
**profile likelihood** [24] for $k_{\mathrm{on}}$, $k_{\mathrm{off}}$ and $d_x$. We
already know analytically that $k_{\mathrm{on}}$ and $k_{\mathrm{off}}$ are
non-identifiable at steady state; *demonstrating* that with flat profiles is a
methodological result in its own right, and more honest than a comparison whose
winner is known in advance.

---

## 7. Parameter provenance

Every parameter is traced to a **primary measurement**, not to a secondary
description. The source hierarchy for the oxidative module:

| Source | What it is used for | What it must **not** be used for |
|---|---|---|
| Delaunay et al. 2000 [7] | Direct Yap1 oxidation state; calibration of $K_{\mathrm{act}}$, $n$; justification of $(1-A)$ and $-k_{\mathrm{off}}AX$ | — |
| Dacquay & McMillen 2021 [10] | Accumulated reporter at 3 h; $EC_{50}$ target; $\beta$ ratio 23.6×; toxicity bound ~600 µM | Timing information |
| Gasch et al. 2000 [11] | TRX2 mRNA kinetics; verification that dose is constant (0.32 ± 0.03 mM, assayed every 3 min) | Parameter fitting — prediction checking only |
| Kuge & Jones 1994 [8] | TRX2 as direct Yap1 target; post-translational activation | **Time-course claims** — all measurements are 1 h endpoints |
| Nevozhay et al. 2009 [5] | Full TetR package ($a$, $d$, $\theta$, $n$, $P_{\min}$) as a **coherent set** | Mixing individual values with other sources |
| Khmelinskii et al. 2012 [6] | sfGFP and mCherry maturation constants | — |
| Ma et al. 2009 [9] | Topological justification of the adaptation motif | Parameter values |

Two coherence rules learned the hard way. **Parameters from one system are
adopted as a package**: mixing $K_R$ from one source with a production rate
from another makes TetR saturate the Hill term at every sweep point.
$k_{\deg,\mathrm{TetR}}$ from Nevozhay is already a combined
degradation-plus-dilution rate, so adding a separate $\mu$ double-counts.

---

## 8. Reproducing this study

Simulation uses libRoadRunner [25] via Tellurium [26] on SBML models [27]
manipulated with libSBML [28].

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python run_all.py
```

Order is enforced: snapshots → checks → arms. If a check fails the arms do not
run, because their results would not mean what this page says they mean. Arm B1
halts cleanly until its data file exists. Runtime is a few minutes; all CSV and
PNG outputs carry the upstream commit hash in the filename.

---

## 9. Summary of what we learned

**Supported.** The decoder works as a population instrument: at N = 100 cells
it dates events 0.5–3 h old to better than 30 minutes and recovers dose to
~10 % beyond 1 h. The TetR autoregulation measurably linearises the dose response
in both modules ($n_H$ 3.01 → 1.97 in ox, 5.56 → 2.79 in er; linearity $R^2$
0.41 → 0.91 and 0.60 → 0.94). The tandem-timer split is necessary for timing:
matching red to green costs a factor of 89 in $\sigma_t$. Both observables are
required — ablating any circuit element drives the $t$/dose correlation to −1,
making them degenerate. $K_d^{\mathrm{TIP\text{-}TetR}} \le 1$ nM is a
design requirement.

**Reported against our own interest.** At matched gain the static sensor is the
better dose meter; the adaptive sensor is chosen for biological fidelity and
that choice has a measured cost. Part of the apparent linearisation pre-exists
the feedback circuit. The memory window must be narrowed from 0–4 h to
"0.5–3 h to better than 30 min". Single-cell decoding is not feasible at
CV = 25 %; the device is a population instrument. If systematic imaging error
reaches 5 %, no number of cells recovers the target resolution. The first
version of this study used a metric that could not distinguish a Hill $n = 1$
from a straight line, and an earlier version of this section misread the
readout-time sweep as a protocol recommendation when it is a characterisation
of the device.

**Module-specific, not general.** The ox module is a partial adaptor and its
lesioned circuit saturates; the er module is a perfect adaptor, does not
saturate, and retains dose information roughly thirty times longer at 12 h. er
parameters are provisional and er results are descriptive.

**Open.** Arm B1 awaits the digitised Delaunay time course. θ awaits the
image-pipeline noise floor. Both are named blockers with named owners rather
than silent gaps.

---

## 10. References

**Methods**

- **[1]** Meyes, R., Lu, M., de Puiseau, C. W., & Meisen, T. (2019). Ablation studies in artificial neural networks. *arXiv*:1901.08644.
- **[2]** Akaike, H. (1974). A new look at the statistical model identification. *IEEE Transactions on Automatic Control*, 19(6), 716–723.
- **[3]** Hurvich, C. M., & Tsai, C.-L. (1989). Regression and time series model selection in small samples. *Biometrika*, 76(2), 297–307.
- **[4]** Burnham, K. P., & Anderson, D. R. (2002). *Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach* (2nd ed.). Springer.
- **[14]** Hill, A. V. (1910). The possible effects of the aggregation of the molecules of haemoglobin on its dissociation curves. *Journal of Physiology*, 40(Suppl.), iv–vii.
- **[15]** Goutelle, S., Maurin, M., Rougier, F., Barbaut, X., Bourguignon, L., Ducher, M., & Maire, P. (2008). The Hill equation: a review of its capabilities in pharmacological modelling. *Fundamental & Clinical Pharmacology*, 22(6), 633–648.
- **[16]** Branch, M. A., Coleman, T. F., & Li, Y. (1999). A subspace, interior, and conjugate gradient method for large-scale bound-constrained minimization problems. *SIAM Journal on Scientific Computing*, 21(1), 1–23.
- **[17]** Virtanen, P., et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261–272.
- **[18]** Rao, C. R. (1945). Information and the accuracy attainable in the estimation of statistical parameters. *Bulletin of the Calcutta Mathematical Society*, 37, 81–91.
- **[19]** Cramér, H. (1946). *Mathematical Methods of Statistics*. Princeton University Press.
- **[20]** Kay, S. M. (1993). *Fundamentals of Statistical Signal Processing, Volume I: Estimation Theory*. Prentice Hall.
- **[21]** Komorowski, M., Costa, M. J., Rand, D. A., & Stumpf, M. P. H. (2011). Sensitivity, robustness, and identifiability in stochastic chemical kinetics models. *PNAS*, 108(21), 8645–8650.
- **[22]** Cheong, R., Rhee, A., Wang, C. J., Nemenman, I., & Levchenko, A. (2011). Information transduction capacity of noisy biochemical signaling networks. *Science*, 334(6054), 354–358.
- **[23]** Saltelli, A., et al. (2008). *Global Sensitivity Analysis: The Primer*. Wiley.
- **[24]** Raue, A., Kreutz, C., Maiwald, T., Bachmann, J., Schilling, M., Klingmüller, U., & Timmer, J. (2009). Structural and practical identifiability analysis of partially observed dynamical models by exploiting the profile likelihood. *Bioinformatics*, 25(15), 1923–1929.

**Software**

- **[25]** Somogyi, E. T., et al. (2015). libRoadRunner: a high performance SBML simulation and analysis library. *Bioinformatics*, 31(20), 3315–3321.
- **[26]** Choi, K., et al. (2018). Tellurium: an extensible Python-based modeling environment for systems and synthetic biology. *Biosystems*, 171, 74–79.
- **[27]** Hucka, M., et al. (2003). The systems biology markup language (SBML): a medium for representation and exchange of biochemical network models. *Bioinformatics*, 19(4), 524–531.
- **[28]** Bornstein, B. J., Keating, S. M., Jouraku, A., & Hucka, M. (2008). LibSBML: an API library for SBML. *Bioinformatics*, 24(6), 880–881.

**Biology and parameter sources**

- **[5]** Nevozhay, D., Adams, R. M., Murphy, K. F., Josić, K., & Balázsi, G. (2009). Negative autoregulation linearizes the dose–response and suppresses the heterogeneity of gene expression. *PNAS*, 106(13), 5123–5128.
- **[6]** Khmelinskii, A., et al. (2012). Tandem fluorescent protein timers for in vivo analysis of protein dynamics. *Nature Biotechnology*, 30(7), 708–714.
- **[7]** Delaunay, A., Isnard, A.-D., & Toledano, M. B. (2000). H₂O₂ sensing through oxidation of the Yap1 transcription factor. *The EMBO Journal*, 19(19), 5157–5166.
- **[8]** Kuge, S., & Jones, N. (1994). YAP1 dependent activation of TRX2 is essential for the response of *Saccharomyces cerevisiae* to oxidative stress by hydroperoxides. *The EMBO Journal*, 13(3), 655–664.
- **[9]** Ma, W., Trusina, A., El-Samad, H., Lim, W. A., & Tang, C. (2009). Defining network topologies that can achieve biochemical adaptation. *Cell*, 138(4), 760–773.
- **[10]** Dacquay, L. C., & McMillen, D. R. (2021). Improving the design of an oxidative stress sensing biosensor in yeast. *FEMS Yeast Research*.
- **[11]** Gasch, A. P., Spellman, P. T., Kao, C. M., Carmel-Harel, O., Eisen, M. B., Storz, G., Botstein, D., & Brown, P. O. (2000). Genomic expression programs in the response of yeast cells to environmental changes. *Molecular Biology of the Cell*, 11(12), 4241–4257.
- **[12]** Cox, J. S., Shamu, C. E., & Walter, P. (1993). Transcriptional induction of genes encoding endoplasmic reticulum resident proteins requires a transmembrane protein kinase. *Cell*, 73(6), 1197–1206.
- **[13]** Pincus, D., Chevalier, M. W., Aragón, T., van Anken, E., Vidal, S. E., El-Samad, H., & Walter, P. (2010). BiP binding to the ER-stress sensor Ire1 tunes the homeostatic behavior of the unfolded protein response. *PLoS Biology*, 8(7), e1000415.
