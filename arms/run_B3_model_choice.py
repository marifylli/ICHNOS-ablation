"""
ARM B3 -- which sensor model should the DECODER use: FINAL or PULSE?

WHY THIS ARM EXISTS
  B2 compares the two structures on dose-response metrics and finds, honestly,
  that at matched gain the STATIC sensor is the better dose meter. That is the
  wrong criterion for choosing the decoder's forward model: the deliverable
  estimates (dose, time since onset), so the two candidates have to be
  compared on the bound for that joint estimate -- in minutes and percent, at
  the ages the device is actually used at.

  Without this arm, "we chose PULSE because it is closer to the biology" is a
  preference. With it, the choice is a measurement, the price of the choice is
  stated, and a reviewer can check both.

WHAT IS HELD FIXED
  The reporting layer is byte-identical between the two repos (asserted by
  checks/test_downstream_identical.py), so the ONLY difference is the sensing
  module. The gain confound is handled the same way as in B2: a third arm
  re-runs FINAL at PULSE's beta_max, because beta_max_ox is 25 in FINAL and
  118 in PULSE and otherwise "structure" and "gain" are confounded.

WHAT THIS ARM DOES NOT DO
  It does not say which model is TRUE -- that is B1, on data. It says which
  one the decoder should invert, given that both are approximations.

  python arms/run_B3_model_choice.py            # ox
  python arms/run_B3_model_choice.py er
"""
import sys, csv
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import protocol as P, ichnos_ablation as A
from models import build_snapshots as S

RESULTS = Path(__file__).resolve().parents[1] / "results"
REF_DOSE = {"ox": 200, "er": 2000}


def main(variant="ox"):
    mf = S.verify()
    sb = {tag: S.path(variant, tag).read_text() for tag in ("FINAL", "PULSE")}
    dose = REF_DOSE[variant]
    N, cv = P.N_REPORT, P.NOISE_CV["green"]
    bmax = A.pname("beta_max", variant)

    arms = [("FINAL", sb["FINAL"], None), ("PULSE", sb["PULSE"], None)]
    b_pulse = float(A.te.loadSBMLModel(sb["PULSE"])[A.ids(sb["PULSE"])[bmax]])
    try:
        A.assert_override_bites(sb["FINAL"], bmax, b_pulse, f"B3/{variant} gain control")
        arms.append((f"FINAL@beta{b_pulse:g}", sb["FINAL"], {bmax: b_pulse}))
    except RuntimeError as e:
        print(f"  ! gain control SKIPPED: {e}")

    print(f"\n=== B3: which model should the decoder invert? -- {variant}")
    print(f"  dose {dose}, N = {N} cells, per-cell CV = {cv:.0%}, no shared error.")
    print(f"  sigma_t in MINUTES, sigma_dose relative.\n")
    hdr = f"{'age (h)':>8}"
    for name, _, _ in arms:
        hdr += f"{name:>26}"
    print(hdr)

    rows = []
    curves = {name: [] for name, _, _ in arms}
    for age in P.AGE_GRID:
        line = f"{age:8.2f}"
        for name, model, ov in arms:
            J = A.decoding_jacobian(model, variant, dose, t_star=float(age),
                                    overrides=ov)
            v = A.population_variance(cv, 0.0, N)
            c = A.crlb_from_jacobian(J, v, v)
            curves[name].append((60 * c["sigma_t"], c["sigma_lnD"]))
            line += f"{60 * c['sigma_t']:16.1f} min {c['sigma_lnD']:7.0%}"
            rows.append(dict(arm=name, age_h=age, sigma_t_min=60 * c["sigma_t"],
                             sigma_dose_rel=c["sigma_lnD"]))
        print(line)

    # ---- the comparison, stated in the units of the decision ---------------
    ages = np.array(P.AGE_GRID, float)
    win = (ages >= 0.5) & (ages <= 2.0)     # where timing is the job
    print("\n  TIMING, over the 0.5-2 h window where the timer carries the "
          "information:")
    base = np.array([c[0] for c in curves["PULSE"]])[win]
    for name, _, _ in arms:
        if name == "PULSE":
            continue
        other = np.array([c[0] for c in curves[name]])[win]
        print(f"    PULSE is {np.mean(other / base):.1f}x better than {name} "
              f"(range {np.min(other / base):.1f}-{np.max(other / base):.1f}x)")

    print("\n  DOSE, over the whole grid -- reported because it does NOT favour "
          "PULSE:")
    for name, _, _ in arms:
        d = np.array([c[1] for c in curves[name]])
        fin = d[np.isfinite(d) & (d < 3)]
        print(f"    {name:>18}: median {np.median(fin):.0%}, best {fin.min():.0%}")

    late = ages >= 4.0
    if late.any():
        lp = np.array([c[0] for c in curves["PULSE"]])[late]
        lf = np.array([c[0] for c in curves["FINAL"]])[late]
        if np.mean(lf) < np.mean(lp):
            print(f"\n  AND THE PART THAT ARGUES AGAINST US: beyond 4 h the "
                  f"STATIC sensor\n  times better ({np.mean(lf):.0f} vs "
                  f"{np.mean(lp):.0f} min on average). An adaptive sensor "
                  f"switches\n  itself off, and the information goes with it.")

    print("\n  Choosing PULSE is therefore a trade, not a free win: it buys "
          "timing\n  resolution in the operating window and biological fidelity "
          "(B1: a static\n  Hill cannot produce the measured transient at all), "
          "and it costs dose\n  accuracy and late-window timing. State the price "
          "when quoting the choice.")

    # ---- figure ------------------------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for name, _, _ in arms:
        ax[0].plot(ages, [c[0] for c in curves[name]], "-o", ms=3, label=name)
        ax[1].plot(ages, [100 * c[1] for c in curves[name]], "-o", ms=3, label=name)
    ax[0].axvspan(0.5, 2.0, color="grey", alpha=0.12)
    ax[0].set_yscale("log"); ax[0].set_xlabel("age of the stress event (h)")
    ax[0].set_ylabel("timing error (min)")
    ax[0].set_title(f"Timing -- {variant}\n(shaded: window where timing is the job)")
    ax[0].legend(fontsize=8)
    ax[1].set_yscale("log"); ax[1].set_ylim(1, 300)
    ax[1].set_xlabel("age of the stress event (h)")
    ax[1].set_ylabel("dose error (%)")
    ax[1].set_title(f"Dose -- {variant}\n(the axis where PULSE does not win)")
    ax[1].legend(fontsize=8)
    plt.tight_layout()

    RESULTS.mkdir(exist_ok=True)
    commit = mf["repos"]["PULSE"]["commit"][:8]
    stem = f"B3_model_choice_{variant}_{commit}"
    plt.savefig(RESULTS / f"{stem}.png", dpi=170)
    with open(RESULTS / f"{stem}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\n  -> {stem}.png / .csv\n  caption: {P.caption(variant)}")
    return rows


if __name__ == "__main__":
    main(*(sys.argv[1:] or ["ox"]))