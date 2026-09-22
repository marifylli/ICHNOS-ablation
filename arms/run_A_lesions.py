"""
ARM A -- lesion ladder WITHIN a single model (design counterfactual).

QUESTION: what does the biosensor lose if circuit X is removed?
No fitting to data -> AIC/AICc DOES NOT APPLY here.

  python arms/run_A_lesions.py            # PULSE, ox
  python arms/run_A_lesions.py PULSE er
"""
import sys, csv
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import protocol as P, ichnos_ablation as A
from models import build_snapshots as S

RESULTS = Path(__file__).resolve().parents[1] / "results"


def main(tag="PULSE", variant="ox"):
    mf = S.verify()
    commit = mf["repos"][tag]["commit"][:8]
    sbml = S.path(variant, tag).read_text()

    rows = []
    for label, transform, ov, note in A.lesion_set(sbml, variant):
        model = transform(sbml, variant) if transform else sbml
        res = A.sweep(model, variant, ov)
        ok, why = A.gate(res)
        rows.append(dict(
            lesion=label, gate=why, note=note,
            fold_range=res["fold_range"], fold_vs_zero=res["fold_vs_zero"],
            # SHAPE METRICS ARE WRITTEN ONLY IF THE GATE WAS PASSED.
            # n_hill and lin_r2 come from a real fit on the WIDE grid and
            # carry an SE, so they are reported for EVERY lesion, gated or
            # not: a saturated circuit has a perfectly meaningful n_hill
            # (that is precisely the finding). Previously the gate blanked
            # them out and no_TetR_feedback -- the headline result -- showed
            # as "-".
            n_hill=res["n_hill"], se_n_hill=res["se_n_hill"],
            ec50=res["ec50"], lin_r2=res["lin_r2"],
            # The SHAPE metrics stay gated: those are the ones a dead
            # circuit can win.
            ratio_spread_3h=res["ratio_spread_at_t"] if ok else np.nan,
            ratio_spread_max=res["ratio_spread_max"] if ok else np.nan))

    print(f"\n=== A: lesion ladder -- {tag}/{variant} @ {commit}, theta={P.THETA_DEFAULT}")
    print(f"{'lesion':>20} {'gate':>18} {'fold':>7} {'vs0':>8} "
          f"{'n_hill':>14} {'EC50':>8} {'linR2':>7} {'spr3h':>8}")
    f = lambda v: "   -   " if np.isnan(v) else f"{v:7.3f}"
    for r in rows:
        nh = ("      -       " if np.isnan(r["n_hill"])
              else f"{r['n_hill']:6.2f}+/-{r['se_n_hill']:.2f}")
        print(f"{r['lesion']:>20} {r['gate']:>18} {r['fold_range']:7.3f} "
              f"{r['fold_vs_zero']:8.3f} {nh:>14} {r['ec50']:8.0f} "
              f"{f(r['lin_r2'])} {f(r['ratio_spread_3h'])}")

    print("\n  threshold sensitivity (which lesions pass the gate):")
    for th in P.THETA_SWEEP:
        keep = []
        for label, transform, ov, _ in A.lesion_set(sbml, variant):
            model = transform(sbml, variant) if transform else sbml
            if A.gate(A.sweep(model, variant, ov), th)[0]:
                keep.append(label)
        print(f"    theta={th}: " + (", ".join(keep) if keep else "(none)"))

    # --- figure: lesion ladder ---
    # Gated lesions get no n_eff bar by design; the fold panel shows WHY they
    # were gated, which is the informative part (saturated vs dead).
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    lab = [r["lesion"] for r in rows]
    x = np.arange(len(lab))
    ax[0].bar(x - 0.2, [r["fold_range"] for r in rows], 0.4, label="fold_range")
    ax[0].bar(x + 0.2, [r["fold_vs_zero"] for r in rows], 0.4, label="fold_vs_zero")
    ax[0].axhline(P.THETA_DEFAULT, ls=":", c="r", lw=1, label=f"theta={P.THETA_DEFAULT}")
    ax[0].set_yscale("log"); ax[0].set_ylabel("fold-change")
    ax[0].set_title("Viability gate")
    ax[1].bar(x, [0 if np.isnan(r["n_hill"]) else r["n_hill"] for r in rows], 0.5,
              yerr=[0 if np.isnan(r["se_n_hill"]) else r["se_n_hill"] for r in rows],
              capsize=3,
              color=["tab:blue" if r["gate"] == "viable" else "tab:orange" for r in rows])
    ax[1].axhline(rows[0]["n_hill"], ls="--", c="k", lw=1, label="intact circuit")
    ax[1].axhline(1.0, ls="-.", c="grey", lw=.8, label="n=1 (non-cooperative)")
    ax[1].set_ylabel("Hill coefficient n_H (fitted, +/-SE)")
    ax[1].set_title("Dose-response cooperativity\n(orange = gated, n_H still valid)")
    for a in ax:
        a.set_xticks(x); a.set_xticklabels(lab, rotation=30, ha="right", fontsize=8)
        a.legend(fontsize=8)
    fig.suptitle(f"Arm A lesion ladder -- {tag}/{variant}", fontsize=11)
    plt.tight_layout()

    RESULTS.mkdir(exist_ok=True)
    plt.savefig(RESULTS / f"A_lesions_{tag}_{variant}_{commit}.png", dpi=170)
    out = RESULTS / f"A_lesions_{tag}_{variant}_{commit}.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\n  -> {out.name}\n  caption: {P.caption(variant)}")
    return rows


if __name__ == "__main__":
    main(*(sys.argv[1:] or ["PULSE", "ox"]))
