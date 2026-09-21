"""
ARM A, step 4 -- DESIGN SPECIFICATION from the TIP-TetR affinity sweep.

WHY THIS IS A SEPARATE ARM
  Kd_TIP_TetR has no literature baseline: TIP is a novel molecule, so this
  parameter is exploratory rather than measured. Sweeping it does two jobs at
  once:
    (a) it is the cleanest demonstration that shape metrics need a gate --
        as binding weakens, ratio spread falls monotonically toward zero
        while the circuit dies, so "minimise the spread" would select a dead
        device;
    (b) it yields an actual specification for the wet lab: the affinity
        window in which the biosensor is functional.

  python arms/run_A4_kd_spec.py            # PULSE, ox
"""
import sys, csv
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import protocol as P, ichnos_ablation as A
from models import build_snapshots as S

RESULTS = Path(__file__).resolve().parents[1] / "results"
KD_GRID = (0.05, 0.25, 1, 5, 25, 100, 500, 2000, 5000)   # nM


def main(tag="PULSE", variant="ox"):
    mf = S.verify(); commit = mf["repos"][tag]["commit"][:8]
    sbml = S.path(variant, tag).read_text()

    rows = []
    print(f"\n=== A4: Kd_TIP_TetR specification -- {tag}/{variant} @ {commit}")
    print(f"{'Kd (nM)':>9} {'gate':>12} {'fold':>8} {'vs0':>8} {'spread@3h':>11} {'n_eff':>8}")
    for kd in KD_GRID:
        res = A.sweep(sbml, variant, {"Kd_TIP_TetR": kd})
        ok, why = A.gate(res)
        rows.append(dict(kd_nM=kd, gate=why, fold_range=res["fold_range"],
                         fold_vs_zero=res["fold_vs_zero"],
                         ratio_spread_3h=res["ratio_spread_at_t"],
                         n_eff=res["n_eff"]))
        print(f"{kd:>9} {why:>12} {res['fold_range']:8.3f} {res['fold_vs_zero']:8.3f} "
              f"{res['ratio_spread_at_t']:11.5f} {res['n_eff']:8.3f}")

    viable = [r["kd_nM"] for r in rows if r["gate"] == "viable"]
    if viable:
        print(f"\n  SPECIFICATION: functional for Kd <= {max(viable):g} nM "
              f"(theta = {P.THETA_DEFAULT}); non-responsive above.")
    print("  NOTE: ratio spread falls monotonically as the circuit dies -- "
          "this is why shape metrics are gated, never optimised directly.")

    kd = [r["kd_nM"] for r in rows]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(kd, [r["fold_range"] for r in rows], "-o", label="fold_range (signal)")
    ax.plot(kd, [r["ratio_spread_3h"] for r in rows], "-s",
            label="ratio spread @3h (apparent 'timer quality')")
    ax.axhline(P.THETA_DEFAULT, ls=":", c="r", lw=1,
               label=f"viability threshold theta = {P.THETA_DEFAULT}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Kd_TIP_TetR (nM)"); ax.set_ylabel("metric value")
    ax.set_title("Shape metrics improve as the circuit dies\n"
                 "(why a viability gate is required)")
    ax.legend(fontsize=8); plt.tight_layout()

    RESULTS.mkdir(exist_ok=True)
    stem = f"A4_kd_spec_{tag}_{variant}_{commit}"
    plt.savefig(RESULTS / f"{stem}.png", dpi=170)
    with open(RESULTS / f"{stem}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\n  -> {stem}.png / .csv\n  caption: {P.caption(variant)}")
    return rows


if __name__ == "__main__":
    main(*(sys.argv[1:] or ["PULSE", "ox"]))
