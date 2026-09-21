"""
ARM A, step 3 -- SENSITIVITY COMPARISON (iGEM Engineering Success framework).

Percent change in the primary output per ablated part, relative to the intact
circuit, at a reference dose. This is what tells the Wet Lab which parts must
be tuned precisely and which are tolerant.

Reported at the readout time t* (not at true steady state): the decoder
operates at 3 h, and the circuit is still transient there. Steady-state
numbers would describe a regime the device is never read in.

  python arms/run_A3_sensitivity.py            # PULSE, ox, dose 200 uM
"""
import sys, csv
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import protocol as P, ichnos_ablation as A
from models import build_snapshots as S

RESULTS = Path(__file__).resolve().parents[1] / "results"
REF_DOSE = {"ox": 200, "er": 2000}   # near K_act: most informative operating point


def extra_lesions(sbml, variant):
    """Parameter ablations in the sense of the iGEM framework (set to zero)."""
    return [("no_leaky_expression", None, {A.pname("beta_basal", variant): 0.0},
             "alpha_0 = 0: removes basal TIP production"),
            ("no_FRET_correction", None, {"E": 0.0},
             "E = 0: removes the FRET bleed-through term")]


def main(tag="PULSE", variant="ox"):
    mf = S.verify(); commit = mf["repos"][tag]["commit"][:8]
    sbml = S.path(variant, tag).read_text()
    dose = REF_DOSE[variant]

    def outputs(model, ov):
        t, og, ra = A.run(model, variant, dose, ov)
        return A.at(t, og, P.PRIMARY_T), A.at(t, ra, P.PRIMARY_T)

    g0, r0 = outputs(sbml, {})
    rows = []
    for label, tf, ov, note in list(A.lesion_set(sbml, variant)) + extra_lesions(sbml, variant):
        if label == "full":
            continue
        model = tf(sbml, variant) if tf else sbml
        g, r = outputs(model, ov)
        res = A.sweep(model, variant, ov)
        ok, why = A.gate(res)
        rows.append(dict(part=label, note=note, gate=why,
                         d_green_pct=100 * (g - g0) / g0,
                         d_ratio_pct=100 * (r - r0) / r0,
                         n_eff=res["n_eff"] if ok else np.nan,
                         fold_range=res["fold_range"]))
    rows.sort(key=lambda x: -abs(x["d_green_pct"]))

    print(f"\n=== A3: sensitivity -- {tag}/{variant} @ {commit}, "
          f"reference dose {dose} uM, t*={P.PRIMARY_T} h")
    print(f"{'ablated part':>22} {'gate':>12} {'dGreen%':>9} {'dRatio%':>9} {'n_eff':>7}")
    for r in rows:
        ne = "   -   " if np.isnan(r["n_eff"]) else f"{r['n_eff']:7.3f}"
        print(f"{r['part']:>22} {r['gate']:>12} {r['d_green_pct']:+9.1f} "
              f"{r['d_ratio_pct']:+9.1f} {ne}")

    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / f"A3_sensitivity_{tag}_{variant}_{commit}.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\n  -> {out.name}\n  caption: {P.caption(variant)}")
    return rows


if __name__ == "__main__":
    main(*(sys.argv[1:] or ["PULSE", "ox"]))
