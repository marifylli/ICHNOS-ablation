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
             "alpha_0 = 0: removes basal TIP production")]


def observation_model_lesions(sbml, variant):
    """NOT circuit parts -- parameters of how we OBSERVE the circuit.

    E (FRET bleed-through) was previously listed next to no_leaky_expression,
    which conflates "what the device does" with "how we measure it". A judge
    reading one table cannot tell which rows are biology. Separate section.
    """
    return [("no_FRET_correction", None, {"E": 0.0},
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
    groups = ([("circuit", l) for l in list(A.lesion_set(sbml, variant))
               + extra_lesions(sbml, variant)]
              + [("observation", l) for l in observation_model_lesions(sbml, variant)])
    for kind, (label, tf, ov, note) in groups:
        if label == "full":
            continue
        model = tf(sbml, variant) if tf else sbml
        g, r = outputs(model, ov)
        res = A.sweep(model, variant, ov)
        ok, why = A.gate(res)
        rows.append(dict(part=label, kind=kind, note=note, gate=why,
                         d_green_pct=100 * (g - g0) / g0,
                         d_ratio_pct=100 * (r - r0) / r0,
                         n_hill=res["n_hill"], se_n_hill=res["se_n_hill"],
                         lin_r2=res["lin_r2"],
                         fold_range=res["fold_range"]))
    rows.sort(key=lambda x: (x["kind"] != "circuit", -abs(x["d_green_pct"])))

    print(f"\n=== A3: sensitivity -- {tag}/{variant} @ {commit}, "
          f"reference dose {dose} uM, t*={P.PRIMARY_T} h")
    print(f"{'ablated part':>22} {'kind':>11} {'gate':>12} {'dGreen%':>9} "
          f"{'dRatio%':>9} {'n_H':>7} {'linR2':>7}")
    for r in rows:
        nh = "   -   " if np.isnan(r["n_hill"]) else f"{r['n_hill']:7.2f}"
        lr = "   -   " if np.isnan(r["lin_r2"]) else f"{r['lin_r2']:7.3f}"
        print(f"{r['part']:>22} {r['kind']:>11} {r['gate']:>12} "
              f"{r['d_green_pct']:+9.1f} {r['d_ratio_pct']:+9.1f} {nh} {lr}")

    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / f"A3_sensitivity_{tag}_{variant}_{commit}.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\n  -> {out.name}\n  caption: {P.caption(variant)}")
    return rows


if __name__ == "__main__":
    main(*(sys.argv[1:] or ["PULSE", "ox"]))
