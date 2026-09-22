"""
ARM B2 -- static Hill sensor (FINAL) vs adaptive sensor (PULSE).

QUESTION: what does the sensor's STRUCTURE change in the READOUT?
No fitting to data -> AIC/AICc DOES NOT APPLY here. ("Which structure fits
the data" is arm B1, a separate script.)

WHY THIS IS CONTROLLED
  reporter_module_v2.sbml and TIP_TetR_binding.sbml are byte-identical
  across the two repos -- build_snapshots checks this and records it in the
  MANIFEST. If that stops being true, this arm is INVALID and the script
  warns.

THE CONFOUND
  beta_max_ox = 25 (FINAL) vs 118 (PULSE, after the 5 -> 23.6 correction).
  So "structural difference" is also a gain difference. The third arm re-runs
  FINAL with the PULSE beta_max. DO NOT omit it in a presentation.

ABOUT er
  ERModule is NOT identical across repos: structure, n_er (2.0 -> 3.97) and
  beta_basal (5 -> 0.77) all change at once. The er comparison is therefore
  NOT controlled and is reported as descriptive only.

  python arms/run_B2_structure.py            # ox
  python arms/run_B2_structure.py er
"""
import sys, csv
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import protocol as P, ichnos_ablation as A
from models import build_snapshots as S

RESULTS = Path(__file__).resolve().parents[1] / "results"


def main(variant="ox"):
    mf = S.verify()
    if not mf["downstream_identical"]:
        print("  ! WARNING: downstream differs across repos -- B2 is no longer "
              "a controlled comparison.")
    sb = {t: S.path(variant, t).read_text() for t in ("FINAL", "PULSE")}
    bmax = A.pname("beta_max", variant)
    b_pulse = __import__("tellurium").loadSBMLModel(sb["PULSE"])[A.ids(sb["PULSE"])[bmax]]

    # Confound control: re-run FINAL at PULSE's gain. GUARDED -- if the two
    # values already coincide the "control" is a no-op and must not be shown
    # (this silently happened for er: three identical rows presented as a
    # controlled comparison).
    arms = [("FINAL", "FINAL", {}), ("PULSE", "PULSE", {})]
    try:
        A.assert_override_bites(sb["FINAL"], bmax, b_pulse, f"B2/{variant} gain control")
        arms.append((f"FINAL_beta{b_pulse:g}", "FINAL", {bmax: b_pulse}))
    except RuntimeError as e:
        print(f"  ! gain-confound control SKIPPED: {e}")

    R, rows = {}, []
    print(f"\n=== B2: sensor structure -- variant {variant}")
    print(f"{'arm':>18} {'fold':>7} {'vs0':>7} {'n_H':>7} {'spr3h':>8} {'spr6h':>8} {'spr12h':>8}")
    for label, tag, ov in arms:
        res = A.sweep(sb[tag], variant, ov); R[label] = res
        t, s = res["ratio_spread_t"]
        sp = {q: float(np.interp(q, t, s)) for q in (3, 6, 12)}
        ok, why = A.gate(res)
        print(f"{label:>18} {res['fold_range']:7.3f} {res['fold_vs_zero']:7.3f} "
              f"{res['n_hill']:7.3f} " + " ".join(f"{sp[q]:8.4f}" for q in (3, 6, 12)))
        rows.append(dict(arm=label, gate=why, fold_range=res["fold_range"],
                         fold_vs_zero=res["fold_vs_zero"], n_hill=res["n_hill"],
                         se_n_hill=res["se_n_hill"], lin_r2=res["lin_r2"],
                         spread_3h=sp[3], spread_6h=sp[6], spread_12h=sp[12]))

    # --- figure ---
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))
    doses = P.DOSES[variant]
    cmap = plt.cm.viridis(np.linspace(0, .9, len(doses)))
    for a, tag in [(ax[0, 0], "FINAL"), (ax[0, 1], "PULSE")]:
        for c, d in zip(cmap, doses):
            t, og, ra = R[tag]["traj"][d]
            a.plot(t, ra, color=c, lw=1.6, label=f"{d}")
        a.axvline(P.PRIMARY_T, ls=":", c="k", lw=1)
        a.set_title(f"{tag}\nMeasured_Ratio_RG"); a.set_xlabel("t (h)"); a.set_ylabel("R/G")
    ax[0, 1].legend(fontsize=7, ncol=2, title="dose (uM)")
    a = ax[1, 0]
    for label, st in zip(R, ("-o", "-s", "--^")):
        a.plot(doses, [R[label]["green_at_t"][d] for d in doses], st, label=label)
    a.set_xscale("log"); a.set_xlabel("dose (uM)")
    a.set_ylabel(f"Observed_Green @ {P.PRIMARY_T} h")
    a.set_title("Dose-response"); a.legend(fontsize=8)
    a = ax[1, 1]
    for label, st in zip(R, ("-", "-", "--")):
        t, s = R[label]["ratio_spread_t"]; a.plot(t, s * 100, st, label=label)
    a.set_yscale("log"); a.axvline(P.PRIMARY_T, ls=":", c="k", lw=1)
    a.set_xlabel("t (h)"); a.set_ylabel("max relative ratio spread (%)")
    a.set_title("Dose-invariance of the ratio\n(asymptotic, not operational)")
    a.legend(fontsize=8)
    plt.tight_layout()

    RESULTS.mkdir(exist_ok=True)
    stem = f"B2_structure_{variant}_{mf['repos']['PULSE']['commit'][:8]}"
    plt.savefig(RESULTS / f"{stem}.png", dpi=170)
    with open(RESULTS / f"{stem}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\n  -> {stem}.png / .csv\n  caption: {P.caption(variant)}")
    return rows


if __name__ == "__main__":
    main(*(sys.argv[1:] or ["ox"]))
