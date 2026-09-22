"""
ARM A, step 5 -- DECODING BOUND per lesion, and justification of t*.

WHY THIS ARM EXISTS
  Every other metric in this study is a proxy. The deliverable is a decoder
  returning (dose, time-since-onset) per cell, so the question a judge will
  ask is "how much worse does the estimate get without part X?", in hours
  and in percent -- not in units of ratio spread.

  decoding_crlb() answers exactly that: the Cramer-Rao lower bound on the
  joint estimate of (t, ln dose) from the two observables under
  multiplicative measurement noise. It is a LOWER bound: no decoder, however
  clever, does better. So a lesion that doubles sigma_t has halved the
  device's timing resolution, full stop.

  It also removes the need for the viability gate as a *scoring* device: a
  dead circuit has a singular Fisher matrix and its bound is infinite, so it
  cannot win by being dead -- unlike ratio_spread, which it minimises.

THREE OUTPUTS
  1. sigma_t / sigma_dose per lesion at the reference dose and t*.
  2. TIMING RESOLUTION vs AGE OF THE STRESS EVENT, which is NOT a choice of
     when to image. In use the age is the UNKNOWN the decoder estimates, so
     this curve characterises the device: how precisely can an event be
     dated, given how long ago it happened. It yields the honest version of
     the "memory window" claim, and it is also what tells the wet lab which
     time points a validation experiment must cover.
  3. POPULATION decoding. The decoder reads cultures, not single cells, so
     the relevant bound is for a mean over N cells, with a noise floor from
     errors shared by the whole field. Output: sigma_t(N) per shared-noise
     level, the floor it cannot beat, and the minimum N that reaches the
     target resolution -- a number the imaging protocol can use.

SCOPE -- what this arm deliberately does NOT model
  Measurement physics (Poisson photon statistics, autofluorescence,
  photobleaching, per-channel quantum yield, detection threshold) belongs to
  the image-pipeline Fisher analysis, not here. The two connect at exactly
  one point: the per-channel variances passed to crlb_from_jacobian(). When
  that model produces a CV, it replaces P.NOISE_CV and nothing else changes.

  python arms/run_A5_decoding.py            # PULSE, ox
  python arms/run_A5_decoding.py PULSE er
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


def main(tag="PULSE", variant="ox"):
    mf = S.verify(); commit = mf["repos"][tag]["commit"][:8]
    sbml = S.path(variant, tag).read_text()
    dose = REF_DOSE[variant]

    # ---- 1. per-lesion decoding bound -----------------------------------
    rows = []
    print(f"\n=== A5: decoding bound -- {tag}/{variant} @ {commit}, "
          f"dose {dose}, t*={P.PRIMARY_T} h, CV={P.NOISE_CV['green']:.0%}")
    print(f"{'lesion':>22} {'sigma_t (h)':>12} {'sigma_dose':>11} "
          f"{'corr':>7} {'vs intact':>10}")
    base = None
    for label, tf, ov, _ in A.lesion_set(sbml, variant):
        model = tf(sbml, variant) if tf else sbml
        c = A.decoding_crlb(model, variant, dose, overrides=ov)
        if label == "full":
            base = c["sigma_t"]
        rel = c["sigma_t"] / base if base else np.nan
        rows.append(dict(lesion=label, sigma_t_h=c["sigma_t"],
                         sigma_dose_rel=c["sigma_lnD"], corr=c["corr"],
                         sigma_t_vs_intact=rel))
        st = "   inf   " if not np.isfinite(c["sigma_t"]) else f"{c['sigma_t']:12.3f}"
        print(f"{label:>22} {st} {c['sigma_lnD']:10.1%} "
              f"{c['corr']:+7.3f} {rel:9.2f}x")

    print("\n  corr near +/-1 means the two channels are nearly redundant at\n"
          "  this operating point: dose and time trade off against each other\n"
          "  and the pair is only weakly identifiable. That is the quantitative\n"
          "  form of 'both observables are structurally necessary'.")

    # ---- 2. timing resolution vs AGE of the event -------------------------
    # The project previously fixed t* = 3 h and treated it as "the readout
    # time". That framing only applies to a calibration experiment, where the
    # onset is controlled. In deployment the age IS the unknown, so the right
    # object is sigma_t(age): the precision with which an event of that age
    # can be dated. Reported for a population, because the decoder reads
    # cultures.
    N = P.N_REPORT
    ages, st, sd = [], [], []
    print(f"\n  timing resolution vs age of the event "
          f"(N = {N} cells, per-cell CV = {P.NOISE_CV['green']:.0%}):")
    print(f"    {'age (h)':>8} {'sigma_t (min)':>14} {'sigma_dose':>11}  "
          f"{'shared 2%':>10} {'shared 5%':>10}")
    for age in P.AGE_GRID:
        Ja = A.decoding_jacobian(sbml, variant, dose, t_star=float(age))
        def sig(cs):
            v = A.population_variance(P.NOISE_CV["green"], cs, N)
            return A.crlb_from_jacobian(Ja, v, v)
        c0, c2, c5 = sig(0.0), sig(0.02), sig(0.05)
        ages.append(age); st.append(c0["sigma_t"]); sd.append(c0["sigma_lnD"])
        print(f"    {age:8.2f} {60 * c0['sigma_t']:14.1f} {c0['sigma_lnD']:10.0%}  "
              f"{60 * c2['sigma_t']:9.1f}' {60 * c5['sigma_t']:9.1f}'")

    st = np.array(st); sd = np.array(sd); ages_a = np.array(ages)
    ok_t = ages_a[st <= P.TARGET_SIGMA_T_H]
    ok_d = ages_a[sd <= 0.30]
    print(f"\n  USABLE WINDOW at N={N}, no shared error:")
    if ok_t.size:
        print(f"    timing to better than {60 * P.TARGET_SIGMA_T_H:.0f} min: "
              f"events {ok_t.min():g}-{ok_t.max():g} h old")
    if ok_d.size:
        print(f"    dose to better than 30%:            "
              f"events {ok_d.min():g}-{ok_d.max():g} h old")
    print("    The two are COMPLEMENTARY: timing degrades with age (the ratio\n"
          "    flattens as the population equilibrates) while dose improves\n"
          "    (signal accumulates). This is the quantitative form of the\n"
          "    'sensor says HOW MUCH, timer says WHEN' design claim.")
    print(f"  NOTE: t* = {P.PRIMARY_T} h is the CALIBRATION readout (Dacquay,\n"
          "  accumulated reporter at 3 h) and is kept for every comparison in\n"
          "  this study. It is not a claim that the device is read at 3 h --\n"
          "  in use the age is unknown. A validation experiment should span\n"
          "  the grid above, not sit at a single time point.")

    # ---- 3. noise sensitivity -------------------------------------------
    print("\n  noise sweep (intact circuit, sigma_t in h):")
    for cv in P.NOISE_CV_SWEEP:
        c = A.decoding_crlb(sbml, variant, dose,
                            cv={"green": cv, "red": cv})
        print(f"    CV={cv:.0%}: sigma_t = {c['sigma_t']:.3f} h, "
              f"sigma_dose = {c['sigma_lnD']:.1%}")
    print("  sigma scales linearly with CV (F ~ 1/CV^2), so the RANKING of\n"
          "  lesions is CV-independent -- only the absolute hours move. That\n"
          "  is why the conclusions survive the image pipeline not being\n"
          "  calibrated yet.")

    # ---- 4. population decoding ------------------------------------------
    # J is computed ONCE; every (N, shared CV) point is then exact algebra on
    # the same Fisher matrix -- no further simulation.
    J = A.decoding_jacobian(sbml, variant, dose)
    cv_cell = P.NOISE_CV["green"]
    target = P.TARGET_SIGMA_T_H
    pop_rows = []
    print(f"\n  population decoding (intact circuit, per-cell CV = {cv_cell:.0%}, "
          f"target sigma_t <= {target} h):")
    head = "    shared CV " + "".join(f"{'N=' + str(n):>9}" for n in P.N_CELLS_SWEEP) \
           + f"{'floor':>9}{'N needed':>10}"
    print(head)
    for cs in P.SHARED_CV_SWEEP:
        cells = []
        for n in P.N_CELLS_SWEEP:
            v = A.population_variance(cv_cell, cs, n)
            c = A.crlb_from_jacobian(J, v, v)
            cells.append(c["sigma_t"])
            pop_rows.append(dict(shared_cv=cs, n_cells=n,
                                 sigma_t_h=c["sigma_t"],
                                 sigma_dose_rel=c["sigma_lnD"]))
        floor = (A.crlb_from_jacobian(J, cs ** 2, cs ** 2)["sigma_t"]
                 if cs > 0 else 0.0)
        # smallest N reaching the target, solved exactly rather than read
        # off the grid: sigma_t scales as sqrt(var), var = cv^2/N + cs^2
        s1 = A.crlb_from_jacobian(J, 1.0, 1.0)["sigma_t"]   # sigma_t at var = 1
        need_var = (target / s1) ** 2
        if need_var <= cs ** 2:
            n_req = "never"            # the floor alone is above the target
        else:
            n_req = f"{int(np.ceil(cv_cell ** 2 / (need_var - cs ** 2)))}"
        print(f"    {cs:>8.0%} " + "".join(f"{x:9.3f}" for x in cells)
              + f"{floor:9.3f}{n_req:>10}")
    print("  sigma_t in hours. 'floor' = N -> infinity: the part no amount of\n"
          "  cell counting removes. Once the N=... columns approach it, image\n"
          "  quality (calibration), not cell number, is the limit.")
    xs = [f"{(cv_cell / cs) ** 2:.0f} cells at {cs:.0%}"
          for cs in P.SHARED_CV_SWEEP if cs > 0]
    print("  crossover N* = (CV_cell / CV_shared)^2, beyond which more cells\n"
          "  buy almost nothing: " + "; ".join(xs) + ".")

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.2))
    lab = [r["lesion"] for r in rows]
    v = [min(r["sigma_t_h"], 1e3) for r in rows]
    ax[0].bar(np.arange(len(lab)), v, 0.55,
              color=["tab:green" if r["lesion"] == "full" else "tab:blue" for r in rows])
    ax[0].axhline(base, ls="--", c="k", lw=1, label="intact circuit")
    ax[0].set_xticks(np.arange(len(lab)))
    ax[0].set_xticklabels(lab, rotation=30, ha="right", fontsize=8)
    ax[0].set_ylabel("sigma_t (h), Cramer-Rao lower bound")
    ax[0].set_title("Timing resolution lost per ablation"); ax[0].legend(fontsize=8)
    ax[1].plot(ages, 60 * st, "-o", ms=3, label=f"sigma_t (N={P.N_REPORT})")
    ax[1].axhline(60 * P.TARGET_SIGMA_T_H, ls=":", c="r", lw=1,
                  label=f"target {60 * P.TARGET_SIGMA_T_H:.0f} min")
    a1 = ax[1].twinx()
    a1.plot(ages, 100 * sd, "-s", ms=3, c="tab:green", label="sigma_dose (%)")
    a1.set_ylabel("dose error (%)", color="tab:green"); a1.set_ylim(0, 100)
    ax[1].set_yscale("log"); ax[1].set_xlabel("age of the stress event (h)")
    ax[1].set_ylabel("timing error (min)")
    ax[1].set_title("Resolution vs age of the event\n(not a choice of readout time)")
    ax[1].legend(fontsize=7, loc="upper left")
    ns = np.array(P.N_CELLS_SWEEP, float)
    for cs in P.SHARED_CV_SWEEP:
        y = [r["sigma_t_h"] for r in pop_rows if r["shared_cv"] == cs]
        ax[2].plot(ns, y, "-o", ms=3, label=f"shared CV {cs:.0%}")
    ax[2].axhline(target, ls=":", c="r", lw=1, label=f"target {target} h")
    ax[2].set_xscale("log"); ax[2].set_yscale("log")
    ax[2].set_xlabel("cells averaged per measurement, N")
    ax[2].set_ylabel("sigma_t (h)")
    ax[2].set_title("Population decoding\n(floors where shared noise dominates)")
    ax[2].legend(fontsize=7)
    plt.tight_layout()

    RESULTS.mkdir(exist_ok=True)
    stem = f"A5_decoding_{tag}_{variant}_{commit}"
    plt.savefig(RESULTS / f"{stem}.png", dpi=170)
    with open(RESULTS / f"{stem}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with open(RESULTS / f"{stem}_population.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pop_rows[0])); w.writeheader(); w.writerows(pop_rows)
    print(f"\n  -> {stem}.png / .csv / _population.csv\n  caption: {P.caption(variant)}")
    return rows


if __name__ == "__main__":
    main(*(sys.argv[1:] or ["PULSE", "ox"]))