"""
ICHNOS ablation — shared core (simulation, metrics, gate, lesions).

THE THREE ARMS AND WHY THEY ARE SEPARATE
  A  (arms/run_A_lesions.py)    design counterfactual WITHIN one model.
                                No fitting -> AIC/AICc DOES NOT APPLY.
  B2 (arms/run_B2_structure.py) same protocol, two sensor STRUCTURES.
                                No fitting -> AIC/AICc DOES NOT APPLY.
  B1 (arms/run_B1_aicc.py)      fit to DATA (Delaunay 2000).
                                Only here is there a likelihood -> only here AICc.
  aicc() lives exclusively in run_B1_aicc.py and is deliberately NOT imported
  here, so it cannot be called by accident from A or B2.

THE VIABILITY GATE
  ratio_spread is MINIMISED by a dead circuit (measured: at Kd = 5000 nM,
  spread = 0.0004 while fold = 1.003), so shape metrics without a gate
  reward the dead circuit. The gate is therefore BINARY and COMES FIRST.

  The deeper fix is decoding_crlb(): a dead circuit has a singular Fisher
  matrix and its bound diverges, so it cannot win by being dead. The gate is
  kept because it is what SEPARATES dead from saturated -- a classification,
  not a score.

  TWO FAILURE MODES -- do not conflate them in the presentation:
    dead      : fold_range AND fold_vs_zero both low -> no response at all
    saturated : fold_vs_zero HIGH, fold_range low -> responds strongly but
                cannot discriminate doses (the linear range is gone)
  The second is a RESULT, not a rejection.
"""
import numpy as np, libsbml, tellurium as te
from scipy.optimize import curve_fit
import protocol as P


def pname(base, variant):
    """Parameter names carry an _ox / _er suffix in the SBML."""
    return f"{base}_{variant}"


def ids(sbml_str):
    doc = libsbml.readSBMLFromString(sbml_str)
    m = doc.getModel()
    d = {}
    for lst in (m.getListOfSpecies(), m.getListOfParameters()):
        for e in lst:
            d[e.getName() or e.getId()] = e.getId()
    return d


def has(sbml_str, name):
    return name in ids(sbml_str)


def run(sbml_str, variant, dose, overrides=None):
    """One simulation with pre-equilibration.

    Overrides are applied BEFORE pre-equilibration so that the lesion also
    holds at baseline -- otherwise you would be comparing circuits that
    start from different states.
    """
    I = ids(sbml_str)
    S = I[pname("S", variant)]
    r = te.loadSBMLModel(sbml_str); r.reset()
    for k, v in (overrides or {}).items():
        r[I[k]] = v
    r[S] = 0.0
    if P.PRE_T > 0:
        r.simulate(0, P.PRE_T, 200)
    r[S] = dose
    sel = ["time", f"[{I['Reporter_green']}]", f"[{I['Reporter_red']}]",
           f"[{I['Reporter_dark_red']}]", f"[{I['Reporter_intermediate']}]"]
    res = r.simulate(0, P.T_END, P.NPTS, selections=sel)
    t   = np.array(res[:, 0])
    g   = np.array(res[:, 1])
    red = np.array(res[:, 2])
    # COMMON BUG: Observed_Green and Measured_Ratio_RG are assignment-rule
    # PARAMETERS. roadrunner does not return them in the default selections --
    # you would silently read back the initial value 1.0. Recomputed here
    # from the species.
    tot = red + np.array(res[:, 3]) + np.array(res[:, 4])
    eps, E = r[I["eps"]], r[I["E"]]
    bf = red / (tot + eps)
    og = g * (1 - bf * E)
    return t, og, red / (og + eps)


def at(t, y, tq):
    """Linear interpolation of a trajectory at a query time."""
    return float(np.interp(tq, t, y))


def hill_fit(doses, resp):
    """4-parameter Hill fit:  y = base + (top - base) * d^n / (d^n + EC50^n).

    REPLACES the old logit-slope n_eff, which was NOT a measurement of
    cooperativity. That estimator min-max normalised the response on the
    6-point grid and regressed logit(y) on log10(d). Verified failure modes:

        true Hill n=1         -> returned 2.03
        a PURELY LINEAR y = d -> returned 2.03
        true Hill n=4         -> returned 4.29 or 8.85 depending only on EC50

    i.e. "n_eff ~ 2.0 = linearised" was the estimator's fixed point for any
    gentle curve, not a property of the circuit. Reported numbers must come
    from an actual fit, with an uncertainty, or they are not measurements.

    Returns base, top, ec50, n_hill, their standard errors, and the fit R^2.
    Requires the WIDE grid (P.FIT_DOSES): on a grid that samples neither
    plateau the four parameters are not jointly identifiable and the SEs
    correctly blow up -- which is the point of reporting them.
    """
    d = np.asarray(doses, float)
    y = np.asarray(resp, float)
    if y.max() - y.min() <= 0 or d.size < 5:
        return dict(base=np.nan, top=np.nan, ec50=np.nan, n_hill=np.nan,
                    se_n_hill=np.nan, se_ec50=np.nan, r2=np.nan)
    p0 = [y.min(), y.max(), float(np.exp(np.mean(np.log(d)))), 2.0]
    try:
        p, cov = curve_fit(_hill, d, y, p0=p0,
                           bounds=([0, 0, d.min() / 10, 0.2],
                                   [np.inf, np.inf, d.max() * 10, 20.0]),
                           maxfev=40000)
        se = np.sqrt(np.diag(cov))
    except Exception:
        return dict(base=np.nan, top=np.nan, ec50=np.nan, n_hill=np.nan,
                    se_n_hill=np.nan, se_ec50=np.nan, r2=np.nan)
    ss = ((y - _hill(d, *p)) ** 2).sum()
    return dict(base=p[0], top=p[1], ec50=p[2], n_hill=p[3],
                se_ec50=se[2], se_n_hill=se[3],
                r2=1 - ss / ((y - y.mean()) ** 2).sum())


def _hill(d, base, top, ec50, n):
    return base + (top - base) * d ** n / (d ** n + ec50 ** n)


def linearity_r2(doses, resp, variant):
    """R^2 of a STRAIGHT LINE inside the operating window only.

    The design claim is "the feedback linearises the response", and this is
    the direct test of it: no normalisation, no transform, nothing that can
    manufacture a number on a dead signal. n_hill says how cooperative the
    curve is; this says how usable it is where the device is actually read.
    """
    lo, hi = P.OPERATING_WINDOW[variant]
    d = np.asarray(doses, float); y = np.asarray(resp, float)
    k = (d >= lo) & (d <= hi)
    if k.sum() < 3 or y[k].max() - y[k].min() <= 0:
        return np.nan
    res = y[k] - np.polyval(np.polyfit(d[k], y[k], 1), d[k])
    return float(1 - (res ** 2).sum() / ((y[k] - y[k].mean()) ** 2).sum())


# --------------------------------------------------------------------------
# DECODING BOUND (Cramer-Rao) -- the metric the DEVICE is actually judged on
# --------------------------------------------------------------------------
def decoding_crlb(sbml_str, variant, dose, t_star=None, cv=None,
                  overrides=None):
    """Lower bound on the decoder's error in (t_since_onset, log-dose).

    WHY THIS EXISTS
      ratio_spread is MINIMISED BY A DEAD CIRCUIT, which is exactly why the
      viability gate had to be bolted on in front of it. That is a patch over
      the wrong metric. The deliverable is a decoder returning (dose, t), so
      the honest figure of merit is HOW WELL (dose, t) CAN BE RECOVERED.

      Here the gate is not needed as a concept: a dead circuit has a singular
      Fisher matrix, so the bound goes to infinity automatically and the row
      excuses itself. (The gate is KEPT, because it is what classifies dead
      vs saturated, and "saturated" is a result.)

    MODEL
      Observables: Observed_Green and Reporter_red, each with multiplicative
      lognormal noise of coefficient of variation cv -> on a log scale the
      noise is additive with constant sigma, so the Fisher matrix is
        F = J^T Sigma^-1 J,   J = d(log observable) / d(t, log dose).
      Sigma^-1 is diagonal, which is an ASSUMPTION: it treats the two
      channels as independently noisy. Shared extrinsic noise (cell size,
      expression capacity) is correlated and would make the true bound
      slightly worse. State this rather than hide it.

      The time derivative is analytic-free (finite difference on the
      trajectory); the dose derivative is a central difference in ln(dose),
      which is why dose error comes out as a RELATIVE error -- the natural
      unit for a dose-response device.

    RETURNS
      sigma_t     : hours -- irreducible timing error at that (dose, t*)
      sigma_lnD   : relative dose error (0.20 = +/-20%)
      corr        : t/dose correlation. Near +/-1 means the two observables
                    are nearly redundant and the pair is weakly identifiable
                    -- the quantitative version of "we need both channels".
    """
    t_star = P.PRIMARY_T if t_star is None else t_star
    cv = P.NOISE_CV if cv is None else cv
    h = P.FD_REL_STEP

    def obs(d):
        t, og, _ = run(sbml_str, variant, d, overrides)
        red = _red_trace(sbml_str, variant, d, overrides)
        return t, np.log(np.maximum(og, 1e-12)), np.log(np.maximum(red, 1e-12))

    t, lg0, lr0 = obs(dose)
    _, lg_p, lr_p = obs(dose * np.exp(h))
    _, lg_m, lr_m = obs(dose * np.exp(-h))

    # d/dt at t*, from the trajectory itself
    dg_dt = float(np.interp(t_star, t, np.gradient(lg0, t)))
    dr_dt = float(np.interp(t_star, t, np.gradient(lr0, t)))
    # d/dln(dose), central difference
    dg_dD = (at(t, lg_p, t_star) - at(t, lg_m, t_star)) / (2 * h)
    dr_dD = (at(t, lr_p, t_star) - at(t, lr_m, t_star)) / (2 * h)

    J = np.array([[dg_dt, dg_dD], [dr_dt, dr_dD]])
    return crlb_from_jacobian(J, cv["green"] ** 2, cv["red"] ** 2)


def decoding_jacobian(sbml_str, variant, dose, t_star=None, overrides=None):
    """J = d(ln Green, ln Red) / d(t, ln dose) at (dose, t*).

    Exposed separately because J depends only on the MODEL, while the noise
    covariance depends only on the MEASUREMENT. Computing J once and then
    varying the noise analytically is exact (F = J^T Sigma^-1 J is linear in
    Sigma^-1) and needs no further simulation -- which is what makes the
    population sweep below essentially free.

    This split is also the interface to the image-pipeline noise model:
    whatever that model produces (Poisson photon counts, autofluorescence,
    detection floor) enters ONLY through the variances passed to
    crlb_from_jacobian(). Nothing here has to change when it arrives.
    """
    t_star = P.PRIMARY_T if t_star is None else t_star
    h = P.FD_REL_STEP

    def obs(d):
        t, og, _ = run(sbml_str, variant, d, overrides)
        red = _red_trace(sbml_str, variant, d, overrides)
        return t, np.log(np.maximum(og, 1e-12)), np.log(np.maximum(red, 1e-12))

    t, lg0, lr0 = obs(dose)
    _, lg_p, lr_p = obs(dose * np.exp(h))
    _, lg_m, lr_m = obs(dose * np.exp(-h))
    return np.array([
        [float(np.interp(t_star, t, np.gradient(lg0, t))),
         (at(t, lg_p, t_star) - at(t, lg_m, t_star)) / (2 * h)],
        [float(np.interp(t_star, t, np.gradient(lr0, t))),
         (at(t, lr_p, t_star) - at(t, lr_m, t_star)) / (2 * h)]])


def crlb_from_jacobian(J, var_green, var_red):
    """Cramer-Rao bound from a Jacobian and per-channel variances (on the
    log scale, i.e. squared CVs for small noise)."""
    Sinv = np.diag([1 / var_green, 1 / var_red])
    F = J.T @ Sinv @ J
    if not np.all(np.isfinite(F)) or abs(np.linalg.det(F)) < 1e-14:
        return dict(sigma_t=np.inf, sigma_lnD=np.inf, corr=np.nan)
    C = np.linalg.inv(F)
    return dict(sigma_t=float(np.sqrt(C[0, 0])),
                sigma_lnD=float(np.sqrt(C[1, 1])),
                corr=float(C[0, 1] / np.sqrt(C[0, 0] * C[1, 1])))


def population_variance(cv_cell, cv_shared, n_cells):
    """Variance of a POPULATION-AVERAGED measurement, per channel.

    The decoder reads cultures, not single cells, so the observable is a mean
    over N cells. Its noise has two parts that behave completely differently:

        var = cv_cell^2 / N  +  cv_shared^2

      cv_cell   -- independent per cell (expression noise, per-cell shot
                   noise, segmentation error). Averages away as 1/N.
      cv_shared -- common to every cell in the field (illumination drift,
                   exposure, flat-field, background subtraction, crosstalk
                   coefficient). Identical for all N cells, so averaging
                   does NOTHING to it.

    Consequence: sigma_t does not go to zero with N -- it floors at the value
    set by cv_shared alone. Past the point where cv_cell^2/N ~ cv_shared^2,
    measuring more cells buys nothing and only calibration helps. That
    crossover, N* = (cv_cell / cv_shared)^2, is the actionable number for
    the imaging protocol.

    ASSUMPTION stated rather than hidden: the shared term is taken as
    independent between the two channels. Illumination drift that hits both
    channels equally partly cancels in the ratio, so the true floor for a
    ratio-dominated estimate may be lower than this.
    """
    return cv_cell ** 2 / n_cells + cv_shared ** 2


def _red_trace(sbml_str, variant, dose, overrides=None):
    """Reporter_red on the same time grid as run(). Second channel of the
    decoder: the ratio alone discards the absolute level, which carries dose
    information -- the decoder inverts BOTH observables, so the bound must
    be computed on both."""
    I = ids(sbml_str)
    r = te.loadSBMLModel(sbml_str); r.reset()
    for k, v in (overrides or {}).items():
        r[I[k]] = v
    S = I[pname("S", variant)]
    r[S] = 0.0
    if P.PRE_T > 0:
        r.simulate(0, P.PRE_T, 200)
    r[S] = dose
    res = r.simulate(0, P.T_END, P.NPTS,
                     selections=["time", f"[{I['Reporter_red']}]"])
    return np.array(res[:, 1])


def sweep(sbml_str, variant, overrides=None):
    """Full dose grid plus the S=0 control. Returns trajectories and metrics."""
    doses = P.DOSES[variant]
    tr = {d: run(sbml_str, variant, d, overrides) for d in doses + (0,)}
    t0 = tr[doses[0]][0]
    M  = np.vstack([np.interp(t0, tr[d][0], tr[d][2]) for d in doses])
    s  = (M.max(0) - M.min(0)) / (np.abs(M.mean(0)) + 1e-12)
    k  = t0 > P.WARMUP
    g  = {d: at(tr[d][0], tr[d][1], P.PRIMARY_T) for d in doses + (0,)}

    out = {
        "variant": variant, "doses": doses, "traj": tr, "green_at_t": g,
        # TWO folds, because they diagnose different things:
        "fold_range":   g[doses[-1]] / g[doses[0]],   # does it discriminate doses?
        "fold_vs_zero": g[doses[-1]] / g[0],          # does it respond at all?
        "ratio_spread_t":  (t0[k], s[k]),
        "ratio_spread_at_t": float(np.interp(P.PRIMARY_T, t0[k], s[k])),
        "ratio_spread_max":  float(s[k].max()),
        "t_at_max_spread":   float(t0[k][s[k].argmax()]),
    }

    # --- dose-response fit on the WIDE grid ---------------------------------
    # The operating grid (P.DOSES) is what the device is READ on; it does not
    # sample either plateau, so fitting on it alone is unidentifiable. The
    # extra doses only pin baseline/top -- they are not claims about the cell.
    fd = P.FIT_DOSES[variant]
    gf = {d: at(*run(sbml_str, variant, d, overrides)[:2], P.PRIMARY_T)
          for d in fd}
    out["fit"] = hill_fit(fd, [gf[d] for d in fd])
    out["n_hill"] = out["fit"]["n_hill"]
    out["se_n_hill"] = out["fit"]["se_n_hill"]
    out["ec50"] = out["fit"]["ec50"]
    out["lin_r2"] = linearity_r2(fd, [gf[d] for d in fd], variant)
    out["fit_doses"] = fd
    out["green_fit_grid"] = gf
    return out


def assert_override_bites(sbml_str, name, value, label=""):
    """Raise if an override is a no-op.

    B2's confound control arm "FINAL_beta25" (er) printed numbers IDENTICAL
    to the uncontrolled arm, because the override value already equalled the
    model's value. A control that controls nothing is worse than no control:
    it looks like evidence. Nothing in the codebase caught it, so this does.
    """
    I = ids(sbml_str)
    if name not in I:
        raise KeyError(f"{label}: parameter {name} not in the model")
    cur = float(te.loadSBMLModel(sbml_str)[I[name]])
    if abs(cur - value) <= 1e-12 * max(1.0, abs(cur)):
        raise RuntimeError(
            f"{label}: override {name}={value:g} is a NO-OP (current value is "
            f"already {cur:g}). This arm is not a control -- remove it or fix "
            "the value, do not present it as one.")


def gate(res, theta=None):
    """Returns (viable, reason). The reason separates dead from saturated."""
    theta = P.THETA_DEFAULT if theta is None else theta
    fr, fz = res["fold_range"], res["fold_vs_zero"]
    if fr >= theta and fz >= theta:
        return True, "viable"
    if fz >= theta > fr:
        return False, "saturated"   # responds, but cannot discriminate doses
    if fr >= theta > fz:
        return False, "no-basal-contrast"
    return False, "dead"


# --------------------------------------------------------------------------
# LESIONS
# --------------------------------------------------------------------------
def _replace_token(formula, old, new):
    """Replace a WHOLE identifier (not a substring)."""
    out, buf = [], ""
    for ch in formula + " ":
        if ch.isalnum() or ch == "_":
            buf += ch
        else:
            out.append(new if buf == old else buf); out.append(ch); buf = ""
    return "".join(out)[:-1]


def ablate_tetr_feedback(sbml_str, variant):
    """Freeze P ONLY in TetR production (the autoregulation), leaving the
    reporter's P untouched (the signal path).

    TWO DESIGN CHOICES YOU MUST BE ABLE TO DEFEND:
    1. The reaction is located TOPOLOGICALLY (product = TetR_active AND the
       rate law contains P), not by name: names change on merge, topology
       does not. If two such reactions are found it raises rather than guess.
    2. The frozen value is P at ZERO stress, so both circuits start from the
       same baseline TetR level. Otherwise part of the measured difference
       would simply be a different amount of repressor, not the feedback.
    """
    I = ids(sbml_str)
    r = te.loadSBMLModel(sbml_str); r.reset()
    r[I[pname("S", variant)]] = 0.0
    r.simulate(0, P.PRE_T, 200)
    p_open = float(r[I["P"]])

    doc = libsbml.readSBMLFromString(sbml_str); m = doc.getModel()
    tetr, pid = I["TetR_active"], I["P"]
    target = None
    for rx in m.getListOfReactions():
        if tetr not in [x.getSpecies() for x in rx.getListOfProducts()]:
            continue
        f = libsbml.formulaToL3String(rx.getKineticLaw().getMath())
        if pid in f:
            if target is not None:
                raise RuntimeError("Two TetR-producing reactions contain P -- check the merge.")
            target = rx
    if target is None:
        raise RuntimeError("TetR autoregulation reaction not found.")

    prm = m.createParameter()
    prm.setId("P_open_control"); prm.setName("P_open_control")
    prm.setValue(p_open); prm.setConstant(True)
    f  = libsbml.formulaToL3String(target.getKineticLaw().getMath())
    nf = _replace_token(f, pid, "P_open_control")
    if nf == f:
        raise RuntimeError("Substitution of P failed.")
    target.getKineticLaw().setMath(libsbml.parseL3Formula(nf))
    return libsbml.writeSBMLToString(doc)


def lesion_set(sbml_str, variant):
    """(label, SBML transform | None, parameter overrides, note).

    Lesions that are meaningless for a given model are skipped automatically
    (e.g. removing adaptation from a static Hill sensor).
    """
    I = ids(sbml_str)
    m_green = te.loadSBMLModel(sbml_str)[I["m"]]
    L = [("full", None, {}, "reference"),
         ("no_TetR_feedback", ablate_tetr_feedback, {},
          "removes the linearisation; saturation expected")]
    r0 = te.loadSBMLModel(sbml_str)
    dx, koff = pname("d_x", variant), pname("k_off", variant)
    # A lesion whose target ALREADY has the lesioned value is not a null
    # result, it is a non-applicable row. d_x_er is already 0 (the ER sensor
    # is a perfect adaptor), so "perfect_adaptation" there produced an exact
    # 0.0% line that reads like a finding. Skip it instead.
    if dx in I and abs(r0[I[dx]] - 0.0) > 1e-12:
        L.append(("perfect_adaptation", None, {dx: 0.0},
                  "perfect adaptor -- adaptation-depth axis"))
    if koff in I and abs(r0[I[koff]] - 0.0) > 1e-12:
        L.append(("no_adaptation", None, {koff: 0.0},
                  "kills the sensor's own feedback"))
    # POSITIVE CONTROL -- red matures on the SAME TOTAL TIMESCALE as green.
    # The old version set m1 = m2 = 2*m, which gives a total mean maturation
    # time of 1/m1 + 1/m2 = 1/m only by coincidence of the factor 2, but also
    # silently rescales the red pool, so the +63% it produced in d_ratio_pct
    # was a units artefact, not an ablation effect. Stated explicitly: the
    # two-step red chain has mean time 1/m1 + 1/m2; matching it to the
    # single-step green (1/m) requires m1 = m2 = 2*m -- which is what this
    # is, now with the reason written down and the ratio reported relative
    # to its own t=0 value (see arms) so the scale change cannot masquerade
    # as a result.
    L.append(("no_timer_split", None, {"m1": 2 * m_green, "m2": 2 * m_green},
              "positive control: red total maturation time matched to green"))
    return L