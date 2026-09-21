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
  ratio_spread is MINIMISED by a dead circuit, and n_eff stays well defined
  on top of negligible signal (measured: at Kd = 5000 nM, spread = 0.0004
  while fold = 1.003). So shape metrics without a gate reward the dead
  circuit. The gate is therefore BINARY and COMES FIRST.

  TWO FAILURE MODES -- do not conflate them in the presentation:
    dead      : fold_range AND fold_vs_zero both low -> no response at all
    saturated : fold_vs_zero HIGH, fold_range low -> responds strongly but
                cannot discriminate doses (the linear range is gone)
  The second is a RESULT, not a rejection.
"""
import numpy as np, libsbml, tellurium as te
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
    return float(np.interp(tq, t, y))


def hill_neff(doses, resp):
    """n_eff = slope of logit(response) vs log10(dose).

    This is a SHAPE property: it stays well defined on top of negligible
    signal. That is exactly why it is never reported without the gate.
    """
    d = np.log10(np.asarray(doses, float))
    y = np.asarray(resp, float)
    if y.max() - y.min() <= 0:
        return np.nan
    frac = (y - y.min()) / (y.max() - y.min())
    keep = (frac > .05) & (frac < .95)
    if keep.sum() < 3:
        keep = np.ones_like(frac, bool)
    p = np.clip(frac[keep], 1e-6, 1 - 1e-6)
    return float(np.polyfit(d[keep], np.log10(p / (1 - p)), 1)[0])


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
        "n_eff": hill_neff(doses, [g[d] for d in doses]),
    }
    out["n_eff_by_t"] = {tq: hill_neff(doses, [at(tr[d][0], tr[d][1], tq) for d in doses])
                         for tq in P.READOUT}
    return out


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
       same baseline TetR level. Otherwise part of the n_eff difference would
       simply be a different amount of repressor, not the feedback itself.
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
    dx, koff = pname("d_x", variant), pname("k_off", variant)
    if dx in I:
        L.append(("perfect_adaptation", None, {dx: 0.0},
                  "perfect adaptor -- adaptation-depth axis"))
    if koff in I:
        L.append(("no_adaptation", None, {koff: 0.0},
                  "kills the sensor's own feedback"))
    L.append(("no_timer_split", None, {"m1": 2 * m_green, "m2": 2 * m_green},
              "positive control: red matures as fast as green"))
    return L
