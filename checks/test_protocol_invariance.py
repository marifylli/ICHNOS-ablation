"""
CHECK 2 -- the conclusions do not hinge on arbitrary protocol choices.

PRE_T  : already converged at 20 h. WITHOUT pre-equilibration the dose
         information is underestimated ~5x (measured: 2.5% instead of 13%) --
         i.e. you would conclude the ratio is nearly dose-invariant, which is
         the WRONG conclusion and leads to a wrong decoder design. This check
         documents that.
WARMUP : no effect for 0-0.5 h (maximum spread occurs around 0.7 h).
n_eff  : must be invariant under both.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import protocol as P, ichnos_ablation as A
from models import build_snapshots as S

TOL = 0.02   # relative tolerance on n_eff


def _spread3(sbml, variant):
    res = A.sweep(sbml, variant)
    return res["ratio_spread_at_t"], res["n_eff"]


def test_pre_t_converged(variant="ox", tag="PULSE"):
    S.verify(); sbml = S.path(variant, tag).read_text()
    old = P.PRE_T; got = {}
    try:
        for pre in (0.0, 20.0, 50.0, 200.0):
            P.PRE_T = pre; got[pre] = _spread3(sbml, variant)
    finally:
        P.PRE_T = old
    for pre in (20.0, 50.0, 200.0):
        assert abs(got[pre][0] - got[50.0][0]) < 1e-3, f"PRE_T={pre} has not converged"
    assert got[0.0][0] < 0.5 * got[50.0][0], (
        "Expected the cold start to UNDERESTIMATE the spread; it does not -- "
        "check the merged model's initial conditions.")
    for v in got.values():
        assert abs(v[1] - got[50.0][1]) / got[50.0][1] < TOL, "n_eff is sensitive to PRE_T"
    return got


def test_warmup_irrelevant(variant="ox", tag="PULSE"):
    S.verify(); sbml = S.path(variant, tag).read_text()
    old = P.WARMUP; got = {}
    try:
        for w in (0.0, 0.1, 0.25, 0.5):
            P.WARMUP = w; got[w] = _spread3(sbml, variant)
    finally:
        P.WARMUP = old
    for w, v in got.items():
        assert abs(v[0] - got[0.25][0]) < 1e-6, f"WARMUP={w} changes the result"
    return got


if __name__ == "__main__":
    g = test_pre_t_converged()
    print("PRE_T  :", {k: (round(v[0], 4), round(v[1], 3)) for k, v in g.items()})
    g = test_warmup_irrelevant()
    print("WARMUP :", {k: (round(v[0], 4), round(v[1], 3)) for k, v in g.items()})
    print("ok: conclusions are invariant to the protocol choices")
