"""
ARM B1 -- which sensor STRUCTURE describes the DATA?

THE ONLY ARM WITH A LIKELIHOOD, THEREFORE THE ONLY ONE WITH AICc.
aicc() is defined HERE and not in the shared core, so it cannot be called by
accident from A or B2, where it is meaningless (no fitting happens there).

STATUS: SKELETON. Requires data/delaunay2000_fig2b.csv (digitised).
The script HALTS cleanly if it is missing -- it never produces a fake result.

RULES THAT MUST BE FOLLOWED (otherwise the analysis invites criticism)
  * AICc, NOT AIC: required when n/k < 40, which holds with few points.
  * Akaike weights, NOT declaring a winner. dAIC < 2 = indistinguishable,
    dAIC > 10 = strong preference.
  * SAME data and SAME error model for both structures -- otherwise AICc is
    comparing error models, not biology.
  * A static Hill CANNOT structurally produce a transient response. If it
    loses, report it as STRUCTURAL INADEQUACY, not as a statistical victory:
    the outcome is close to predetermined.
  * k counts ONLY the free parameters of this fit (+1 for sigma).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np

DATA = Path(__file__).resolve().parents[1] / "data" / "delaunay2000_fig2b.csv"


def aicc(rss, n, k):
    """k = free parameters OF THE FIT, including sigma."""
    if n - k - 1 <= 0:
        return np.inf
    return n * np.log(rss / n) + 2 * k + (2 * k * (k + 1)) / (n - k - 1)


def akaike_weights(aiccs):
    a = np.asarray(aiccs, float)
    w = np.exp(-(a - a.min()) / 2)
    return w / w.sum()


def report(fits):
    """fits: [(name, rss, n, k)] -> table with dAICc and Akaike weights."""
    vals = [aicc(rss, n, k) for _, rss, n, k in fits]
    w = akaike_weights(vals)
    d = np.array(vals) - min(vals)
    print(f"{'model':>22} {'k':>3} {'AICc':>10} {'dAICc':>8} {'weight':>8}")
    for (name, _, _, k), a, dd, ww in zip(fits, vals, d, w):
        print(f"{name:>22} {k:3d} {a:10.2f} {dd:8.2f} {ww:8.3f}")
    if d.max() < 2:
        print("\n  dAICc < 2 -> the models are NOT distinguishable from this data.")
    elif d.max() > 10:
        print("\n  dAICc > 10 -> strong preference. Check whether it is driven by "
              "structural inadequacy (e.g. a static Hill cannot produce a "
              "transient) and report it AS SUCH.")
    return dict(zip([f[0] for f in fits], zip(vals, d, w)))


def main():
    if not DATA.exists():
        raise SystemExit(
            f"Missing {DATA}.\n"
            "Requires a digitisation of Delaunay 2000 Fig. 2B/2C with a\n"
            "provenance header: which figure, which edition, who digitised it,\n"
            "with which tool. (Source-mixing has already occurred 4 times in\n"
            "the POC.)\n"
            "Then: fit A_ox(t) for static vs adaptive using the SAME error\n"
            "model, and call report([('static_hill', rss, n, k), ('adaptive', ...)]).")
    raise SystemExit("TODO: fit routine -- see the rules in the docstring.")


if __name__ == "__main__":
    main()
