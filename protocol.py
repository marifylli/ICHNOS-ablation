"""
ICHNOS ablation — SINGLE SOURCE OF TRUTH FOR THE PROTOCOL.

No other file re-declares these values. Figure captions are generated from
here (caption()), so code and captions cannot drift apart.
"""

# --- doses -----------------------------------------------------------------
# ox: valid window 20-600 uM (above ~600 the signal collapses due to
#     toxicity, Dacquay Fig. 3A). K_act_ox = 208 uM.
# er: K_act_er ~ 2345, so the grid shifts accordingly.
DOSES = {"ox": (50, 100, 200, 300, 400, 600),
         "er": (500, 1000, 2000, 3000, 4000, 6000)}

# --- timing ----------------------------------------------------------------
T_END, NPTS = 12.0, 3000
READOUT     = (1.0, 3.0, 6.0)
PRIMARY_T   = 3.0          # Dacquay: accumulated reporter at 3 h

# --- initial conditions ----------------------------------------------------
# All species start at 0 in the merged SBML. That is not a cell, it is an
# empty vessel. Without pre-equilibration the R/G ratio during the first
# hours is dominated by the FILLING of the reporter pools rather than by the
# stress response -- same timescale, therefore not separable.
# Measured: without pre-eq the dose information is underestimated ~5x.
# Already converged by 20 h; 50 h is a safe margin.
PRE_T  = 50.0
WARMUP = 0.25   # the ratio is 0/0-ish at t~0; verified: no effect on results

# --- viability gate --------------------------------------------------------
# PROVISIONAL. Must be tied to the image-pipeline noise floor (cell-to-cell
# variability) -- do not leave it arbitrary in a presentation.
# Measured: theta = 3.0 rejects even the intact circuit on the current grid.
THETA_DEFAULT = 2.0
THETA_SWEEP   = (1.5, 2.0, 3.0)


def caption(variant):
    d = DOSES[variant]
    return (f"Pre-equilibration {PRE_T:.0f} h at S_{variant}=0 "
            f"(results identical for 20-200 h); doses {d[0]}-{d[-1]} uM, "
            f"n={len(d)}; readout t*={PRIMARY_T} h; spread computed for "
            f"t > {WARMUP} h.")
