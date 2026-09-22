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

# Grid for the DOSE-RESPONSE FIT (in silico only). Must span BOTH plateaus,
# otherwise baseline/top/EC50/n_H are not jointly identifiable. Geometric,
# 14 points, 3 decades, centred on K_act. Doses above the operating window
# are NOT claims about the cell (toxicity) -- they only pin the top plateau.
FIT_DOSES = {"ox": (5, 8.5, 14.5, 24.8, 42.3, 72.2, 123, 211, 360, 614,
                    1048, 1789, 3054, 5000),
             "er": (50, 85, 145, 248, 423, 722, 1233, 2105, 3594, 6135,
                    10474, 17882, 30529, 50000)}

# Where the device is actually USED. Linearity is judged only here.
# ox: 20-600 uM (Dacquay Fig. 3A toxicity). er: PROVISIONAL = DOSES range,
# no documented upper bound yet -- replace when one exists.
OPERATING_WINDOW = {"ox": (20, 600), "er": (500, 6000)}

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


# --- measurement noise (decoding bound) -------------------------------------
# PROVISIONAL. Coefficient of variation of each fluorescence channel per cell.
# Must be replaced by the image-pipeline estimate (cell-to-cell CV of the
# single-fluorophore control strain). Until then every decoding-bound result
# is reported for the whole sweep, never for one value.
NOISE_CV       = {"green": 0.25, "red": 0.25}
NOISE_CV_SWEEP = (0.10, 0.25, 0.40)
FD_REL_STEP    = 0.05   # relative step in ln(dose) for the dose derivative

# --- population decoding ----------------------------------------------------
# The decoder reads CULTURES (confirmed design target), so the measurement is
# a mean over N cells. See ichnos_ablation.population_variance().
N_CELLS_SWEEP   = (1, 10, 30, 100, 300, 1000)
# PROVISIONAL, like NOISE_CV: the part of the noise common to every cell in a
# field (illumination, exposure, flat-field, background, crosstalk
# coefficient). Swept, never fixed, until the image pipeline measures it.
SHARED_CV_SWEEP = (0.0, 0.02, 0.05, 0.10)
# Target timing resolution used to turn the sweep into a protocol number
# ("image at least N cells"). 0.5 h = one eighth of the 0-4 h memory window.
TARGET_SIGMA_T_H = 0.5
# Cells averaged per measurement in the headline population numbers.
N_REPORT = 100
# Ages of the stress event at which the decoder is characterised. NOT a
# choice of when to image: in use, the age is the UNKNOWN being estimated,
# so this grid is the range over which the device is characterised.
AGE_GRID = (0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0)


def caption(variant):
    d = DOSES[variant]
    return (f"Pre-equilibration {PRE_T:.0f} h at S_{variant}=0 "
            f"(results identical for 20-200 h); doses {d[0]}-{d[-1]} uM, "
            f"n={len(d)}; readout t*={PRIMARY_T} h; spread computed for "
            f"t > {WARMUP} h.")