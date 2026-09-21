# data/

## Required file: `delaunay2000_fig2b.csv`

Digitised time course of oxidised Yap1 from Delaunay et al. 2000.

**Mandatory provenance header** (as `#` comments at the top):

```
# source: Delaunay A, Isnard AD, Toledano MB. EMBO J 2000
# figure: 2B   (panel, axis units, dose)
# digitised_by: <name>   date: <YYYY-MM-DD>   tool: <WebPlotDigitizer x.y>
# columns: time_h, oxidised_fraction
```

**Why:** four separate source-mixing incidents were found in the POC. Every
parameter must trace back to the primary measurement, not to a secondary
description. Source hierarchy for the oxidative module:

| Source | Use |
|---|---|
| Delaunay 2000 | direct Yap1 oxidation state — calibration of `A_ox` |
| Dacquay & McMillen 2021 | accumulated reporter at 3 h — EC50 target |
| Gasch 2000 | TRX2 mRNA kinetics — **prediction checking only** |
| Kuge & Jones 1994 | **single endpoint only** (1 h, 1 mM) — consistency check |
