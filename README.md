# dipsim

Interactive Ashenfelter-dip simulator for the `budgetsc` synthetic-control
stage, comparing plain synthetic control (SC) with ridge-augmented synthetic
control (ASCM).

```r
remotes::install_github("andrebonfrer/dipsim")   # needs budgetsc >= 0.6.2, augMultiSynth >= 0.3.9
dipsim::dip_app()
```

**Single run.** Set the data-generating process, the estimators and the design,
run, and read the diagnostics: the latent shock around adoption, the gap by
event week for each estimator (fitted on the same simulated panel), the bias
table (absolute, relative to the pre-adoption level, and relative to the true
effect), the augmentation diagnostics (penalty, negative weight mass,
pre-period fit), the remedies, and a generated interpretation.

**Compare designs.** Sweep L, match_end and the matching set; every cell fits
each selected estimator on the same draws. Post window fixed at 40 weeks
(adoption week + 39).

**ASCM penalty.** Tuned by in-time placebo forecasts (one-SE rule or
minimum), or fixed. With shock-triggered adoption the tuning's held-out weeks
contain the dip; see the Notes tab.

**Remedies** (observable-based, applied to each estimator's gaps): drop
shock-triggered adopters; reversion adjustment (shock persistence estimated
from donors); placebo bias estimate (donors with trigger-rule pseudo-onsets,
second fit).

**Experiment script.** `inst/scripts/ascm_dip_experiment.R` reproduces the
SC-vs-ASCM comparison summarised in the Notes tab (`dip_findings()`).
