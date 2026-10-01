# dipsim

Interactive Ashenfelter-dip simulator for the `budgetsc` synthetic-control stage.

```r
install.packages("dipsim_0.2.0.tar.gz", repos = NULL, type = "source")   # needs budgetsc >= 0.4.2
library(Matrix)                                                           # until augMultiSynth 0.3.6
dipsim::dip_app()
```

Single run: set the data-generating process and the design, run, read the four
diagnostics, the bias table (absolute, relative to the pre-adoption level, and
relative to the true effect), the remedies, and a generated interpretation.
Compare designs: sweep L, match_end and the matching set. Post window fixed at
40 weeks (adoption week + 39).

Remedies (all observable-based): drop shock-triggered adopters; reversion
adjustment (shock persistence estimated from donors); placebo bias estimate
(donors with trigger-rule pseudo-onsets, second SC fit). See the Notes tab for
which of them work and when.
