# dipsim 0.4.2

- True effects are now computed exactly for the units and window the
  estimator uses. `dip_run()` builds a counterfactual twin of the simulated
  panel (treated post-adoption spend divided by `exp(spend_log)`, arrears
  count by `exp(arrears_log)`, arrears length reduced by `arrears_length`),
  passes both panels through budgetsc's `load_panel()` and
  `derive_outcomes()`, and takes the per-adopter post-minus-pre difference,
  the same contrast as the gap estimator. Previously the spend effect was
  scaled by the mean spend of all untreated rows (donors and treated
  pre-periods pooled), and derived outcomes such as
  `liquidity_deficit_rate` were given a true effect of zero even when a
  spend effect changed them.
- Bias columns are computed after all estimators are fitted
  (internal `.add_bias()`); results are otherwise unchanged.

# dipsim 0.4.1

- App: after each run the gap plot and interpretation switch to the first
  outcome with a non-zero true effect (previously they stayed on
  `numarrears`, so a spend effect looked as if it had not been applied).
  The plot title states the true effect, the true-effect line is drawn over
  the post window, kept inside the y-range and added to the legend.
- App: the spend effect input is labelled as a log change (-0.3 is about -26%).

# dipsim 0.4.0

- Ridge-augmented synthetic control (ASCM). `dip_run()` fits every estimator in
  `p$estimators` (`"sc"`, `"ascm"`) on the same simulated panel, through
  budgetsc's `sc.augment = "ridge"` (budgetsc >= 0.6.2, augMultiSynth >= 0.3.9),
  so the comparison is paired. Results gain an `estimator` column and an
  `augment` table (penalty, mean and max negative weight mass, pre-period RMSE
  with SC and ASCM weights).
- ASCM penalty: tuned by in-time placebo forecasts with the one-SE rule
  (default) or the minimum, or fixed (`augment_lambda`; `Inf` = no
  augmentation).
- App: estimator and penalty controls; gap plot overlays the estimators;
  augmentation diagnostics table; interpretation compares SC and ASCM; the
  design grid draws SC solid and ASCM dashed; Notes explain augmentation and
  report the findings of `inst/scripts/ascm_dip_experiment.R` (`dip_findings()`).
- Remedies are applied to each estimator's gaps.
- `dip_grid()` returns one row per cell x replication x estimator x outcome,
  with the augmentation diagnostics; `augment_lambda` can be a grid factor.
- Requires budgetsc >= 0.6.2.
- Findings (`dip_findings()`, 80 paired panels): the one-SE tuning never augments in this design (tuned ASCM equals
  SC); a fixed penalty shifts every estimate down by 25 to 112 on total_spend with heavy extrapolation (negative
  weight mass 0.8 to 2.9); augmentation does not remove the dip; the reversion adjustment helps SC with or without
  held-out weeks, which revises the earlier note that it double-counts with the spike inside the window.
