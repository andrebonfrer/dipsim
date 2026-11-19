# findings.R ------------------------------------------------------------------
# Plain-language findings shown in the app's Notes tab. Update when the
# experiment in inst/scripts/ascm_dip_experiment.R is re-run.

#' Findings from the SC vs ridge-ASCM dip experiment
#'
#' From `inst/scripts/ascm_dip_experiment.R` (dipsim 0.4.0, budgetsc 0.6.2,
#' augMultiSynth 0.3.9): 80 simulated panels (60 adopters, 120 later adopters,
#' phi = 0.9, lean matching, per-unit eligibility, true effects zero), each
#' fitted with SC and with ASCM at three penalty settings, 10 replications per
#' cell of dip selection (0, 3) x baseline L (23, 70) x held-out weeks (0, 4).
#' Numbers are for total_spend unless stated (pre-adoption level about 1,000).
#' The cell-level summary is in
#' `system.file("extdata/ascm_dip_experiment_summary.csv", package = "dipsim")`.
#'
#' @return Character vector of paragraphs (shown in the app's Notes tab).
#' @export
dip_findings <- function() {
  c(
    paste("Setup: 80 simulated panels, each fitted with SC and with ridge ASCM at three penalty settings (tuned by the one-SE rule, fixed 1, fixed 0.1),",
          "so every comparison is on the same draw. 10 replications per cell of dip selection (0 or 3) x baseline L (23 or 70) x held-out weeks (0 or 4);",
          "persistence 0.9, lean matching, true effects zero. Figures are for total_spend (pre-adoption level about 1,000) unless stated.",
          "Reproduce with inst/scripts/ascm_dip_experiment.R."),
    paste("1. With the default tuning, augmentation never switches on. The one-SE rule chose no augmentation (penalty Inf) in all 80 panels, with or without",
          "a dip, so 'tuned ASCM' reproduced SC exactly. The in-time placebo forecasts found no gain from augmenting in this design."),
    paste("2. A fixed penalty shifts every estimate down, dip or no dip: by 25 to 75 with penalty 1 and by 52 to 112 with penalty 0.1, in every cell",
          "(standard errors 6 to 23). SC's bias is positive throughout (5 to 89), so the shift often lands closer to zero: at penalty 1 the mean absolute bias",
          "falls from 8 to 10% of the pre-adoption level to 7 to 8%, and for liquidity_deficit_rate from 17 to 22% to 14 to 19%. numarrears gets worse",
          "in all four cells (for example 13.7% to 16.2% with a dip and no held-out weeks; 11.38% to 11.44% in the closest one). The shift overshoots",
          "where SC is already close (L = 70 with held-out weeks and no dip: 5 becomes -32)."),
    paste("3. Augmentation does not remove the dip. With the spike held out (match_end = 4), the first-8-week bias with a dip is 214 and 208 for SC and",
          "still 171 and 174 with ASCM at penalty 1: the early bulge that marks a dip survives. With the spike inside the matching window, ASCM",
          "flattens the early bulge (first-8-week bias 63 and 152 for SC, -17 and 42 for ASCM) but turns the 40-week bias negative at L = 23 (21 becomes",
          "-45): it fits the spike and carries it into the counterfactual rather than learning that it unwinds."),
    paste("4. The augmented counterfactual leans heavily on extrapolation. Mean negative weight mass is 0.8 to 1.3 at penalty 1 and 2.3 to 2.9 at",
          "penalty 0.1 (largest single adopter 5.8 and 16). Each adopter has few eligible donors relative to the 5 outcomes x 23 to 70 weeks being",
          "matched, so the ridge correction extrapolates far beyond the donors. Check the augmentation diagnostics table before trusting an ASCM estimate."),
    paste("5. Remedies (with a dip). The reversion adjustment helps SC with or without held-out weeks: without them the bias goes from 21 to -8",
          "(L = 23) and from 66 to 24 (L = 70); with four held-out weeks from 64 to 27 and from 77 to 36. This differs from the earlier dipsim note,",
          "based on 3 replications and the pre-0.3.8 solver, which found the adjustment double-counted when the spike was inside the window.",
          "Adding ASCM at penalty 1 to held-out weeks plus reversion gave -20 (L = 23) and -1 (L = 70); the -1 is the smallest dip bias in the",
          "experiment at L = 70, while at L = 23 SC with reversion and no held-out weeks (-8) did better. Treat these as indications, not settled:",
          "10 replications, one outcome, and a penalty chosen by hand rather than by the data."),
    paste("Reading for the paper: report SC as the headline. ASCM with its default tuning is SC in this setting. A fixed-penalty ASCM is a robustness",
          "check whose value depends on extrapolation, so report its negative weight mass alongside it. The reversion adjustment is the remedy that",
          "helped consistently; whether to combine it with held-out weeks depends on L and deserves more replications."))
}
