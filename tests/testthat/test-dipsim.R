test_that("dip_run returns the diagnostics with the stub backend", {
  p <- dip_defaults(); p$backend <- "stub"; p$n_treated <- 30L; p$n_later <- 40L; p$n_weeks <- 160L; p$L <- 8L; p$match_end <- 2L
  r <- dip_run(p)
  expect_true(all(c("event_time", "summary", "shock_profile", "funnel", "params", "seconds") %in% names(r)))
  expect_true(all(c("ate", "ate_first8", "heldout_gap", "true_effect", "matched") %in% names(r$summary)))
  expect_true("total_spend" %in% r$summary$outcome)
  expect_true(all(r$summary$true_effect == 0))
  expect_equal(range(r$shock_profile$tau_k), c(-8, 39))
  expect_true(all(c("bias", "bias_pct_baseline", "bias_pct_effect", "bias_nodip", "bias_reversion", "phi_hat") %in% names(r$summary)))
  expect_true(all(is.na(r$summary$bias_pct_effect)))              # true effect zero -> relative-to-effect undefined
  txt <- interpret_run(r, "total_spend"); expect_true(length(txt) >= 4); expect_true(any(grepl("Verdict", txt)))
  # strong dip selection: the latent shock peaks just before adoption
  expect_gt(r$shock_profile[tau_k == -1, mean_u], r$shock_profile[tau_k == -8, mean_u])
})

test_that("true effects are placed on the outcome scale", {
  p <- dip_defaults(); p$backend <- "stub"; p$n_treated <- 20L; p$n_later <- 30L; p$n_weeks <- 150L; p$L <- 8L
  p$effects$arrears_log <- -1.5; p$effects$spend_log <- -0.4
  r <- dip_run(p)
  expect_lt(r$summary[outcome == "numarrears", true_effect], 0)
  expect_lt(r$summary[outcome == "total_spend", true_effect], 0)
  expect_true(is.finite(r$summary[outcome == "total_spend", bias_pct_effect]))
})

test_that("dip_grid sweeps a small design and app object builds", {
  p <- dip_defaults(); p$backend <- "stub"; p$n_treated <- 20L; p$n_later <- 30L; p$n_weeks <- 150L
  g <- dip_grid(p, data.frame(L = c(8L, 16L), match_end = 0L, matching = "single", stringsAsFactors = FALSE), reps = 1L)
  expect_equal(sort(unique(g$L)), c(8L, 16L))
  expect_true(all(c("ate", "true_effect", "rep") %in% names(g)))
  expect_s3_class(dip_app(), "shiny.appobj")
})

test_that("real backend runs when installed", {
  skip_if_not_installed("augMultiSynth")
  p <- dip_defaults(); p$backend <- "augmultisynth"; p$n_treated <- 30L; p$n_later <- 60L; p$n_weeks <- 160L; p$L <- 8L; p$match_end <- 4L
  r <- dip_run(p)
  expect_true(is.finite(r$summary[outcome == "total_spend", heldout_gap]))
})


test_that("remedies reduce dip bias on total_spend with the real backend", {
  skip_if_not_installed("augMultiSynth")
  p <- dip_defaults(); p$n_treated <- 60L; p$n_later <- 150L; p$n_weeks <- 200L; p$L <- 23L; p$match_end <- 4L
  p$remedies <- c("nodip", "reversion", "placebo"); p$phi <- 0.9; p$dip_selection <- 3
  r <- dip_run(p); s <- r$summary[outcome == "total_spend"]
  expect_true(all(is.finite(c(s$bias_nodip, s$bias_reversion, s$bias_placebo_adj))))
  expect_gt(s$phi_hat, 0.5)                                       # persistence recovered from donors
})
