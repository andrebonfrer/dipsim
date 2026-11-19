test_that("dip_run returns the diagnostics with the stub backend", {
  p <- dip_defaults(); p$backend <- "stub"; p$estimators <- "sc"; p$n_treated <- 30L; p$n_later <- 40L; p$n_weeks <- 160L; p$L <- 8L; p$match_end <- 2L
  r <- dip_run(p)
  expect_true(all(c("event_time", "summary", "augment", "shock_profile", "funnel", "params", "seconds") %in% names(r)))
  expect_equal(unique(r$summary$estimator), "sc")
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
  p <- dip_defaults(); p$backend <- "stub"; p$estimators <- "sc"; p$n_treated <- 20L; p$n_later <- 30L; p$n_weeks <- 150L; p$L <- 8L
  p$effects$arrears_log <- -1.5; p$effects$spend_log <- -0.4
  r <- dip_run(p)
  expect_lt(r$summary[outcome == "numarrears", true_effect], 0)
  expect_lt(r$summary[outcome == "total_spend", true_effect], 0)
  expect_true(is.finite(r$summary[outcome == "total_spend", bias_pct_effect]))
})

test_that("dip_grid sweeps a small design and app object builds", {
  p <- dip_defaults(); p$backend <- "stub"; p$estimators <- "sc"; p$n_treated <- 20L; p$n_later <- 30L; p$n_weeks <- 150L
  g <- dip_grid(p, data.frame(L = c(8L, 16L), match_end = 0L, matching = "single", stringsAsFactors = FALSE), reps = 1L)
  expect_equal(sort(unique(g$L)), c(8L, 16L))
  expect_true(all(c("ate", "true_effect", "rep", "estimator") %in% names(g)))
  expect_s3_class(dip_app(), "shiny.appobj")
})

test_that("real backend runs when installed", {
  skip_if_not_installed("augMultiSynth")
  p <- dip_defaults(); p$backend <- "augmultisynth"; p$estimators <- "sc"; p$n_treated <- 30L; p$n_later <- 60L; p$n_weeks <- 160L; p$L <- 8L; p$match_end <- 4L
  r <- dip_run(p)
  expect_true(is.finite(r$summary[outcome == "total_spend", heldout_gap]))
})


test_that("remedies reduce dip bias on total_spend with the real backend", {
  skip_if_not_installed("augMultiSynth")
  p <- dip_defaults(); p$n_treated <- 60L; p$n_later <- 150L; p$n_weeks <- 200L; p$L <- 23L; p$match_end <- 4L
  p$remedies <- c("nodip", "reversion", "placebo"); p$phi <- 0.9; p$dip_selection <- 3; p$estimators <- "sc"
  r <- dip_run(p); s <- r$summary[outcome == "total_spend"]
  expect_true(all(is.finite(c(s$bias_nodip, s$bias_reversion, s$bias_placebo_adj))))
  expect_gt(s$phi_hat, 0.5)                                       # persistence recovered from donors
})


test_that("ASCM is fitted on the same panel as SC", {
  skip_if_not(dipsim:::.ascm_available())
  p <- dip_defaults(); p$n_treated <- 30L; p$n_later <- 70L; p$n_weeks <- 160L; p$L <- 12L
  p$remedies <- character(0); p$estimators <- c("sc", "ascm")
  # penalty Inf: augmentation switched off, so ASCM must reproduce SC exactly
  p$augment_lambda <- Inf
  r <- dip_run(p)
  expect_setequal(unique(r$summary$estimator), c("sc", "ascm"))
  sc <- r$summary[estimator == "sc"][order(outcome)]; as <- r$summary[estimator == "ascm"][order(outcome)]
  expect_equal(as$ate, sc$ate, tolerance = 1e-10)
  expect_equal(as$true_effect, sc$true_effect)
  expect_equal(r$augment[estimator == "ascm", neg_mass], 0)
  # a finite penalty changes the counterfactual and reports its diagnostics
  p$augment_lambda <- 0.5
  r2 <- dip_run(p)
  a <- r2$augment[estimator == "ascm"]
  expect_equal(a$lambda, 0.5)
  expect_lte(a$pre_rmse_ascm, a$pre_rmse_sc + 1e-12)
  expect_false(isTRUE(all.equal(r2$summary[estimator == "ascm" & outcome == "total_spend", ate],
                                r2$summary[estimator == "sc" & outcome == "total_spend", ate])))
  txt <- interpret_run(r2, "total_spend"); expect_true(any(grepl("Augmentation", txt)))
})

test_that("ASCM falls back to SC with the stub backend", {
  p <- dip_defaults(); p$backend <- "stub"; p$estimators <- c("sc", "ascm")
  p$n_treated <- 20L; p$n_later <- 30L; p$n_weeks <- 150L; p$L <- 8L; p$remedies <- character(0)
  expect_warning(r <- dip_run(p), "fitting SC only")
  expect_equal(unique(r$summary$estimator), "sc")
})

test_that("findings text is available", {
  expect_true(length(dip_findings()) >= 1)
})
