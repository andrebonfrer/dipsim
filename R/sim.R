# sim.R -----------------------------------------------------------------------
# Simulation engine. One `dip_run()` = one simulated panel, fitted with each
# requested estimator (plain synthetic control and/or ridge-augmented SC) on
# the SAME panel, plus the statistics that reveal a dip: the event-time gap
# profile, the 40-week and first-8-week ATEs, the held-out gap (when
# match_end > 0) and the latent shock profile around adoption. Fitting both
# estimators on one draw makes the comparison paired: differences come from
# the estimator, not the simulation. `dip_grid()` repeats it over design
# choices. Everything runs through budgetsc's own functions, so what the app
# shows is what the pipeline does.
#' @importFrom data.table data.table := .N .SD fread rbindlist
NULL

#' Default simulation parameters
#'
#' `estimators` chooses what is fitted on each simulated panel: `"sc"` (simplex
#' synthetic control) and/or `"ascm"` (ridge-augmented synthetic control,
#' budgetsc spec `sc.augment = "ridge"`, augMultiSynth >= 0.3.9). For ASCM,
#' `augment_lambda = NULL` tunes the relative ridge penalty by in-time placebo
#' forecasts on the fit's own donor pools (`augMultiSynth::tune_augment()`),
#' with `augment_one_se` choosing the one-standard-error rule (default) or the
#' minimum; a number fixes the penalty (`Inf` = no augmentation).
#' @return Named list.
#' @export
dip_defaults <- function() list(
  n_treated = 60L, n_later = 120L, n_weeks = 200L, phi = 0.9, sigma_u = 0.25,
  dip_selection = 3, type_selection = 1, covid = TRUE,
  effects = list(spend_log = 0, arrears_log = 0, arrears_length = 0, share_FastFood = 0, share_Alcohol = 0, share_Gambling = 0),
  L = 23L, H = 39L, match_end = 0L, matching = "lean", donor_eligibility = "per_unit", custom_outcomes = NULL, max_donors = 1000L,
  backend = if (requireNamespace("augMultiSynth", quietly = TRUE)) "augmultisynth" else "stub",
  estimators = if (.ascm_available()) c("sc", "ascm") else "sc",
  augment_lambda = NULL, augment_one_se = TRUE, augment_tune_units = 40L,
  outcomes = c("numarrears", "total_spend", "liquidity_deficit_rate"), seed = 1L,
  remedies = c("nodip", "reversion"), z_cut = 1, k_pre = 4L)

# ridge augmentation needs augMultiSynth >= 0.3.9 and budgetsc >= 0.6.2
.ascm_available <- function() {
  requireNamespace("augMultiSynth", quietly = TRUE) &&
    utils::packageVersion("augMultiSynth") >= "0.3.9" &&
    utils::packageVersion("budgetsc") >= "0.6.2"
}

.estimator_label <- function(e) unname(c(sc = "SC", ascm = "ridge ASCM")[e])

.matching_set <- function(matching, custom = NULL) switch(matching,
  single = "total_spend", lean = budgetsc::match_outcomes_lean(), wide = budgetsc::match_outcomes_wide(),
  custom = unique(c(custom, "total_spend")), stop("unknown matching set"))

# budgetsc spec overrides for one estimator
.estimator_spec <- function(p, estimator, mo) {
  s <- budgetsc::spec_main("sample.n_lags" = as.integer(p$L), "sample.n_leads" = as.integer(p$H),
                           "sc.match_outcomes" = mo, "sc.match_end" = as.integer(p$match_end),
                           "sc.max_donors" = as.integer(p$max_donors), "sc.backend" = p$backend,
                           "sc.donor_eligibility" = p$donor_eligibility,
                           "sc.parallel" = FALSE, "post.outcomes" = p$outcomes,
                           "sc.holdout_outcomes" = c("weekly_income", "monthly_loan_payment_w"))
  if (identical(estimator, "ascm")) {
    s$sc$augment <- "ridge"
    s$sc$augment_tune <- list(one_se = isTRUE(p$augment_one_se), tune_units = as.integer(p$augment_tune_units))
    if (!is.null(p$augment_lambda) && length(p$augment_lambda) == 1L && !is.na(p$augment_lambda)) {
      s$sc$augment_lambda <- as.numeric(p$augment_lambda)
    }
  }
  s$name <- paste0("dipsim_", estimator, "_", format(Sys.time(), "%H%M%S"), "_", sample.int(1e6, 1))
  s
}

# Fit one estimator on the saved panel and collect its statistics.
.fit_estimator <- function(p, estimator, mo, root, panel) {
  est <- estimator
  spec <- .estimator_spec(p, est, mo)
  run <- budgetsc::run_register(spec, root)
  on.exit(unlink(run$dir, recursive = TRUE), add = TRUE)
  budgetsc::run_sample(run, root); budgetsc::fit_sc(run); d <- budgetsc::sc_diagnostics(run)

  rem <- lapply(p$outcomes, function(o) {
    ug <- unit_gaps(run, o); out <- list(outcome = o, pre_mean = ug$pre_mean)
    if ("nodip" %in% p$remedies) { r <- remedy_nodip(ug, p$k_pre, p$z_cut); out$ate_nodip <- r$ate; out$n_dropped_nodip <- r$n_dropped }
    if ("reversion" %in% p$remedies) { r <- remedy_reversion(ug); out$ate_reversion <- r$ate; out$phi_hat <- r$phi_hat }
    if ("placebo" %in% p$remedies) { r <- tryCatch(remedy_placebo(run, o, k = p$dip_selection, seed = p$seed), error = function(e) NULL)
      if (!is.null(r)) out$pseudo_ate <- r$pseudo_ate }
    data.table::as.data.table(out) })
  rem <- rbindlist(rem, fill = TRUE)
  et <- data.table::copy(d$pretrend)
  ate <- et[window == "post", list(ate = mean(mean_gap)), by = outcome]
  first8 <- et[window == "post" & tau_k < 8, list(ate_first8 = mean(mean_gap)), by = outcome]
  held <- d$fit_summary[window == "heldout_pre", list(outcome, heldout_gap = mean_gap)]
  summ <- merge(merge(ate, first8, by = "outcome"), held, by = "outcome", all.x = TRUE)
  summ[, matched := outcome %in% mo]
  summ <- merge(summ, rem, by = "outcome", all.x = TRUE)
  if ("pseudo_ate" %in% names(summ)) summ[, ate_placebo_adj := ate - pseudo_ate]

  # augmentation diagnostics (penalty, extrapolation, pre-period fit before/after)
  scfit <- readRDS(file.path(run$dir, "sc_fit.rds"))
  aug <- scfit$fit$augment
  augd <- data.table::data.table(estimator = est, lambda = NA_real_, neg_mass = 0, max_neg_mass = 0,
                                 pre_rmse_sc = NA_real_, pre_rmse_ascm = NA_real_)
  if (!is.null(aug)) {
    dg <- aug$diagnostics
    if (is.data.frame(dg)) {
      augd <- data.table::data.table(estimator = est, lambda = as.numeric(aug$lambda)[1],
                                     neg_mass = mean(dg$neg_mass), max_neg_mass = max(dg$neg_mass),
                                     pre_rmse_sc = mean(dg$pre_rmse_scm), pre_rmse_ascm = mean(dg$pre_rmse_aug))
    } else {
      augd$lambda <- as.numeric(aug$lambda)[1]
    }
  }
  tr <- readRDS(file.path(run$dir, "sample.rds"))$ids[role == "treated"]
  funnel <- fread(file.path(run$dir, "tables", "sample_funnel.csv"))
  et[, estimator := est]; summ[, estimator := est]
  list(event_time = et, summary = summ, augment = augd, treated = tr, funnel = funnel)
}

#' One simulated panel, fitted with each requested estimator
#'
#' @param p Parameter list as from [dip_defaults()] (missing entries filled).
#' @param root Scratch directory for the budgetsc runs (deleted afterwards).
#' @param keep Ignored (kept for compatibility); run directories are removed.
#' @return A list: `event_time` (gap by event week, matched and holdout
#'   outcomes, one block per estimator), `summary` (per estimator and outcome:
#'   ate, ate_first8, heldout_gap, true effect, bias, remedies), `augment`
#'   (per estimator: ridge penalty, mean and max negative weight mass,
#'   pre-period RMSE before and after augmentation), `shock_profile` (mean
#'   latent shock u by event week, treated units), `funnel`, `params`,
#'   `seconds`.
#' @export
dip_run <- function(p = dip_defaults(), root = file.path(tempdir(), "dipsim"), keep = FALSE) {
  p <- utils::modifyList(dip_defaults(), p, keep.null = TRUE)
  est <- unique(intersect(p$estimators, c("sc", "ascm")))
  if (!length(est)) stop("estimators must include 'sc' and/or 'ascm'.", call. = FALSE)
  if ("ascm" %in% est && (!identical(p$backend, "augmultisynth") || !.ascm_available())) {
    warning("ridge ASCM needs the augmultisynth backend, augMultiSynth >= 0.3.9 and budgetsc >= 0.6.2; fitting SC only.",
            call. = FALSE)
    est <- setdiff(est, "ascm"); if (!length(est)) est <- "sc"
  }
  p$estimators <- est
  t0 <- Sys.time()
  dir.create(file.path(root, "Processed"), recursive = TRUE, showWarnings = FALSE)
  old <- options(budgetsc.root = root); on.exit(options(old), add = TRUE)
  sim <- budgetsc::sim_panel(n_treated = p$n_treated, n_later = p$n_later, n_never = 0L, n_weeks = p$n_weeks,
                             phi = p$phi, sigma_u = p$sigma_u, dip_selection = p$dip_selection,
                             type_selection = p$type_selection, effects = p$effects, covid = p$covid,
                             keep_latent = TRUE, seed = p$seed)
  panel <- sim$panel
  saveRDS(panel[, !c("u", "lat", "D")], file.path(root, "Processed", "analysis_panel.rds"))
  mo <- .matching_set(p$matching, p$custom_outcomes)

  fits <- lapply(est, function(e) .fit_estimator(p, e, mo, root, panel))
  names(fits) <- est
  tr <- fits[[1]]$treated
  te <- .true_effect(p, panel, tr, root, mo)
  for (k in names(fits)) fits[[k]]$summary <- .add_bias(fits[[k]]$summary, te)
  sp <- panel[tr, on = "customer_id"][, list(mean_u = mean(u)), by = list(tau_k = wID - treat_wID)][tau_k >= -p$L & tau_k <= p$H][order(tau_k)]
  list(event_time = rbindlist(lapply(fits, `[[`, "event_time"), fill = TRUE),
       summary = rbindlist(lapply(fits, `[[`, "summary"), fill = TRUE),
       augment = rbindlist(lapply(fits, `[[`, "augment"), fill = TRUE),
       shock_profile = sp, funnel = fits[[1]]$funnel, params = p,
       seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

# True effect per outcome: the effect on the treated units the estimator uses,
# averaged over the same post window and net of the same pre-window intercept
# as the gap estimator. Built from a counterfactual twin of the simulated
# panel: the treated rows' spend columns divided by exp(spend_log) (spend is
# multiplicative in sim_panel, so this is Y(0) up to rounding to cents), the
# arrears count divided by exp(arrears_log) (its expectation under no
# effect) and the arrears length reduced by arrears_length. Both panels then
# go through budgetsc's own load_panel() and derive_outcomes(), so winsorising,
# rolling rates and liquidity indicators are applied exactly as for the fit.
.true_effect <- function(p, panel, tr, root, mo) {
  out <- stats::setNames(numeric(length(p$outcomes)), p$outcomes)
  eff <- p$effects
  if (all(unlist(eff[c("spend_log", "arrears_log", "arrears_length")]) == 0)) return(out)
  y0 <- data.table::copy(panel)
  sc <- grep("^(spend_|irs_spend_)", names(y0), value = TRUE)
  if (eff$spend_log != 0) y0[D == 1L, (sc) := lapply(.SD, function(x) x * exp(-eff$spend_log)), .SDcols = sc]
  if (eff$arrears_log != 0) { y0[, numarrears := as.numeric(numarrears)]; y0[D == 1L, numarrears := numarrears * exp(-eff$arrears_log)] }
  if (eff$arrears_length != 0) y0[D == 1L & numarrears > 0, lengtharrears := pmax(7, lengtharrears - eff$arrears_length)]
  f0 <- file.path(root, "Processed", "analysis_panel_y0.rds")
  saveRDS(y0[, !c("u", "lat", "D")], f0); on.exit(unlink(f0), add = TRUE)
  spec <- .estimator_spec(p, "sc", mo)
  der <- function(file) {
    s <- spec; s$data$panel_file <- file
    x <- budgetsc::load_panel(s, root)
    budgetsc::derive_outcomes(x, focal = s$derive$focal_spend, winsor_probs = s$derive$winsor_probs,
                              roll_n = s$derive$roll_n, cv_basis = s$derive$cv_basis)
    x[tr[, list(customer_id, treat_wID)], on = "customer_id", nomatch = NULL]
  }
  a <- der(file.path(root, "Processed", "analysis_panel.rds")); b <- der(f0)
  data.table::setkey(a, customer_id, wID); data.table::setkey(b, customer_id, wID)
  tau <- a$wID - a$treat_wID
  post <- tau >= 0 & tau <= p$H
  pre <- tau >= -p$L & tau <= -(p$match_end + 1L)
  for (o in intersect(p$outcomes, names(a))) {          # per adopter, then averaged, as the gap estimator does
    d <- data.table::data.table(id = a$customer_id, d = a[[o]] - b[[o]], post, pre)
    u <- d[, list(v = mean(d[post]) - mean(d[pre])), by = id]
    out[[o]] <- mean(u$v, na.rm = TRUE)
  }
  out
}

# Bias of each estimate (and remedy) against the true effect.
.add_bias <- function(summ, te) {
  summ[, true_effect := unname(te[outcome])]
  summ[is.na(true_effect), true_effect := 0]
  summ[, bias := ate - true_effect]
  summ[, bias_first8 := ate_first8 - true_effect]
  summ[, bias_pct_baseline := 100 * bias / abs(pre_mean)]
  summ[, bias_pct_effect := data.table::fifelse(true_effect != 0, 100 * bias / abs(true_effect), NA_real_)]
  for (v in intersect(c("ate_nodip", "ate_reversion", "ate_placebo_adj"), names(summ)))
    summ[, (sub("^ate_", "bias_", v)) := get(v) - true_effect]
  summ
}

#' Repeat `dip_run()` over a design grid
#'
#' Each cell simulates one panel per replication and fits every estimator in
#' `p$estimators` on it, so SC and ridge ASCM are compared on the same draws.
#' @param p Base parameters.
#' @param grid data.frame of columns to vary (any of L, match_end, matching,
#'   n_later, phi, dip_selection, max_donors, augment_lambda).
#' @param reps Replications per cell.
#' @param progress Optional function(i, n) called after each run.
#' @return data.table: one row per cell x rep x estimator x outcome with the
#'   summary statistics, the augmentation diagnostics and the cell's factors.
#' @export
dip_grid <- function(p = dip_defaults(), grid, reps = 1L, progress = NULL) {
  grid <- data.table::as.data.table(grid); out <- list(); k <- 0L; n <- nrow(grid) * reps
  for (i in seq_len(nrow(grid))) for (r in seq_len(reps)) {
    pp <- utils::modifyList(p, as.list(grid[i]), keep.null = TRUE); pp$seed <- p$seed + 1000L * r + i
    res <- tryCatch(dip_run(pp), error = function(e) { message("dip_grid: cell ", i, " rep ", r, " failed: ", conditionMessage(e)); NULL })
    k <- k + 1L; if (!is.null(progress)) progress(k, n)
    if (is.null(res)) next
    s <- merge(res$summary, res$augment, by = "estimator", all.x = TRUE)
    out[[length(out) + 1L]] <- cbind(grid[i], rep = r, s)
  }
  rbindlist(out, fill = TRUE)
}

utils::globalVariables(c("window", "mean_gap", "tau_k", "outcome", "true_effect", "matched", "role", "pseudo_ate",
                         "ate_placebo_adj", "bias", "pre_mean", "bias_pct_baseline", "bias_pct_effect", "ate",
                         "treat_wID", "u", "wID", "D", "numarrears", "customer_id", "estimator", "bias_first8",
                         "ate_first8", "lambda", "neg_mass", "lengtharrears", "id", "v"))


#' Plain-language interpretation of one run
#' @param r Result of [dip_run()].
#' @param outcome Outcome to describe. Default `"total_spend"`.
#' @return Character vector of paragraphs.
#' @export
interpret_run <- function(r, outcome = "total_spend") {
  p <- r$params; oc <- outcome
  S <- r$summary[outcome == oc]
  if (!nrow(S) || all(is.na(S$ate))) return("No result for this outcome.")
  dgp <- sprintf(paste("The simulation created %d adopters and %d later adopters over %d weeks. Each customer's weekly outcomes carry a transitory shock",
    "with persistence %.2f (%s). Adoption was %s; a persistent stress type %s adoption. True effects were %s."),
    p$n_treated, p$n_later, p$n_weeks, p$phi,
    if (p$phi == 0) "a one-week blip" else if (p$phi < 0.6) "fading within a few weeks" else "unwinding over months",
    if (p$dip_selection == 0) "unrelated to the shock" else sprintf("%s likely in weeks following a positive shock (selection %.1f)", if (p$dip_selection < 2) "somewhat more" else "much more", p$dip_selection),
    if (p$type_selection == 0) "did not affect" else "also raised",
    if (all(unlist(p$effects) == 0)) "set to zero, so everything estimated is bias" else "non-zero, so bias is estimate minus truth")
  dsg <- sprintf(paste("The weights were fitted on %d outcome(s) over %d pre-adoption weeks%s, and the effect is averaged over the 40-week post window",
    "(adoption week plus 39). Estimators: %s."), length(.matching_set(p$matching, p$custom_outcomes)), p$L - p$match_end,
    if (p$match_end > 0) sprintf(", holding out the last %d weeks before adoption", p$match_end) else "",
    paste(.estimator_label(S$estimator), collapse = " and "))

  one <- function(s) {
    rel <- if (is.finite(s$bias_pct_effect)) sprintf("%.0f%% of the true effect", abs(s$bias_pct_effect)) else if (is.finite(s$bias_pct_baseline)) sprintf("%.1f%% of the outcome's pre-adoption level", abs(s$bias_pct_baseline)) else "relative size undefined"
    where <- if (abs(s$bias_first8) > 2 * abs(s$bias)) "concentrated in the first eight weeks and fading afterwards, the signature of a dip" else "spread over the post window, which does not look like a dip"
    held <- if (p$match_end > 0 && is.finite(s$heldout_gap)) sprintf(" Held-out pre-weeks gap: %.3f.", s$heldout_gap) else ""
    sprintf("%s: 40-week bias %.3f (%s), first-8-week bias %.3f, %s.%s", .estimator_label(s$estimator), s$bias, rel, s$bias_first8, where, held)
  }
  res <- vapply(seq_len(nrow(S)), function(i) one(S[i]), character(1))

  cmp <- character(0)
  if (all(c("sc", "ascm") %in% S$estimator)) {
    a <- S[estimator == "sc"]; b <- S[estimator == "ascm"]; ag <- r$augment[estimator == "ascm"]
    pen <- if (nrow(ag) && is.finite(ag$lambda)) sprintf("relative penalty %.3g", ag$lambda) else if (nrow(ag) && isTRUE(is.infinite(ag$lambda))) "penalty Inf (tuning chose no augmentation)" else "penalty n/a"
    cmp <- sprintf(paste("Augmentation (%s; mean negative weight mass %.2f; pre-period fit RMSE %.3g -> %.3g) changes the 40-week bias from %.3f to %.3f",
      "and the first-8-week bias from %.3f to %.3f: it %s."),
      pen, ag$neg_mass, ag$pre_rmse_sc, ag$pre_rmse_ascm, a$bias, b$bias, a$bias_first8, b$bias_first8,
      if (abs(b$bias) < 0.8 * abs(a$bias)) "reduces the dip bias" else if (abs(b$bias) > 1.25 * abs(a$bias)) "makes the dip bias worse" else "leaves the dip bias roughly unchanged")
  }

  rem <- character(0)
  for (i in seq_len(nrow(S))) {
    s <- S[i]; lab <- .estimator_label(s$estimator); rr <- character(0)
    if (!is.null(s$bias_nodip) && is.finite(s$bias_nodip)) rr <- c(rr, sprintf("dropping the %d shock-triggered adopters leaves %.3f", s$n_dropped_nodip, s$bias_nodip))
    if (!is.null(s$bias_reversion) && is.finite(s$bias_reversion)) rr <- c(rr, sprintf("the reversion adjustment (persistence %.2f from donors) leaves %.3f", s$phi_hat, s$bias_reversion))
    if (!is.null(s$bias_placebo_adj) && is.finite(s$bias_placebo_adj)) rr <- c(rr, sprintf("the placebo adjustment leaves %.3f", s$bias_placebo_adj))
    if (length(rr)) rem <- c(rem, sprintf("Remedies with %s: %s.", lab, paste(rr, collapse = "; ")))
  }

  best <- S[which.min(abs(bias))]
  verdict <- if (!is.finite(best$bias_pct_baseline)) "Verdict: relative bias undefined (pre-adoption level is zero)."
             else if (abs(best$bias_pct_baseline) < 2) sprintf("Verdict: %s is close to clean in this setting.", .estimator_label(best$estimator))
             else if (abs(best$bias_first8) > 2 * abs(best$bias)) "Verdict: the dip is visible in the early weeks; the 40-week average contains it."
             else "Verdict: the bias is material and not confined to the early weeks; a longer baseline, held-out weeks or a remedy is needed."
  c(dgp, dsg, res, cmp, rem, verdict)
}
