# sim.R -----------------------------------------------------------------------
# Simulation engine. One `dip_run()` = one simulated panel + one SC fit + the
# statistics that reveal a dip: the event-time gap profile, the 41-week and
# first-8-week ATEs, the held-out gap (when match_end > 0) and the latent
# shock profile around adoption. `dip_grid()` repeats it over design choices.
# Everything runs through budgetsc's own functions, so what the app shows is
# what the pipeline does.
#' @importFrom data.table data.table := .N .SD fread rbindlist
NULL

#' Default simulation parameters
#' @return Named list.
#' @export
dip_defaults <- function() list(
  n_treated = 60L, n_later = 120L, n_weeks = 200L, phi = 0.9, sigma_u = 0.25,
  dip_selection = 3, type_selection = 1, covid = TRUE,
  effects = list(spend_log = 0, arrears_log = 0, arrears_length = 0, share_FastFood = 0, share_Alcohol = 0, share_Gambling = 0),
  L = 23L, H = 39L, match_end = 0L, matching = "lean", custom_outcomes = NULL, max_donors = 1000L,
  backend = if (requireNamespace("augMultiSynth", quietly = TRUE)) "augmultisynth" else "stub",
  outcomes = c("numarrears", "total_spend", "liquidity_deficit_rate"), seed = 1L,
  remedies = c("nodip", "reversion"), z_cut = 1, k_pre = 4L)

.matching_set <- function(matching, custom = NULL) switch(matching,
  single = "total_spend", lean = budgetsc::match_outcomes_lean(), wide = budgetsc::match_outcomes_wide(),
  custom = unique(c(custom, "total_spend")), stop("unknown matching set"))

#' One simulated panel, one SC fit
#'
#' @param p Parameter list as from [dip_defaults()] (missing entries filled).
#' @param root Scratch directory for the budgetsc run (deleted afterwards
#'   unless `keep = TRUE`).
#' @param keep Logical; keep the run directory.
#' @return A list: `event_time` (gap by event week, matched and holdout
#'   outcomes), `summary` (per outcome: ate, ate_first8, heldout_gap, true
#'   effect), `shock_profile` (mean latent shock u by event week, treated
#'   units), `funnel`, `params`, `seconds`.
#' @export
dip_run <- function(p = dip_defaults(), root = file.path(tempdir(), "dipsim"), keep = FALSE) {
  p <- utils::modifyList(dip_defaults(), p)
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
  spec <- budgetsc::spec_main("sample.n_lags" = as.integer(p$L), "sample.n_leads" = as.integer(p$H),
                              "sc.match_outcomes" = mo, "sc.match_end" = as.integer(p$match_end),
                              "sc.max_donors" = as.integer(p$max_donors), "sc.backend" = p$backend,
                              "sc.parallel" = FALSE, "post.outcomes" = p$outcomes,
                              "sc.holdout_outcomes" = c("weekly_income", "monthly_loan_payment_w"))
  spec$name <- paste0("dipsim_", format(Sys.time(), "%H%M%S"), "_", sample.int(1e6, 1))
  run <- budgetsc::run_register(spec, root)
  budgetsc::run_sample(run, root); budgetsc::fit_sc(run); d <- budgetsc::sc_diagnostics(run)
  # per-unit gaps and remedies for the post-estimated outcomes
  rem <- lapply(p$outcomes, function(o) {
    ug <- unit_gaps(run, o); out <- list(outcome = o, pre_mean = ug$pre_mean)
    if ("nodip" %in% p$remedies) { r <- remedy_nodip(ug, p$k_pre, p$z_cut); out$ate_nodip <- r$ate; out$n_dropped_nodip <- r$n_dropped }
    if ("reversion" %in% p$remedies) { r <- remedy_reversion(ug); out$ate_reversion <- r$ate; out$phi_hat <- r$phi_hat }
    if ("placebo" %in% p$remedies) { r <- tryCatch(remedy_placebo(run, o, k = p$dip_selection, seed = p$seed), error = function(e) NULL)
      if (!is.null(r)) out$pseudo_ate <- r$pseudo_ate }
    data.table::as.data.table(out) })
  rem <- rbindlist(rem, fill = TRUE)
  et <- d$pretrend
  ate <- et[window == "post", list(ate = mean(mean_gap)), by = outcome]
  first8 <- et[window == "post" & tau_k < 8, list(ate_first8 = mean(mean_gap)), by = outcome]
  held <- d$fit_summary[window == "heldout_pre", list(outcome, heldout_gap = mean_gap)]
  summ <- merge(merge(ate, first8, by = "outcome"), held, by = "outcome", all.x = TRUE)
  summ[, true_effect := .true_effect(outcome, p$effects, panel)]
  summ[, matched := outcome %in% mo]
  summ <- merge(summ, rem, by = "outcome", all.x = TRUE)
  if ("pseudo_ate" %in% names(summ)) summ[, ate_placebo_adj := ate - pseudo_ate]
  summ[, bias := ate - true_effect]
  summ[, bias_pct_baseline := 100 * bias / abs(pre_mean)]                      # relative to the outcome's pre-adoption level
  summ[, bias_pct_effect := data.table::fifelse(true_effect != 0, 100 * bias / abs(true_effect), NA_real_)]  # relative to the true effect
  for (v in intersect(c("ate_nodip", "ate_reversion", "ate_placebo_adj"), names(summ)))
    summ[, (sub("^ate_", "bias_", v)) := get(v) - true_effect]
  tr <- readRDS(file.path(run$dir, "sample.rds"))$ids[role == "treated"]
  sp <- panel[tr, on = "customer_id"][, list(mean_u = mean(u)), by = list(tau_k = wID - treat_wID)][tau_k >= -p$L & tau_k <= p$H][order(tau_k)]
  funnel <- fread(file.path(run$dir, "tables", "sample_funnel.csv"))
  if (!keep) unlink(run$dir, recursive = TRUE)
  list(event_time = et, summary = summ, shock_profile = sp, funnel = funnel, params = p,
       seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

# True effect per outcome, on the outcome's scale, from treated units'
# pre-adoption means (approximate for the multiplicative effects).
.true_effect <- function(outcomes, eff, panel) {
  pre <- panel[D == 0]
  spend_cols <- grep("^spend_", names(pre), value = TRUE)
  base_spend <- mean(rowSums(pre[, spend_cols, with = FALSE]))
  vapply(outcomes, function(o) switch(o,
    total_spend = (exp(eff$spend_log) - 1) * base_spend,
    numarrears  = (exp(eff$arrears_log) - 1) * mean(pre$numarrears),
    lengtharrears = eff$arrears_length * mean(pre$numarrears > 0),
    0), numeric(1))
}

#' Repeat `dip_run()` over a design grid
#' @param p Base parameters.
#' @param grid data.frame of columns to vary (any of L, match_end, matching,
#'   n_later, phi, dip_selection, max_donors).
#' @param reps Replications per cell.
#' @param progress Optional function(i, n) called after each run.
#' @return data.table: one row per cell x rep x outcome with the summary
#'   statistics and the cell's factors.
#' @export
dip_grid <- function(p = dip_defaults(), grid, reps = 1L, progress = NULL) {
  grid <- data.table::as.data.table(grid); out <- list(); k <- 0L; n <- nrow(grid) * reps
  for (i in seq_len(nrow(grid))) for (r in seq_len(reps)) {
    pp <- utils::modifyList(p, as.list(grid[i])); pp$seed <- p$seed + 1000L * r + i
    res <- tryCatch(dip_run(pp), error = function(e) NULL)
    k <- k + 1L; if (!is.null(progress)) progress(k, n)
    if (is.null(res)) next
    out[[length(out) + 1L]] <- cbind(grid[i], rep = r, res$summary)
  }
  rbindlist(out, fill = TRUE)
}

utils::globalVariables(c("window", "mean_gap", "tau_k", "outcome", "true_effect", "matched", "role", "pseudo_ate",
                         "ate_placebo_adj", "bias", "pre_mean", "bias_pct_baseline", "bias_pct_effect", "ate",
                         "treat_wID", "u", "wID", "D", "numarrears", "customer_id"))


#' Plain-language interpretation of one run
#' @param r Result of [dip_run()].
#' @param outcome Outcome to describe. Default `"total_spend"`.
#' @return Character vector of paragraphs.
#' @export
interpret_run <- function(r, outcome = "total_spend") {
  p <- r$params; oc <- outcome; s <- r$summary[outcome == oc][1]
  if (!nrow(s) || is.na(s$ate)) return("No result for this outcome.")
  dgp <- sprintf(paste("The simulation created %d adopters and %d later adopters over %d weeks. Each customer's weekly outcomes carry a transitory shock",
    "with persistence %.2f (%s). Adoption was %s; a persistent stress type %s adoption. True effects were %s."),
    p$n_treated, p$n_later, p$n_weeks, p$phi,
    if (p$phi == 0) "a one-week blip" else if (p$phi < 0.6) "fading within a few weeks" else "unwinding over months",
    if (p$dip_selection == 0) "unrelated to the shock" else sprintf("%s likely in weeks following a positive shock (selection %.1f)", if (p$dip_selection < 2) "somewhat more" else "much more", p$dip_selection),
    if (p$type_selection == 0) "did not affect" else "also raised",
    if (all(unlist(p$effects) == 0)) "set to zero, so everything estimated is bias" else "non-zero, so bias is estimate minus truth")
  dsg <- sprintf(paste("The synthetic control matched %d outcome(s) over %d pre-adoption weeks%s, and the effect is averaged over the 40-week post window",
    "(adoption week plus 39)."), length(.matching_set(p$matching, p$custom_outcomes)), p$L - p$match_end,
    if (p$match_end > 0) sprintf(", holding out the last %d weeks before adoption", p$match_end) else "")
  rel <- if (is.finite(s$bias_pct_effect)) sprintf("%.0f%% of the true effect", abs(s$bias_pct_effect)) else if (is.finite(s$bias_pct_baseline)) sprintf("%.1f%% of the outcome's pre-adoption level", abs(s$bias_pct_baseline)) else "relative size undefined"
  where <- if (abs(s$ate_first8 - s$true_effect) > 2 * abs(s$bias)) "concentrated in the first eight weeks after adoption and fading afterwards, the signature of a dip"
           else "spread evenly over the post window, which does not look like a dip"
  held <- if (p$match_end > 0 && is.finite(s$heldout_gap)) sprintf(" The held-out pre-weeks show a gap of %.3f, %s.", s$heldout_gap,
            if (abs(s$heldout_gap) > abs(s$bias)) "larger than the full-window bias and pointing the same way the early post weeks do" else "small relative to the full-window bias") else ""
  res <- sprintf("For %s the 40-week bias is %.3f, %s, with the bias %s.%s", outcome, s$bias, rel, where, held)
  rem <- character(0)
  if (!is.null(s$bias_nodip) && is.finite(s$bias_nodip)) rem <- c(rem, sprintf("dropping the %d adopters who looked shock-triggered leaves a bias of %.3f", s$n_dropped_nodip, s$bias_nodip))
  if (!is.null(s$bias_reversion) && is.finite(s$bias_reversion)) rem <- c(rem, sprintf("subtracting the expected unwinding of each adopter's pre-adoption deviation (persistence estimated from donors as %.2f) leaves %.3f", s$phi_hat, s$bias_reversion))
  if (!is.null(s$bias_placebo_adj) && is.finite(s$bias_placebo_adj)) rem <- c(rem, sprintf("subtracting the pseudo-effect from donors given shock-triggered pseudo-onsets leaves %.3f", s$bias_placebo_adj))
  remtxt <- if (length(rem)) paste0("Remedies: ", paste(rem, collapse = "; "), ". A remedy works when its bias is closer to zero than the raw one.") else ""
  verdict <- if (!is.finite(s$bias_pct_baseline)) "Verdict: relative bias undefined (pre-adoption level is zero)."
             else if (abs(s$bias_pct_baseline) < 2) "Verdict: the headline estimate is close to clean in this setting."
             else if (abs(s$ate_first8 - s$true_effect) > 2 * abs(s$bias)) "Verdict: the dip is visible in the early weeks and the diagnostics, but the 40-week average contains it."
             else "Verdict: the bias is material and not confined to the early weeks; a longer baseline or a remedy is needed."
  c(dgp, dsg, res, remtxt, verdict)
}
