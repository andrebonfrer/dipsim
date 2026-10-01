# remedies.R ------------------------------------------------------------------
# Observable-based corrections for a dip, applied after the SC fit. Each takes
# the per-unit gap matrix (treated x event week, pre-window demeaned) and the
# outcome matrix, and returns an adjusted post-window ATE. All use only what
# is observed in real data: outcomes and the weights.
#
#   nodip     : drop treated units whose standardised deviation from their own
#               matched-window mean, averaged over the last `k_pre` pre-weeks,
#               exceeds `z_cut` (the adopters who look shock-triggered)
#   reversion : estimate the shock persistence phi from donors' within-unit
#               residuals, then subtract phi^(tau+1) * (deviation at tau = -1)
#               from each treated unit's post gaps (Ashenfelter & Card logic)
#   placebo   : give half the donors pseudo-onsets drawn with the same trigger
#               rule (hazard ~ exp(k * z_{t-1}) on the observed deviation),
#               refit the SC with them as treated and the other half as donors,
#               and subtract the pseudo-ATE from the real one

# Per-unit event-time gaps for one outcome from a saved budgetsc run.
unit_gaps <- function(run, outcome) {
  scfit <- readRDS(file.path(run$dir, "sc_fit.rds")); p <- readRDS(file.path(run$dir, "panel.rds")); data.table::setDT(p)
  Yp <- budgetsc::build_Ylist(p, outcome); Y <- Yp$Y_list[[1]]
  w <- scfit$fit$weights_mat; tt <- scfit$treat_time; tr <- which(is.finite(tt))
  L <- run$spec$sample$n_lags; H <- run$spec$sample$n_leads; me <- scfit$fit$match_end %||% 0L; Tn <- ncol(Y)
  synth <- w %*% Y
  g <- matrix(NA_real_, length(tr), L + H + 1); dev <- matrix(NA_real_, length(tr), L + H + 1)
  for (j in seq_along(tr)) {
    t0 <- tt[tr[j]]; idx <- (t0 - L):(t0 + H); ok <- idx >= 1 & idx <= Tn
    gap <- Y[tr[j], idx[ok]] - synth[j, idx[ok]]; pre_fit <- idx[ok] < t0 - me
    g[j, which(ok)] <- gap - mean(gap[pre_fit], na.rm = TRUE)
    y <- Y[tr[j], idx[ok]]; mu <- mean(y[pre_fit], na.rm = TRUE); s <- stats::sd(y[pre_fit], na.rm = TRUE)
    dev[j, which(ok)] <- (y - mu) / max(s, 1e-8)
  }
  # pre-adoption level of the outcome among treated units (for relative bias)
  pre_mean <- mean(vapply(seq_along(tr), function(j) { t0 <- tt[tr[j]]; idx <- (t0 - L):(t0 - 1); idx <- idx[idx >= 1]
    mean(Y[tr[j], idx], na.rm = TRUE) }, numeric(1)), na.rm = TRUE)
  list(gap = g, dev = dev, ev = seq(-L, H), ids = Yp$unit_vals[tr], Y = Y, treat_time = tt, tr = tr, L = L, H = H, me = me, pre_mean = pre_mean)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

remedy_nodip <- function(ug, k_pre = 4L, z_cut = 1) {
  pre_cols <- which(ug$ev >= -k_pre & ug$ev <= -1); post <- which(ug$ev >= 0)
  z <- rowMeans(ug$dev[, pre_cols, drop = FALSE], na.rm = TRUE)
  keep <- is.finite(z) & z <= z_cut
  list(ate = mean(ug$gap[keep, post], na.rm = TRUE), n_kept = sum(keep), n_dropped = sum(!keep))
}

remedy_reversion <- function(ug) {
  Y <- ug$Y; dn <- which(!is.finite(ug$treat_time)); post <- which(ug$ev >= 0); m1 <- which(ug$ev == -1)
  # phi from donors: lag-1 autocorrelation of within-unit residuals
  R <- Y[dn, , drop = FALSE] - rowMeans(Y[dn, , drop = FALSE], na.rm = TRUE)
  x <- as.vector(R[, -1]); l <- as.vector(R[, -ncol(R)]); ok <- is.finite(x) & is.finite(l)
  phi <- stats::cov(x[ok], l[ok]) / stats::var(l[ok]); phi <- max(min(phi, 0.99), 0)
  # raw deviation at tau = -1 (outcome units): sd * standardised deviation
  pre_fit <- which(ug$ev < -ug$me)
  dev_raw <- vapply(seq_len(nrow(ug$dev)), function(j) {
    y <- Y[ug$tr[j], ]; t0 <- ug$treat_time[ug$tr[j]]; idx <- (t0 - ug$L):(t0 - ug$me - 1)
    idx <- idx[idx >= 1]; y[t0 - 1] - mean(y[idx], na.rm = TRUE) }, numeric(1))
  adj <- ug$gap[, post, drop = FALSE] - outer(dev_raw, phi^(ug$ev[post] + 1))
  list(ate = mean(adj, na.rm = TRUE), phi_hat = phi)
}

remedy_placebo <- function(run, outcome, k, seed = 1L) {
  # pseudo-onsets for half the donors with the same trigger rule on the observed deviation
  scfit <- readRDS(file.path(run$dir, "sc_fit.rds")); smp <- readRDS(file.path(run$dir, "sample.rds"))
  p <- readRDS(file.path(run$dir, "panel.rds")); data.table::setDT(p)
  tr_w <- smp$ids[role == "treated", treat_wID]; dn <- smp$ids[role == "donor", customer_id]
  set.seed(seed); pseudo <- sample(dn, floor(length(dn) / 2)); rest <- setdiff(dn, pseudo)
  Y <- budgetsc::build_Ylist(p, outcome)$Y_list[[1]]; wk <- seq(min(tr_w), max(tr_w))
  z <- function(id) { y <- Y[as.character(id), ]; (y - mean(y[1:(min(wk) - 1)], na.rm = TRUE)) / max(stats::sd(y[1:(min(wk) - 1)], na.rm = TRUE), 1e-8) }
  onset <- vapply(pseudo, function(id) { zz <- z(id)[wk - 1]; pr <- exp(k * zz); pr[!is.finite(pr)] <- 0
    if (sum(pr) <= 0) sample(wk, 1) else sample(wk, 1, prob = pr / sum(pr)) }, integer(1))
  tl <- budgetsc::bsc_timeline()
  p2 <- p[customer_id %in% c(pseudo, rest)]
  p2[, minBudgetDate := as.Date(NA)][, donor := 0L]
  p2[data.table::data.table(customer_id = pseudo, d = budgetsc::week_start(onset, tl$origin)), on = "customer_id", minBudgetDate := i.d]
  p2[customer_id %in% rest, minBudgetDate := budgetsc::week_start(max(wk) + run$spec$sample$n_leads + 2L, tl$origin)]  # later adopters
  root2 <- file.path(tempdir(), "dipsim_placebo"); dir.create(file.path(root2, "Processed"), recursive = TRUE, showWarnings = FALSE)
  for (v in intersect(c("total_spend", "net_cashflow", "cf_vol", "income_share", "income_share_vol", "liquidity_deficit", "liquidity_deficit_rate",
                        "total_spend_roll", "liquidity_buffer", "liquidity_buffer_w", "cf_vol_w", "income_share_vol_w", "income_share_w",
                        "cv_2020_spend", "homeloan_mean", "br_under", "br_target", "br_over", "avg_transaction_size", "treat_wID", "budgetdummy"),
                      names(p2))) p2[, (v) := NULL]
  data.table::setnames(p2, grep("^Spend", names(p2), value = TRUE), sub("^Spend", "spend_", grep("^Spend", names(p2), value = TRUE)))
  p2[, c("spend_Alcohol", "spend_Tobacco") := NULL]
  p2 <- p2[, !grep("^walletshare_|^budgetever_", names(p2)), with = FALSE]
  saveRDS(p2, file.path(root2, "Processed", "analysis_panel.rds"))
  spec2 <- run$spec; spec2$name <- paste0(run$name, "_placebo"); spec2$post$outcomes <- outcome
  old <- options(budgetsc.root = root2); on.exit(options(old), add = TRUE)
  r2 <- budgetsc::run_register(spec2, root2); budgetsc::run_sample(r2, root2); budgetsc::fit_sc(r2); d2 <- budgetsc::sc_diagnostics(r2)
  oc <- outcome; pseudo_ate <- d2$pretrend[outcome == oc & window == "post", mean(mean_gap)]
  unlink(r2$dir, recursive = TRUE)
  list(pseudo_ate = pseudo_ate, n_pseudo = length(pseudo))
}

utils::globalVariables(c("treat_wID", "minBudgetDate", "donor", "i.d", "role", "customer_id"))
