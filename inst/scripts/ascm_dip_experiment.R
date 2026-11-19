# inst/scripts/ascm_dip_experiment.R
#
# How does ridge augmentation (ASCM) behave under an Ashenfelter dip?
#
# Paired design: every panel is fitted with SC and with ASCM at three penalty
# settings (tuned by the one-SE rule, fixed 1, fixed 0.1), so differences come
# from the estimator, not the draw. True effects are zero, so every estimate is
# bias. Factors:
#   dip_selection : 0 (no dip) or 3 (adoption follows a positive shock)
#   L             : 23 and 70 (the paper's pilot and public-launch baselines)
#   match_end     : 0 (spike inside the matching window) or 4 (held out)
# Reported per cell: mean bias of SC and of each ASCM variant, and the mean
# paired change (ASCM - SC) with its standard error over replications.
#
# Run (about 30 minutes on two cores):
#   Rscript -e 'source(system.file("scripts/ascm_dip_experiment.R", package = "dipsim"))'

library(dipsim); library(data.table)
out_dir <- getwd()
reps <- 10L
outcome_main <- "total_spend"

variants <- list(`ASCM tuned (one-SE)` = list(augment_lambda = NULL, augment_one_se = TRUE),
                 `ASCM fixed 1`        = list(augment_lambda = 1),
                 `ASCM fixed 0.1`      = list(augment_lambda = 0.1))
cells <- CJ(dip_selection = c(0, 3), L = c(23L, 70L), match_end = c(0L, 4L))

base <- dip_defaults()
base$remedies <- "reversion"; base$estimators <- c("sc", "ascm"); base$n_weeks <- 220L

rows <- list(); k <- 0L; n <- nrow(cells) * reps
for (i in seq_len(nrow(cells))) for (r in seq_len(reps)) {
  p <- utils::modifyList(base, as.list(cells[i]), keep.null = TRUE)
  p$seed <- 50000L + 100L * i + r                          # same panel for every variant
  for (v in names(variants)) {
    pv <- utils::modifyList(p, variants[[v]], keep.null = TRUE)
    if (is.null(variants[[v]]$augment_lambda)) pv["augment_lambda"] <- list(NULL)
    res <- tryCatch(dip_run(pv), error = function(e) { message("cell ", i, " rep ", r, " ", v, ": ", conditionMessage(e)); NULL })
    if (is.null(res)) next
    s <- merge(res$summary, res$augment, by = "estimator", all.x = TRUE)
    s[, variant := ifelse(estimator == "sc", "SC", v)]
    rows[[length(rows) + 1L]] <- cbind(cells[i], rep = r, s)
  }
  k <- k + 1L; if (k %% 10 == 0) message(k, " / ", n, " panels")
}
res <- rbindlist(rows, fill = TRUE)
res <- unique(res, by = c("dip_selection", "L", "match_end", "rep", "outcome", "variant"))   # SC fitted once per variant run
fwrite(res, file.path(out_dir, "ascm_dip_experiment.csv"))

# ---- paired summary -----------------------------------------------------------
m <- res[outcome == outcome_main]
sc <- m[variant == "SC", .(dip_selection, L, match_end, rep, sc_bias = bias, sc_first8 = bias_first8)]
pr <- merge(m[variant != "SC"], sc, by = c("dip_selection", "L", "match_end", "rep"))
summ <- pr[, .(sc_bias = mean(sc_bias), ascm_bias = mean(bias),
               change = mean(bias - sc_bias), change_se = sd(bias - sc_bias) / sqrt(.N),
               abs_smaller = mean(abs(bias) < abs(sc_bias)),
               sc_first8 = mean(sc_first8), ascm_first8 = mean(bias_first8),
               lambda = paste(unique(signif(lambda, 2)), collapse = ","), neg_mass = mean(neg_mass), n = .N),
           by = .(dip_selection, L, match_end, variant)]
setorder(summ, dip_selection, match_end, L, variant)
fwrite(summ, file.path(out_dir, "ascm_dip_experiment_summary.csv"))
print(summ, digits = 3, nrows = 200)
