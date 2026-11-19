# app.R -----------------------------------------------------------------------
# Two modes. "Single run": choose a data-generating process and a design, run
# once, and look at the diagnostics that reveal a dip, for plain synthetic
# control (SC) and ridge-augmented SC (ASCM) fitted on the same simulated
# panel. "Compare designs": hold the DGP fixed and sweep baseline length,
# held-out weeks and the matching set, with both estimators in every cell.
# Both call dip_run().
#
# To extend: add a parameter to dip_defaults() and an input here; add a remedy
# by adding a spec option in budgetsc and exposing it in .params_from_input().

`%||%` <- function(a, b) if (is.null(a)) b else a

.est_cols <- c(sc = "#222222", ascm = "#1f5fa8")

.params_from_input <- function(input) {
  p <- dip_defaults()
  p$n_treated <- as.integer(input$n_treated); p$n_later <- as.integer(input$n_later); p$n_weeks <- as.integer(input$n_weeks)
  p$phi <- input$phi; p$sigma_u <- input$sigma_u; p$dip_selection <- input$dip_selection; p$type_selection <- input$type_selection
  p$covid <- isTRUE(input$covid)
  p$effects <- list(spend_log = input$eff_spend, arrears_log = input$eff_arrears, arrears_length = input$eff_length,
                    share_FastFood = 0, share_Alcohol = 0, share_Gambling = 0)
  p$L <- as.integer(input$L); p$H <- 39L; p$match_end <- as.integer(input$match_end)
  p$remedies <- input$remedies %||% character(0); p$z_cut <- input$z_cut; p$k_pre <- as.integer(input$k_pre)
  p$donor_eligibility <- input$donor_eligibility %||% "per_unit"; p$matching <- input$matching; p$custom_outcomes <- input$custom_outcomes; p$max_donors <- as.integer(input$max_donors)
  p$backend <- input$backend; p$seed <- as.integer(input$seed)
  p$estimators <- input$estimators %||% "sc"
  pm <- input$penalty_mode %||% "tuned_1se"
  p$augment_one_se <- !identical(pm, "tuned_min")
  p$augment_lambda <- if (identical(pm, "fixed")) input$augment_lambda else NULL
  p$augment_tune_units <- as.integer(input$augment_tune_units %||% 40L)
  p
}

#' Shiny UI
#' @return shiny UI.
#' @export
app_ui <- function() {
  d <- dip_defaults()
  shiny::fluidPage(
    shiny::titlePanel("Ashenfelter-dip simulator: synthetic control and ridge-augmented synthetic control (budgetsc)"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(width = 3,
        shiny::h4("Data-generating process"),
        shiny::sliderInput("phi", "Shock persistence phi (0 = white noise)", 0, 0.95, d$phi, step = 0.05),
        shiny::sliderInput("dip_selection", "Selection on last week's shock", 0, 5, d$dip_selection, step = 0.5),
        shiny::sliderInput("type_selection", "Selection on persistent stress type", 0, 3, d$type_selection, step = 0.5),
        shiny::numericInput("sigma_u", "Shock SD", d$sigma_u, step = 0.05),
        shiny::checkboxInput("covid", "Common lockdown-like shocks in 2020", d$covid),
        shiny::fluidRow(shiny::column(6, shiny::numericInput("n_treated", "treated", d$n_treated)), shiny::column(6, shiny::numericInput("n_later", "later adopters (donors)", d$n_later))),
        shiny::numericInput("n_weeks", "weeks in panel", d$n_weeks),
        shiny::h5("True effects (0 = pure bias)"),
        shiny::fluidRow(shiny::column(4, shiny::numericInput("eff_spend", "spend (log)", 0, step = 0.05)), shiny::column(4, shiny::numericInput("eff_arrears", "arrears (log)", 0, step = 0.1)),
                        shiny::column(4, shiny::numericInput("eff_length", "length (wks)", 0, step = 1))),
        shiny::hr(), shiny::h4("Estimators"),
        shiny::checkboxGroupInput("estimators", NULL,
          c("Synthetic control (SC)" = "sc", "Ridge-augmented SC (ASCM)" = "ascm"), d$estimators),
        shiny::conditionalPanel("input.estimators.indexOf('ascm') > -1",
          shiny::radioButtons("penalty_mode", "ASCM penalty",
            c("Tuned, one-SE rule" = "tuned_1se", "Tuned, minimum" = "tuned_min", "Fixed" = "fixed"), "tuned_1se"),
          shiny::conditionalPanel("input.penalty_mode == 'fixed'",
            shiny::numericInput("augment_lambda", "relative penalty (Inf = none)", 1, min = 0, step = 0.5)),
          shiny::numericInput("augment_tune_units", "units used for tuning", d$augment_tune_units, min = 5),
          shiny::helpText("Tuning moves adoption back a few weeks and scores forecasts of the held-out pre-weeks. With shock-triggered adoption those weeks contain the dip, so the tuning itself sees it. See Notes.")),
        shiny::hr(), shiny::h4("Design"),
        shiny::numericInput("L", "L: baseline weeks", d$L, min = 4),
        shiny::p(shiny::em("Post window fixed at 40 weeks (adoption week + 39).")),
        shiny::numericInput("match_end", "Held-out pre-weeks (match_end)", d$match_end, min = 0),
        shiny::helpText("The last match_end weeks before adoption are excluded when the weights are fitted, then used as a check: if the synthetic control tracks them, the weights generalise beyond the weeks they were fitted on. See Notes."),
        shiny::selectInput("matching", "Matching set", c("single", "lean", "wide", "custom"), d$matching),
        shiny::conditionalPanel("input.matching == 'custom'",
          shiny::checkboxGroupInput("custom_outcomes", NULL, budgetsc::match_outcomes_wide(), budgetsc::match_outcomes_lean())),
        shiny::numericInput("max_donors", "max donors (screening)", d$max_donors),
        shiny::selectInput("donor_eligibility", "Donor eligibility", c("per_unit", "global"), d$donor_eligibility),
        shiny::helpText("per_unit: a donor must adopt more than H weeks after the treated unit (the rule the paper states). global: every later adopter is a donor for every treated unit, including those that adopt inside a treated unit's post window, which biases effects toward zero."),
        shiny::selectInput("backend", "SC backend", c("augmultisynth", "stub"), d$backend),
        shiny::numericInput("seed", "seed", d$seed),
        shiny::h4("Remedies"),
        shiny::checkboxGroupInput("remedies", NULL, c("Drop shock-triggered adopters" = "nodip", "Reversion adjustment (phi from donors)" = "reversion",
                                                     "Placebo bias estimate (second fit)" = "placebo"), d$remedies),
        shiny::fluidRow(shiny::column(6, shiny::numericInput("z_cut", "no-dip z cut", d$z_cut, step = 0.25)), shiny::column(6, shiny::numericInput("k_pre", "no-dip weeks", d$k_pre))),
        shiny::actionButton("run", "Run once", class = "btn-primary"), shiny::br(), shiny::br(),
        shiny::textOutput("timing")),
      shiny::mainPanel(width = 9, shiny::tabsetPanel(
        shiny::tabPanel("Single run",
          shiny::p(shiny::em("A dip shows as: (1) a latent shock that peaks just before adoption, (2) an early post-adoption bulge in the gap that fades, (3) a non-zero held-out gap. The 40-week ATE minus the true effect is the bias that reaches the paper. SC and ASCM are fitted on the same simulated panel.")),
          shiny::fluidRow(shiny::column(6, shiny::h5("Latent shock around adoption (what the DGP did)"), shiny::plotOutput("p_shock", height = 280)),
                          shiny::column(6, shiny::h5("Estimated gap by event week"), shiny::selectInput("oc", NULL, d$outcomes), shiny::plotOutput("p_gap", height = 280))),
          shiny::h5("Bias summary (estimate minus true effect; relative to the outcome's pre-adoption level, and to the true effect when non-zero)"),
          shiny::tableOutput("t_summary"),
          shiny::h5("Augmentation diagnostics"),
          shiny::tableOutput("t_augment"),
          shiny::h5("Interpretation"), shiny::uiOutput("interp"),
          shiny::h5("Sample funnel"), shiny::tableOutput("t_funnel")),
        shiny::tabPanel("Compare designs",
          shiny::p("Same DGP and estimators as in the sidebar; the design factors below are swept. Each cell is a full simulate + fit, with every selected estimator fitted on the same panel."),
          shiny::fluidRow(shiny::column(3, shiny::textInput("grid_L", "L values", "8, 23, 70")),
                          shiny::column(3, shiny::textInput("grid_me", "match_end values", "0, 4")),
                          shiny::column(3, shiny::checkboxGroupInput("grid_matching", "matching sets", c("single", "lean", "wide"), c("lean"))),
                          shiny::column(3, shiny::numericInput("reps", "replications", 1, min = 1, max = 10))),
          shiny::actionButton("run_grid", "Run grid", class = "btn-primary"), shiny::span("  "), shiny::downloadButton("dl_grid", "Download results"),
          shiny::br(), shiny::br(), shiny::selectInput("grid_oc", "Outcome", d$outcomes),
          shiny::fluidRow(shiny::column(6, shiny::h5("40-week bias vs baseline length"), shiny::plotOutput("p_grid_ate", height = 320)),
                          shiny::column(6, shiny::h5("First-8-week bias vs baseline length"), shiny::plotOutput("p_grid_first8", height = 320))),
          shiny::p(shiny::em("Solid lines: SC. Dashed lines: ridge ASCM. Colours: design (held-out weeks, matching set).")),
          shiny::h5("Reading"), shiny::uiOutput("grid_reading"),
          shiny::h5("Cell means"), shiny::tableOutput("t_grid")),
        shiny::tabPanel("Notes", .notes_ui())
      ))
    )
  )
}

.notes_ui <- function() {
  shiny::tagList(
    shiny::h4("Glossary"),
    shiny::p("The simulator reproduces, on invented customers, the estimation used in the paper, and then checks it against a truth that is known because it was simulated."),
    shiny::tags$dl(
      shiny::tags$dt("Synthetic control (SC)"), shiny::tags$dd("For each adopter, a weighted average of later adopters (the donors), with non-negative weights that sum to one, chosen so that before adoption the average tracks the adopter's own outcomes as closely as possible. After adoption it stands in for what the adopter would have done without the budget. The weights are fitted on several outcomes at once (the matching set)."),
      shiny::tags$dt("Ridge-augmented SC (ASCM)"), shiny::tags$dd(paste("SC plus a correction from a ridge regression fitted on the donors: how a customer's pre-adoption outcomes predict their later outcomes.",
        "Where the SC weights leave the adopter's pre-adoption history imperfectly matched, the correction adds what the regression predicts that leftover difference to mean later on.",
        "The result is again a weighted average of donors, but some weights can be negative (extrapolation beyond the donors). It is the method of Ben-Michael, Feller and Rothstein, run through budgetsc's sc.augment option.")),
      shiny::tags$dt("Penalty"), shiny::tags$dd(paste("How strongly the ridge regression is shrunk. Large values give back plain SC; small values fit the pre-adoption history more closely and extrapolate more.",
        "'Tuned' chooses it by moving adoption back a few weeks, forecasting the held-out pre-weeks, and picking the penalty with the smallest error (the one-SE rule picks the largest penalty within one standard error of that minimum).")),
      shiny::tags$dt("Negative weight mass"), shiny::tags$dd("The total of the negative ASCM weights for an adopter: how much the counterfactual relies on extrapolation. Zero for SC."),
      shiny::tags$dt("Gap"), shiny::tags$dd("The adopter's outcome minus the synthetic control's, week by week, after removing the average difference over the matching window (an intercept). Averaged over adopters it is the event-time plot; averaged over the post window it is the estimated effect."),
      shiny::tags$dt("Baseline weeks (L)"), shiny::tags$dd("How many pre-adoption weeks the weights are fitted on. The paper uses 70 for public-launch adopters and 23 for pilot adopters."),
      shiny::tags$dt("Post window"), shiny::tags$dd("The 40 weeks from adoption (week 0) to week 39 over which the effect is averaged."),
      shiny::tags$dt("Donor eligibility"), shiny::tags$dd("Donors are later adopters. With global eligibility every later adopter can be a donor for every adopter, including donors that start budgeting inside an adopter's post-adoption window; their outcomes then carry a treatment effect and the estimate is pulled toward zero. Per-unit eligibility allows a donor only if it adopts more than H weeks after that adopter."),
      shiny::tags$dt("Held-out pre-weeks (match_end)"), shiny::tags$dd(paste("Weeks immediately before adoption that are deliberately NOT used when fitting the weights (or, for ASCM, the ridge correction).",
        "If the synthetic control still tracks the adopter through those weeks, it is capturing the adopter's underlying pattern rather than fitting noise. The average gap over them is the held-out gap.",
        "Holding out weeks also changes where a pre-adoption spike sits: inside the matching window the weights partly reproduce it, outside they do not.")),
      shiny::tags$dt("Transitory shock, persistence (phi)"), shiny::tags$dd("Each customer's outcomes carry a temporary disturbance that decays by the factor phi each week. phi = 0 means a spike is gone next week; phi = 0.9 means it unwinds over months."),
      shiny::tags$dt("Selection on the shock (the Ashenfelter dip)"), shiny::tags$dd("Customers are made more likely to adopt in the weeks after a positive shock, as a real customer might set a budget after an expensive month. The post-adoption improvement that follows is partly the shock unwinding, not the budget. This is the mechanism named after Ashenfelter's observation that trainees' earnings dipped just before they enrolled."),
      shiny::tags$dt("Selection on stress type"), shiny::tags$dd("Adoption also depends on a permanent customer trait that raises arrears. The synthetic control handles this kind of selection; it is included so the two kinds can be compared."),
      shiny::tags$dt("True effect"), shiny::tags$dd("What adoption actually did in the simulation. Zero by default, so every estimated effect is error."),
      shiny::tags$dt("Bias"), shiny::tags$dd("Estimated effect minus true effect, over the 40-week window. Reported in the outcome's units, as a percentage of the adopters' pre-adoption level of that outcome, and (when the true effect is not zero) as a percentage of the true effect."),
      shiny::tags$dt("First-8-week bias"), shiny::tags$dd("The same quantity over weeks 0 to 7 only. A dip shows up here first: a shock unwinding after adoption makes the early gap large while the 40-week average dilutes it."),
      shiny::tags$dt("Latent shock profile"), shiny::tags$dd("The simulated shock itself, averaged over adopters by week relative to adoption. Visible only because the data are simulated. A peak just before week 0 is the dip.")),
    shiny::h5("Why augmentation might or might not help with a dip"),
    shiny::tags$ul(
      shiny::tags$li("The ridge correction is learned from donors. Donors carry the same kind of transitory shocks, which also unwind, so in principle the regression could learn that a high recent value predicts a lower later one and predict part of the adopter's unwinding, as the reversion remedy does by hand. The experiment below did not find this: with the spike inside the window, ASCM carried the spike into the counterfactual instead."),
      shiny::tags$li("With the spike held out (match_end > 0), the correction never sees it and cannot help; ASCM then mainly improves the match on the earlier weeks."),
      shiny::tags$li("The tuning scores forecasts of the last pre-weeks. For shock-triggered adopters those weeks contain the dip, so the tuning is itself distorted by it and may choose no augmentation even when a fixed penalty would help. Compare 'Tuned' with 'Fixed' in the sidebar."),
      shiny::tags$li("Extrapolation (negative weights) lets the counterfactual follow an adopter whose recent history is unusual. With a dip, that unusual history is the shock, so close fitting can carry the spike into the counterfactual and move the bias the other way.")),
    shiny::h5("Findings"),
    shiny::uiOutput("notes_findings"),
    shiny::h5("What each control does"),
    shiny::tags$ul(
      shiny::tags$li(shiny::strong("phi"), ": persistence of the idiosyncratic shock u."),
      shiny::tags$li(shiny::strong("Selection on last week's shock"), ": the hazard of adopting in a week is proportional to exp(k * u_{t-1}). This is the Ashenfelter mechanism; 0 switches it off."),
      shiny::tags$li(shiny::strong("Selection on stress type"), ": adoption depends on a persistent customer trait that also raises arrears."),
      shiny::tags$li(shiny::strong("Estimators / ASCM penalty"), ": which estimators are fitted on each panel, and how the ridge penalty is chosen."),
      shiny::tags$li(shiny::strong("L"), ": weeks matched before adoption. Dilutes the spike's leverage on the weights and the outcome intercept."),
      shiny::tags$li(shiny::strong("match_end"), ": the last k pre-weeks are held out of matching; their gap is the held-out diagnostic."),
      shiny::tags$li(shiny::strong("Matching set"), ": how many outcomes the weights match on. Pins down latent loadings; does not dilute an idiosyncratic shock."),
      shiny::tags$li(shiny::strong("Remedies"), ": observable-based corrections, applied to each estimator's gaps and judged by whether their bias is closer to zero than the raw one.")),
    shiny::h5("Extending"),
    shiny::p("Add a DGP parameter in dipsim::dip_defaults() and budgetsc::sim_panel(); add a remedy as a budgetsc spec option and expose it in .params_from_input(). The grid tab sweeps whatever dip_grid() is given.")
  )
}

#' Shiny server
#' @param input,output,session Shiny objects.
#' @export
app_server <- function(input, output, session) {
  res <- shiny::eventReactive(input$run, {
    shiny::req(length(input$estimators) > 0)
    shiny::withProgress(message = "simulating and fitting", value = 0.3, {
      r <- dip_run(.params_from_input(input)); shiny::incProgress(0.7); r })
  })
  output$timing <- shiny::renderText({ r <- res(); sprintf("last run: %.1f s (%s backend; %s)", r$seconds, r$params$backend,
                                                           paste(.estimator_label(r$params$estimators), collapse = " + ")) })
  output$p_shock <- shiny::renderPlot({ r <- res(); sp <- r$shock_profile
    graphics::plot(sp$tau_k, sp$mean_u, type = "l", lwd = 2, xlab = "weeks since adoption", ylab = "mean latent shock u (treated)")
    graphics::abline(h = 0, lty = 2); graphics::abline(v = 0, col = "red") })
  output$p_gap <- shiny::renderPlot({ r <- res(); e <- r$event_time[outcome == input$oc]
    tv <- r$summary[outcome == input$oc, true_effect][1]
    ests <- unique(e$estimator)
    graphics::plot(range(e$tau_k), range(c(e$mean_gap, 0), na.rm = TRUE), type = "n", xlab = "weeks since adoption",
                   ylab = "treated - synthetic", main = input$oc)
    me <- r$params$match_end; if (me > 0) graphics::rect(-me - 0.5, -1e9, -0.5, 1e9, col = grDevices::adjustcolor("orange", 0.15), border = NA)
    graphics::abline(h = 0, lty = 2); graphics::abline(v = 0, col = "red")
    if (length(tv) && is.finite(tv) && tv != 0) graphics::abline(h = tv, col = "darkgreen", lty = 3)
    for (k in ests) { s <- e[estimator == k][order(tau_k)]
      graphics::lines(s$tau_k, s$mean_gap, type = "b", pch = 16, cex = 0.6, col = .est_cols[[k]], lty = if (k == "sc") 1 else 2) }
    graphics::legend("topleft", .estimator_label(ests), col = .est_cols[ests], lty = ifelse(ests == "sc", 1, 2), pch = 16, bty = "n", cex = 0.85) })
  output$t_summary <- shiny::renderTable({ s <- data.table::copy(res()$summary)
    s[, estimator := .estimator_label(estimator)]
    keep <- intersect(c("estimator", "outcome", "matched", "pre_mean", "true_effect", "ate", "bias", "bias_pct_baseline", "bias_pct_effect", "bias_first8", "heldout_gap",
                        "bias_nodip", "n_dropped_nodip", "bias_reversion", "phi_hat", "bias_placebo_adj"), names(s))
    s[outcome %in% res()$params$outcomes, keep, with = FALSE][order(outcome, estimator)] }, digits = 3)
  output$t_augment <- shiny::renderTable({ a <- data.table::copy(res()$augment)[estimator == "ascm"]
    shiny::validate(shiny::need(nrow(a) > 0, "ASCM not fitted in this run."))
    data.table::data.table(`relative penalty` = format(a$lambda), `mean negative weight mass` = a$neg_mass,
                           `max negative weight mass` = a$max_neg_mass,
                           `pre-period RMSE, SC weights` = a$pre_rmse_sc, `pre-period RMSE, ASCM weights` = a$pre_rmse_ascm) }, digits = 3)
  output$interp <- shiny::renderUI({ txt <- interpret_run(res(), input$oc); shiny::tagList(lapply(txt[nzchar(txt)], shiny::p)) })
  output$t_funnel <- shiny::renderTable(res()$funnel)
  output$notes_findings <- shiny::renderUI(shiny::tagList(lapply(dip_findings(), shiny::p)))

  grid_res <- shiny::eventReactive(input$run_grid, {
    Ls <- as.integer(trimws(strsplit(input$grid_L, ",")[[1]])); mes <- as.integer(trimws(strsplit(input$grid_me, ",")[[1]]))
    shiny::req(length(input$grid_matching) > 0, length(input$estimators) > 0)
    g <- expand.grid(L = Ls, match_end = mes, matching = input$grid_matching, stringsAsFactors = FALSE)
    g <- g[g$match_end < g$L, ]
    shiny::withProgress(message = "running grid", value = 0, {
      dip_grid(.params_from_input(input), g, reps = as.integer(input$reps),
               progress = function(i, n) shiny::setProgress(i / n, detail = sprintf("%d / %d", i, n))) })
  })
  grid_cells <- shiny::reactive({ g <- grid_res(); shiny::req(nrow(g) > 0)
    cols <- intersect(c("bias_nodip", "bias_reversion", "bias_placebo_adj"), names(g))
    g[outcome == input$grid_oc, c(list(bias = mean(bias), bias_pct_baseline = mean(bias_pct_baseline), bias_first8 = mean(bias_first8),
                                       heldout_gap = mean(heldout_gap, na.rm = TRUE), neg_mass = mean(neg_mass, na.rm = TRUE), n = .N),
                                  lapply(.SD, mean, na.rm = TRUE)),
      by = list(L, match_end, matching, estimator), .SDcols = cols][order(matching, match_end, estimator, L)] })
  output$grid_reading <- shiny::renderUI({ g <- grid_cells(); shiny::req(nrow(g) > 1)
    best <- g[which.min(abs(bias))]; worst <- g[which.max(abs(bias))]
    txt <- sprintf(paste("Across %d estimator x design cells for %s, the 40-week bias is smallest with %s, L = %d, match_end = %d, %s matching (%.3f, %.1f%% of the pre-adoption level)",
      "and largest with %s, L = %d, match_end = %d, %s matching (%.3f, %.1f%%)."),
      nrow(g), input$grid_oc, .estimator_label(best$estimator), best$L, best$match_end, best$matching, best$bias, best$bias_pct_baseline,
      .estimator_label(worst$estimator), worst$L, worst$match_end, worst$matching, worst$bias, worst$bias_pct_baseline)
    if (all(c("sc", "ascm") %in% g$estimator)) {
      w <- data.table::dcast(g, L + match_end + matching ~ estimator, value.var = "bias")
      better <- sum(abs(w$ascm) < abs(w$sc)); txt <- paste(txt, sprintf("ASCM has the smaller absolute bias in %d of %d designs.", better, nrow(w)))
    }
    shiny::p(txt) })
  .grid_plot <- function(col, ylab) {
    g <- grid_cells(); keys <- unique(g[, list(match_end, matching)]); cols <- grDevices::hcl.colors(max(2, nrow(keys)), "Dark 3")[seq_len(nrow(keys))]
    graphics::plot(range(g$L), range(c(0, g[[col]])), type = "n", xlab = "L (baseline weeks)", ylab = ylab); graphics::abline(h = 0, lty = 2)
    for (k in seq_len(nrow(keys))) for (e in unique(g$estimator)) {
      s <- g[match_end == keys$match_end[k] & matching == keys$matching[k] & estimator == e][order(L)]
      graphics::lines(s$L, s[[col]], col = cols[k], lwd = 2, type = "b", pch = if (e == "sc") 16 else 1, lty = if (e == "sc") 1 else 2) }
    ests <- unique(g$estimator)
    graphics::legend("topright", legend = c(paste0(keys$matching, ", match_end=", keys$match_end), .estimator_label(ests)),
                     col = c(cols, rep("grey30", length(ests))), lwd = 2, lty = c(rep(1, nrow(keys)), ifelse(ests == "sc", 1, 2)),
                     pch = c(rep(NA, nrow(keys)), ifelse(ests == "sc", 16, 1)), cex = 0.8, bg = "white", box.col = "grey80")
  }
  output$p_grid_ate <- shiny::renderPlot(.grid_plot("bias", "40-week bias"))
  output$p_grid_first8 <- shiny::renderPlot(.grid_plot("bias_first8", "first-8-week bias"))
  output$t_grid <- shiny::renderTable({ g <- data.table::copy(grid_cells()); g[, estimator := .estimator_label(estimator)]; g }, digits = 4)
  output$dl_grid <- shiny::downloadHandler(filename = function() "dip_grid_results.csv",
                                           content = function(f) data.table::fwrite(grid_res(), f))
}

#' Launch the simulator
#' @param ... Passed to `shiny::shinyApp()`.
#' @export
dip_app <- function(...) shiny::shinyApp(app_ui(), app_server, ...)

utils::globalVariables(c("bias", "bias_first8", "ate", "ate_first8", "heldout_gap", "L", "match_end", "matching", "matched",
                         "bias_pct_baseline", "estimator", "neg_mass"))
