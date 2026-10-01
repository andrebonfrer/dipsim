# app.R -----------------------------------------------------------------------
# Two modes. "Single run": choose a data-generating process and a design, run
# once, and look at the four diagnostics that reveal a dip. "Compare designs":
# hold the DGP fixed and sweep baseline length, held-out weeks and the
# matching set, to see which choices contain the bias. Both call dip_run().
#
# To extend: add a parameter to dip_defaults() and an input here; add a remedy
# by adding a spec option in budgetsc and exposing it in .params_from_input().

`%||%` <- function(a, b) if (is.null(a)) b else a

.params_from_input <- function(input) {
  p <- dip_defaults()
  p$n_treated <- as.integer(input$n_treated); p$n_later <- as.integer(input$n_later); p$n_weeks <- as.integer(input$n_weeks)
  p$phi <- input$phi; p$sigma_u <- input$sigma_u; p$dip_selection <- input$dip_selection; p$type_selection <- input$type_selection
  p$covid <- isTRUE(input$covid)
  p$effects <- list(spend_log = input$eff_spend, arrears_log = input$eff_arrears, arrears_length = input$eff_length,
                    share_FastFood = 0, share_Alcohol = 0, share_Gambling = 0)
  p$L <- as.integer(input$L); p$H <- 39L; p$match_end <- as.integer(input$match_end)
  p$remedies <- input$remedies %||% character(0); p$z_cut <- input$z_cut; p$k_pre <- as.integer(input$k_pre)
  p$matching <- input$matching; p$custom_outcomes <- input$custom_outcomes; p$max_donors <- as.integer(input$max_donors)
  p$backend <- input$backend; p$seed <- as.integer(input$seed)
  p
}

#' Shiny UI
#' @return shiny UI.
#' @export
app_ui <- function() {
  d <- dip_defaults()
  shiny::fluidPage(
    shiny::titlePanel("Ashenfelter-dip simulator (budgetsc synthetic-control stage)"),
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
        shiny::hr(), shiny::h4("Design"),
        shiny::numericInput("L", "L: baseline weeks", d$L, min = 4),
        shiny::p(shiny::em("Post window fixed at 40 weeks (adoption week + 39).")),
        shiny::numericInput("match_end", "Held-out pre-weeks (match_end)", d$match_end, min = 0),
        shiny::helpText("The last match_end weeks before adoption are excluded when the weights are fitted, then used as a check: if the synthetic control tracks them, the weights generalise beyond the weeks they were fitted on. See Notes."),
        shiny::selectInput("matching", "Matching set", c("single", "lean", "wide", "custom"), d$matching),
        shiny::conditionalPanel("input.matching == 'custom'",
          shiny::checkboxGroupInput("custom_outcomes", NULL, budgetsc::match_outcomes_wide(), budgetsc::match_outcomes_lean())),
        shiny::numericInput("max_donors", "max donors (screening)", d$max_donors),
        shiny::selectInput("backend", "SC backend", c("augmultisynth", "stub"), d$backend),
        shiny::numericInput("seed", "seed", d$seed),
        shiny::h4("Remedies"),
        shiny::checkboxGroupInput("remedies", NULL, c("Drop shock-triggered adopters" = "nodip", "Reversion adjustment (phi from donors)" = "reversion",
                                                     "Placebo bias estimate (second SC fit)" = "placebo"), d$remedies),
        shiny::fluidRow(shiny::column(6, shiny::numericInput("z_cut", "no-dip z cut", d$z_cut, step = 0.25)), shiny::column(6, shiny::numericInput("k_pre", "no-dip weeks", d$k_pre))),
        shiny::actionButton("run", "Run once", class = "btn-primary"), shiny::br(), shiny::br(),
        shiny::textOutput("timing")),
      shiny::mainPanel(width = 9, shiny::tabsetPanel(
        shiny::tabPanel("Single run",
          shiny::p(shiny::em("A dip shows as: (1) a latent shock that peaks just before adoption, (2) an early post-adoption bulge in the gap that fades, (3) a non-zero held-out gap. The 41-week ATE minus the true effect is the bias that reaches the paper.")),
          shiny::fluidRow(shiny::column(6, shiny::h5("Latent shock around adoption (what the DGP did)"), shiny::plotOutput("p_shock", height = 260)),
                          shiny::column(6, shiny::h5("Estimated gap by event week"), shiny::selectInput("oc", NULL, d$outcomes), shiny::plotOutput("p_gap", height = 260))),
          shiny::h5("Bias summary (estimate minus true effect; relative to the outcome's pre-adoption level, and to the true effect when non-zero)"),
          shiny::tableOutput("t_summary"),
          shiny::h5("Interpretation"), shiny::uiOutput("interp"),
          shiny::h5("Sample funnel"), shiny::tableOutput("t_funnel")),
        shiny::tabPanel("Compare designs",
          shiny::p("Same DGP as in the sidebar; the design factors below are swept. Each cell is a full simulate + fit."),
          shiny::fluidRow(shiny::column(3, shiny::textInput("grid_L", "L values", "8, 23, 70")),
                          shiny::column(3, shiny::textInput("grid_me", "match_end values", "0, 4")),
                          shiny::column(3, shiny::checkboxGroupInput("grid_matching", "matching sets", c("single", "lean", "wide"), c("lean", "wide"))),
                          shiny::column(3, shiny::numericInput("reps", "replications", 1, min = 1, max = 10))),
          shiny::actionButton("run_grid", "Run grid", class = "btn-primary"), shiny::span("  "), shiny::downloadButton("dl_grid", "Download results"),
          shiny::br(), shiny::br(), shiny::selectInput("grid_oc", "Outcome", d$outcomes),
          shiny::fluidRow(shiny::column(6, shiny::h5("41-week bias vs baseline length"), shiny::plotOutput("p_grid_ate", height = 300)),
                          shiny::column(6, shiny::h5("First-8-week bias vs baseline length"), shiny::plotOutput("p_grid_first8", height = 300))),
          shiny::h5("Reading"), shiny::uiOutput("grid_reading"),
          shiny::h5("Cell means"), shiny::tableOutput("t_grid")),
        shiny::tabPanel("Notes",
          shiny::h4("Glossary"),
          shiny::p("The simulator reproduces, on invented customers, the estimation used in the paper, and then checks it against a truth that is known because it was simulated."),
          shiny::tags$dl(
            shiny::tags$dt("Synthetic control"), shiny::tags$dd("For each adopter, a weighted average of later adopters (the donors) chosen so that, before adoption, the average tracks the adopter's own outcomes as closely as possible. After adoption it stands in for what the adopter would have done without the budget. The weights are fitted on several outcomes at once (the matching set)."),
            shiny::tags$dt("Gap"), shiny::tags$dd("The adopter's outcome minus the synthetic control's, week by week, after removing the average difference over the matching window (an intercept). Averaged over adopters it is the event-time plot; averaged over the post window it is the estimated effect."),
            shiny::tags$dt("Baseline weeks (L)"), shiny::tags$dd("How many pre-adoption weeks the weights are fitted on. The paper uses 70 for public-launch adopters and 23 for pilot adopters."),
            shiny::tags$dt("Post window"), shiny::tags$dd("The 40 weeks from adoption (week 0) to week 39 over which the effect is averaged."),
            shiny::tags$dt("Held-out pre-weeks (match_end)"), shiny::tags$dd(paste("Weeks immediately before adoption that are deliberately NOT used when fitting the weights. They are set aside so they can act as a test:",
              "the weights were chosen without seeing them, so if the synthetic control still tracks the adopter through those weeks, it is capturing the adopter's underlying pattern rather than fitting noise.",
              "The average gap over those weeks is the held-out gap; near zero is reassuring, a clear positive or negative value means the counterfactual already diverges before adoption.",
              "Holding out weeks also changes where a pre-adoption spike sits: inside the matching window the weights partly reproduce it, outside they do not, and that changes the sign of the bias (see below).")),
            shiny::tags$dt("Transitory shock, persistence (phi)"), shiny::tags$dd("Each customer's outcomes carry a temporary disturbance that decays by the factor phi each week. phi = 0 means a spike is gone next week; phi = 0.9 means it unwinds over months."),
            shiny::tags$dt("Selection on the shock (the Ashenfelter dip)"), shiny::tags$dd("Customers are made more likely to adopt in the weeks after a positive shock, as a real customer might set a budget after an expensive month. The post-adoption improvement that follows is partly the shock unwinding, not the budget. This is the mechanism named after Ashenfelter's observation that trainees' earnings dipped just before they enrolled."),
            shiny::tags$dt("Selection on stress type"), shiny::tags$dd("Adoption also depends on a permanent customer trait that raises arrears. The synthetic control handles this kind of selection; it is included so the two kinds can be compared."),
            shiny::tags$dt("True effect"), shiny::tags$dd("What adoption actually did in the simulation. Zero by default, so every estimated effect is error."),
            shiny::tags$dt("Bias"), shiny::tags$dd("Estimated effect minus true effect, over the 40-week window. Reported in the outcome's units, as a percentage of the adopters' pre-adoption level of that outcome, and (when the true effect is not zero) as a percentage of the true effect."),
            shiny::tags$dt("First-8-week bias"), shiny::tags$dd("The same quantity over weeks 0 to 7 only. A dip shows up here first: a shock unwinding after adoption makes the early gap large while the 40-week average dilutes it."),
            shiny::tags$dt("Latent shock profile"), shiny::tags$dd("The simulated shock itself, averaged over adopters by week relative to adoption. Visible only because the data are simulated; it shows what the estimator is up against. A peak just before week 0 is the dip.")),
          shiny::h5("What each control does"),
          shiny::tags$ul(
            shiny::tags$li(shiny::strong("phi"), ": persistence of the idiosyncratic shock u. With phi = 0 a spike is gone next week and cannot bias post-adoption outcomes; with phi = 0.9 it unwinds over months."),
            shiny::tags$li(shiny::strong("Selection on last week's shock"), ": the hazard of adopting in a week is proportional to exp(k * u_{t-1}). This is the Ashenfelter mechanism; 0 switches it off."),
            shiny::tags$li(shiny::strong("Selection on stress type"), ": adoption depends on a persistent customer trait that also raises arrears. Handled by the synthetic control; kept here so the two channels can be told apart."),
            shiny::tags$li(shiny::strong("L"), ": weeks matched before adoption. Dilutes the spike's leverage on the weights and the outcome intercept."),
            shiny::tags$li(shiny::strong("match_end"), ": the last k pre-weeks are held out of matching; their gap is the held-out diagnostic. With the spike outside the matching window the early bias is positive (the spike unwinds); with it inside, and L short, the counterfactual carries the spike forward and the bias turns negative."),
            shiny::tags$li(shiny::strong("Matching set"), ": how many outcomes the weights match on. Pins down latent loadings; does not dilute an idiosyncratic shock."),
            shiny::tags$li(shiny::strong("True effects"), ": set to zero, everything estimated is bias; set to the paper's magnitudes to see attenuation or inflation."),
            shiny::tags$li(shiny::strong("Remedies"), ": observable-based corrections that could be applied to real data, judged by whether their bias is closer to zero than the raw one. Findings so far (phi 0.9, strong selection, true effect 0, three replications): the reversion adjustment removes most of the bias when the spike is held out of the matching window (match_end > 0: L=23 bias 39 -> 4; L=70 30 -> -11), because then the unwinding is in the gap and can be subtracted; with the spike inside the window it double-counts and makes the bias worse. The placebo bias estimate works at L=23 with hold-out (39 -> -1) and is unreliable otherwise. Dropping shock-triggered adopters makes the bias worse in every design tested: the spike is in the weights, not only in the dropped units. So the usable combination is hold-out matching plus the reversion adjustment.")),
          shiny::h5("Extending"),
          shiny::p("Add a DGP parameter in dipsim::dip_defaults() and budgetsc::sim_panel(); add a remedy as a budgetsc spec option and expose it in .params_from_input(). The grid tab sweeps whatever dip_grid() is given."))
      ))
    )
  )
}

#' Shiny server
#' @param input,output,session Shiny objects.
#' @export
app_server <- function(input, output, session) {
  res <- shiny::eventReactive(input$run, {
    shiny::withProgress(message = "simulating and fitting", value = 0.3, {
      r <- dip_run(.params_from_input(input)); shiny::incProgress(0.7); r })
  })
  output$timing <- shiny::renderText({ r <- res(); sprintf("last run: %.1f s (%s backend)", r$seconds, r$params$backend) })
  output$p_shock <- shiny::renderPlot({ r <- res(); sp <- r$shock_profile
    graphics::plot(sp$tau_k, sp$mean_u, type = "l", lwd = 2, xlab = "weeks since adoption", ylab = "mean latent shock u (treated)")
    graphics::abline(h = 0, lty = 2); graphics::abline(v = 0, col = "red") })
  output$p_gap <- shiny::renderPlot({ r <- res(); e <- r$event_time[outcome == input$oc]; tv <- r$summary[outcome == input$oc, true_effect]
    graphics::plot(e$tau_k, e$mean_gap, type = "b", pch = 16, cex = 0.7, xlab = "weeks since adoption", ylab = "treated - synthetic", main = input$oc)
    graphics::abline(h = 0, lty = 2); graphics::abline(v = 0, col = "red")
    if (length(tv) && is.finite(tv) && tv != 0) graphics::abline(h = tv, col = "blue", lty = 3)
    me <- r$params$match_end; if (me > 0) graphics::rect(-me - 0.5, -1e9, -0.5, 1e9, col = grDevices::adjustcolor("orange", 0.15), border = NA) })
  output$t_summary <- shiny::renderTable({ s <- data.table::copy(res()$summary)
    s[, bias_first8 := ate_first8 - true_effect]
    keep <- intersect(c("outcome", "matched", "pre_mean", "true_effect", "ate", "bias", "bias_pct_baseline", "bias_pct_effect", "bias_first8", "heldout_gap",
                        "bias_nodip", "n_dropped_nodip", "bias_reversion", "phi_hat", "bias_placebo_adj"), names(s))
    s[, keep, with = FALSE] }, digits = 3)
  output$interp <- shiny::renderUI({ txt <- interpret_run(res(), input$oc); shiny::tagList(lapply(txt[nzchar(txt)], shiny::p)) })
  output$t_funnel <- shiny::renderTable(res()$funnel)

  grid_res <- shiny::eventReactive(input$run_grid, {
    Ls <- as.integer(trimws(strsplit(input$grid_L, ",")[[1]])); mes <- as.integer(trimws(strsplit(input$grid_me, ",")[[1]]))
    g <- expand.grid(L = Ls, match_end = mes, matching = input$grid_matching, stringsAsFactors = FALSE)
    g <- g[g$match_end < g$L, ]
    shiny::withProgress(message = "running grid", value = 0, {
      dip_grid(.params_from_input(input), g, reps = as.integer(input$reps),
               progress = function(i, n) shiny::setProgress(i / n, detail = sprintf("%d / %d", i, n))) })
  })
  grid_cells <- shiny::reactive({ g <- grid_res(); shiny::req(nrow(g) > 0)
    g[, bias := ate - true_effect][, bias_first8 := ate_first8 - true_effect]
    cols <- intersect(c("bias_nodip", "bias_reversion", "bias_placebo_adj"), names(g))
    g[outcome == input$grid_oc, c(list(bias = mean(bias), bias_pct_baseline = mean(bias_pct_baseline), bias_first8 = mean(bias_first8),
                                       heldout_gap = mean(heldout_gap, na.rm = TRUE), n = .N),
                                  lapply(.SD, mean, na.rm = TRUE)),
      by = list(L, match_end, matching), .SDcols = cols][order(matching, match_end, L)] })
  output$grid_reading <- shiny::renderUI({ g <- grid_cells(); shiny::req(nrow(g) > 1)
    best <- g[which.min(abs(bias))]; worst <- g[which.max(abs(bias))]
    shiny::p(sprintf(paste("Across %d design cells for %s, the 40-week bias is smallest with L = %d, match_end = %d, %s matching (%.3f, %.1f%% of the pre-adoption level)",
      "and largest with L = %d, match_end = %d, %s matching (%.3f, %.1f%%). Bias in the first eight weeks ranges from %.3f to %.3f: where it exceeds the",
      "40-week bias by a wide margin the dip is being averaged down rather than removed."),
      nrow(g), input$grid_oc, best$L, best$match_end, best$matching, best$bias, best$bias_pct_baseline, worst$L, worst$match_end, worst$matching, worst$bias, worst$bias_pct_baseline,
      min(g$bias_first8), max(g$bias_first8))) })
  .grid_plot <- function(col, ylab) {
    g <- grid_cells(); keys <- unique(g[, list(match_end, matching)]); cols <- grDevices::rainbow(nrow(keys))
    graphics::plot(range(g$L), range(c(0, g[[col]])), type = "n", xlab = "L (baseline weeks)", ylab = ylab); graphics::abline(h = 0, lty = 2)
    for (k in seq_len(nrow(keys))) { s <- g[match_end == keys$match_end[k] & matching == keys$matching[k]][order(L)]
      graphics::lines(s$L, s[[col]], col = cols[k], lwd = 2, type = "b", pch = 16) }
    graphics::legend("topright", legend = paste0(keys$matching, ", match_end=", keys$match_end), col = cols, lwd = 2, cex = 0.8, bty = "n")
  }
  output$p_grid_ate <- shiny::renderPlot(.grid_plot("bias", "41-week bias"))
  output$p_grid_first8 <- shiny::renderPlot(.grid_plot("bias_first8", "first-8-week bias"))
  output$t_grid <- shiny::renderTable(grid_cells(), digits = 4)
  output$dl_grid <- shiny::downloadHandler(filename = function() "dip_grid_results.csv",
                                           content = function(f) data.table::fwrite(grid_res(), f))
}

#' Launch the simulator
#' @param ... Passed to `shiny::shinyApp()`.
#' @export
dip_app <- function(...) shiny::shinyApp(app_ui(), app_server, ...)

utils::globalVariables(c("bias", "bias_first8", "ate", "ate_first8", "heldout_gap", "L", "match_end", "matching", "matched", "bias_pct_baseline"))
