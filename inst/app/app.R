# Deployable launcher: rsconnect::deployApp(system.file("app", package = "dipsim"))
# Matrix is attached because augMultiSynth < 0.3.6 needs it on the search path
# (fit_weights_fw crossprod issue); harmless once 0.3.6 is installed.
suppressPackageStartupMessages({ library(Matrix); library(dipsim) })
dip_app()
