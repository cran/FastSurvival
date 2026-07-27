## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>"
)

## ----load-pkg-----------------------------------------------------------------
library(FastSurvival)

## ----quick--------------------------------------------------------------------
library(survival)

# Single-time-point Kaplan-Meier estimate
ord <- order(ovarian$futime)
survfit_fast(ovarian$futime[ord], ovarian$fustat[ord],
             t_eval = 500, conf.type = "log")

# Log-rank test
survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
              control = 1, side = 1)

# Hazard ratio via the Pike-Halley Estimator
coxph_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
           control = 1, side = 1)

