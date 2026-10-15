## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>"
)

## ----load-pkg-----------------------------------------------------------------
library(FastSurvival)

## ----quick, eval = requireNamespace("survival", quietly = TRUE)---------------
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

## ----formula, eval = requireNamespace("survival", quietly = TRUE)-------------
fit <- coxph_fast(Surv(futime, fustat) ~ rx, data = ovarian, side = 1)
as.data.frame(fit)
exp(confint(fit))

as.data.frame(rmst_fast(Surv(futime, fustat) ~ rx, data = ovarian,
                        tau = 500))

## ----km-plot, eval = requireNamespace("survival", quietly = TRUE), fig.width = 7, fig.height = 5----
km <- kmcurve_fast(Surv(futime, fustat) ~ rx, data = ovarian)
plot(km, risk_table = TRUE, conf.style = "band")

