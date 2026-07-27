## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  eval     = FALSE
)

## ----load---------------------------------------------------------------------
# library(FastSurvival)
# library(survival)
# library(microbenchmark)

## ----data---------------------------------------------------------------------
# dataset <- simdata_fast(
#   nsim     = 1,
#   n        = 500,
#   a.time   = c(0, 12.5),
#   a.rate   = 40,
#   e.median = list(5.811, 4.3),
#   seed     = 1
# )
# 
# # Sort once and reuse, the intended pattern for the pre-sorted fast path.
# ord <- order(dataset$tte)
# t_s <- dataset$tte[ord]
# e_s <- dataset$event[ord]
# g_s <- dataset$group[ord]
# 
# # Control is group 1, treatment is group 2.
# arm <- as.integer(dataset$group == 2)
# 
# # Restriction horizon within both arms' follow-up.
# tau <- floor(min(tapply(t_s, g_s, max)))
# 
# # Factor arm for the nphRCT reference used in the rmw_fast benchmark.
# df_rmw <- data.frame(
#   tte   = dataset$tte,
#   event = dataset$event,
#   arm   = factor(ifelse(dataset$group == 1, "control", "treatment"),
#                  levels = c("control", "treatment"))
# )

## ----bench-survfit------------------------------------------------------------
# microbenchmark(
#   fast = survfit_fast(t_s, e_s, t_eval = tau, presorted = TRUE),
#   ref  = summary(survfit(Surv(tte, event) ~ 1, data = dataset), times = tau),
#   times = 1000
# )

## ----bench-survdiff-----------------------------------------------------------
# microbenchmark(
#   fast = survdiff_fast(t_s, e_s, g_s, control = 1, side = 1, presorted = TRUE),
#   ref  = survdiff(Surv(tte, event) ~ group, data = dataset),
#   times = 1000
# )

## ----bench-coxph--------------------------------------------------------------
# microbenchmark(
#   fast = coxph_fast(t_s, e_s, g_s, control = 1, side = 1, presorted = TRUE),
#   ref  = coxph(Surv(tte, event) ~ I(group == 2), data = dataset),
#   times = 1000
# )

## ----bench-rmst---------------------------------------------------------------
# microbenchmark(
#   fast = rmst_fast(t_s, e_s, g_s, control = 1, tau = tau, side = 1,
#                    presorted = TRUE),
#   ref  = survRM2::rmst2(time = dataset$tte, status = dataset$event,
#                         arm = arm, tau = tau),
#   times = 1000
# )

## ----bench-wlr----------------------------------------------------------------
# microbenchmark(
#   fast = survdiff_fast(t_s, e_s, g_s, control = 1, side = 1,
#                        weight = "fh", rho = 0, gamma = 1, presorted = TRUE),
#   ref  = nph::logrank.test(dataset$tte, dataset$event, dataset$group,
#                            rho = 0, gamma = 1),
#   times = 1000
# )

## ----bench-milestone----------------------------------------------------------
# microbenchmark(
#   fast = milestone_fast(t_s, e_s, g_s, control = 1, tau = tau, side = 1,
#                         presorted = TRUE),
#   ref  = summary(survfit(Surv(tte, event) ~ group, data = dataset),
#                  times = tau),
#   times = 1000
# )

## ----bench-medsurv------------------------------------------------------------
# microbenchmark(
#   fast = medsurv_fast(t_s, e_s, g_s, control = 1, side = 1,
#                       method = "nph", presorted = TRUE),
#   ref  = nph::nphparams(time = dataset$tte, event = dataset$event,
#                         group = as.integer(dataset$group == 2),
#                         param_type = "Q", param_par = 0.5),
#   times = 1000
# )

## ----bench-maxcombo-----------------------------------------------------------
# microbenchmark(
#   fast = maxcombo_fast(t_s, e_s, g_s, control = 1, side = 1,
#                        rho = c(0, 0, 1), gamma = c(0, 1, 0), presorted = TRUE),
#   ref  = nph::logrank.maxtest(dataset$tte, dataset$event,
#                               as.integer(dataset$group == 2)),
#   times = 1000
# )

## ----bench-rmw----------------------------------------------------------------
# microbenchmark(
#   fast = rmw_fast(t_s, e_s, g_s, control = 1, side = 1, s_star = 0.5,
#                   presorted = TRUE),
#   ref  = {
#     nphRCT::wlrt(Surv(tte, event) ~ arm, data = df_rmw,
#                  method = "mw", s_star = 1)
#     nphRCT::wlrt(Surv(tte, event) ~ arm, data = df_rmw,
#                  method = "mw", s_star = 0.5)
#   },
#   times = 1000
# )

## ----bench-ahsw---------------------------------------------------------------
# microbenchmark(
#   fast = ahsw_fast(t_s, e_s, g_s, control = 1, tau = tau, side = 1,
#                    presorted = TRUE),
#   ref  = survAH::ah2(time = dataset$tte, status = dataset$event,
#                      arm = arm, tau = tau),
#   times = 1000
# )

