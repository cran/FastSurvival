## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  eval     = FALSE
)

## ----load---------------------------------------------------------------------
# library(FastSurvival)

## ----simulate-----------------------------------------------------------------
# df <- simdata_fast(
#   nsim     = 10000,
#   n        = c(241, 241),
#   a.time   = c(0, 23),
#   a.prop   = 1,
#   e.median = list(12.9, 9.0),
#   d.hazard = list(-log(1 - 0.05) / 12, -log(1 - 0.05) / 12),
#   seed     = 1
# )

## ----analyze------------------------------------------------------------------
# res <- analysis_fast(
#   df, control = 2,
#   event.looks = c(252, 336),
#   stat = c("logrank", "coxph"), side = 2
# )

## ----boundaries---------------------------------------------------------------
# library(gsDesign)
# 
# gsd <- gsDesign(
#   k      = 2,
#   timing = c(252, 336) / 336,
#   alpha  = 0.025,
#   beta   = 0.1,
#   sfu    = sfLDOF,
#   test.type = 1
# )
# 
# spend_alpha <- 2 * pnorm(gsd$upper$bound, lower.tail = FALSE)
# spend_alpha

## ----summary------------------------------------------------------------------
# simsummary_fast(
#   res,
#   p.col     = "logrank.p",
#   alpha     = spend_alpha,
#   direction = "lower"
# )

## ----rpact--------------------------------------------------------------------
# library(rpact)
# 
# design <- getDesignGroupSequential(
#   kMax = 2,
#   alpha = 0.05,
#   beta  = 0.1,
#   sided = 2,
#   typeOfDesign = "asOF",
#   informationRates = c(252, 336) / 336
# )
# 
# results <- getPowerSurvival(
#   design,
#   maxNumberOfEvents   = 336,
#   median1             = 12.9,
#   median2             = 9.0,
#   maxNumberOfSubjects = 482,
#   accrualTime         = c(0, 23),
#   dropoutRate1        = -log(1 - 0.05),
#   dropoutRate2        = -log(1 - 0.05),
#   allocationRatioPlanned = 1
# )

## ----nph----------------------------------------------------------------------
# df_delay <- simdata_fast(
#   nsim     = 10000,
#   n        = c(241, 241),
#   a.time   = c(0, 23),
#   a.prop   = 1,
#   e.hazard = list(c(0.077, 0.045), c(0.077, 0.077)),
#   e.time   = c(0, 3, Inf),
#   d.hazard = list(-log(1 - 0.05) / 12, -log(1 - 0.05) / 12),
#   seed     = 1
# )
# 
# res_delay <- analysis_fast(
#   df_delay, control = 2,
#   event.looks = c(252, 336),
#   stat = "maxcombo", side = 2
# )

