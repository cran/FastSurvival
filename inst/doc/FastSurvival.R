## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.align = "center"
)

## ----load-pkg-----------------------------------------------------------------
library(FastSurvival)
library(survival)

## ----data---------------------------------------------------------------------
ord <- order(ovarian$futime)
t_s <- ovarian$futime[ord]
e_s <- ovarian$fustat[ord]
g_s <- ovarian$rx[ord]

## ----survfit-agreement--------------------------------------------------------
# FastSurvival
fit_fast <- survfit_fast(t_s, e_s, t_eval = 500, conf.type = "log")
fit_fast

# survival
fit_ref <- survfit(Surv(futime, fustat) ~ 1, data = ovarian, conf.type = "log")
summary(fit_ref, times = 500)

## ----survdiff-agreement-------------------------------------------------------
# FastSurvival (two-sided chi-square statistic)
sd_fast <- survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
                         control = 1, side = 2)
sd_fast

# survival
sd_ref <- survdiff(Surv(futime, fustat) ~ rx, data = ovarian)
sd_ref

cat("survdiff_fast chi-square:", as.numeric(sd_fast), "\n")
cat("survdiff      chi-square:", sd_ref$chisq,       "\n")

## ----coxph-agreement----------------------------------------------------------
# FastSurvival
cx_fast <- coxph_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1)
cx_fast

# survival
cx_ref <- summary(coxph(Surv(futime, fustat) ~ rx, data = ovarian))
cx_ref$coefficients
cx_ref$conf.int

cat("coxph_fast HR :", cx_fast["exp(coef)"], "\n")
cat("coxph      HR :", cx_ref$coefficients[, "exp(coef)"], "\n")

## ----benchmark, message = FALSE-----------------------------------------------
library(microbenchmark)

ord <- order(ovarian$futime)
t_s <- ovarian$futime[ord]
e_s <- ovarian$fustat[ord]
g_s <- ovarian$rx[ord]

bm <- microbenchmark(
  survfit_fast  = survfit_fast(t_s, e_s, t_eval = 500),
  survfit       = summary(survfit(Surv(futime, fustat) ~ 1, data = ovarian),
                          times = 500),
  survdiff_fast = survdiff_fast(t_s, e_s, g_s, control = 1, side = 2,
                                presorted = TRUE),
  survdiff      = survdiff(Surv(futime, fustat) ~ rx, data = ovarian),
  coxph_fast    = coxph_fast(t_s, e_s, g_s, control = 1, presorted = TRUE),
  coxph         = coxph(Surv(futime, fustat) ~ rx, data = ovarian),
  times = 100
)
print(bm)

## ----simulate-----------------------------------------------------------------
set.seed(42)
df <- simdata_fast(
  nsim     = 100,
  n        = c(100, 100),
  a.time   = c(0, 12),
  a.rate   = 200 / 12,
  e.median = list(18, 24),
  d.hazard = list(0.01, 0.01),
  seed     = 42
)
head(df)

## ----analyse------------------------------------------------------------------
nsim   <- 100L
t_land <- 18

km_mat <- matrix(NA_real_, nrow = nsim, ncol = 4L,
                 dimnames = list(NULL,
                                 c("surv", "std.err", "lower", "upper")))
ld_vec <- numeric(nsim)   # log-rank chi-square
hr_mat <- matrix(NA_real_, nrow = nsim, ncol = 5L,
                 dimnames = list(NULL,
                                 c("coef", "exp(coef)", "se(coef)",
                                   "lower .95", "upper .95")))

for (s in seq_len(nsim)) {
  d <- df[df$sim == s, ]

  # Treatment-group KM at t = 18 (sort once)
  d_trt <- d[d$group == 2L, ]
  ord_t <- order(d_trt$tte)
  km_mat[s, ] <- unclass(
    survfit_fast(d_trt$tte[ord_t], d_trt$event[ord_t], t_eval = t_land)
  )

  # Pooled sort for log-rank and Cox
  ord_p <- order(d$tte)
  t_p   <- d$tte[ord_p]
  e_p   <- d$event[ord_p]
  g_p   <- d$group[ord_p]

  ld_vec[s] <- as.numeric(
    survdiff_fast(t_p, e_p, g_p, control = 1L, side = 2L, presorted = TRUE)
  )
  hr_mat[s, ] <- unclass(
    coxph_fast(t_p, e_p, g_p, control = 1L, presorted = TRUE)
  )
}

## ----summarize----------------------------------------------------------------
# Landmark survival in the treatment group at t = 18
km_summary <- c(
  mean_surv    = mean(km_mat[, "surv"]),
  sd_surv      = sd(km_mat[, "surv"]),
  mean_std_err = mean(km_mat[, "std.err"])
)
km_summary

# Log-rank test: empirical power at alpha = 0.05 (two-sided)
p_logrank   <- pchisq(ld_vec, df = 1L, lower.tail = FALSE)
power_empir <- mean(p_logrank < 0.05)
cat(sprintf("Empirical power (log-rank, two-sided, alpha = 0.05): %.2f\n",
            power_empir))

# Hazard ratio summary
hr_summary <- c(
  mean_HR     = mean(hr_mat[, "exp(coef)"]),
  sd_HR       = sd(hr_mat[, "exp(coef)"]),
  mean_log_HR = mean(hr_mat[, "coef"]),
  sd_log_HR   = sd(hr_mat[, "coef"]),
  cover_95    = mean(hr_mat[, "lower .95"] <= (18 / 24) &
                     hr_mat[, "upper .95"] >= (18 / 24))
)
hr_summary

