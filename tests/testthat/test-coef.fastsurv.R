test_that("coef.fastsurv returns the estimates named by term", {
  set.seed(91)
  n <- 120
  d <- data.frame(tt = round(rexp(n, 0.1), 1), ee = rbinom(n, 1, 0.8),
                  g = rep(c(0, 1), n / 2))

  cx <- coxph_fast(d$tt, d$ee, d$g, control = 0)
  expect_identical(coef(cx), c(coef = cx[["coef"]]))

  rm <- rmst_fast(d$tt, d$ee, d$g, control = 0, tau = 5)
  cf <- coef(rm)
  expect_identical(names(cf), c("rmst.ctrl", "rmst.trt", "diff", "log.ratio"))
  expect_equal(unname(cf[1:3]),
               c(rm[["rmst.ctrl"]], rm[["rmst.trt"]], rm[["diff"]]))
  expect_equal(unname(cf["log.ratio"]), log(rm[["ratio"]]))

  ms <- milestone_fast(d$tt, d$ee, d$g, control = 0, tau = 5)
  expect_identical(coef(ms), c(surv.control = ms$surv[["control"]],
                               surv.treatment = ms$surv[["treatment"]],
                               diff = ms$diff))

  ar <- ahr_fast(d$tt, d$ee, d$g, control = 0, tau = 5)
  expect_equal(coef(ar), c(log.ahr = log(ar$ahr),
                           theta = unname(ar$theta[2L])))

  # An estimate that is not available keeps its place
  na_cx <- coxph_fast(c(1, 2, 3, 4), c(0, 0, 0, 0), c(0, 0, 1, 1),
                      control = 0)
  expect_identical(coef(na_cx), c(coef = NA_real_))
})

test_that("coef.fastsurv gives an error for a test without an estimate", {
  set.seed(92)
  n <- 80
  d <- data.frame(tt = rexp(n, 0.1), ee = 1, g = rep(c(0, 1), n / 2))
  msg <- "test without an estimate"
  expect_error(coef(survdiff_fast(d$tt, d$ee, d$g, control = 0)), msg)
  expect_error(coef(maxcombo_fast(d$tt, d$ee, d$g, control = 0,
                                  rho = c(0, 0), gamma = c(0, 1))), msg)
  expect_error(coef(rmw_fast(d$tt, d$ee, d$g, control = 0)), msg)
})
