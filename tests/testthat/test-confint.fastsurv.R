test_that("confint.fastsurv returns the limits in the layout of stats", {
  set.seed(93)
  n <- 120
  d <- data.frame(tt = round(rexp(n, 0.1), 1), ee = rbinom(n, 1, 0.8),
                  g = rep(c(0, 1), n / 2))

  cx <- coxph_fast(d$tt, d$ee, d$g, control = 0)
  ci <- confint(cx)
  # The column labels of stats::confint() for an lm fit
  lm_fit <- lm(dist ~ speed, data = cars)
  expect_identical(colnames(ci), colnames(confint(lm_fit)))
  expect_identical(rownames(ci), "coef")
  expect_equal(exp(unname(ci[1L, ])),
               unname(c(cx[["lower .95"]], cx[["upper .95"]])))

  # The level is that of the result
  cx90 <- coxph_fast(d$tt, d$ee, d$g, control = 0, conf.level = 0.9)
  expect_identical(colnames(confint(cx90)),
                   colnames(confint(lm_fit, level = 0.9)))
  expect_identical(confint(cx90, level = 0.9), confint(cx90))
  expect_error(confint(cx90, level = 0.95), "conf.level = 0.9")
  expect_error(confint(cx90, level = "a"), "conf.level = 0.9")
})

test_that("confint.fastsurv selects estimates with parm", {
  set.seed(94)
  n <- 120
  d <- data.frame(tt = round(rexp(n, 0.1), 1), ee = rbinom(n, 1, 0.8),
                  g = rep(c(0, 1), n / 2))
  rm <- rmst_fast(d$tt, d$ee, d$g, control = 0, tau = 5)
  ci <- confint(rm)
  expect_identical(rownames(ci),
                   c("rmst.ctrl", "rmst.trt", "diff", "log.ratio"))
  # rmst_fast() gives no interval for the RMST of each group
  expect_true(all(is.na(ci[1:2, ])))
  expect_equal(unname(ci["diff", ]),
               c(rm[["diff.lower"]], rm[["diff.upper"]]))
  expect_equal(exp(unname(ci["log.ratio", ])),
               c(rm[["ratio.lower"]], rm[["ratio.upper"]]))
  expect_identical(confint(rm, parm = "diff"), ci["diff", , drop = FALSE])
  expect_identical(confint(rm, parm = 3:4), ci[3:4, , drop = FALSE])
  expect_error(confint(rm, parm = c("diff", "nope")), "unknown 'parm': nope")
  expect_error(confint(rm, parm = 5), "between 1 and 4")
  expect_error(confint(rm, parm = 2.5), "whole numbers")
  expect_error(confint(rm, parm = TRUE), "character or numeric")
  expect_error(confint(survdiff_fast(d$tt, d$ee, d$g, control = 0)),
               "test without an estimate")
})
