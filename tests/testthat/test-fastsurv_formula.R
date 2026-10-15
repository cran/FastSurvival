test_that("fastsurv_formula reads Surv(time, event) ~ group from data", {
  d <- data.frame(tt = c(3, 1, 2, 5, 4, 6), ee = c(1, 0, 1, 1, 0, 1),
                  arm = factor(c("B", "A", "B", "A", "B", "A"),
                               levels = c("B", "A")),
                  num = c(2, 1, 2, 1, 2, 1))
  f <- fastsurv_formula(Surv(tt, ee) ~ arm, d, group_mode = "required")
  expect_identical(f$time, d$tt)
  expect_identical(f$event, d$ee)
  expect_identical(f$group, d$arm)
  expect_null(f$strata)
  # The control group is the first level of a factor, otherwise the smallest
  # value, and FALSE for a logical group
  expect_identical(f$control, "B")
  expect_identical(fastsurv_formula(Surv(tt, ee) ~ num, d,
                                    group_mode = "required")$control, 1)
  expect_identical(fastsurv_formula(Surv(tt, ee) ~ I(num == 2), d,
                                    group_mode = "required")$control, FALSE)
  # A single group
  f1 <- fastsurv_formula(Surv(tt, ee) ~ 1, d, group_mode = "optional")
  expect_null(f1$group)
  expect_null(f1$control)
})

test_that("fastsurv_formula accepts the forms of Surv() for right censoring", {
  d <- data.frame(tt = c(3, 1, 2, 5), st = c(2, 1, 2, 2), g = c(0, 1, 0, 1))
  f_named <- fastsurv_formula(Surv(event = st == 2, time = tt) ~ g, d,
                              group_mode = "required")
  expect_identical(f_named$time, d$tt)
  expect_identical(f_named$event, d$st == 2)
  f_ns <- fastsurv_formula(survival::Surv(tt, st == 2) ~ g, d,
                           group_mode = "required")
  expect_identical(f_ns$event, c(TRUE, FALSE, TRUE, TRUE))
  # The second argument is the data frame, and without data the variables
  # come from the environment of the formula
  expect_identical(fastsurv_formula(Surv(tt, st == 2) ~ g, event = d,
                                    group_mode = "required")$time, d$tt)
  tt2 <- c(9, 8)
  ee2 <- c(1, 1)
  expect_identical(fastsurv_formula(Surv(tt2, ee2) ~ 1)$time, tt2)
})

test_that("fastsurv_formula combines several strata without merging them", {
  # Pasting the labels with "." would merge ("a.b", "c") and ("a", "b.c");
  # the four distinct pairs must stay four strata (hand count)
  n <- 80
  d <- data.frame(tt = seq_len(n), ee = 1, g = rep(0:1, n / 2),
                  s1 = rep(c("a.b", "a"), each = n / 2),
                  s2 = rep(c("c", "b.c"), n / 2))
  f2 <- fastsurv_formula(Surv(tt, ee) ~ g + strata(s1, s2), d,
                         group_mode = "required", strata_ok = TRUE)
  expect_length(unique(f2$strata), 4L)
  expect_length(unique(paste(d$s1, d$s2, sep = ".")), 3L)
  f2b <- fastsurv_formula(Surv(tt, ee) ~ strata(s1) + g + strata(s2), d,
                          group_mode = "required", strata_ok = TRUE)
  expect_identical(f2b$strata, f2$strata)
  expect_identical(f2b$group, d$g)
  # One stratification variable is returned as it is
  f1 <- fastsurv_formula(Surv(tt, ee) ~ g + strata(s1), d,
                         group_mode = "required", strata_ok = TRUE)
  expect_identical(f1$strata, d$s1)
})

test_that("fastsurv_formula rejects formulas outside the supported forms", {
  d <- data.frame(tt = 1:4, ee = c(1, 0, 1, 1), g = c(0, 1, 0, 1),
                  h = c(1, 1, 2, 2))
  req <- function(f, ...) {
    fastsurv_formula(f, d, group_mode = "required", ...)
  }
  expect_error(fastsurv_formula(~ g, d), "response")
  expect_error(req(tt ~ g), "Surv\\(time, event\\)")
  expect_error(req(Surv(tt, tt, ee) ~ g), "right-censored")
  expect_error(req(Surv(tt, ee, type = "right") ~ g), "right-censored")
  expect_error(req(Surv(tt, ee) ~ g:h), "interactions")
  expect_error(req(Surv(tt, ee) ~ g * h), "interactions")
  expect_error(req(Surv(tt, ee) ~ g - 1), "interactions")
  expect_error(req(Surv(tt, ee) ~ 0 + g), "constants")
  expect_error(req(Surv(tt, ee) ~ g + h), "more than one")
  expect_error(req(Surv(tt, ee) ~ .), "not supported")
  expect_error(req(Surv(tt, ee) ~ g + strata(h)), "only by survdiff_fast")
  expect_error(req(Surv(tt, ee) ~ g + strata(h, na.group = TRUE),
                   strata_ok = TRUE), "no named")
  expect_error(req(Surv(tt, ee) ~ 1), "must name the group")
  expect_error(fastsurv_formula(Surv(tt, ee) ~ g, d, group_mode = "none",
                                fun = "survfit_fast()"),
               "survfit_fast\\(\\) is for a single group")
})

test_that("fastsurv_formula checks the data and the other arguments", {
  d <- data.frame(tt = 1:4, ee = c(1, 0, 1, 1), g = c(0, 1, 0, 1))
  expect_error(fastsurv_formula(Surv(tt, ee) ~ g, d, event = d$ee,
                                group_mode = "required"),
               "second argument")
  expect_error(fastsurv_formula(Surv(tt, ee) ~ g, d, group = d$g,
                                group_mode = "required"),
               "right-hand side")
  expect_error(fastsurv_formula(Surv(tt, ee) ~ g, d, strata = d$g,
                                group_mode = "required", strata_ok = TRUE),
               "strata\\(\\) in the formula")
  expect_error(fastsurv_formula(Surv(tt, ee) ~ g, 1:4,
                                group_mode = "required"),
               "data frame or a list")
  d_na <- d
  d_na$g[2] <- NA
  expect_error(fastsurv_formula(Surv(tt, ee) ~ g, d_na,
                                group_mode = "required"),
               "missing values")
  short <- c(0, 1)
  expect_error(fastsurv_formula(Surv(tt, ee) ~ short, d,
                                group_mode = "required"),
               "same length")
})

test_that("fastsurv_formula points a status coded 1 and 2 to status == 2", {
  d <- data.frame(tt = 1:6, st = c(2, 1, 2, 2, 1, 2), g = rep(0:1, 3))
  expect_error(fastsurv_formula(Surv(tt, st) ~ g, d,
                                group_mode = "required"),
               "Surv\\(time, status == 2\\)")
  expect_error(coxph_fast(Surv(tt, st) ~ g, d), "status == 2")
  expect_identical(fastsurv_formula(Surv(tt, st == 2) ~ g, d,
                                    group_mode = "required")$event,
                   d$st == 2)
  # A status of all 1 is a valid 0/1 indicator (all events)
  expect_identical(fastsurv_formula(Surv(tt, rep(1, 6)) ~ g, d,
                                    group_mode = "required")$event,
                   rep(1, 6))
})
