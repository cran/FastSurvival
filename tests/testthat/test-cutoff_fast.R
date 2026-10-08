make_cut_data <- function(nsim = 30, seed = 11) {
  simdata_fast(nsim = nsim, n = c(100, 100), a.time = c(0, 12),
               a.rate = 200 / 12, e.median = list(10, 14),
               d.hazard = 0.01, seed = seed)
}

# Reference: calendar time of the d-th counted event of one simulation,
# computed directly in base R.
ref_event_time <- function(d_sim, d, subset = NULL) {
  keep <- d_sim$event == 1
  if (!is.null(subset)) keep <- keep & subset
  cal <- sort(d_sim$accrual_time[keep] + d_sim$tte[keep])
  if (d <= length(cal)) cal[d] else NA_real_
}

test_that("cutoff_fast: event looks equal the calendar time of the target event", {
  df  <- make_cut_data()
  cut <- cutoff_fast(df, event.looks = c(60, 120))
  expect_equal(dim(cut), c(30L, 2L))
  expect_equal(rownames(cut), as.character(1:30))
  expect_equal(attr(cut, "look.value"), c(60, 120))
  for (s in 1:30) {
    d_s <- df[df$sim == s, ]
    expect_equal(unname(cut[s, ]),
                 c(ref_event_time(d_s, 60), ref_event_time(d_s, 120)))
  }
})

test_that("cutoff_fast: event looks reproduce analysis_fast event.looks", {
  df    <- make_cut_data()
  res_e <- analysis_fast(df, control = 1, event.looks = c(60, 120),
                         stat = c("logrank", "coxph", "rmst"), tau = 12)
  cut   <- cutoff_fast(df, event.looks = c(60, 120))
  res_c <- analysis_fast(df, control = 1, cutoff.looks = cut,
                         stat = c("logrank", "coxph", "rmst"), tau = 12)
  expect_equal(res_c, res_e)
})

test_that("cutoff_fast: an unmet target is NA unless capped by max.time", {
  df  <- make_cut_data(nsim = 10)
  cut <- cutoff_fast(df, event.looks = 10000)
  expect_true(all(is.na(cut)))
  cap <- cutoff_fast(df, event.looks = 10000, max.time = 40)
  expect_equal(unname(cap[, 1]), rep(40, 10))

  res <- analysis_fast(df, control = 1, cutoff.looks = cut)
  ref <- analysis_fast(df, control = 1, event.looks = 10000)
  expect_false(any(res$reached))
  expect_equal(res$reached, ref$reached)
  expect_equal(res$logrank.z, ref$logrank.z)
})

test_that("cutoff_fast: combined rules follow min(max(...), max.time)", {
  df <- make_cut_data()
  te <- unname(cutoff_fast(df, event.looks = 120)[, 1])

  # 120 events or month 20, whichever comes first.
  first <- cutoff_fast(df, event.looks = 120, max.time = 20)
  expect_equal(unname(first[, 1]), pmin(te, 20))
  # 120 events, but not before month 30.
  later <- cutoff_fast(df, event.looks = 120, time.looks = 30)
  expect_equal(unname(later[, 1]), pmax(te, 30))
  expect_equal(attr(later, "look.value"), 120)

  # A minimum gap after the previous look.
  two <- cutoff_fast(df, event.looks = c(60, 70), min.gap = c(NA, 6))
  t60 <- unname(cutoff_fast(df, event.looks = 60)[, 1])
  t70 <- unname(cutoff_fast(df, event.looks = 70)[, 1])
  expect_equal(unname(two[, 1]), t60)
  expect_equal(unname(two[, 2]), pmax(t70, t60 + 6))

  # A minimum follow-up after the 150th enrolled subject.
  enr <- cutoff_fast(df, min.enrolled = 150, min.followup = 12)
  ref <- vapply(split(df$accrual_time, df$sim),
                function(a) sort(a)[150] + 12, numeric(1))
  expect_equal(unname(enr[, 1]), unname(ref))
  expect_true(is.na(attr(enr, "look.value")))

  # A calendar look alone is constant across simulations.
  cal <- cutoff_fast(df, time.looks = c(15, 30))
  expect_equal(unname(cal[, 2]), rep(30, 30))
  expect_equal(attr(cal, "look.value"), c(15, 30))
})

test_that("cutoff_fast: event.subset counts only the selected rows", {
  df  <- make_cut_data(nsim = 10)
  cut <- cutoff_fast(df, event.looks = 50, event.subset = df$group == 1)
  for (s in 1:10) {
    d_s <- df[df$sim == s, ]
    expect_equal(unname(cut[s, 1]),
                 ref_event_time(d_s, 50, d_s$group == 1))
  }
})

test_that("cutoff_fast: a matrix of event targets is applied per simulation", {
  df  <- make_cut_data(nsim = 10)
  tg  <- cbind(40 + 0:9, 80 + 2 * (0:9))
  cut <- cutoff_fast(df, event.looks = tg)
  for (s in 1:10) {
    d_s <- df[df$sim == s, ]
    expect_equal(unname(cut[s, ]),
                 c(ref_event_time(d_s, tg[s, 1]),
                   ref_event_time(d_s, tg[s, 2])))
  }
  expect_true(all(is.na(attr(cut, "look.value"))))
})

test_that("cutoff_fast: events can be counted on another endpoint", {
  df <- simdata_fast(nsim = 10, n = c(150, 150), a.time = c(0, 12),
                     a.rate = 300 / 12,
                     h01.hazard = list(0.08, 0.05),
                     h02.hazard = list(0.04, 0.03), seed = 3)
  cut <- cutoff_fast(df, event.looks = 100, tte.col = "e1_tte",
                     event.col = "e1_event")
  for (s in 1:10) {
    d_s <- df[df$sim == s & df$e1_event == 1, ]
    cal <- sort(d_s$accrual_time + d_s$e1_tte)
    expect_equal(unname(cut[s, 1]), cal[100])
  }
})

test_that("cutoff_fast: agrees with simtrial::get_analysis_date", {
  skip_if_not_installed("simtrial")
  df  <- make_cut_data(nsim = 20)
  cut <- cutoff_fast(df, event.looks = c(80, 140), time.looks = c(14, NA),
                     max.time = c(NA, 40), min.gap = c(NA, 8),
                     min.enrolled = c(180, NA), min.followup = c(3, NA))
  for (s in 1:20) {
    d_s <- df[df$sim == s, ]
    x <- data.frame(enroll_time = d_s$accrual_time,
                    cte = d_s$accrual_time + d_s$tte,
                    fail = d_s$event, stratum = "All")
    r1 <- tryCatch(
      simtrial::get_analysis_date(x, planned_calendar_time = 14,
                                  target_event_overall = 80,
                                  min_n_overall = 180, min_followup = 3),
      error = function(e) NA_real_)
    r2 <- tryCatch(
      simtrial::get_analysis_date(x, target_event_overall = 140,
                                  max_extension_for_target_event = 40,
                                  previous_analysis_date = r1,
                                  min_time_after_previous_analysis = 8),
      error = function(e) NA_real_)
    expect_equal(unname(cut[s, ]), c(r1, r2))
  }
})

test_that("cutoff_fast: input validation", {
  df <- make_cut_data(nsim = 3)
  expect_error(cutoff_fast(df), "at least one")
  expect_error(cutoff_fast(df, event.looks = 1.5), "whole numbers")
  expect_error(cutoff_fast(df, event.looks = c(10, 20), time.looks = c(1, 2, 3)),
               "length 1 or")
  expect_error(cutoff_fast(df, event.looks = matrix(10, 2, 1)), "one row per")
  expect_error(cutoff_fast(df, min.followup = 3), "requires 'min.enrolled'")
  expect_error(cutoff_fast(df, time.looks = -1), "positive")
  expect_error(cutoff_fast(df, event.looks = 10, event.subset = TRUE),
               "event.subset")
  expect_error(cutoff_fast(df, event.looks = 10, tte.col = "nope"),
               "must be a data frame")
  expect_error(cutoff_fast(df[0, ], event.looks = 10), "no rows")
})

test_that("cutoff_fast: a later look before the previous one gives a warning", {
  df <- simdata_fast(nsim = 10, n = c(150, 150), a.time = c(0, 12),
                     a.rate = 25, e.median = list(12, 18), d.hazard = 0.01,
                     seed = 1)
  # About 100 events are expected by month 15, so the capped second look
  # precedes the 150-event first look.
  expect_warning(cut <- cutoff_fast(df, event.looks = c(150, 220),
                                    max.time = c(NA, 15)),
                 "earlier than that of the previous look")
  expect_true(any(cut[, 2] < cut[, 1], na.rm = TRUE))
  expect_warning(cutoff_fast(df, event.looks = c(150, 220)), NA)
})
