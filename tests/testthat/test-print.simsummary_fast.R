# ------------------------------------------------------------------ #
#  print.simsummary_fast
# ------------------------------------------------------------------ #

# Two looks, four simulated trials, Z-mode efficacy boundaries.
make_summary <- function() {
  df <- data.frame(sim = rep(1:4, each = 2), look = rep(1:2, 4),
                   logrank.z = c(-2.5, -3.1, -1.0, -1.5, -3.0, -2.2, 0.5, -0.4),
                   n.event = rep(c(50, 100), 4), cutoff = rep(c(12, 24), 4))
  simsummary_fast(df, eff.col = "logrank.z", efficacy = c(-2.8, -1.96))
}

test_that("print.simsummary_fast returns its input invisibly", {
  df <- data.frame(sim = 1:4, look = 1L, logrank.z = c(-2.5, -1.0, -3.0, 0.5),
                   n.event = c(50, 52, 48, 55), cutoff = rep(24, 4))
  res <- simsummary_fast(df, eff.col = "logrank.z", efficacy = -1.96)
  expect_output(print(res), "Group-Sequential Operating Characteristics")
  # capture.output keeps the printed report out of the test log
  utils::capture.output(vis <- withVisible(print(res)))
  expect_false(vis$visible)
  expect_identical(vis$value, res)
})

test_that("print.simsummary_fast prints a selection of rows as a report", {
  res <- make_summary()
  sub <- res[res$look != "overall", ]
  expect_s3_class(sub, "simsummary_fast")
  expect_false(is.null(attr(sub, "boundary")))
  expect_output(print(sub), "Stopping Boundaries: Look by Look")
})

test_that("print.simsummary_fast prints a selection of columns as a data frame", {
  res <- make_summary()
  # Columns of the report are missing.
  sub1 <- res[, c("look", "cum.reject")]
  expect_s3_class(sub1, "simsummary_fast")
  out1 <- utils::capture.output(v <- withVisible(print(sub1)))
  expect_false(any(grepl("Group-Sequential", out1, fixed = TRUE)))
  expect_true(any(grepl("cum.reject", out1, fixed = TRUE)))
  expect_false(v$visible)
  expect_identical(v$value, sub1)
  # All the columns of the report are kept, but not the boundary settings.
  sub2 <- res[, c("population", "look", "cum.reject", "prob.stop.efficacy")]
  expect_null(attr(sub2, "boundary"))
  out2 <- utils::capture.output(print(sub2))
  expect_false(any(grepl("Group-Sequential", out2, fixed = TRUE)))
  expect_true(any(grepl("prob.stop.efficacy", out2, fixed = TRUE)))
})

test_that("print.simsummary_fast prints a selection without look rows as a data frame", {
  res <- make_summary()
  sub <- res[res$look == "overall", ]
  out <- utils::capture.output(print(sub))
  expect_false(any(grepl("Group-Sequential", out, fixed = TRUE)))
  expect_true(any(grepl("overall", out, fixed = TRUE)))
})

test_that("print.simsummary_fast prints a selection of some looks as a data frame", {
  res <- make_summary()
  # Look 2 only: the boundaries refer to both looks of the design.
  out1 <- utils::capture.output(print(res[res$look %in% c("2", "overall"), ]))
  expect_false(any(grepl("Stopping Boundaries", out1, fixed = TRUE)))
  expect_true(any(grepl("cum.reject", out1, fixed = TRUE)))
  # Reordered rows.
  out2 <- utils::capture.output(print(res[c(2, 1, 3), ]))
  expect_false(any(grepl("Stopping Boundaries", out2, fixed = TRUE)))
  # Two of three looks: printed without an error.
  df3 <- data.frame(sim = rep(1:4, each = 3), look = rep(1:3, 4),
                    logrank.z = rep(c(-1, -2, -3), 4),
                    n.event = rep(c(30, 60, 90), 4),
                    cutoff = rep(c(6, 12, 18), 4))
  res3 <- simsummary_fast(df3, eff.col = "logrank.z",
                          efficacy = c(-3.5, -2.8, -1.96))
  out3 <- utils::capture.output(
    v <- withVisible(print(res3[res3$look %in% c("1", "3", "overall"), ])))
  expect_false(v$visible)
  expect_false(any(grepl("Stopping Boundaries", out3, fixed = TRUE)))
})

test_that("print.simsummary_fast pairs each look with its own boundary", {
  res <- make_summary()
  out <- utils::capture.output(print(res))
  # Header of the boundary columns (the statistic need not be a Z-score).
  expect_true(any(grepl("Efficacy Bound", out, fixed = TRUE)))
  tab_line <- function(lk) {
    out[grepl(paste0("^ +", lk, " +[0-9.]+ +[0-9.]+ +-"), out)][1L]
  }
  expect_match(tab_line(1), "-2.8000", fixed = TRUE)
  expect_match(tab_line(2), "-1.9600", fixed = TRUE)
  # A selection of whole blocks (here the single population) keeps the report.
  sub <- res[res$population == res$population[1L], ]
  expect_output(print(sub), "Stopping Boundaries: Look by Look")
})

test_that("print.simsummary_fast prints results combined by rbind() as a data frame", {
  # Two results with different boundaries: rbind() keeps the attributes of the
  # first, so the second block would be printed with the wrong boundaries.
  df <- data.frame(sim = rep(1:4, each = 2), look = rep(1:2, 4),
                   logrank.z = c(-2.5, -3.1, -1.0, -1.5, -3.0, -2.2, 0.5, -0.4),
                   n.event = rep(c(50, 100), 4), cutoff = rep(c(12, 24), 4))
  dA <- df
  dA$population <- "overall"
  dB <- df
  dB$population <- "subgroup"
  resA <- simsummary_fast(dA, eff.col = "logrank.z", efficacy = c(-2.8, -1.96))
  resB <- simsummary_fast(dB, eff.col = "logrank.z", efficacy = c(-3.5, -2.5))
  out1 <- utils::capture.output(print(rbind(resA, resB)))
  expect_false(any(grepl("Group-Sequential", out1, fixed = TRUE)))
  # The same result twice duplicates the looks of its block.
  out2 <- utils::capture.output(print(rbind(resA, resA)))
  expect_false(any(grepl("Group-Sequential", out2, fixed = TRUE)))
  # One result with both populations is still printed as a report.
  res2 <- simsummary_fast(rbind(dA, dB), eff.col = "logrank.z",
                          efficacy = c(-2.8, -1.96))
  out3 <- utils::capture.output(print(res2))
  expect_true(any(grepl("Population: subgroup", out3, fixed = TRUE)))
})

test_that("print.simsummary_fast checks the rows of each block against the result", {
  # One data set with two populations summarized with two sets of boundaries:
  # both results have blocks of the same names.
  df <- data.frame(sim = rep(1:4, each = 2), look = rep(1:2, 4),
                   logrank.z = c(-2.5, -3.1, -1.0, -1.5, -3.0, -2.2, 0.5, -0.4),
                   n.event = rep(c(50, 100), 4), cutoff = rep(c(12, 24), 4))
  d2 <- rbind(transform(df, population = "overall"),
              transform(df, population = "subgroup"))
  resA <- simsummary_fast(d2, eff.col = "logrank.z", efficacy = c(-2.8, -1.96))
  resB <- simsummary_fast(d2, eff.col = "logrank.z", efficacy = c(-3.5, -2.5))
  is_report <- function(x) {
    any(grepl("Group-Sequential", utils::capture.output(print(x)),
              fixed = TRUE))
  }
  # Rows of the subgroup from the other boundaries.
  combo <- rbind(resA[resA$population == "overall", ],
                 resB[resB$population == "subgroup", ])
  expect_false(is_report(combo))
  expect_true(any(grepl("cum.reject", utils::capture.output(print(combo)),
                        fixed = TRUE)))
  # One population of either result stays a report with its own boundaries.
  expect_true(is_report(resB[resB$population == "subgroup", ]))
  out3 <- utils::capture.output(print(resA[resA$population == "subgroup", ]))
  expect_true(any(grepl("-2.8000", out3, fixed = TRUE)))
  # Without the overall rows the looks are still checked.
  expect_true(is_report(resA[resA$look != "overall", ]))
  # The rows of the first result, combined again in the same order.
  expect_true(is_report(rbind(resA[resA$population == "overall", ],
                              resA[resA$population == "subgroup", ])))
})

test_that("print.simsummary_fast prints arm results and old attributes as reports", {
  df <- data.frame(sim = rep(1:4, each = 2), look = rep(1:2, 4),
                   logrank.z = c(-2.5, -3.1, -1.0, -1.5, -3.0, -2.2, 0.5, -0.4),
                   n.event = rep(c(50, 100), 4), cutoff = rep(c(12, 24), 4))
  da <- rbind(transform(df, arm = 2),
              transform(df, arm = 3, logrank.z = -logrank.z))
  resC <- simsummary_fast(da, eff.col = "logrank.z", efficacy = c(-2.8, -1.96))
  out4 <- utils::capture.output(print(resC))
  expect_true(any(grepl("Arm: 3", out4, fixed = TRUE)))
  out4b <- utils::capture.output(print(resC[resC$arm == 3, ]))
  expect_true(any(grepl("Group-Sequential", out4b, fixed = TRUE)))
  # A boundary attribute without 'looks' and 'rows', as before 1.2.0.
  old <- make_summary()
  bd <- attr(old, "boundary")
  bd$looks <- NULL
  bd$rows <- NULL
  attr(old, "boundary") <- bd
  expect_output(print(old), "Group-Sequential Operating Characteristics")
})

test_that("print.simsummary_fast checks the overall rows against the result", {
  # z1 and z2 stop different trials at look 1 with the same probabilities, so
  # the look rows agree and only the overall rows differ.
  df <- data.frame(sim = rep(1:4, each = 2), look = rep(1:2, 4),
                   z1 = c(-3.0, -1.0, -1.0, -1.0, -1.0, -2.5, 0.5, 0.5),
                   z2 = c(-1.0, -1.0, -3.0, -1.0, -1.0, -2.5, 0.5, 0.5),
                   n.event = rep(c(50, 100), 4),
                   cutoff = c(10, 24, 14, 24, 12, 24, 12, 24))
  r1 <- simsummary_fast(df, eff.col = "z1", efficacy = c(-2.8, -1.96))
  r2 <- simsummary_fast(df, eff.col = "z2", efficacy = c(-2.8, -1.96))
  expect_equal(r1$cum.reject, r2$cum.reject)
  expect_equal(r1$cutoff.mean, c(12, 24, 20.5))
  expect_equal(r2$cutoff.mean, c(12, 24, 21.5))
  out <- utils::capture.output(
    print(rbind(r1[r1$look != "overall", ], r2[r2$look == "overall", ])))
  expect_false(any(grepl("Group-Sequential", out, fixed = TRUE)))
  expect_true(any(grepl("cum.reject", out, fixed = TRUE)))
})
