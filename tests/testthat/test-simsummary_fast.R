# ------------------------------------------------------------------ #
#  Helper: pull a single cell from a simsummary_fast block
# ------------------------------------------------------------------ #
cell <- function(res, look, col, pop = "overall") {
  res[res$population == pop & res$look == as.character(look), col]
}

# ------------------------------------------------------------------ #
#  Fixed design (single look): Z mode and p mode
# ------------------------------------------------------------------ #

test_that("simsummary_fast fixed design Z mode counts lower-tail crossings", {
  df <- data.frame(
    sim        = 1:6,
    look       = 1L,
    logrank.z  = c(-2.5, -1.0, -3.0, 0.5, -2.1, -1.9),
    n.enrolled = rep(100, 6),
    n.event    = c(50, 52, 48, 55, 51, 49),
    cutoff     = rep(24, 6)
  )
  res <- simsummary_fast(df, eff.col = "logrank.z", efficacy = -1.96,
                         direction = "lower")
  expect_s3_class(res, "simsummary_fast")
  # z <= -1.96 for -2.5, -3.0, -2.1 -> 3 of 6
  expect_equal(cell(res, 1, "prob.stop.efficacy"), 3 / 6)
  expect_equal(cell(res, "overall", "prob.stop.efficacy"), 3 / 6)
  expect_equal(cell(res, "overall", "prob.stop.futility"), 0)
  expect_equal(cell(res, "overall", "n.event.mean"), mean(c(50, 52, 48, 55, 51, 49)))
})

test_that("simsummary_fast fixed design p mode counts p <= alpha", {
  df <- data.frame(
    sim       = 1:5,
    look      = 1L,
    logrank.p = c(0.01, 0.20, 0.04, 0.06, 0.001),
    n.event   = c(40, 41, 42, 43, 44),
    cutoff    = rep(18, 5)
  )
  res <- simsummary_fast(df, p.col = "logrank.p", alpha = 0.025)
  expect_equal(cell(res, "overall", "prob.stop.efficacy"), 2 / 5)
})

test_that("simsummary_fast upper direction matches positive-Z benefit", {
  df <- data.frame(
    sim     = 1:4,
    look    = 1L,
    rmst.z  = c(2.5, 1.0, -0.5, 3.0),
    n.event = c(30, 31, 32, 33),
    cutoff  = rep(12, 4)
  )
  res <- simsummary_fast(df, eff.col = "rmst.z", efficacy = 1.96,
                         direction = "upper")
  expect_equal(cell(res, "overall", "prob.stop.efficacy"), 2 / 4)
})

# ------------------------------------------------------------------ #
#  Group-sequential design: sequential crossing, futility, timing
# ------------------------------------------------------------------ #

test_that("simsummary_fast applies boundaries in look order with futility", {
  df <- data.frame(
    sim        = rep(1:4, each = 2),
    look       = rep(1:2, 4),
    logrank.z  = c(-3.0, -2.0,   # sim1: efficacy at look 1
                   -1.0, -2.5,   # sim2: efficacy at look 2
                   -0.1, -0.2,   # sim3: never crosses -> continue
                   1.0, -2.5),  # sim4: futility at look 1
    n.enrolled = rep(c(100, 200), 4),
    n.event    = c(40, 80, 41, 82, 39, 79, 42, 83),
    cutoff     = rep(c(12, 24), 4)
  )
  res <- simsummary_fast(df, eff.col = "logrank.z",
                         efficacy = c(-2.8, -1.96),
                         futility = c(0.5, NA),
                         direction = "lower")

  expect_equal(cell(res, "overall", "prob.stop.efficacy"), 2 / 4)  # sim1, sim2
  expect_equal(cell(res, "overall", "prob.stop.futility"), 1 / 4)  # sim4
  # events at stop: sim1=40, sim2=82, sim3 final=79, sim4=42
  expect_equal(cell(res, "overall", "n.event.mean"), mean(c(40, 82, 79, 42)))

  expect_equal(cell(res, 1, "prob.stop.efficacy"), 1 / 4)  # sim1
  expect_equal(cell(res, 2, "prob.stop.efficacy"), 1 / 4)  # sim2
  expect_equal(cell(res, 1, "prob.stop.futility"), 1 / 4)  # sim4
  expect_equal(cell(res, 2, "cum.reject"), 2 / 4)
})

test_that("simsummary_fast NA statistic at a look triggers no boundary", {
  df <- data.frame(
    sim       = rep(1:2, each = 2),
    look      = rep(1:2, 2),
    logrank.z = c(NA, -2.5,    # sim1: look 1 not computable, efficacy at look 2
                  NA, -1.0),   # sim2: never crosses
    n.event   = c(0, 70, 0, 71),
    cutoff    = c(6, 24, 6, 24)
  )
  res <- simsummary_fast(df, eff.col = "logrank.z",
                         efficacy = c(-2.8, -1.96), direction = "lower")
  expect_equal(cell(res, "overall", "prob.stop.efficacy"), 1 / 2)
  expect_equal(cell(res, 1, "prob.stop.efficacy"), 0)
  expect_equal(cell(res, 2, "prob.stop.efficacy"), 1 / 2)
})

test_that("simsummary_fast efficacy precedes futility at the same look", {
  df <- data.frame(
    sim       = c(1, 1),
    look      = c(1, 2),
    logrank.z = c(-3.0, -3.0),
    n.event   = c(40, 80),
    cutoff    = c(12, 24)
  )
  res <- simsummary_fast(df, eff.col = "logrank.z",
                         efficacy = c(-2.8, -1.96),
                         futility = c(-2.9, NA),  # also "crossed" at look 1
                         direction = "lower")
  expect_equal(cell(res, "overall", "prob.stop.efficacy"), 1)
  expect_equal(cell(res, "overall", "prob.stop.futility"), 0)
})

# ------------------------------------------------------------------ #
#  Separate efficacy / futility columns and scales
# ------------------------------------------------------------------ #

test_that("simsummary_fast judges efficacy and futility on different columns", {
  # Efficacy on Z (logrank.z), futility on log HR (cox.coef).
  df <- data.frame(
    sim       = rep(1:3, each = 2),
    look      = rep(1:2, 3),
    logrank.z = c(-0.5, -3.0,   # sim1: efficacy at look 2
                  -0.2, -1.0,   # sim2: no efficacy
                  -0.1, -0.3),  # sim3: no efficacy
    cox.coef  = c(0.25, -0.5,   # sim1: log(1.2)=0.182, look1 cox>=0.182 -> futility?
                  0.30, -0.1,   # sim2: futility at look 1
                  0.05, -0.2),  # sim3: no futility
    n.event   = c(30, 60, 31, 62, 29, 58),
    cutoff    = rep(c(12, 24), 3)
  )
  # Futility only at look 1 (cox.coef >= log(1.2)); efficacy only at look 2.
  res <- simsummary_fast(df, eff.col = "logrank.z",
                         efficacy = c(NA, -1.96),
                         fut.col = "cox.coef",
                         futility = c(log(1.2), NA),
                         direction = "lower")
  # sim1: look1 cox.coef 0.25 >= 0.182 -> futility stop at look 1 (efficacy NA)
  # sim2: look1 cox.coef 0.30 >= 0.182 -> futility stop at look 1
  # sim3: look1 cox.coef 0.05 < 0.182 -> continue; look2 z -0.3 > -1.96 -> no rej
  expect_equal(cell(res, "overall", "prob.stop.futility"), 2 / 3)
  expect_equal(cell(res, "overall", "prob.stop.efficacy"), 0)
  expect_equal(cell(res, 1, "prob.stop.futility"), 2 / 3)
})

test_that("simsummary_fast efficacy-only and futility-only timing via NA", {
  # First look futility-only, later looks efficacy-only.
  df <- data.frame(
    sim       = rep(1:2, each = 3),
    look      = rep(1:3, 2),
    logrank.z = c(-0.1, -3.0, -2.0,   # sim1: efficacy at look 2
                  0.5, -0.2, -0.3),  # sim2: continue past look1 futility? z not used at 1
    cox.coef  = c(-0.1, 0.0, 0.0,     # sim1: no futility at look 1
                  0.3, 0.0, 0.0),    # sim2: futility at look 1 (0.3 >= log(1.2))
    n.event   = c(20, 40, 60, 21, 41, 61),
    cutoff    = rep(c(6, 12, 18), 2)
  )
  res <- simsummary_fast(df, eff.col = "logrank.z",
                         efficacy = c(NA, -2.96, -1.97),
                         fut.col = "cox.coef",
                         futility = c(log(1.2), NA, NA),
                         direction = "lower")
  # sim1: look1 no fut (cox -0.1 < 0.182), look2 z -3.0 <= -2.96 -> efficacy
  # sim2: look1 cox 0.3 >= 0.182 -> futility at look 1
  expect_equal(cell(res, "overall", "prob.stop.efficacy"), 1 / 2)
  expect_equal(cell(res, "overall", "prob.stop.futility"), 1 / 2)
  expect_equal(cell(res, 1, "prob.stop.futility"), 1 / 2)
  expect_equal(cell(res, 2, "prob.stop.efficacy"), 1 / 2)
})

# ------------------------------------------------------------------ #
#  Output structure: overall row appended per population
# ------------------------------------------------------------------ #

test_that("simsummary_fast appends an overall row after each population", {
  df <- data.frame(
    sim       = rep(1:2, each = 2),
    look      = rep(1:2, 2),
    logrank.z = c(-3.0, -2.0, -0.1, -0.2),
    n.event   = c(40, 80, 38, 78),
    cutoff    = rep(c(12, 24), 2)
  )
  res <- simsummary_fast(df, eff.col = "logrank.z",
                         efficacy = c(-2.8, -1.96), direction = "lower")
  expect_equal(res$look, c("1", "2", "overall"))
  # cum.reject at the last look equals the overall rejection rate
  expect_equal(cell(res, 2, "cum.reject"),
               cell(res, "overall", "prob.stop.efficacy"))
})

test_that("simsummary_fast aggregates each population separately", {
  df <- data.frame(
    population = rep(c("overall", "subgroup_1"), each = 4),
    sim        = rep(rep(1:2, each = 2), 2),
    look       = rep(1:2, 4),
    logrank.z  = c(-3.0, -2.0,  -0.1, -0.2,
                   -0.5, -2.5,  -0.4, -0.3),
    n.event    = c(40, 80, 38, 78, 20, 40, 19, 39),
    cutoff     = rep(c(12, 24), 4)
  )
  res <- simsummary_fast(df, eff.col = "logrank.z",
                         efficacy = c(-2.8, -1.96), direction = "lower")
  expect_setequal(unique(res$population), c("overall", "subgroup_1"))
  expect_equal(cell(res, "overall", "prob.stop.efficacy", pop = "overall"), 1 / 2)
  expect_equal(cell(res, "overall", "prob.stop.efficacy", pop = "subgroup_1"), 1 / 2)
  # Each population has 2 look rows + 1 overall row
  expect_equal(nrow(res), (2L + 1L) * 2L)
})

# ------------------------------------------------------------------ #
#  Input validation
# ------------------------------------------------------------------ #

test_that("simsummary_fast errors when both or neither boundary mode supplied", {
  df <- data.frame(sim = 1:3, look = 1L, logrank.z = c(-1, -2, -3),
                   logrank.p = c(0.3, 0.02, 0.001))
  expect_error(
    simsummary_fast(df, eff.col = "logrank.z", efficacy = -1.96,
                    p.col = "logrank.p", alpha = 0.025),
    "exactly one"
  )
  expect_error(simsummary_fast(df), "exactly one")
})

test_that("simsummary_fast errors on boundary length mismatch", {
  df <- data.frame(sim = rep(1:2, each = 2), look = rep(1:2, 2),
                   logrank.z = c(-1, -2, -3, -4))
  expect_error(
    simsummary_fast(df, eff.col = "logrank.z", efficacy = -1.96),
    "one boundary per look"
  )
})

test_that("simsummary_fast errors on a missing statistic column", {
  df <- data.frame(sim = 1:3, look = 1L, logrank.z = c(-1, -2, -3))
  expect_error(
    simsummary_fast(df, eff.col = "not_a_column", efficacy = -1.96),
    "eff.col"
  )
})

test_that("simsummary_fast errors when fut.col is missing", {
  df <- data.frame(sim = rep(1:2, each = 2), look = rep(1:2, 2),
                   logrank.z = c(-1, -2, -3, -4))
  expect_error(
    simsummary_fast(df, eff.col = "logrank.z", efficacy = c(-2.8, -1.96),
                    fut.col = "no_such_col", futility = c(0.5, NA)),
    "fut.col"
  )
})

test_that("simsummary_fast errors when futility supplied in p mode", {
  df <- data.frame(sim = 1:3, look = 1L, logrank.p = c(0.3, 0.02, 0.001))
  expect_error(
    simsummary_fast(df, p.col = "logrank.p", alpha = 0.025, futility = 0.5),
    "only with the Z mode"
  )
})

# ------------------------------------------------------------------ #
#  print method
# ------------------------------------------------------------------ #

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

# ------------------------------------------------------------------ #
#  End-to-end with analysis_fast output
# ------------------------------------------------------------------ #

test_that("simsummary_fast runs on analysis_fast output", {
  dat <- simdata_fast(nsim = 40, n = c(120, 120), a.time = c(0, 12),
                      a.rate = 240 / 12, e.hazard = list(0.05, 0.03),
                      seed = 7)
  res_a <- analysis_fast(dat, control = 1, event.looks = c(60, 120),
                         stat = c("logrank", "coxph"), side = 1)
  ss <- simsummary_fast(res_a, eff.col = "logrank.z",
                        efficacy = c(-2.96, -1.97),
                        fut.col = "cox.coef", futility = c(log(1.3), NA),
                        direction = "lower")
  expect_s3_class(ss, "simsummary_fast")
  rr <- cell(ss, "overall", "prob.stop.efficacy")
  expect_true(rr >= 0 && rr <= 1)
  expect_true(all(ss$cum.reject >= 0 & ss$cum.reject <= 1))
})

test_that("simsummary_fast handles by.subgroup analysis_fast output", {
  dat <- simdata_fast(nsim = 30, n = c(120, 120), a.time = c(0, 12),
                      a.rate = 240 / 12, e.hazard = list(list(0.05, 0.035), 0.03),
                      prevalence = c(0.5, 0.5), seed = 8)
  res_a <- analysis_fast(dat, control = 1, time.looks = c(18, 30),
                         stat = "logrank", side = 1, by.subgroup = TRUE)
  ss <- simsummary_fast(res_a, p.col = "logrank.p",
                        alpha = c(0.0015, 0.0245))
  expect_setequal(unique(ss$population),
                  c("overall", "subgroup_1", "subgroup_2"))
})

# ------------------------------------------------------------------ #
#  Consume gsDesign boundaries and check the stage-wise stopping logic
#  against an independent manual first-crossing count.
# ------------------------------------------------------------------ #
#  This does NOT compare with gsDesign's analytic crossing probabilities:
#  matching those would require the simulated information (events, effect
#  size, timing) to reproduce the design gsDesign was built for, which is a
#  separate task. Instead gsDesign is used only as a boundary supplier, and
#  the test verifies that simsummary_fast applies those boundaries with the
#  same sequential first-crossing rule as a direct R-level computation.

test_that("simsummary_fast consumes gsDesign boundaries with correct logic", {
  skip_if_not_installed("gsDesign")
  skip_on_cran()

  d_looks <- c(120, 200)

  gs <- gsDesign::gsDesign(k = 2, test.type = 1, sfu = gsDesign::sfLDOF,
                           alpha = 0.025, beta = 0.2,
                           timing = d_looks / max(d_looks))
  # gsDesign returns an upper-Z boundary; map to the lower-tail sign of
  # logrank.z (treatment benefit is a negative Z).
  eff_z <- -gs$upper$bound

  set.seed(1)
  dat <- simdata_fast(nsim = 500, n = c(250, 250), a.time = c(0, 18),
                      a.rate = 500 / 18,
                      e.hazard = list(0.03, 0.03 * 0.7),
                      seed = 99)
  res_a <- analysis_fast(dat, control = 1, event.looks = d_looks,
                         stat = "logrank", side = 1)
  ss <- simsummary_fast(res_a, eff.col = "logrank.z", efficacy = eff_z,
                        direction = "lower")

  # Independent manual sequential first-crossing computation on the same
  # logrank.z values, reshaped to a (nsim x 2) matrix.
  sims <- sort(unique(res_a$sim))
  zm   <- matrix(NA_real_, length(sims), 2L)
  zm[cbind(match(res_a$sim, sims), match(res_a$look, c(1, 2)))] <-
    res_a$logrank.z
  cross <- sweep(zm, 2L, eff_z, FUN = "<=")
  cross[is.na(cross)] <- FALSE
  first <- apply(cross, 1L, function(r) {
    w <- which(r); if (length(w)) w[1L] else NA_integer_
  })
  ref_stage1 <- mean(!is.na(first) & first == 1L)
  ref_stage2 <- mean(!is.na(first) & first == 2L)
  ref_cum    <- cumsum(c(ref_stage1, ref_stage2))

  sim_cum <- ss$cum.reject[ss$look %in% c("1", "2")]
  expect_equal(sim_cum, ref_cum, tolerance = 1e-12)
  # The overall rejection rate equals the final cumulative power.
  expect_equal(cell(ss, "overall", "prob.stop.efficacy"), ref_cum[2L],
               tolerance = 1e-12)
})

test_that("simsummary_fast summarizes pairwise_fast output arm by arm", {
  dfk <- simdata_fast(nsim = 20, n = c(100, 100, 100), a.time = c(0, 12),
                      a.rate = 300 / 12, e.median = list(12, 16, 20),
                      seed = 606)
  pw <- pairwise_fast(dfk, control = 1, time.looks = 30, stat = "logrank",
                      side = 1)
  ss <- simsummary_fast(pw, p.col = "logrank.p", alpha = 0.025)
  expect_equal(unique(ss$arm), c(2, 3))
  # Each arm's block equals the summary of that arm's rows alone.
  for (a in c(2, 3)) {
    ref <- simsummary_fast(pw[pw$arm == a, setdiff(names(pw), "arm")],
                           p.col = "logrank.p", alpha = 0.025)
    blk <- ss[ss$arm == a, setdiff(names(ss), "arm")]
    rownames(blk) <- NULL
    expect_equal(as.data.frame(blk), as.data.frame(ref),
                 ignore_attr = TRUE)
  }
  expect_output(print(ss), "Arm: 3")
})

test_that("simsummary_fast rejects duplicate (sim, look) rows", {
  df <- data.frame(sim = c(1, 1, 2, 2), look = 1L,
                   logrank.p = c(0.01, 0.20, 0.03, 0.40))
  expect_error(simsummary_fast(df, p.col = "logrank.p", alpha = 0.025),
               "more than one row")
})
