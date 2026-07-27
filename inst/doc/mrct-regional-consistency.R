## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>"
)

## ----load, message = FALSE----------------------------------------------------
library(FastSurvival)

## ----scenario-----------------------------------------------------------------
NSIM    <- 10000       # number of simulated trials
SEED    <- 1

mst.T   <- 5.811       # treatment median survival
mst.C   <- 4.3         # control median survival
N.T     <- c(25, 112, 113)   # treatment sample size per region (region 1 = Japan)
N.C     <- c(25, 112, 113)   # control sample size per region
a       <- 12.5        # accrual duration
a_delay <- 3           # accrual delay for region 1 (Japan starts later than others)
PI      <- 0.5         # effect-retention fraction for Method 1

E     <- c(142, 248, 354)  # target total events at the three looks

## ----boundaries---------------------------------------------------------------
eff.bound <- c(NA, -2.437, -2)   # lower efficacy boundary on logrank.z
fut.bound <- c(0.381, NA, -2)    # upper (non-binding) futility boundary

## ----simulate-----------------------------------------------------------------
K      <- length(N.T)
blocks <- lapply(seq_len(K), function(k) {
  a.time_k <- if (k == 1L) c(a_delay, a) else c(0, a)
  dat <- simdata_fast(
    nsim     = NSIM,
    n        = c(N.C[k], N.T[k]),
    a.time   = a.time_k,
    a.prop   = 1,
    e.median = list(mst.C, mst.T),
    seed     = SEED + k
  )
  dat$subgroup <- k
  dat
})
data_all <- do.call(rbind, blocks)
head(data_all)

## ----analyze------------------------------------------------------------------
res <- analysis_fast(
  data        = data_all,
  control     = 1,
  event.looks = E,
  stat        = c("logrank", "coxph"),
  side        = 1,
  by.subgroup = TRUE
)
head(res[, c("sim", "look", "population", "n.event", "logrank.z", "cox.hr")])

## ----indicators---------------------------------------------------------------
ov  <- res[res$population == "overall", ]
key <- function(d) paste(d$sim, d$look)

# Region-1 hazard ratio aligned to the overall rows
r1  <- res[res$population == "subgroup_1", ]
ov$HR.1 <- r1$cox.hr[match(key(ov), key(r1))]

# Method 2: all regional hazard ratios below 1 at this (sim, look)
reg <- res[res$population != "overall", ]
all_below_1 <- tapply(reg$cox.hr, key(reg), function(h) all(h < 1))
ov$M2 <- as.integer(all_below_1[key(ov)])

# Method 1: region 1 retains at least a fraction PI of the overall effect
ov$M1 <- as.integer((1 - ov$HR.1) > PI * (1 - ov$cox.hr))

# Undemonstrated criteria (NA hazard ratio) count as not met
ov$M1[is.na(ov$M1)] <- 0L
ov$M2[is.na(ov$M2)] <- 0L

## ----gsd----------------------------------------------------------------------
sims  <- sort(unique(ov$sim));  nsim <- length(sims)
looks <- sort(unique(ov$look)); nlk  <- length(looks)
cell  <- cbind(match(ov$sim, sims), match(ov$look, looks))
mk    <- function(v) { m <- matrix(NA_real_, nsim, nlk); m[cell] <- v; m }

Zm   <- mk(ov$logrank.z)
M1m  <- mk(ov$M1)
M2m  <- mk(ov$M2)
Cutm <- mk(ov$cutoff)

effB <- matrix(eff.bound, nsim, nlk, byrow = TRUE)
futB <- matrix(fut.bound, nsim, nlk, byrow = TRUE)
eff_cross <- !is.na(effB) & Zm <= effB
fut_cross <- !is.na(futB) & Zm >= futB
stop_any  <- eff_cross | fut_cross

first_stop <- apply(stop_any, 1L, function(r) {
  w <- which(r); if (length(w)) w[1L] else nlk
})
reject <- eff_cross[cbind(seq_len(nsim), first_stop)]

## ----consistency--------------------------------------------------------------
con1 <- con2 <- joi1 <- joi2 <- eff <- atime <- rep(NA_real_, nlk)
for (k in seq_len(nlk)) {
  atime[k] <- mean(Cutm[, k], na.rm = TRUE)   # calendar time of look k (all trials)
  sel      <- reject & first_stop == k        # trials stopping for efficacy at look k
  eff[k]   <- mean(sel)
  joi1[k]  <- sum(M1m[sel, k]) / nsim          # joint prob (0 when no efficacy stop)
  joi2[k]  <- sum(M2m[sel, k]) / nsim
  if (any(sel)) {                              # conditional prob (undefined otherwise)
    con1[k] <- mean(M1m[sel, k])
    con2[k] <- mean(M2m[sel, k])
  }
}

consistency <- data.frame(
  look          = looks,
  Info.Fraction = E / max(E),
  Events        = E,
  Analysis_time = atime,
  Efficacy      = eff,
  cum.power     = cumsum(eff),
  CON_M1        = con1,
  JOI_M1        = joi1,
  CON_M2        = con2,
  JOI_M2        = joi2
)
consistency

## ----summary------------------------------------------------------------------
gs <- simsummary_fast(
  res[res$population == "overall", ],
  eff.col   = "logrank.z", efficacy = eff.bound,
  fut.col   = "logrank.z", futility = fut.bound,
  direction = "lower"
)
gs

