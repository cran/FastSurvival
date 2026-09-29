#' Print Method for Sequential Analysis Summaries
#'
#' @description
#' Formats an object returned by \code{\link{simsummary_fast}} in the style of a
#' group-sequential design report. After a header with the simulation count and
#' the boundary settings, two look-by-look tables are shown: a stopping-boundary
#' table (information fraction, events, sample size, the efficacy and futility
#' boundaries, and the cumulative efficacy crossing probability) and an
#' analysis-timing table (sample size, events, dropouts, pipeline, analysis time,
#' and the per-look efficacy and futility crossing probabilities). An overall
#' block reports the rejection rate and the expected counts and timing at the
#' stopping look.
#'
#' @details
#' Column labels follow the convention of group-sequential design software:
#' \code{Events (s)}, \code{Sample (n)}, \code{Dropouts (d)}, \code{Pipeline}
#' (the enrolled count minus events minus dropouts), \code{Analysis Time} (mean
#' calendar time), and \code{Info. Frac.} (the information fraction, computed as
#' the mean events at a look divided by the mean events at the final look, or
#' from \code{look.value} when no event count is available). Probabilities are
#' printed to \code{digits} decimal places and counts and times to fewer. Because
#' the summary is a Monte Carlo estimate under a single data-generating truth, it
#' does not carry the separate null and alternative columns or the alpha and beta
#' spending of an analytic design report. When the summary has several blocks
#' (arms of \code{\link{pairwise_fast}} output or subgroup populations), each
#' block is printed under its own heading. The underlying object is an ordinary
#' data frame, so the unrounded values remain available by subsetting it directly.
#'
#' @param x An object of class \code{"simsummary_fast"} from
#'   \code{\link{simsummary_fast}}.
#' @param digits A single positive integer, the number of decimal places used for
#'   the printed probabilities. Defaults to 4.
#' @param ... Further arguments, currently ignored.
#'
#' @return The object \code{x}, invisibly.
#'
#' @seealso \code{\link{simsummary_fast}}.
#'
#' @export
print.simsummary_fast <- function(x, digits = 4, ...) {
  bd   <- attr(x, "boundary")
  nsim <- attr(x, "nsim")
  dat  <- as.data.frame(x, stringsAsFactors = FALSE)

  is_z    <- !is.null(bd) && identical(bd$mode, "z")
  has_fut <- is_z && !is.null(bd$futility)
  has_col <- function(nm) nm %in% names(dat)

  fmt <- function(v, d) {
    out <- formatC(as.numeric(v), format = "f", digits = d)
    out[is.na(v)] <- "NA"
    out
  }

  cat("Group-Sequential Operating Characteristics (simsummary_fast)\n")
  if (!is.null(nsim)) {
    lab <- if (length(nsim) == 1L) {
      as.character(nsim)
    } else {
      paste0(paste(nsim, collapse = ", "), " (per block)")
    }
    cat("  Simulations: ", lab, "\n", sep = "")
  }
  if (is_z) {
    cat("  Boundaries: efficacy on '", bd$eff.col, "' (direction = ",
        bd$direction, ")", sep = "")
    if (has_fut) cat(", futility on '", bd$fut.col, "'", sep = "")
    cat("\n")
  } else if (!is.null(bd)) {
    cat("  Boundaries: nominal p-value on '", bd$p.col, "'\n", sep = "")
  }
  cat("\n")

  has_arm    <- has_col("arm")
  pops       <- unique(dat$population)
  multi_pop  <- length(pops) > 1L
  keys       <- unique(dat[, c(if (has_arm) "arm", "population"), drop = FALSE])
  single_blk <- nrow(keys) == 1L

  print_block <- function(block) {
    lk    <- block[block$look != "overall", , drop = FALSE]
    ov    <- block[block$look == "overall", , drop = FALSE]
    nlook <- nrow(lk)

    info_frac <- rep(NA_real_, nlook)
    if (has_col("n.event.mean")) {
      ev  <- lk$n.event.mean
      fin <- ev[nlook]
      if (is.finite(fin) && fin > 0) info_frac <- ev / fin
    } else if (has_col("look.value")) {
      lv <- lk$look.value
      if (any(is.finite(lv)) && max(lv, na.rm = TRUE) > 0) {
        info_frac <- lv / max(lv, na.rm = TRUE)
      }
    }

    eff_b <- if (is_z) bd$efficacy else bd$alpha

    # Table 1: stopping boundaries
    t1 <- list("Look"        = as.character(lk$look),
               "Info. Frac." = fmt(info_frac, 2))
    if (has_col("n.event.mean"))    t1[["Events (s)"]] <- fmt(lk$n.event.mean, 1)
    if (has_col("n.enrolled.mean")) t1[["Sample (n)"]] <- fmt(lk$n.enrolled.mean, 1)
    if (is_z) {
      t1[["Efficacy Z"]] <- fmt(eff_b, digits)
      if (has_fut) t1[["Futility Z"]] <- fmt(bd$futility, digits)
    } else {
      t1[["Nominal p"]] <- fmt(eff_b, digits)
    }
    t1[["Cum. Cross. Eff."]] <- fmt(lk$cum.reject, digits)
    t1 <- data.frame(t1, check.names = FALSE, stringsAsFactors = FALSE)

    cat("Stopping Boundaries: Look by Look\n")
    print(t1, row.names = FALSE, right = TRUE)
    cat("\n")

    # Table 2: timing and crossing probabilities
    t2 <- list("Look"        = as.character(lk$look),
               "Info. Frac." = fmt(info_frac, 2))
    if (has_col("n.enrolled.mean")) t2[["Sample (n)"]]   <- fmt(lk$n.enrolled.mean, 1)
    if (has_col("n.event.mean"))    t2[["Events (s)"]]   <- fmt(lk$n.event.mean, 1)
    if (has_col("n.dropout.mean"))  t2[["Dropouts (d)"]] <- fmt(lk$n.dropout.mean, 1)
    if (has_col("n.pipeline.mean")) t2[["Pipeline"]]     <- fmt(lk$n.pipeline.mean, 1)
    if (has_col("cutoff.mean"))     t2[["Analysis Time"]] <- fmt(lk$cutoff.mean, 2)
    t2[["Cross. Eff."]] <- fmt(lk$prob.stop.efficacy, digits)
    if (has_fut) t2[["Cross. Fut."]] <- fmt(lk$prob.stop.futility, digits)
    t2 <- data.frame(t2, check.names = FALSE, stringsAsFactors = FALSE)

    cat("Events, Sample Size, Dropouts, Pipeline and Analysis Times: Look by Look\n")
    print(t2, row.names = FALSE, right = TRUE)
    cat("\n")

    # Overall summary
    if (nrow(ov) == 1L) {
      kv <- function(label, value) {
        cat("  ", formatC(label, width = -32L), value, "\n", sep = "")
      }
      cat("Overall\n")
      kv("Rejection rate (efficacy):", fmt(ov$cum.reject, digits))
      if (has_fut) {
        kv("Futility-stop rate:", fmt(ov$prob.stop.futility, digits))
      }
      if (has_col("n.event.mean")) {
        kv("Expected events at stop:", fmt(ov$n.event.mean, 1))
      }
      if (has_col("n.enrolled.mean")) {
        kv("Expected sample size at stop:", fmt(ov$n.enrolled.mean, 1))
      }
      if (has_col("cutoff.mean")) {
        kv("Expected analysis time at stop:", fmt(ov$cutoff.mean, 2))
      }
    }
  }

  for (b in seq_len(nrow(keys))) {
    p   <- keys$population[b]
    sel <- dat$population == p
    if (has_arm) sel <- sel & dat$arm == keys$arm[b]
    if (!single_blk) {
      hdr <- c(if (has_arm) paste0("Arm: ", keys$arm[b]),
               if (multi_pop) paste0("Population: ", p))
      cat(paste(hdr, collapse = ", "), "\n", sep = "")
    }
    print_block(dat[sel, , drop = FALSE])
    if (!single_blk) cat("\n")
  }

  invisible(x)
}
