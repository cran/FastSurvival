#' Look rows of each block of a simsummary_fast result
#'
#' Internal helper of \code{plot.simsummary_fast()}. It splits a result of
#' \code{simsummary_fast()} into its blocks (arms and populations) and returns
#' the look rows of each block, provided that every block keeps all the looks
#' of the design in their original order and that its rows are those of the
#' result stored with the boundary settings. This is the condition under which
#' the print method prints the report, except that the overall rows, which the
#' plot does not use, are not compared.
#'
#' @param x An object of class \code{"simsummary_fast"}.
#'
#' @return A list with one data frame of look rows per block, each with the
#'   attribute \code{label} (the arm and population, or \code{""} for a
#'   single block), or \code{NULL} when the condition does not hold.
#'
#' @keywords internal
#' @noRd
simsummary_blocks <- function(x) {
  bd  <- attr(x, "boundary")
  dat <- as.data.frame(x, stringsAsFactors = FALSE)
  need <- c("look", "population", "cum.reject", "prob.stop.efficacy",
            "prob.stop.futility")
  if (is.null(bd) || is.null(bd$looks) || !all(need %in% names(dat))) {
    return(NULL)
  }

  blk_cols <- c(if ("arm" %in% names(dat)) "arm", "population")
  keys     <- unique(dat[, blk_cols, drop = FALSE])
  multi_pop <- length(unique(dat$population)) > 1L
  out <- vector("list", nrow(keys))
  for (b in seq_len(nrow(keys))) {
    in_blk <- rep(TRUE, nrow(dat))
    for (nm in blk_cols) {
      in_blk <- in_blk & as.character(dat[[nm]]) == as.character(keys[[nm]][b])
    }
    lk <- dat[in_blk & dat$look != "overall", , drop = FALSE]
    if (!identical(as.character(lk$look), bd$looks)) return(NULL)

    # The rows must be those of the stored result (not, for example, rows of
    # another result combined by rbind())
    if (!is.null(bd$rows)) {
      ref <- bd$rows
      for (nm in blk_cols) {
        if (!(nm %in% names(ref))) return(NULL)
        ref <- ref[as.character(ref[[nm]]) == as.character(keys[[nm]][b]), ,
                   drop = FALSE]
      }
      ref  <- ref[ref$look != "overall", , drop = FALSE]
      cols <- intersect(names(lk), names(ref))
      if (nrow(ref) != nrow(lk) ||
          !isTRUE(all.equal(ref[, cols, drop = FALSE], lk[, cols, drop = FALSE],
                            check.attributes = FALSE))) {
        return(NULL)
      }
    }

    lab <- c(if ("arm" %in% blk_cols) paste0("Arm: ", keys$arm[b]),
             if (multi_pop) paste0("Population: ", keys$population[b]))
    rownames(lk) <- NULL
    attr(lk, "label") <- paste(lab, collapse = ", ")
    out[[b]] <- lk
  }
  out
}
