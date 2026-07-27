#' Pairwise comparisons of each experimental arm against a shared control
#'
#' Runs the two-group analysis of \code{\link{analysis_fast}} once for each
#' experimental arm against a common control arm, on multi-arm simulated data
#' such as the output of \code{\link{simdata_fast}} with a length-greater-than-two
#' \code{n}. The contrasts are analyzed at a shared set of looks, and the results
#' are stacked into one long data frame with an \code{arm} column, optionally with
#' a Bonferroni-adjusted p-value across the contrasts.
#'
#' @details
#' Two timing regimes are supported, exactly one of which must be requested.
#'
#' With \code{time.looks}, every contrast is analyzed at the same fixed calendar
#' time or times, so all contrasts share the same data cutoff by construction.
#'
#' With \code{event.looks}, the analysis is event-driven and \code{primary} names
#' the experimental arm whose control-versus-arm comparison defines the data
#' cutoff. The primary contrast is analyzed at the requested cumulative event
#' counts, its per-simulation calendar cutoffs are recorded, and then every
#' contrast, including the primary one, is analyzed at those same cutoffs. This
#' reproduces the standard design in which the primary event-driven analysis
#' fixes a single data cutoff at which all comparisons are performed. For a
#' simulation in which the primary event target is not reached at a look, that
#' look is marked with \code{reached = FALSE} and \code{NA} statistics for every
#' contrast.
#'
#' The Bonferroni option multiplies each p-value by the number of contrasts and
#' caps it at one, controlling the family-wise error rate across the
#' control-versus-arm comparisons at each look. Multiplicity across looks is a
#' separate matter handled by group-sequential boundaries in
#' \code{\link{simsummary_fast}}, not by this adjustment.
#'
#' This is a single-endpoint helper: it reads the \code{tte} and \code{event}
#' columns and does not support subgroups. Comparisons for a second endpoint are
#' obtained by calling \code{pairwise_fast} again on that endpoint's columns.
#'
#' @param data A data frame of simulated trial data with columns \code{sim},
#'   \code{group}, \code{accrual_time}, \code{tte}, and \code{event}, as produced
#'   by \code{\link{simdata_fast}}.
#' @param control The group label of the control arm.
#' @param event.looks A numeric vector of target cumulative event counts for the
#'   primary contrast, one per look. Mutually exclusive with \code{time.looks};
#'   requires \code{primary}.
#' @param time.looks A numeric vector of calendar times, one per look. Mutually
#'   exclusive with \code{event.looks}.
#' @param primary The group label of the experimental arm whose control-versus-arm
#'   comparison defines the shared calendar cutoff. Required with
#'   \code{event.looks} and ignored with \code{time.looks}.
#' @param arms A vector of experimental arm labels to compare against the control.
#'   Defaults to every group other than \code{control}.
#' @param stat The test statistic passed to \code{\link{analysis_fast}} (for
#'   example \code{"logrank"}, \code{"coxph"}, or \code{"rmst"}).
#' @param adjust Multiplicity adjustment across contrasts. Either \code{"none"}
#'   (default) or \code{"bonferroni"}, which adds a \code{p.adj} column.
#' @param p.col The name of the p-value column to adjust. By default it is the
#'   column matching the chosen \code{stat} (for example \code{"logrank.p"}).
#' @param ... Further arguments passed to \code{\link{analysis_fast}}, such as
#'   \code{side}, \code{conf.level}, \code{tau}, \code{rho}, and \code{gamma}.
#'   The argument \code{by.subgroup} is not supported.
#'
#' @return A data frame with one row per contrast, look, and simulation. The
#'   first column is \code{arm} (the experimental arm compared against
#'   \code{control}), followed by the columns returned by
#'   \code{\link{analysis_fast}}: \code{sim}, \code{look}, \code{look.value},
#'   \code{cutoff}, \code{reached}, the enrollment and event counts, and the
#'   statistic columns for the chosen \code{stat}. When \code{adjust =
#'   "bonferroni"} a \code{p.adj} column is appended.
#'
#' @examples
#' # Three-arm trial: control (group 1) and two experimental arms.
#' dfk <- simdata_fast(
#'   nsim     = 100,
#'   n        = c(120, 120, 120),
#'   a.time   = c(0, 12),
#'   a.rate   = 360 / 12,
#'   e.median = list(12, 16, 20),
#'   seed     = 8
#' )
#'
#' # Fixed calendar look at month 30, Bonferroni across the two contrasts.
#' pw <- pairwise_fast(dfk, control = 1, time.looks = 30,
#'                     stat = "logrank", side = 1, adjust = "bonferroni")
#' head(pw)
#'
#' # Event-driven: arm 3 is the primary contrast and its 200th control-plus-arm-3
#' # event fixes the cutoff at which both contrasts are analyzed.
#' pw2 <- pairwise_fast(dfk, control = 1, event.looks = 200, primary = 3,
#'                      stat = "logrank", side = 1, adjust = "bonferroni")
#' head(pw2)
#'
#' @seealso \code{\link{simdata_fast}}, \code{\link{analysis_fast}}
#'
#' @export
pairwise_fast <- function(data, control,
                          event.looks = NULL, time.looks = NULL,
                          primary = NULL, arms = NULL,
                          stat = "logrank",
                          adjust = c("none", "bonferroni"),
                          p.col = NULL, ...) {
  adjust <- match.arg(adjust)
  dots <- list(...)
  if (isTRUE(dots$by.subgroup)) {
    stop("pairwise_fast() does not support 'by.subgroup'; it compares each arm ",
         "against the control on the pooled population.")
  }
  if (is.null(data[["sim"]]) || is.null(data[["group"]]) ||
      is.null(data[["tte"]]) || is.null(data[["event"]]) ||
      is.null(data[["accrual_time"]])) {
    stop("'data' must have 'sim', 'group', 'accrual_time', 'tte', and 'event' ",
         "columns, as produced by simdata_fast().")
  }

  has_event <- !is.null(event.looks)
  has_time  <- !is.null(time.looks)
  if (has_event == has_time) {
    stop("Supply exactly one of 'event.looks' or 'time.looks'.")
  }

  groups <- sort(unique(data$group))
  if (!control %in% groups) {
    stop("'control' (", control, ") is not among the group labels in 'data'.")
  }
  if (is.null(arms)) arms <- setdiff(groups, control)
  arms <- sort(unique(arms))
  if (length(arms) < 1L) {
    stop("No experimental arms to compare against the control.")
  }
  if (any(arms == control)) {
    stop("'arms' must not include the control group.")
  }
  if (!all(arms %in% groups)) {
    stop("'arms' contains labels not present in 'data': ",
         paste(setdiff(arms, groups), collapse = ", "))
  }

  # Resolve the p-value column to adjust from the analysis output names.
  p_col_of <- function(nms) {
    if (!is.null(p.col)) {
      if (!p.col %in% nms) {
        stop("'p.col' (", p.col, ") is not a column of the analysis output.")
      }
      return(p.col)
    }
    cand <- grep("\\.p$", nms, value = TRUE)
    pref <- paste0(stat, ".p")
    if (pref %in% cand) return(pref)
    if (length(cand) == 1L) return(cand)
    stop("Could not determine the p-value column automatically; supply 'p.col'. ",
         "Candidates: ", paste(cand, collapse = ", "))
  }

  add_bonferroni <- function(out) {
    if (adjust == "none") return(out)
    pc <- p_col_of(names(out))
    out[["p.adj"]] <- pmin(1, out[[pc]] * length(arms))
    out
  }

  # ---- Fixed calendar looks: all contrasts share the same time.looks --------
  if (has_time) {
    blocks <- lapply(arms, function(j) {
      sub <- data[data$group %in% c(control, j), , drop = FALSE]
      res <- analysis_fast(sub, control = control, time.looks = time.looks,
                           stat = stat, ...)
      cbind(arm = j, res)
    })
    out <- do.call(rbind, blocks)
    out <- out[order(out$look, out$arm, out$sim), , drop = FALSE]
    rownames(out) <- NULL
    return(add_bonferroni(out))
  }

  # ---- Event-driven: the primary contrast defines the shared cutoff ---------
  if (is.null(primary)) {
    stop("With 'event.looks', 'primary' must name the experimental arm whose ",
         "event-driven analysis defines the shared calendar cutoff.")
  }
  if (length(primary) != 1L || !primary %in% arms) {
    stop("'primary' must be a single experimental arm listed in 'arms'.")
  }

  sims_all <- sort(unique(data$sim))
  nsim     <- length(sims_all)
  L        <- length(event.looks)

  prim_sub <- data[data$group %in% c(control, primary), , drop = FALSE]
  prim_res <- analysis_fast(prim_sub, control = control,
                            event.looks = event.looks, stat = stat, ...)
  # analysis_fast orders rows (sim outer, look inner), so a byrow fill recovers
  # the per-simulation-by-look cutoff and reached matrices.
  A_mat <- matrix(prim_res$cutoff,  nrow = nsim, ncol = L, byrow = TRUE)
  R_mat <- matrix(prim_res$reached, nrow = nsim, ncol = L, byrow = TRUE)
  finiteA <- A_mat[is.finite(A_mat)]
  big     <- if (length(finiteA)) max(finiteA) + 1 else 1

  id_cols <- c("arm", "sim", "look", "look.value", "cutoff", "reached")

  blocks <- vector("list", length(arms) * L)
  b <- 0L
  for (j in arms) {
    sub     <- data[data$group %in% c(control, j), , drop = FALSE]
    row_sim <- match(sub$sim, sims_all)
    for (l in seq_len(L)) {
      A_sim   <- A_mat[, l]
      reached <- R_mat[, l]
      # Simulations that did not reach the target keep all subjects (cutoff set
      # to 'big') so every simulation stays present for a clean row alignment;
      # their statistics are blanked afterward.
      A_use <- ifelse(is.finite(A_sim), A_sim, big)
      A_row <- A_use[row_sim]

      enrolled <- sub$accrual_time <= A_row
      obs_t    <- pmin(sub$tte, A_row - sub$accrual_time)
      obs_e    <- sub$event * as.integer(sub$accrual_time + sub$tte <= A_row)
      cut_dat  <- data.frame(
        sim          = sub$sim,
        group        = sub$group,
        accrual_time = sub$accrual_time,
        tte          = obs_t,
        event        = obs_e
      )[enrolled, , drop = FALSE]

      res <- analysis_fast(cut_dat, control = control, time.looks = big,
                           stat = stat, ...)
      # Align the one-row-per-simulation result to the full simulation set.
      res <- res[match(sims_all, res$sim), , drop = FALSE]

      res$look       <- l
      res$look.value <- event.looks[l]
      res$cutoff     <- A_sim
      res$reached    <- reached

      out_j <- cbind(arm = j, res)
      if (any(!reached)) {
        blank_cols <- setdiff(names(out_j), id_cols)
        out_j[!reached, blank_cols] <- NA
      }

      b <- b + 1L
      blocks[[b]] <- out_j
    }
  }
  out <- do.call(rbind, blocks)
  out <- out[order(out$look, out$arm, out$sim), , drop = FALSE]
  rownames(out) <- NULL
  add_bonferroni(out)
}
