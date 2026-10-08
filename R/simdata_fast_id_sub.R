# ------------------------------------------------------------------ #
#  Internal wrapper: illness-death simulation with subgroups
# ------------------------------------------------------------------ #
# Called by simdata_fast_id() when 'prevalence' is supplied. The hazards have
# already been converted from medians. Every transition, switching, and dropout
# specification follows the rules of 'e.hazard' with subgroups: in a two-group
# simulation it is shared (not a list) or a list with one element per group,
# and each group's element is shared or a list with one element per subgroup
# cell; in a one-group simulation a list holds one element per cell. All
# subjects share one accrual process, and the cells are assigned after accrual
# (see simdata_core_id_sub()). The seed is already set by the caller.
simdata_fast_id_sub <- function(nsim, n, alloc, alloc_given,
                                a.time, a.rate, a.prop,
                                h01.hazard, h01.time,
                                h02.hazard, h02.time,
                                h12.hazard, h12.time,
                                switch.prop,
                                h12.switch.hazard, h12.switch.time,
                                has_dropout, d.hazard, d.time,
                                prevalence, fixed.alloc) {
  group_specific_prev <- is_group_specific_prev(prevalence)
  if (group_specific_prev) {
    spec_ctrl <- build_prevalence_spec(prevalence[["control"]])
    spec_trt  <- build_prevalence_spec(prevalence[["treatment"]])
    lev_c <- apply(spec_ctrl$level_table, 2L, max)
    lev_t <- apply(spec_trt$level_table, 2L, max)
    if (spec_ctrl$n_cell != spec_trt$n_cell || length(lev_c) != length(lev_t) ||
        any(lev_c != lev_t)) {
      stop("'control' and 'treatment' prevalence must use the same ",
           "number of factors and the same number of levels per factor")
    }
  } else {
    spec_ctrl <- build_prevalence_spec(prevalence)
    spec_trt  <- spec_ctrl
  }
  n_cell <- spec_ctrl$n_cell

  # Number of groups, as for 'e.hazard' with subgroups: a length-two 'n', an
  # explicit 'alloc', a group-specific prevalence, or a length-two list with a
  # list element (per-cell values within a group) gives two groups.
  spec_args <- list(h01.hazard = h01.hazard, h02.hazard = h02.hazard,
                    h12.hazard = h12.hazard,
                    h12.switch.hazard = h12.switch.hazard,
                    switch.prop = switch.prop, d.hazard = d.hazard)
  any_nested <- any(vapply(spec_args, function(x) {
    is.list(x) && length(x) == 2L && any(vapply(x, is.list, logical(1L)))
  }, logical(1L)))
  if (length(n) == 2L) {
    two <- TRUE
  } else if (length(n) == 1L) {
    two <- group_specific_prev || alloc_given || any_nested
  } else {
    stop("'n' must be a scalar (total N) or a vector of length 2 (per-group)")
  }
  n_grp <- if (length(n) == 2L) n else if (two) split_total(n, alloc) else n
  n_grp_int <- if (two) as.integer(n_grp[1:2]) else as.integer(n_grp[1L])
  check_output_size(nsim, n_grp_int)
  if (two) {
    for (nm in names(spec_args)) {
      x <- spec_args[[nm]]
      if (is.list(x) && length(x) != 2L) {
        stop("For a two-group simulation, '", nm, "' must be a single ",
             "(shared) specification or a list of length 2 (one element per ",
             "group); it is a list of length ", length(x), ".")
      }
    }
  }

  # Group element (two-group lists are per group), then cell element (lists
  # within a group, or in a one-group simulation, are per cell).
  pick_g <- function(x, g) if (two && is.list(x)) x[[g]] else x
  pick_c <- function(x, cc, nm) {
    if (!is.list(x)) return(x)
    if (length(x) != n_cell) {
      stop("A per-cell list for '", nm, "' must have one element per ",
           "subgroup cell (", n_cell, " here) but has ", length(x), ".")
    }
    x[[cc]]
  }
  cell_spec <- function(x, g, cc, nm) pick_c(pick_g(x, g), cc, nm)

  empty_spec <- list(hazard = 1, fin_time = numeric(0), cum_haz = numeric(0))
  build_group <- function(g) {
    lapply(seq_len(n_cell), function(cc) {
      h01 <- piecewise_precompute(cell_spec(h01.hazard, g, cc, "h01.hazard"),
                                  cell_spec(h01.time, g, cc, "h01.time"))
      h02 <- piecewise_precompute(cell_spec(h02.hazard, g, cc, "h02.hazard"),
                                  cell_spec(h02.time, g, cc, "h02.time"))
      # The post-event no-switch hazard defaults to the direct terminal hazard.
      h12 <- if (is.null(h12.hazard)) h02 else
        piecewise_precompute(cell_spec(h12.hazard, g, cc, "h12.hazard"),
                             cell_spec(h12.time, g, cc, "h12.time"))
      sw <- if (is.null(switch.prop)) 0 else
        cell_spec(switch.prop, g, cc, "switch.prop")
      if (is.null(sw)) sw <- 0
      if (!is.numeric(sw) || length(sw) != 1L || !is.finite(sw) || sw < 0 ||
          sw > 1) {
        stop("Each 'switch.prop' value must be a single probability in [0, 1]")
      }
      if (sw > 0 && is.null(h12.switch.hazard)) {
        stop("'h12.switch.hazard' (or 'h12.switch.median') must be supplied ",
             "when any 'switch.prop' is positive")
      }
      h12s <- if (sw > 0) {
        piecewise_precompute(
          cell_spec(h12.switch.hazard, g, cc, "h12.switch.hazard"),
          cell_spec(h12.switch.time, g, cc, "h12.switch.time"))
      } else {
        empty_spec
      }
      d <- if (has_dropout) {
        piecewise_precompute(cell_spec(d.hazard, g, cc, "d.hazard"),
                             cell_spec(d.time, g, cc, "d.time"))
      } else {
        empty_spec
      }
      list(h01_haz = h01$hazard, h01_fin = h01$fin_time, h01_cum = h01$cum_haz,
           h02_haz = h02$hazard, h02_fin = h02$fin_time, h02_cum = h02$cum_haz,
           h12_haz = h12$hazard, h12_fin = h12$fin_time, h12_cum = h12$cum_haz,
           h12s_haz = h12s$hazard, h12s_fin = h12s$fin_time,
           h12s_cum = h12s$cum_haz,
           d_haz = d$hazard, d_fin = d$fin_time, d_cum = d$cum_haz,
           sw = as.numeric(sw))
    })
  }
  spec_c <- build_group(1L)
  spec_t <- if (two) build_group(2L) else list()

  acc <- resolve_accrual_counts(n_grp_int, a.time, a.rate, a.prop)
  fixed_c <- if (fixed.alloc) {
    fixed_cell_counts(n_grp_int[1L], spec_ctrl$cell_prob)
  } else {
    integer(n_cell)
  }
  fixed_t <- if (fixed.alloc && two) {
    fixed_cell_counts(n_grp_int[2L], spec_trt$cell_prob)
  } else {
    integer(n_cell)
  }
  level_tab <- function(spec) {
    matrix(as.integer(spec$level_table), nrow = nrow(spec$level_table))
  }

  simdata_core_id_sub(
    as.integer(nsim), n_grp_int, as.numeric(acc$a.time_full),
    acc$acc_counts_c, acc$acc_counts_t,
    spec_c, spec_t, has_dropout,
    as.numeric(spec_ctrl$cum_prev), as.numeric(spec_trt$cum_prev),
    level_tab(spec_ctrl), level_tab(spec_trt),
    as.character(spec_ctrl$sub_names),
    isTRUE(fixed.alloc),
    as.integer(fixed_c), as.integer(fixed_t)
  )
}
