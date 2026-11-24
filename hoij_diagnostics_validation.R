# Diagnostics and targeted bootstrap-refit validation for hoij_lavaan()
# Compatible with the HOIJ-lavaan source interface examined on 2026-11-23.
# Source after hoij_core.R and hoij_lavaan.R.
# Requires lavaan only for hoij_validate(). Uses base R otherwise.
# These diagnostics do not change or trim the primary HOIJ-2 interval.
# Updated 2026-11-23: matched-weight random-subset percentile interval
# comparisons, method-specific boundary errors and ij diagnostics.
# The selected random subset is kept separate from targeted selections.

.hd_check <- function(x, label = "res", expected_order = NULL) {
  if (!is.list(x) || is.null(x$replicates) || is.null(x$weights))
    stop(label, " must come from hoij_lavaan(..., details = TRUE).",
         call. = FALSE)
  if (!is.matrix(x$replicates) || !is.matrix(x$weights) ||
      nrow(x$replicates) != nrow(x$weights) ||
      is.null(colnames(x$replicates)))
    stop(label, " contains inconsistent replicates, weights or parameter names.",
         call. = FALSE)
  if (!is.null(expected_order) &&
      !identical(as.integer(x$diagnostics$order), as.integer(expected_order)))
    stop(label, " must have order = ", expected_order, ".", call. = FALSE)
  invisible(TRUE)
}

.hd_check_ij <- function(res, ij) {
  if (is.null(ij)) return(invisible(TRUE))
  .hd_check(ij, "ij", expected_order = 1L)
  if (!identical(res$weights, ij$weights))
    stop("res and ij must have identical weight matrices. Run both with ",
         "the same B and seed.", call. = FALSE)
  if (!identical(colnames(res$replicates), colnames(ij$replicates)))
    stop("The parameter orders in res and ij differ.", call. = FALSE)
  invisible(TRUE)
}

.hd_fns <- function(functional, th_names) {
  if (!exists(".hoij_make_functionals", mode = "function"))
    stop("First source hoij_core.R and hoij_lavaan.R.", call. = FALSE)
  fns <- .hoij_make_functionals(functional, th_names)
  if (anyDuplicated(names(fns)))
    stop("Functional names must be unique.", call. = FALSE)
  fns
}

.hd_eval <- function(x, fns) {
  B <- nrow(x$replicates)
  ans <- vapply(fns, function(fun) {
    vapply(seq_len(B), function(i) {
      z <- tryCatch(suppressWarnings(fun(x$replicates[i, ])),
                    error = function(e) NA_real_)
      if (!is.numeric(z) || length(z) != 1L) return(NA_real_)
      as.numeric(z)
    }, numeric(1))
  }, numeric(B))
  if (is.null(dim(ans))) ans <- matrix(ans, ncol = length(fns))
  colnames(ans) <- names(fns)
  ans
}

.hd_inadmissible <- function(x) {
  nm <- colnames(x$replicates)
  parts <- strsplit(nm, "~~", fixed = TRUE)
  vi <- which(vapply(parts, function(z) length(z) == 2L && z[1L] == z[2L],
                     logical(1)))
  if (!length(vi)) return(rep(FALSE, nrow(x$replicates)))
  apply(x$replicates[, vi, drop = FALSE] < 0, 1L,
        function(z) any(z, na.rm = TRUE))
}

.hd_keep <- function(x) {
  if (identical(x$diagnostics$admissibility, "drop"))
    !.hd_inadmissible(x)
  else rep(TRUE, nrow(x$replicates))
}

.hd_quantiles <- function(x, level = 0.95) {
  z <- x[is.finite(x)]
  if (length(z) < 2L) return(c(lower = NA_real_, upper = NA_real_))
  p <- c((1 - level) / 2, 1 - (1 - level) / 2)
  setNames(as.numeric(stats::quantile(z, p, names = FALSE, type = 7)),
           c("lower", "upper"))
}

.hd_bounds <- function(bounds, fns) {
  if (!is.list(bounds)) stop("bounds must be a named list.", call. = FALSE)
  if (length(bounds) && (is.null(names(bounds)) ||
      any(!nzchar(names(bounds))) ||
      !all(names(bounds) %in% names(fns))))
    stop("bounds must be named for functionals in 'functional'.", call. = FALSE)
  for (nm in names(bounds)) {
    z <- bounds[[nm]]
    if (!is.numeric(z) || length(z) != 2L || anyNA(z) || z[1] > z[2])
      stop("Each bound must be c(lower, upper) (possibly infinite).", call. = FALSE)
  }
  bounds
}

# Summary of HOIJ-2 replicates, optional first-order comparator and boundaries.
# 'n_below' and 'n_above' use *all finite HOIJ-2 values*, even when an
# optional admissibility='drop' rule was used for the primary intervals.
hoij_diagnostics <- function(res, functional = NULL, ij = NULL,
                             bounds = list(), fit = NULL, level = NULL) {
  .hd_check(res, expected_order = 2L)
  .hd_check_ij(res, ij)
  fns <- .hd_fns(functional, colnames(res$replicates))
  bounds <- .hd_bounds(bounds, fns)
  if (is.null(level)) level <- res$diagnostics$level
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must lie strictly between zero and one.", call. = FALSE)

  val <- .hd_eval(res, fns)
  val_ij <- if (!is.null(ij)) .hd_eval(ij, fns) else NULL
  keep <- .hd_keep(res)
  keep_ij <- if (!is.null(ij)) .hd_keep(ij) else NULL
  inad <- .hd_inadmissible(res)
  B <- nrow(val)

  tbl <- do.call(rbind, lapply(seq_along(fns), function(k) {
    nm <- names(fns)[k]
    v <- val[, k]
    valid <- is.finite(v)
    v_used <- v[valid & keep]
    q <- .hd_quantiles(v_used, level)
    has_bounds <- nm %in% names(bounds)
    lim <- if (has_bounds) bounds[[nm]] else c(-Inf, Inf)
    within <- v_used[v_used >= lim[1] & v_used <= lim[2]]
    q_bounded <- if (has_bounds) .hd_quantiles(within, level) else
      c(lower = NA_real_, upper = NA_real_)
    below <- if (is.finite(lim[1])) sum(valid & v < lim[1]) else NA_integer_
    above <- if (is.finite(lim[2])) sum(valid & v > lim[2]) else NA_integer_

    q1 <- c(lower = NA_real_, upper = NA_real_)
    median_delta <- p95_delta <- NA_real_
    if (!is.null(val_ij)) {
      vi <- val_ij[, k]
      q1 <- .hd_quantiles(vi[is.finite(vi) & keep_ij], level)
      paired <- is.finite(v) & is.finite(vi) & keep & keep_ij
      if (any(paired)) {
        abs_diff <- abs(v[paired] - vi[paired])
        median_delta <- stats::median(abs_diff)
        p95_delta <- unname(stats::quantile(abs_diff, 0.95))
      }
    }
    data.frame(
      functional = nm,
      n_generated = B,
      n_finite = sum(valid),
      n_nonfinite = sum(!valid),
      n_used = length(v_used),
      n_inadmissible_parameters = sum(inad),
      n_below = below,
      n_above = above,
      pct_below_of_finite = if (is.na(below) || !any(valid)) NA_real_
                            else 100 * below / sum(valid),
      pct_above_of_finite = if (is.na(above) || !any(valid)) NA_real_
                            else 100 * above / sum(valid),
      hoij_lo = unname(q[1]), hoij_hi = unname(q[2]),
      restricted_lo = unname(q_bounded[1]),
      restricted_hi = unname(q_bounded[2]),
      n_restricted = if (has_bounds) length(within) else NA_integer_,
      ij_lo = unname(q1[1]), ij_hi = unname(q1[2]),
      median_abs_hoij_minus_ij = median_delta,
      p95_abs_hoij_minus_ij = p95_delta,
      stringsAsFactors = FALSE
    )
  }))
  rownames(tbl) <- NULL

  numerical <- list(
    derivative_spread = res$diagnostics$grad_spread,
    observed_information_condition_number = NA_real_
  )
  if (!is.null(fit)) {
    H <- lavaan::lavTech(fit, "information.observed")
    numerical$observed_information_condition_number <-
      tryCatch(kappa(H, exact = TRUE), error = function(e) NA_real_)
  }

  structure(list(summary = tbl,
                 numerical = numerical,
                 values = val,
                 ij_values = val_ij,
                 weights = res$weights,
                 admissible_for_primary_interval = keep,
                 note = paste(
                   "Bounds and ij comparisons are descriptive diagnostics.",
                   "A large correction alone is not evidence of inaccuracy.",
                   "Restricted intervals are sensitivity analyses, not replacements.",
                   "All finite replicates remain in the default interval.")),
            class = "hoij_diagnostics")
}

print.hoij_diagnostics <- function(x, ...) {
  print(x$summary, row.names = FALSE, ...)
  cat("Derivative spread:", x$numerical$derivative_spread, "\n")
  if (is.finite(x$numerical$observed_information_condition_number))
    cat("Condition number (scale-dependent):",
        x$numerical$observed_information_condition_number, "\n")
  invisible(x)
}

.hd_data <- function(fit, data = NULL) {
  X <- fit@Data@X[[1L]]
  nm <- fit@Data@ov.names[[1L]]
  if (length(nm) != ncol(X) || anyDuplicated(nm))
    stop("Cannot determine the internal observed-variable order.", call. = FALSE)
  if (is.null(data)) {
    data <- as.data.frame(X)
    names(data) <- nm
  } else {
    if (!is.data.frame(data)) data <- as.data.frame(data)
    if (!all(nm %in% names(data)))
      stop("data lacks observed variables: ",
           paste(setdiff(nm, names(data)), collapse = ", "), call. = FALSE)
    data <- data[, nm, drop = FALSE]
  }
  if (nrow(data) != nrow(X) || anyNA(data))
    stop("data must contain the original complete cases in their original order.",
         call. = FALSE)
  check <- isTRUE(all.equal(unname(as.matrix(data)), unname(X),
                            tolerance = 1e-10, check.attributes = FALSE))
  if (!check)
    stop("data and the internal model data differ in values or row order.",
         call. = FALSE)
  data
}

.hd_check_count <- function(x, name) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x < 0 || x != as.integer(x))
    stop(name, " must be a nonnegative integer.", call. = FALSE)
  as.integer(x)
}

# Exact refits of a stratified *diagnostic subset* of HOIJ-2 weight vectors.
# Random and targeted subsets are reported separately: do not pool them
# to estimate whole-distribution accuracy or percentile coverage.
# 'n_tail' is the TOTAL number of near-limit indices per functional
# (approximately half near each HOIJ-2 percentile).
# 'n_correction' is per functional and requires an aligned ij object.
hoij_validate <- function(fit, res, functional = NULL, ij = NULL,
                          data = NULL, n_random = 100L, n_tail = 0L,
                          n_correction = 0L, seed = 42L,
                          level = NULL) {
  .hd_check(res, expected_order = 2L)
  .hd_check_ij(res, ij)
  .hoij_check_fit(fit)
  if (is.null(level)) level <- res$diagnostics$level
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must lie strictly between zero and one.", call. = FALSE)
  n_random <- .hd_check_count(n_random, "n_random")
  n_tail <- .hd_check_count(n_tail, "n_tail")
  n_correction <- .hd_check_count(n_correction, "n_correction")
  if (n_correction > 0L && is.null(ij))
    stop("n_correction > 0 requires an aligned ij object.", call. = FALSE)

  fns <- .hd_fns(functional, colnames(res$replicates))
  v2 <- .hd_eval(res, fns)
  v1 <- if (!is.null(ij)) .hd_eval(ij, fns) else NULL
  keep2 <- .hd_keep(res)
  keep1 <- if (!is.null(ij)) .hd_keep(ij) else NULL
  B <- nrow(res$weights)
  N <- ncol(res$weights)
  original <- .hd_data(fit, data)
  if (nrow(original) != N)
    stop("The number of columns in the weight matrix must equal N.",
         call. = FALSE)
  if (!identical(colnames(res$replicates),
                 names(lavaan::coef(fit, type = "free"))))
    stop("Parameter names/order do not match the fitted model.", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)
  pick_random <- if (n_random > 0L)
    sample.int(B, min(n_random, B), replace = FALSE) else integer(0)
  flag_random <- flag_tail <- flag_correction <- rep(FALSE, B)
  flag_random[pick_random] <- TRUE
  tail_names <- correction_names <- rep("", B)
  qpr <- c((1 - level) / 2, 1 - (1 - level) / 2)

  for (k in seq_along(fns)) {
    nm <- names(fns)[k]
    v <- v2[, k]
    use <- which(is.finite(v) & keep2)
    if (length(use) > 1L && n_tail > 0L) {
      q <- stats::quantile(v[use], qpr, names = FALSE, type = 7)
      side_n <- c(ceiling(n_tail / 2), floor(n_tail / 2))
      for (side in seq_len(2L)) {
        if (side_n[side] == 0L) next
        ord <- order(abs(v[use] - q[side]), use)
        take <- head(use[ord], side_n[side])
        flag_tail[take] <- TRUE
        tag <- paste0(nm, if (side == 1L) ":lower" else ":upper")
        tail_names[take] <- ifelse(nzchar(tail_names[take]),
                                   paste(tail_names[take], tag, sep = ";"), tag)
      }
    }
    if (n_correction > 0L) {
      use <- which(is.finite(v) & is.finite(v1[, k]) & keep2 & keep1)
      if (length(use)) {
        ord <- order(abs(v[use] - v1[use, k]), decreasing = TRUE)
        take <- head(use[ord], n_correction)
        flag_correction[take] <- TRUE
        correction_names[take] <- ifelse(nzchar(correction_names[take]),
          paste(correction_names[take], nm, sep = ";"), nm)
      }
    }
  }

  selected <- which(flag_random | flag_tail | flag_correction)
  if (!length(selected)) stop("No weight vectors were selected.", call. = FALSE)
  D <- ncol(res$replicates)
  boot_theta <- matrix(NA_real_, length(selected), D,
                       dimnames = list(as.character(selected),
                                       colnames(res$replicates)))
  refit_status <- rep("not_started", length(selected))
  refit_message <- rep("", length(selected))
  for (j in seq_along(selected)) {
    b <- selected[j]
    w <- as.numeric(res$weights[b, ])
    if (any(!is.finite(w)) || any(w < 0) ||
        any(abs(w - round(w)) > 1e-8) || abs(sum(w) - N) > 1e-8) {
      refit_status[j] <- "invalid_weights"
      refit_message[j] <- "Weights must be nonnegative integer counts summing to N."
      next
    }
    rows <- rep(seq_len(N), times = as.integer(round(w)))
    attempt <- tryCatch(
      suppressWarnings(stats::update(fit, data = original[rows, , drop = FALSE])),
      error = function(e) e
    )
    if (inherits(attempt, "error")) {
      refit_status[j] <- "refit_error"
      refit_message[j] <- conditionMessage(attempt)
      next
    }
    if (!isTRUE(lavaan::lavInspect(attempt, "converged"))) {
      refit_status[j] <- "nonconverged"
      refit_message[j] <- "lavaan reports nonconvergence."
      next
    }
    th <- lavaan::coef(attempt, type = "free")
    if (!identical(names(th), colnames(res$replicates))) {
      refit_status[j] <- "parameter_mismatch"
      refit_message[j] <- "Refit parameter names/order differ."
      next
    }
    if (any(!is.finite(th))) {
      refit_status[j] <- "nonfinite_parameters"
      refit_message[j] <- "Refit includes nonfinite parameter values."
      next
    }
    boot_theta[j, ] <- th
    refit_status[j] <- "converged"
  }

  vb <- matrix(NA_real_, nrow(boot_theta), length(fns),
               dimnames = list(NULL, names(fns)))
  for (j in which(refit_status == "converged")) {
    for (k in seq_along(fns)) {
      z <- tryCatch(suppressWarnings(fns[[k]](boot_theta[j, ])),
                    error = function(e) NA_real_)
      if (is.numeric(z) && length(z) == 1L) vb[j, k] <- as.numeric(z)
    }
  }

  reps <- data.frame(
    weight_index = selected,
    random = flag_random[selected],
    near_limits = flag_tail[selected],
    large_correction = flag_correction[selected],
    near_limits_for = tail_names[selected],
    large_correction_for = correction_names[selected],
    refit_status = refit_status,
    refit_message = refit_message,
    stringsAsFactors = FALSE
  )
  pairs <- do.call(rbind, lapply(seq_along(fns), function(k) {
    h <- v2[selected, k]
    b <- vb[, k]
    i <- if (!is.null(v1)) v1[selected, k] else rep(NA_real_, length(selected))
    data.frame(reps,
      functional = names(fns)[k],
      bootstrap = b, hoij2 = h, ij = i,
      error_hoij2 = ifelse(is.finite(b) & is.finite(h), h - b, NA_real_),
      error_ij = ifelse(is.finite(b) & is.finite(i), i - b, NA_real_),
      stringsAsFactors = FALSE
    )
  }))
  rownames(pairs) <- NULL

  # Overlapping groups are intentional. Never summarize the mixture
  # of random and targeted selections as if it were a random sample.
  groups <- c("random", "near_limits", "large_correction")
  summaries <- list()
  for (nm in names(fns)) for (gr in groups) {
    belongs <- if (gr == "random") {
      pairs$random
    } else if (gr == "near_limits") {
      grepl(paste0(nm, ":"), pairs$near_limits_for, fixed = TRUE)
    } else {
      vapply(strsplit(pairs$large_correction_for, ";", fixed = TRUE),
             function(tags) nm %in% tags, logical(1))
    }
    z <- pairs[pairs$functional == nm & belongs, , drop = FALSE]
    if (!nrow(z)) next
    good <- is.finite(z$error_hoij2)
    err <- z$error_hoij2[good]
    ib <- z$bootstrap[good]
    sd_b <- if (length(ib) >= 2L) stats::sd(ib) else NA_real_
    matched_ij <- good & is.finite(z$error_ij)
    e1 <- z$error_ij[matched_ij]
    e2 <- z$error_hoij2[matched_ij]
    summaries[[length(summaries) + 1L]] <- data.frame(
      functional = nm,
      subset = gr,
      n_selected = nrow(z),
      n_refits_converged = sum(z$refit_status == "converged"),
      n_finite_pairs = sum(good),
      mean_error = if (length(err)) mean(err) else NA_real_,
      median_absolute_error = if (length(err)) stats::median(abs(err)) else NA_real_,
      RMSE = if (length(err)) sqrt(mean(err^2)) else NA_real_,
      max_absolute_error = if (length(err)) max(abs(err)) else NA_real_,
      relative_RMSE = if (is.finite(sd_b) && sd_b > 0)
                        sqrt(mean(err^2)) / sd_b else NA_real_,
      n_finite_three_way = length(e1),
      median_absolute_error_ij = if (length(e1)) stats::median(abs(e1))
                                  else NA_real_,
      median_absolute_error_hoij2_matched =
        if (length(e2)) stats::median(abs(e2)) else NA_real_,
      RMSE_ij = if (length(e1)) sqrt(mean(e1^2)) else NA_real_,
      RMSE_hoij2_matched = if (length(e2)) sqrt(mean(e2^2))
                           else NA_real_,
      proportion_hoij2_better = if (length(e1)) mean(abs(e2) < abs(e1))
                                else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL

  # Confidence limits MUST use the random subset only.  For a fair
  # method-to-method comparison, retain the same weight indices with finite
  # bootstrap and HOIJ-2 values (and ij values if supplied).
  # Targeted weights are never included here unless they also happened to
  # be selected by the independent random draw.
  interval_rows <- list()
  error_rows <- list()
  for (nm in names(fns)) {
    zz <- pairs[pairs$functional == nm & pairs$random, , drop = FALSE]
    if (!nrow(zz)) next
    three_way <- !is.null(ij)
    matched <- is.finite(zz$bootstrap) & is.finite(zz$hoij2) &
               keep2[zz$weight_index]
    if (three_way) matched <- matched & is.finite(zz$ij) &
                              keep1[zz$weight_index]
    use <- zz[matched, , drop = FALSE]
    n_pair <- nrow(use)
    # For n < 40, quantiles are intentionally not calculated.  Even at
    # n = 1000, only approximately 25 observations determine each tail.
    ci <- function(x) if (n_pair >= 40L) .hd_quantiles(x, level) else
                       c(lower = NA_real_, upper = NA_real_)
    bci <- ci(use$bootstrap)
    for (method in c("bootstrap", if (three_way) "ij", "HOIJ-2")) {
      ci_m <- switch(method, bootstrap = bci,
                     ij = ci(use$ij), `HOIJ-2` = ci(use$hoij2))
      interval_rows[[length(interval_rows) + 1L]] <- data.frame(
        functional = nm,
        method = method,
        n_random_selected = nrow(zz),
        n_refits_converged = sum(zz$refit_status == "converged"),
        n_bootstrap_finite = sum(is.finite(zz$bootstrap)),
        n_common = n_pair,
        approx_tail_count = n_pair * (1 - level) / 2,
        lower = unname(ci_m[1L]), upper = unname(ci_m[2L]),
        width = unname(ci_m[2L] - ci_m[1L]),
        stringsAsFactors = FALSE
      )
      if (method != "bootstrap") {
        error_rows[[length(error_rows) + 1L]] <- data.frame(
          functional = nm,
          method = method,
          n_common = n_pair,
          lower_error = unname(ci_m[1L] - bci[1L]),
          upper_error = unname(ci_m[2L] - bci[2L]),
          sum_abs_boundary_error = unname(sum(abs(ci_m - bci))),
          max_abs_boundary_error = unname(max(abs(ci_m - bci))),
          width_error = unname(diff(ci_m) - diff(bci)),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  intervals <- if (length(interval_rows)) do.call(rbind, interval_rows)
               else data.frame()
  boundary_errors <- if (length(error_rows)) do.call(rbind, error_rows)
                     else data.frame()
  rownames(intervals) <- NULL
  rownames(boundary_errors) <- NULL

  structure(list(summary = summary, intervals = intervals,
                 boundary_errors = boundary_errors,
                 selection = reps, pairs = pairs,
                 bootstrap_parameters = boot_theta,
                 note = paste(
                   "Percentile limits in $intervals are based ONLY on",
                   "the random subset and the SAME finite, matched weight vectors",
                   "for each method. The $boundary_errors are relative to the",
                   "bootstrap limits within that random subset.",
                   "If refits fail or method-specific exclusion rules are used,",
                   "the common-weight subset may be selective.",
                   "Near-limit and large-correction subsets are deliberately",
                   "selected and MUST NOT be pooled with random replicates.",
                   "Estimated limits are subject to Monte Carlo error.",
                   "This is neither a coverage study nor a proof of accuracy.")),
            class = "hoij_validation")
}

print.hoij_validation <- function(x, ...) {
  cat("Replicate-level approximation errors by selection rule:\n")
  print(x$summary, row.names = FALSE, ...)
  if (nrow(x$intervals)) {
    cat("\nPercentile intervals, RANDOM subset, matched finite weights:\n")
    print(x$intervals[, c("functional", "method", "n_common",
                          "approx_tail_count", "lower", "upper", "width")],
          row.names = FALSE, ...)
    cat("\nBoundary errors (method minus bootstrap; random subset):\n")
    print(x$boundary_errors, row.names = FALSE, ...)
  }
  cat("\nUnique exact refits:", nrow(x$selection), "\n")
  cat("Refit statuses:\n")
  print(table(x$selection$refit_status))
  cat("\nInterpretation: random-subset quantile comparisons have Monte Carlo",
      "uncertainty; targeted subsets are not representative.\n")
  invisible(x)
}

plot_hoij_validation <- function(x, functional,
                                 subset = c("random", "near_limits",
                                            "large_correction"), ...) {
  subset <- match.arg(subset)
  belongs <- if (subset == "random") {
    x$pairs$random
  } else if (subset == "near_limits") {
    grepl(paste0(functional, ":"), x$pairs$near_limits_for,
          fixed = TRUE)
  } else {
    vapply(strsplit(x$pairs$large_correction_for, ";", fixed = TRUE),
           function(tags) functional %in% tags, logical(1))
  }
  z <- x$pairs[x$pairs$functional == functional & belongs, , drop = FALSE]
  z <- z[is.finite(z$bootstrap) & is.finite(z$hoij2), , drop = FALSE]
  if (!nrow(z)) stop("No finite paired values for this selection.", call. = FALSE)
  plot(z$bootstrap, z$hoij2, xlab = "Exact bootstrap refit",
       ylab = "HOIJ-2", main = paste(functional, "-", subset), pch = 19, ...)
  abline(0, 1, lty = 2)
  invisible(z)
}
