#' Permutation test for weight informativeness in survey regression
#'
#' Non-parametric test that permutes survey weights (optionally within blocks)
#' to generate the null distribution of a chosen statistic. Supports fast closed-form
#' WLS (linear case) via C++ and a pure R engine (with optional future.apply parallelization).
#'
#' @param model An object of class \code{svyglm} (currently supports Gaussian family best).
#' @param stat Statistic to use. Options:
#'   \itemize{
#'     \item \code{"pred_mean"}:
#'       Compares the mean predicted outcome under weighted vs. unweighted regression.
#'       Simple, interpretable, and directly tied to differences in fitted population means.
#'       Sensitive to shifts in overall prediction levels caused by informative weights.
#'
#'     \item \code{"coef_mahal"}:
#'       Computes the Mahalanobis distance between the weighted and unweighted coefficient vectors,
#'       using the unweighted precision matrix (\eqn{X'X}) as the metric. Captures joint shifts in regression coefficients,
#'       not just mean predictions. More powerful when informativeness manifests as changes in slopes
#'       or multiple coefficients simultaneously.
#'   }
#'
#' @param B Number of permutations (e.g., 1000).
#' @param coef_subset Optional character vector of coefficient names to include.
#' @param block Optional factor for blockwise permutations (e.g., strata), permute within levels.
#' @param normalize Logical; if TRUE (default), normalize weights to have mean 1.
#' @param engine \code{"C++"} for fast WLS or \code{"R"} for pure R loop.
#' @param custom_fun Optional function(model, X, y, wts) -> scalar statistic (overrides \code{stat}).
#' @param parallel Logical; if TRUE and engine is "R", uses future.apply for parallel permutations.
#' @param na.action Function to handle missing data.
#'
#' @return An object of class \code{"perm_test"} with fields:
#'   \item{stat_obs}{Observed statistic with actual weights}
#'   \item{stat_null}{Baseline statistic under equal weights (for centering)}
#'   \item{perm_stats}{Vector of permutation statistics}
#'   \item{p.value}{Permutation p-value (two-sided, centered at baseline)}
#'   \item{effect}{Observed minus median of permutation stats}
#'   \item{stat}{Statistic name}
#'   \item{B}{Number of permutations}
#'   \item{call}{Matched call}
#'   \item{method}{Description string}
#'
#' @importFrom Rcpp evalCpp
#'
#' @export
perm_test <- function(model, stat = c("pred_mean", "coef_mahal"), B = 1000,
                      coef_subset = NULL, block = NULL, normalize = TRUE,
                      engine = c("C++", "R"), custom_fun = NULL, parallel = FALSE,
                      na.action = stats::na.omit) {
  # Quick checks
  if (!inherits(model, "svyglm")) stop("Model must be of class 'svyglm'.")
  stat <- match.arg(stat)
  engine <- match.arg(engine)

  # Extract design matrix, response, weights
  wts <- stats::weights(model$survey.design)
  X <- stats::model.matrix(model)
  y <- stats::model.response(stats::model.frame(model))

  # Handle missing data in a unified frame
  dat <- data.frame(y = y, X, wts = wts)
  dat <- na.action(dat)
  y <- dat$y
  X <- as.matrix(dat[, setdiff(names(dat), c("y", "wts"))])
  wts <- dat$wts
  n <- length(y)

  # Optional subset of coefficients (columns of X)
  if (!is.null(coef_subset)) {
    keep_cols <- colnames(X) %in% coef_subset
    if (!any(keep_cols)) stop("No matching coefficients found in model.")
    X <- X[, keep_cols, drop = FALSE]
  }

  # Equal-weight baseline and chosen weights
  w_null <- rep(mean(wts), n)
  w_use <- if (normalize) (wts / mean(wts)) else wts

  # Fast WLS helper in R
  wls_fit <- function(X, y, w) {
    WX <- X * w  # row-wise multiply columns by w
    XtWX <- crossprod(X, WX)
    XtWy <- crossprod(X, w * y)
    beta <- solve(XtWX, XtWy)
    mu <- as.numeric(X %*% beta)
    resid <- y - mu
    # weighted sigma^2 (pseudo-ML): sum(w*e^2)/sum(w)
    sigma2 <- sum(w * resid^2) / sum(w)
    list(beta = beta, mu = mu, resid = resid, sigma2 = sigma2, XtX = crossprod(X))
  }

  # Statistic calculators
  pred_mean_stat <- function(beta) mean(as.numeric(X %*% beta))
  coef_mahal_stat <- function(beta, beta0, XtX) {
    # Mahalanobis distance with respect to unweighted precision X'X:
    d <- beta - beta0 # delta^T (X'X) delta
    as.numeric(t(d) %*% XtX %*% d)
  }

  # Baseline fits
  fit_null <- wls_fit(X, y, if (normalize) (w_null / mean(w_null)) else w_null)
  fit_alt  <- wls_fit(X, y, w_use)

  stat_obs <- switch(stat,
                     "pred_mean" = pred_mean_stat(fit_alt$beta),
                     "coef_mahal" = coef_mahal_stat(fit_alt$beta, fit_null$beta, fit_null$XtX)
  )
  stat_null <- switch(stat,
                      "pred_mean" = pred_mean_stat(fit_null$beta),
                      "coef_mahal" = 0 # baseline difference is 0
  )

  # Generate permutations (within blocks if provided)
  if (!is.null(block)) {
    block <- droplevels(as.factor(block))
    if (length(block) != n) stop("block must have length equal to number of observations.")
    # Create index permutations by shuffling within block
    perm_matrix <- matrix(NA_integer_, nrow = n, ncol = B)
    for (b in seq_len(B)) {
      idx <- seq_len(n)
      for (lev in levels(block)) {
        which_lev <- which(block == lev)
        idx[which_lev] <- sample(which_lev, length(which_lev), replace = FALSE)
      }
      perm_matrix[, b] <- idx
    }
  } else {
    perm_matrix <- replicate(B, sample.int(n, n, replace = FALSE))
  }

  # Compute permutation statistics
  if (!is.null(custom_fun)) {
    # User-supplied statistic evaluated under actual and permuted weights
    perm_stats <- numeric(B)
    for (b in seq_len(B)) {
      w_b <- w_use[perm_matrix[, b]]
      perm_stats[b] <- custom_fun(model, X, y, w_b)
    }
  } else if (engine == "C++") {
    # Use C++ engine for speed (linear WLS and provided stats)
    perm_stats <- perm_stats_cpp(
      X      = X,
      y      = y,
      w      = w_use,
      w_null = if (normalize) (w_null / mean(w_null)) else w_null,
      perm   = perm_matrix,
      stat   = stat
    )
  } else {
    # Pure R engine (optionally parallel via future.apply)
    # Requires: library(future.apply); plan(multisession) or similar set by user.
    calc_one <- function(col_idx) {
      w_b <- w_use[perm_matrix[, col_idx]]
      fit_b <- wls_fit(X, y, w_b)
      switch(stat,
             "pred_mean"  = pred_mean_stat(fit_b$beta),
             "coef_mahal" = coef_mahal_stat(fit_b$beta, fit_null$beta, fit_null$XtX)
      )
    }

    if (parallel) {
      if (!requireNamespace("future.apply", quietly = TRUE)) {
        warning("future.apply not available; falling back to serial execution.")
        perm_stats <- sapply(seq_len(B), calc_one)
      } else {
        perm_stats <- future.apply::future_sapply(seq_len(B), calc_one)
      }
    } else {
      perm_stats <- sapply(seq_len(B), calc_one)
    }
  }

  # Two-sided centered p-value relative to baseline
  diffs <- abs(perm_stats - stat_null)
  obs_diff <- abs(stat_obs - stat_null)
  pval <- (1 + sum(diffs >= obs_diff)) / (B + 1)

  out <- structure(list(
    stat_obs   = stat_obs,
    stat_null  = stat_null,
    perm_stats = perm_stats,
    p.value    = pval,
    effect     = stat_obs - stats::median(perm_stats),
    stat       = stat,
    B          = B,
    method     = paste0("Permutation test for weight informativeness (", stat, ")"),
    call       = match.call()
  ), class = "perm_test")

  out
}

#' @rdname perm_test
#' @method print perm_test
#' @param x An object of class perm_test
#' @param ... Additional arguments passed to methods
#' @export
print.perm_test <- function(x, ...) {
  cat("\n", x$method, "\n", sep = "")
  cat("Observed =", formatC(x$stat_obs, digits = 4, format = "f"),
      " Null =", formatC(x$stat_null, digits = 4, format = "f"),
      " Effect =", formatC(x$effect, digits = 4, format = "f"),
      " p-value =", formatC(x$p.value, digits = 4, format = "f"), "\n")
  invisible(x)
}

#' @rdname perm_test
#' @method summary perm_test
#' @param object An object of class perm_test
#' @param ... Additional arguments passed to methods
#' @export
summary.perm_test <- function(object, ...) {
  cat("\nPermutation test for weight informativeness\n")
  cat("Call:\n")
  print(object$call)
  cat("\nMethod:\n ", object$method, "\n", sep = "")
  cat("\nObserved statistic:\n ",
      formatC(object$stat_obs, digits = 6, format = "f"), "\n", sep = "")
  cat("Baseline (null) statistic:\n ",
      formatC(object$stat_null, digits = 6, format = "f"), "\n", sep = "")
  cat("Effect size (obs - median perm):\n ",
      formatC(object$effect, digits = 6, format = "f"), "\n", sep = "")
  cat("Permutations:\n ", object$B, "\n", sep = "")
  cat("p-value (two-sided, centered):\n ",
      formatC(object$p.value, digits = 6, format = "f"), "\n", sep = "")
  invisible(object)
}

#' @rdname perm_test
#' @method tidy perm_test
#' @param x An object of class perm_test
#' @param ... Additional arguments passed to methods
tidy.perm_test <- function(x, ...) {
  tibble::tibble(
    term      = paste0("perm_", x$stat),
    stat_obs  = x$stat_obs,
    stat_null = x$stat_null,
    effect    = x$effect,
    p.value   = x$p.value,
    B         = x$B,
    method    = x$method
  )
}

#' @rdname perm_test
#' @method glance perm_test
#' @param x An object of class perm_test
#' @param ... Additional arguments passed to methods
glance.perm_test <- function(x, ...) {
  tibble::tibble(
    stat      = x$stat,
    stat_obs  = x$stat_obs,
    stat_null = x$stat_null,
    effect    = x$effect,
    p.value   = x$p.value,
    B         = x$B,
    method    = x$method
  )
}

#' Plot method for permutation test objects
#'
#' Produces a histogram of the permutation distribution with a vertical line
#' indicating the observed statistic.
#'
#' @param x An object of class \code{"perm_test"}.
#' @param bins Number of histogram bins (default 30).
#' @param col Color for histogram bars (default "lightgray").
#' @param line_col Color for observed statistic line (default "red").
#' @param ... Additional arguments passed to \code{hist()}.
#'
#' @export
plot.perm_test <- function(x, bins = 30, col = "lightgray", line_col = "red", ...) {
  if (is.null(x$perm_stats) || is.null(x$statistic)) {
    stop("Object must contain 'perm_stats' and 'statistic'.")
  }

  graphics::hist(x$perm_stats,
                 breaks = bins,
                 col = col,
                 main = "Permutation Test Null Distribution",
                 xlab = "Test statistic",
                 ...)
  graphics::abline(v = x$statistic, col = line_col, lwd = 2)
  graphics::legend("topright",
                   legend = c("Observed statistic"),
                   col = line_col, lwd = 2, bty = "n")
}
