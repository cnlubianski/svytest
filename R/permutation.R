#' Permutation test for weight informativeness in survey regression
#'
#' Non-parametric test that permutes survey weights (optionally within blocks)
#' to generate the null distribution of a chosen statistic. Supports fast closed-form
#' WLS (linear case) via C++ and a pure R engine.
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
#' @param B Number of permutations (e.g., 1000).
#' @param coef_subset Optional character vector of coefficient names to include.
#' @param block Optional factor for blockwise permutations (e.g., strata), permute within levels.
#' @param normalize Logical; if TRUE (default), normalize weights to have mean 1.
#' @param engine \code{"C++"} for fast WLS or \code{"R"} for pure R loop.
#' @param custom_fun Optional function(model, X, y, wts) -> scalar statistic (overrides \code{stat}).
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
#' @details
#' This procedure implements a non‑parametric randomization test for the
#' informativeness of survey weights. The null hypothesis is that, conditional
#' on the covariates \eqn{X}, the survey weights \eqn{w} are \emph{non‑informative}
#' with respect to the outcome \eqn{y}. Under this null, permuting the weights
#' across observations should not change the distribution of any statistic that
#' measures the effect of weighting.
#'
#' The algorithm is:
#' \enumerate{
#'   \item Fit the unweighted regression
#'     \deqn{\hat\beta_{U} = (X^\top X)^{-1} X^\top y}
#'     and the weighted regression
#'     \deqn{\hat\beta_{W} = (X^\top W X)^{-1} X^\top W y,}
#'     where \eqn{W = \mathrm{diag}(w)}.
#'
#'   \item Compute the observed test statistic \eqn{T_{\mathrm{obs}}}:
#'     \itemize{
#'       \item For \code{"pred_mean"}: the difference in mean predicted outcomes
#'             between weighted and unweighted fits.
#'       \item For \code{"coef_mahal"}: the Mahalanobis distance
#'             \deqn{T = (\hat\beta_{W} - \hat\beta_{U})^\top
#'                        (X^\top X)(\hat\beta_{W} - \hat\beta_{U}),}
#'             using the unweighted precision matrix as the metric.
#'       \item For a user‑supplied \code{custom_fun}, any scalar function of
#'             \eqn{(X,y,w)}.
#'     }
#'
#'   \item Generate the null distribution by permuting the weights:
#'     \deqn{w^{*(b)} = P_b w, \quad b=1,\ldots,B,}
#'     where each \eqn{P_b} is a permutation matrix. If a \code{block} factor
#'     is supplied, permutations are restricted within block levels.
#'
#'   \item Recompute the test statistic \eqn{T^{*(b)}} for each permuted weight
#'     vector. The empirical distribution of \eqn{T^{*(b)}} represents the null
#'     distribution under non‑informative weights.
#'
#'   \item The two‑sided permutation p‑value is
#'     \deqn{p = \frac{1 + \sum_{b=1}^B I\{|T^{*(b)} - T_0| \ge |T_{\mathrm{obs}} - T_0|\}}
#'                 {B+1},}
#'     where \eqn{T_0} is the baseline statistic under equal weights.
#' }
#'
#' Intuitively, if the weights are informative, the observed statistic will lie
#' in the tails of the permutation distribution, leading to a small p‑value.
#' If the weights are non‑informative, shuffling them destroys any spurious
#' association with the outcome, and the observed statistic is typical of the
#' permutation distribution.
#'
#' @examples
#' # Load in survey package (required) and load in example data
#' library(survey)
#' data(api, package = "survey")
#'
#' # Create a survey design and fit a weighted regression model
#' des <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat)
#' fit <- svyglm(api00 ~ ell + meals, design = des)
#'
#' # Run permutation diagnostic test; reports permutation statistics with p-value
#' results <- perm_test(fit, stat = "pred_mean", B = 1000, engine = "R")
#' print(results)
#'
#' @useDynLib svytest, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom survey svyglm
#'
#' @export
perm_test <- function(model, stat = c("pred_mean", "coef_mahal"), B = 1000,
                      coef_subset = NULL, block = NULL, normalize = TRUE,
                      engine = c("C++", "R"), custom_fun = NULL,
                      na.action = stats::na.omit) {

  # Argument checks
  if (!inherits(model, "svyglm")) stop("Model must be of class 'svyglm'.")
  valid_stats <- c("pred_mean", "coef_mahal")
  stat <- tolower(if (length(stat) == 1L) stat else stat[1L])
  if (!is.character(stat) || !stat %in% valid_stats) {
    stop("Invalid `stat`: must be one of ", paste(valid_stats, collapse = ", "), ".")
  }
  valid_engines <- c("c++", "r")
  engine <- tolower(if (length(engine) == 1L) engine else engine[1L])
  if (!is.character(engine) || !engine %in% valid_engines) {
    stop("Invalid `engine`: must be one of 'C++', 'R'.")
  }
  if (length(B) != 1L || !is.numeric(B) || is.na(B) || B < 1 || B != as.integer(B)) {
    stop("`B` must be a positive integer (number of permutations).")
  }
  B <- as.integer(B)
  if (!is.null(coef_subset) && !is.character(coef_subset)) {
    stop("`coef_subset` must be a character vector of coefficient names.")
  }
  if (!is.null(block) && !(is.factor(block) || is.character(block) || is.numeric(block))) {
    stop("`block` must be NULL or a vector coercible to factor (e.g., strata labels).")
  }
  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize)) {
    stop("`normalize` must be a single logical value (TRUE/FALSE).")
  }
  if (!is.null(custom_fun) && !is.function(custom_fun)) {
    stop("`custom_fun` must be NULL or a function (model, X, y, wts) -> numeric scalar.")
  }
  if (!is.function(na.action)) {
    stop("`na.action` must be a function (e.g., stats::na.omit).")
  }
  fam_name <- tryCatch(stats::family(model)$family, error = function(e) NA_character_)
  if (!is.na(fam_name) && !grepl("gaussian", fam_name, ignore.case = TRUE)) {
    warning("Non-Gaussian family detected; permutation WLS is calibrated for Gaussian models.")
  }

  # Data extraction
  wts <- tryCatch(stats::weights(model$survey.design),
                  error = function(e) stop("Could not extract weights from `model$survey.design`."))
  X <- tryCatch(stats::model.matrix(model),
                error = function(e) stop("Could not extract model matrix from `model`."))
  y <- tryCatch(stats::model.response(stats::model.frame(model)),
                error = function(e) stop("Could not extract response from `model`."))
  if (length(y) != nrow(X) || length(y) != length(wts)) {
    stop("Mismatch in dimensions of response, design matrix, and weights.")
  }

  # Missing data handling
  dat <- data.frame(y = y, X, wts = wts)
  dat <- na.action(dat)
  y   <- dat$y
  X   <- as.matrix(dat[, setdiff(names(dat), c("y", "wts"))])
  wts <- dat$wts
  n   <- length(y)

  # Subset coefficients
  if (!is.null(coef_subset)) {
    keep_cols <- colnames(X) %in% coef_subset
    if (!any(keep_cols)) {
      stop("None of the names in `coef_subset` matched the model coefficients: ",
           paste(colnames(X), collapse = ", "))
    }
    X <- X[, keep_cols, drop = FALSE]
  }

  # Block validation
  if (!is.null(block)) {
    if (length(block) != n) stop("`block` length must equal the number of observations (", n, ").")
    block <- droplevels(as.factor(block))
    if (nlevels(block) < 1L) stop("`block` must define at least one level.")
  }

  # Weights and baselines
  w_null <- rep(mean(wts), n)
  w_use  <- if (normalize) (wts / mean(wts)) else wts

  fit_null <- .wls_fit(X, y, if (normalize) (w_null / mean(w_null)) else w_null)
  fit_alt  <- .wls_fit(X, y, w_use)

  stat_obs <- if (stat == "pred_mean") {
    .pred_mean_stat(X, fit_alt$beta)
  } else {
    .coef_mahal_stat(fit_alt$beta, fit_null$beta, fit_null$XtX)
  }
  stat_null <- if (stat == "pred_mean") .pred_mean_stat(X, fit_null$beta) else 0

  # Permutation indices
  if (!is.null(block)) {
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
    check_scalar <- function(val) {
      if (!is.numeric(val) || length(val) != 1L || is.na(val)) {
        stop("`custom_fun` must return a single non-missing numeric value.")
      }
      val
    }
    perm_stats <- numeric(B)
    for (b in seq_len(B)) {
      w_b <- w_use[perm_matrix[, b]]
      perm_stats[b] <- check_scalar(custom_fun(model, X, y, w_b))
    }

  } else if (engine == "c++") {
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
      stop("C++ engine requires Rcpp; install Rcpp or use engine = 'R'.")
    }
    cpp_stats <- c("pred_mean", "coef_mahal")
    if (!stat %in% cpp_stats) {
      stop("C++ engine does not support stat = '", stat, "'. Use engine = 'R' or a supported stat.")
    }
    perm_stats <- perm_stats_cpp(
      X      = X,
      y      = y,
      w      = w_use,
      w_null = if (normalize) (w_null / mean(w_null)) else w_null,
      perm   = perm_matrix,
      stat   = stat
    )

  } else {
    # Pure R engine (serial only)
    calc_one <- function(col_idx) {
      idx <- perm_matrix[, col_idx]
      .perm_eval_one(idx, X, y, w_use, fit_null$beta, fit_null$XtX, stat)
    }
    perm_stats <- sapply(seq_len(B), calc_one)
  }

  # P-value
  diffs <- abs(perm_stats - stat_null)
  obs_diff <- abs(stat_obs - stat_null)
  pval <- (1 + sum(diffs >= obs_diff)) / (B + 1)

  structure(list(
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
#' @return A base R side effect plot. Function returns \code{NULL} invisibly.
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


#' Internal helper: Weighted least squares fit
#'
#' Performs a weighted least squares regression using QR decomposition,
#' with a generalized inverse fallback if the system is singular.
#'
#' @param X Numeric design matrix.
#' @param y Numeric response vector.
#' @param w Numeric vector of weights.
#'
#' @return A list with elements \code{beta}, \code{mu}, \code{resid},
#'   \code{sigma2}, and \code{XtX}.
#' @keywords internal
.wls_fit <- function(X, y, w) {
  WX   <- X * w
  XtWX <- crossprod(X, WX)
  XtWy <- crossprod(X, w * y)
  beta <- tryCatch(
    {
      b <- qr.solve(XtWX, XtWy)
      as.matrix(b)
    },
    error = function(e) MASS::ginv(XtWX) %*% XtWy
  )
  mu     <- as.numeric(X %*% beta)
  resid  <- y - mu
  sigma2 <- sum(w * resid^2) / sum(w)
  list(beta = beta, mu = mu, resid = resid, sigma2 = sigma2, XtX = crossprod(X))
}

#' Internal helper: Predicted mean statistic
#'
#' Computes the mean of fitted values given design matrix and coefficients.
#'
#' @param X Numeric design matrix.
#' @param beta Coefficient vector.
#'
#' @return Numeric scalar, mean of predictions.
#' @keywords internal
.pred_mean_stat <- function(X, beta) {
  mean(as.numeric(X %*% beta))
}

#' Internal helper: Mahalanobis coefficient distance
#'
#' Computes the Mahalanobis distance between two coefficient vectors
#' using a supplied precision matrix.
#'
#' @param beta Estimated coefficient vector.
#' @param beta0 Baseline coefficient vector.
#' @param XtX Precision matrix (X'X).
#'
#' @return Numeric scalar distance.
#' @keywords internal
.coef_mahal_stat <- function(beta, beta0, XtX) {
  d <- beta - beta0
  as.numeric(t(d) %*% XtX %*% d)
}

#' Internal helper: Evaluate one permutation
#'
#' Computes the test statistic for a single permutation of weights.
#'
#' @param idx Integer index vector for permuted weights.
#' @param X Numeric design matrix.
#' @param y Numeric response vector.
#' @param w_use Normalized weight vector.
#' @param fit_null_beta Baseline coefficient vector.
#' @param fit_null_XtX Baseline precision matrix.
#' @param stat Character string, statistic type.
#'
#' @return Numeric scalar statistic.
#' @keywords internal
.perm_eval_one <- function(idx, X, y, w_use, fit_null_beta, fit_null_XtX, stat) {
  w_b   <- w_use[idx]
  fit_b <- .wls_fit(X, y, w_b)
  if (stat == "pred_mean") {
    .pred_mean_stat(X, fit_b$beta)
  } else {
    .coef_mahal_stat(fit_b$beta, fit_null_beta, fit_null_XtX)
  }
}


