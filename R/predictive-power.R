#' Pfeffermann-Nathan Predictive Power Test (svyglm, K-fold CV, fold-mean option) (In production)
#'
#' Reimplements the predictive power test following Wang et al. (2023, Sec. 2.2):
#' split observations into estimation and validation sets; fit unweighted and weighted
#' linear regressions on the estimation set; compute validation squared-error differences
#' \eqn{D_i = (y_i - \hat y_{u,i})^2 - (y_i - \hat y_{w,i})^2}; test \eqn{H_0: E[D_i] = 0}
#' with \eqn{Z = \bar D / (s_D / \sqrt{n_V})}. Supports K-fold CV (default) and a
#' "fold-mean" option to reduce dependence among errors by using per-fold means as
#' the test observations.
#'
#' @param model A fitted \code{svyglm} with \code{family = gaussian(identity)}.
#' @param kfold Logical; if TRUE, use K-fold cross-validation (default TRUE).
#' @param K Integer number of folds (default 5).
#' @param est_split Proportion for estimation set if \code{kfold = FALSE} (default 0.5).
#' @param use_fold_means Logical; if TRUE (default), compute one \eqn{D} per fold
#'   as the mean of within-fold \eqn{D_i}, then form \eqn{Z} using the \eqn{K} fold
#'   means. This stabilizes the test by reducing dependence noted in Wang (2023).
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An object of class \code{"pred_power_test"} with fields:
#'   \item{statistic}{Z statistic}
#'   \item{p.value}{Two-sided p-value}
#'   \item{mean_diff}{Mean of \eqn{D} (fold mean if \code{use_fold_means = TRUE})}
#'   \item{n_val}{Count of observations used in Z (\eqn{K} if \code{use_fold_means = TRUE}, else total validation n)}
#'   \item{K}{Number of folds (if \code{kfold = TRUE})}
#'   \item{method}{Description string}
#'   \item{call}{Matched call}
#'
#' @keywords internal
pred_power_test <- function(model, kfold = TRUE, K = 5, est_split = 0.5,
                            use_fold_means = TRUE, seed = NULL) {
  if (!inherits(model, "svyglm")) stop("Model must be of class 'svyglm'.")
  fam <- model$family
  if (!(fam$family == "gaussian" && fam$link == "identity")) {
    stop("pred_power_test supports only gaussian(identity) models.")
  }
  if (!is.null(seed)) set.seed(seed)

  # Extract design components
  w <- as.numeric(weights(model$survey.design))
  X <- stats::model.matrix(model)
  y <- as.numeric(stats::model.response(stats::model.frame(model)))
  n <- nrow(X)

  # Helper to fit unweighted and weighted OLS on estimation set, predict on validation
  predict_pair <- function(est_idx, val_idx) {
    X_est <- X[est_idx, , drop = FALSE]
    y_est <- y[est_idx]
    w_est <- w[est_idx]

    X_val <- X[val_idx, , drop = FALSE]
    y_val <- y[val_idx]

    # Unweighted OLS
    beta_u <- solve(crossprod(X_est), crossprod(X_est, y_est))
    yhat_u <- as.numeric(X_val %*% beta_u)
    v_u <- y_val - yhat_u

    # Weighted OLS (WLS)
    W <- diag(w_est)
    beta_w <- solve(t(X_est) %*% W %*% X_est, t(X_est) %*% W %*% y_est)
    yhat_w <- as.numeric(X_val %*% beta_w)
    v_w <- y_val - yhat_w

    D_i <- (v_u^2) - (v_w^2)
    return(D_i)
  }

  if (kfold) {
    # K-fold CV splits
    folds <- sample(rep(seq_len(K), length.out = n))
    D_all <- vector("list", K)
    for (k in seq_len(K)) {
      est_idx <- which(folds != k)
      val_idx <- which(folds == k)
      D_all[[k]] <- predict_pair(est_idx, val_idx)
    }
    if (use_fold_means) {
      # One D per fold (mean of within-fold D_i)
      D_vec <- vapply(D_all, mean, numeric(1))
      method <- paste0("Pfeffermann-Nathan Predictive Power Test (", K, "-fold CV; fold means)")
    } else {
      # Pool all validation D_i across folds
      D_vec <- unlist(D_all, use.names = FALSE)
      method <- paste0("Pfeffermann-Nathan Predictive Power Test (", K, "-fold CV)")
    }
  } else {
    # Single random split
    est_n <- floor(est_split * n)
    est_idx <- sample(seq_len(n), est_n)
    val_idx <- setdiff(seq_len(n), est_idx)
    D_vec <- predict_pair(est_idx, val_idx)
    method <- "Pfeffermann-Nathan Predictive Power Test (single split)"
  }

  # Z-test on mean difference in squared errors
  Dbar <- mean(D_vec)
  sD <- stats::sd(D_vec)
  n_val <- length(D_vec)
  Z <- Dbar / (sD / sqrt(n_val))
  pval <- 2 * (1 - stats::pnorm(abs(Z)))

  structure(
    list(
      statistic = Z,
      p.value = pval,
      mean_diff = Dbar,
      n_val = n_val,
      K = if (kfold) K else NA_integer_,
      method = method,
      call = match.call()
    ),
    class = "pred_power_test"
  )
}

#' @export
print.pred_power_test <- function(x, ...) {
  cat("\n", x$method, "\n", sep = "")
  cat("Z =", formatC(x$statistic, digits = 4, format = "f"),
      " p-value =", formatC(x$p.value, digits = 4, format = "f"),
      " mean diff =", formatC(x$mean_diff, digits = 4, format = "f"),
      " (n_val =", x$n_val, ")\n")
  invisible(x)
}

#' @export
summary.pred_power_test <- function(object, ...) {
  cat("\nPfeffermann-Nathan Predictive Power Test\n")
  cat("Call:\n"); print(object$call)
  cat("\nMethod:\n ", object$method, "\n", sep = "")
  if (!is.na(object$K)) cat("K-fold cross-validation with K =", object$K, "\n")
  cat("\nObservations used in Z:", object$n_val, "\n")
  cat("\nZ statistic:", formatC(object$statistic, digits = 6, format = "f"),
      " p-value:", formatC(object$p.value, digits = 6, format = "f"),
      " mean diff:", formatC(object$mean_diff, digits = 6, format = "f"), "\n")
  invisible(object)
}

#' @export
tidy.pred_power_test <- function(x, ...) {
  tibble::tibble(
    term      = "predictive_power",
    statistic = x$statistic,
    p.value   = x$p.value,
    mean_diff = x$mean_diff,
    n_val     = x$n_val,
    K         = x$K,
    method    = x$method
  )
}

#' @export
glance.pred_power_test <- function(x, ...) {
  tibble::tibble(
    statistic = x$statistic,
    p.value   = x$p.value,
    mean_diff = x$mean_diff,
    n_val     = x$n_val,
    K         = x$K,
    method    = x$method
  )
}
