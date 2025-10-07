#' Pfeffermann–Nathan Predictive Power Test (with optional K-fold CV)
#'
#' Implements the Pfeffermann–Nathan predictive power test for informativeness
#' of survey weights in linear regression. By default, uses K-fold cross-validation
#' to compare predictive performance of weighted vs. unweighted regressions.
#'
#' @param data A data frame containing outcome, covariates, and weights.
#' @param y Vector of outcome values.
#' @param x Matrix or vector of covariates (will be coerced to matrix).
#' @param wts Vector of survey weights.
#' @param est_split Proportion of the sample to use for estimation (only used if \code{kfold = FALSE}).
#' @param kfold Logical; if TRUE (default), use K-fold cross-validation instead of a single split.
#' @param K Number of folds for cross-validation (default 5).
#' @param seed Optional random seed for reproducibility.
#'
#' @return An object of class \code{"pred_power_test"} containing:
#'   \item{statistic}{Z statistic}
#'   \item{p.value}{Two-sided p-value}
#'   \item{mean_diff}{Mean difference in squared prediction errors}
#'   \item{n_val}{Total validation sample size across folds}
#'   \item{K}{Number of folds (if kfold = TRUE)}
#'   \item{call}{Matched call}
#'   \item{method}{Description string}
#'
#' @export
pred_power_test <- function(data, y, x, wts,
                            est_split = 0.5,
                            kfold = TRUE,
                            K = 5,
                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(data)
  if (length(y) != n || length(wts) != n) stop("y and wts must match nrow(data).")
  X <- as.matrix(cbind(1, x))  # add intercept

  all_D <- c()

  if (kfold) {
    # K-fold CV
    folds <- sample(rep(1:K, length.out = n))
    for (k in 1:K) {
      est_index <- which(folds != k)
      val_index <- which(folds == k)

      X_est <- X[est_index, , drop = FALSE]
      y_est <- y[est_index]
      w_est <- wts[est_index]

      X_val <- X[val_index, , drop = FALSE]
      y_val <- y[val_index]

      # Unweighted OLS
      beta_u <- solve(crossprod(X_est), crossprod(X_est, y_est))
      yhat_u <- X_val %*% beta_u
      v_u <- y_val - yhat_u

      # Weighted OLS
      W <- diag(w_est)
      beta_w <- solve(t(X_est) %*% W %*% X_est, t(X_est) %*% W %*% y_est)
      yhat_w <- X_val %*% beta_w
      v_w <- y_val - yhat_w

      D <- as.numeric(v_u^2 - v_w^2)
      all_D <- c(all_D, D)
    }
    method <- paste0("Pfeffermann–Nathan Predictive Power Test (", K, "-fold CV)")
  } else {
    # Single split
    est_n <- floor(est_split * n)
    est_index <- sample(seq_len(n), est_n)
    val_index <- setdiff(seq_len(n), est_index)

    X_est <- X[est_index, , drop = FALSE]
    y_est <- y[est_index]
    w_est <- wts[est_index]

    X_val <- X[val_index, , drop = FALSE]
    y_val <- y[val_index]

    # Unweighted OLS
    beta_u <- solve(crossprod(X_est), crossprod(X_est, y_est))
    yhat_u <- X_val %*% beta_u
    v_u <- y_val - yhat_u

    # Weighted OLS
    W <- diag(w_est)
    beta_w <- solve(t(X_est) %*% W %*% X_est, t(X_est) %*% W %*% y_est)
    yhat_w <- X_val %*% beta_w
    v_w <- y_val - yhat_w

    all_D <- as.numeric(v_u^2 - v_w^2)
    method <- "Pfeffermann–Nathan Predictive Power Test (single split)"
  }

  Dbar <- mean(all_D)
  sD <- sd(all_D)
  n_val <- length(all_D)

  Z <- Dbar / (sD / sqrt(n_val))
  pval <- 2 * (1 - pnorm(abs(Z)))

  structure(
    list(
      statistic = Z,
      p.value = pval,
      mean_diff = Dbar,
      n_val = n_val,
      est_split = if (kfold) NA else est_split,
      K = if (kfold) K else NA,
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
  cat("\nPfeffermann–Nathan Predictive Power Test\n")
  cat("Call:\n")
  print(object$call)
  cat("\nMethod:\n ", object$method, "\n", sep = "")
  if (!is.na(object$K)) cat("K-fold cross-validation with K =", object$K, "\n")
  if (!is.na(object$est_split)) cat("Estimation split proportion:", object$est_split, "\n")
  cat("\nValidation sample size:", object$n_val, "\n")
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
    est_split = x$est_split,
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
    est_split = x$est_split,
    K         = x$K,
    method    = x$method
  )
}
