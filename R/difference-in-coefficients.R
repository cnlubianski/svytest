#' Difference-in-Coefficients Test for Survey Weights
#'
#' Implements the Hausman-Pfeffermann Difference-in-Coefficients test
#' to assess whether survey weights significantly affect regression estimates.
#'
#' @param model An object of class \code{svyglm}.
#' @param lower.tail Logical; passed to \code{pchisq()}.
#' @param var_equal Logical; assume equal residual variance between models.
#'   If \code{FALSE}, a heteroskedasticity-robust variance estimator is used.
#' @param robust_type Character; type of heteroskedasticity-robust variance
#'   estimator to use if \code{var_equal = FALSE}. Options are
#'   \code{"HC0"}, \code{"HC1"}, \code{"HC2"}, \code{"HC3"} as used in `sandwich` package.
#' @param coef_subset Character vector of coefficient names to include in the test.
#'   Defaults to all coefficients.
#' @param na.action Function to handle missing data before fitting the test.
#'
#' @return An object of class \code{"diff_in_coef_test"} containing:
#'   \item{statistic}{Chi-squared test statistic}
#'   \item{parameter}{Degrees of freedom}
#'   \item{p.value}{P-value for the test}
#'   \item{betas_unweighted}{Unweighted coefficient estimates}
#'   \item{betas_weighted}{Weighted coefficient estimates}
#'   \item{vcov_diff}{Estimated variance-covariance matrix of coefficient differences}
#'   \item{diff_betas}{Vector of coefficient differences}
#'   \item{call}{Function call}
#'
#' @examples
#' if (requireNamespace("survey", quietly = TRUE) &&
#'     exists("svytestCE", where = "package:svytest")) {
#'   library(survey)
#'   data("svytestCE", package = "svytest")
#'   des <- svydesign(ids = ~1, weights = ~FINLWT21, data = svytestCE)
#'   fit_svy <- svyglm(TOTEXPCQ ~ ROOMSQ + BATHRMQ + BEDROOMQ + FAM_SIZE + AGE,
#'                     design = des)
#'   results <- diff_in_coef_test(fit_svy, var_equal = FALSE, robust_type = "HC3")
#'   print(results)
#' }
#'
#' @details
#' This test adapts the general Hausman specification test
#' (Hausman, 1978) to the survey-weighted regression setting.
#' Pfeffermann (1993) showed how comparing weighted and unweighted
#' coefficient estimates provides a diagnostic for whether survey
#' weights materially affect model results. The procedure is often
#' called the Difference-in-Coefficients test.
#'
#' @references
#' Hausman, J. A. (1978). Specification Tests in Econometrics.
#'   *Econometrica*, 46(6), 1251-1271. \doi{10.2307/1913827}
#'
#' Pfeffermann, D. (1993). The Role of Sampling Weights When Modeling Survey Data.
#'   *International Statistical Review*, 61(2), 317-337. \doi{10.2307/1403631}
#'
#' @seealso
#' \code{\link{svytestCE}} for the curated Consumer Expenditure dataset
#'   included in this package, which can be used to demonstrate the
#'   Difference-in-Coefficients test.
#'
#' @export
diff_in_coef_test <- function(model, lower.tail = FALSE, var_equal = TRUE,
                              robust_type = c("HC0", "HC1", "HC2", "HC3"),
                              coef_subset = NULL, na.action = stats::na.omit) {
  if (!inherits(model, "svyglm")) stop("Model must be of class 'svyglm'.")
  robust_type <- match.arg(robust_type)

  # Extract design matrix, response, and weights
  wts <- weights(model$survey.design)
  X <- stats::model.matrix(model)
  y <- stats::model.response(stats::model.frame(model))

  # Handle missing data
  dat <- data.frame(y = y, X, wts = wts)
  dat <- na.action(dat)
  y <- dat$y
  X <- as.matrix(dat[, setdiff(names(dat), c("y", "wts"))])
  wts <- dat$wts

  # Warn if weights are constant (tolerance ~1e-8)
  if (all(abs(wts - mean(wts)) < .Machine$double.eps^0.5)) {
    warning("All weights are constant; the test is not meaningful.")
  }

  # Optionally subset coefficients
  if (!is.null(coef_subset)) {
    keep_cols <- colnames(X) %in% coef_subset
    if (!any(keep_cols)) stop("No matching coefficients found in model.")
    X <- X[, keep_cols, drop = FALSE]
  }

  # Unweighted coefficients
  betas_u <- solve(t(X) %*% X) %*% t(X) %*% y

  # Weighted coefficients
  W <- diag(wts, nrow = length(wts))
  betas_w <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y

  # A matrix (p x n)
  A <- (solve(t(X) %*% W %*% X) %*% t(X) %*% W) -
    (solve(t(X) %*% X) %*% t(X))

  # Variance estimation
  if (var_equal) {
    # Equal variance assumption
    y_hat <- X %*% betas_u
    residuals <- y - y_hat
    SSE <- sum(residuals^2)
    sigma_sq_hat <- SSE / (length(y) - length(betas_u))
    V_hat <- sigma_sq_hat * A %*% t(A)
  } else {
    # Heteroskedasticity-robust variance
    resid_w <- y - X %*% betas_w
    resid_u <- y - X %*% betas_u
    diff_resid <- as.numeric(resid_w - resid_u)

    # Leverage values for HC2/HC3
    H <- X %*% solve(t(X) %*% X) %*% t(X)
    h_ii <- diag(H)

    # Adjust residuals depending on robust_type
    if (robust_type == "HC0") {
      adj_resid <- diff_resid
    } else if (robust_type == "HC1") {
      adj_resid <- diff_resid * sqrt(length(y) / (length(y) - ncol(X)))
    } else if (robust_type == "HC2") {
      adj_resid <- diff_resid / sqrt(1 - h_ii)
    } else if (robust_type == "HC3") {
      adj_resid <- diff_resid / (1 - h_ii)
    }

    # Middle matrix of sandwich
    V_hat <- A %*% (t(A) * adj_resid^2)
  }

  # Test statistic
  diff_betas <- betas_w - betas_u
  Chi_statistic <- as.numeric(t(diff_betas) %*% solve(V_hat) %*% diff_betas)
  df <- length(betas_u)
  p_value <- stats::pchisq(Chi_statistic, df = df, lower.tail = lower.tail)

  structure(
    list(
      statistic = Chi_statistic,
      parameter = df,
      p.value = p_value,
      method = "Hausman-Pfeffermann Difference-in-Coefficients Test",
      betas_unweighted = as.vector(betas_u),
      betas_weighted = as.vector(betas_w),
      vcov_diff = V_hat,
      diff_betas = as.vector(diff_betas),
      call = match.call()
    ),
    class = "diff_in_coef_test"
  )
}

#' @export
print.diff_in_coef_test <- function(x, ...) {
  cat("\n", x$method, "\n", sep = "")
  cat("Chi-squared =", formatC(x$statistic, digits = 4, format = "f"),
      " df =", x$parameter,
      " p-value =", formatC(x$p.value, digits = 4, format = "f"), "\n")
  cat("\nUnweighted coefficients:\n")
  print(x$betas_unweighted)
  cat("\nWeighted coefficients:\n")
  print(x$betas_weighted)
  invisible(x)
}

#' @export
summary.diff_in_coef_test <- function(object, ...) {
  cat("\nDifference-in-Coefficients Test\n")
  cat("Call:\n")
  print(object$call)
  cat("\nMethod:\n ", object$method, "\n", sep = "")
  cat("\nTest Statistic:\n")
  cat(" Chi-sq =", formatC(object$statistic, digits = 4, format = "f"),
      " on", object$parameter, "df",
      ", p-value =", formatC(object$p.value, digits = 4, format = "f"), "\n")
  invisible(object)
}

#' @export
tidy.diff_in_coef_test <- function(x, ...) {
  terms <- names(x$betas_unweighted)
  if (is.null(terms)) {
    terms <- paste0("V", seq_along(x$betas_unweighted))
  }
  tibble::tibble(
    term       = terms,
    unweighted = x$betas_unweighted,
    weighted   = x$betas_weighted,
    diff       = x$diff_betas
  )
}

#' @export
glance.diff_in_coef_test <- function(x, ...) {
  tibble::tibble(
    statistic = x$statistic,
    df        = x$parameter,
    p.value   = x$p.value,
    method    = x$method
  )
}
