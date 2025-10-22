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
#' # Load in survey package (required) and load in example data
#' library(survey)
#' data(api, package = "survey")
#'
#' # Create a survey design and fit a weighted regression model
#' des <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat)
#' fit <- svyglm(api00 ~ ell + meals, design = des)
#'
#' # Run difference-in-coefficients diagnostic test versions with different variance assumptions
#' # and reports Chi-Squared statistic, df, and p-value
#' summary(diff_in_coef_test(fit, var_equal = TRUE))
#' summary(diff_in_coef_test(fit, var_equal = FALSE, robust_type = "HC3"))
#'
#' @details
#' Let \eqn{X} denote the design matrix and \eqn{y} the response vector.
#' Define the unweighted OLS estimator
#' \deqn{\hat\beta_{U} = (X^\top X)^{-1} X^\top y,}
#' and the survey-weighted estimator
#' \deqn{\hat\beta_{W} = (X^\top W X)^{-1} X^\top W y,}
#' where \eqn{W = \mathrm{diag}(w_1, \ldots, w_n)} is the diagonal matrix of survey weights.
#'
#' The test statistic is based on the difference
#' \deqn{d = \hat\beta_{W} - \hat\beta_{U}.}
#'
#' Under the null hypothesis that weights are not informative,
#' \eqn{d} has mean zero and variance \eqn{V_d}.
#' The test statistic is
#' \deqn{T = d^\top V_d^{-1} d,}
#' which is asymptotically \eqn{\chi^2_p} distributed with
#' \eqn{p} equal to the number of coefficients tested.
#'
#' If \code{var_equal = TRUE}, \eqn{V_d} is estimated assuming equal residual variance
#' across weighted and unweighted models. If \code{var_equal = FALSE}, a
#' heteroskedasticity-robust estimator (e.g. HC0–HC3) is used.
#'
#' This test is a survey-weighted adaptation of the Hausman specification test
#' (Hausman, 1978), as proposed by Pfeffermann (1993).
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
#'   @importFrom survey svyglm
#'
#' @export
diff_in_coef_test <- function(model, lower.tail = FALSE, var_equal = TRUE,
                              robust_type = c("HC0", "HC1", "HC2", "HC3"),
                              coef_subset = NULL, na.action = stats::na.omit) {

  # Argument checks
  if (!inherits(model, "svyglm")) {
    stop("`model` must be an object of class 'svyglm' (from the survey package).")
  }

  if (!is.logical(lower.tail) || length(lower.tail) != 1L || is.na(lower.tail)) {
    stop("`lower.tail` must be a single logical value (TRUE/FALSE).")
  }

  if (!is.logical(var_equal) || length(var_equal) != 1L || is.na(var_equal)) {
    stop("`var_equal` must be a single logical value (TRUE/FALSE).")
  }

  valid_types <- c("HC0", "HC1", "HC2", "HC3")
  if (!var_equal && !robust_type %in% valid_types) {
    stop("`robust_type` must be one of: ", paste(valid_types, collapse = ", "), ".")
  }

  # Argument check for coef_subset type
  if (!is.null(coef_subset) && !is.character(coef_subset)) {
    stop("`coef_subset` must be a character vector of coefficient names.")
  }

  if (!is.function(na.action)) {
    stop("`na.action` must be a function (e.g., stats::na.omit).")
  }

  # Extract design matrix, response, and weights
  wts <- tryCatch(stats::weights(model$survey.design),
                  error = function(e) stop("Could not extract weights from `model$survey.design`."))

  X <- tryCatch(stats::model.matrix(model),
                error = function(e) stop("Could not extract model matrix from `model`."))

  y <- tryCatch(stats::model.response(stats::model.frame(model)),
                error = function(e) stop("Could not extract response from `model`."))

  if (length(y) != nrow(X) || length(y) != length(wts)) {
    stop("Mismatch in dimensions of response, design matrix, and weights.")
  }

  # Handle missing and extract data into objects
  dat <- data.frame(y = y, X, wts = wts)
  dat <- na.action(dat)
  y <- dat$y
  X <- as.matrix(dat[, setdiff(names(dat), c("y", "wts"))])
  wts <- dat$wts

  # Warn if weights are constant
  if (all(abs(wts - mean(wts)) < .Machine$double.eps^0.5)) {
    warning("All weights are constant; the test is not meaningful.")
  }

  # Optionally subset coefficients
  if (!is.null(coef_subset)) {
    keep_cols <- colnames(X) %in% coef_subset
    if (!any(keep_cols)) {
      stop("None of the names in `coef_subset` matched the model coefficients: ",
           paste(colnames(X), collapse = ", "))
    }
    X <- X[, keep_cols, drop = FALSE]
  }

  # Unweighted coefficients
  betas_u <- tryCatch(
    solve(t(X) %*% X, t(X) %*% y),
    error = function(e) stop("Unweighted OLS failed: design matrix may be singular.")
  )

  # Weighted coefficients
  W <- diag(wts, nrow = length(wts))
  betas_w <- tryCatch(
    solve(t(X) %*% W %*% X, t(X) %*% W %*% y),
    error = function(e) stop("Weighted regression failed: weighted design matrix may be singular.")
  )

  # A matrix (p x n)
  A <- (solve(t(X) %*% W %*% X) %*% t(X) %*% W) -
    (solve(t(X) %*% X) %*% t(X))

  # Variance estimation
  if (var_equal) {
    y_hat <- X %*% betas_u
    residuals <- y - y_hat
    SSE <- sum(residuals^2)
    sigma_sq_hat <- SSE / (length(y) - length(betas_u))
    V_hat <- sigma_sq_hat * A %*% t(A)
  } else {
    resid_w <- y - X %*% betas_w
    resid_u <- y - X %*% betas_u
    diff_resid <- as.numeric(resid_w - resid_u)

    H <- X %*% solve(t(X) %*% X) %*% t(X)
    h_ii <- diag(H)

    if (robust_type == "HC0") {
      adj_resid <- diff_resid
    } else if (robust_type == "HC1") {
      adj_resid <- diff_resid * sqrt(length(y) / (length(y) - ncol(X)))
    } else if (robust_type == "HC2") {
      adj_resid <- diff_resid / sqrt(1 - h_ii)
    } else if (robust_type == "HC3") {
      adj_resid <- diff_resid / (1 - h_ii)
    }

    V_hat <- A %*% (t(A) * adj_resid^2)
  }

  # Compute test statistic
  diff_betas <- betas_w - betas_u
  Chi_statistic <- as.numeric(t(diff_betas) %*% solve(V_hat, diff_betas))
  df <- length(betas_u)
  p_value <- stats::pchisq(Chi_statistic, df = df, lower.tail = lower.tail)

  # Return structured result
  structure(
    list(
      statistic = Chi_statistic,
      parameter = df,
      p.value = p_value,
      method = "Hausman-Pfeffermann Difference-in-Coefficients Test",
      betas_unweighted = stats::setNames(as.vector(betas_u), names(stats::coef(model))),
      betas_weighted = stats::setNames(as.vector(betas_w), names(stats::coef(model))),
      vcov_diff = V_hat,
      diff_betas = as.vector(diff_betas),
      call = match.call()
    ),
    class = "diff_in_coef_test"
  )
}



#' @rdname diff_in_coef_test
#' @method print diff_in_coef_test
#' @param x An object of class diff_in_coef_test
#' @param ... Additional arguments passed to methods
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

#' @rdname diff_in_coef_test
#' @method summary diff_in_coef_test
#' @param object An object of class diff_in_coef_test
#' @param ... Additional arguments passed to methods
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

#' @rdname diff_in_coef_test
#' @method tidy diff_in_coef_test
#' @param x An object of class diff_in_coef_test
#' @param ... Additional arguments passed to methods
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

#' @rdname diff_in_coef_test
#' @method glance diff_in_coef_test
#' @param x An object of class diff_in_coef_test
#' @param ... Additional arguments passed to methods
glance.diff_in_coef_test <- function(x, ...) {
  tibble::tibble(
    statistic = x$statistic,
    df        = x$parameter,
    p.value   = x$p.value,
    method    = x$method
  )
}
