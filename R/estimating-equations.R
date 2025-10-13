#' Estimating Equations Test for Informative Sampling (Linear Case)
#'
#' Implements the Pfeffermann-Sverchkov estimating equations test for
#' informativeness of survey weights in the linear regression case
#' (Gaussian with identity link). The test compares unweighted estimating
#' equations to adjusted-weight equations using \eqn{q_i = w_i / E_s(w_i \mid x_i)}.
#'
#' @param model An object of class \code{svyglm} with \code{family = gaussian(identity)}.
#' @param coef_subset Optional character vector of coefficient names to include in the test.
#'   Defaults to all coefficients in the model matrix.
#' @param q_method Method for estimating \eqn{E_s(w \mid x)}:
#'   \code{"linear"} (default, OLS regression of w on X) or
#'   \code{"log"} (regress log(w) on X, then exponentiate).
#' @param stabilize Logical; if \code{TRUE} (default) clips extreme q values.
#' @param na.action Function to handle missing data.
#'
#' @return An object of class \code{"estim_eq_test"} containing:
#'   \item{statistic}{Hotelling F statistic}
#'   \item{p.value}{p-value under F distribution}
#'   \item{df1}{Numerator df (# tested coefficients)}
#'   \item{df2}{Denominator df (n - p)}
#'   \item{Rbar}{Mean estimating-equation contrast vector}
#'   \item{S}{Sample covariance of R}
#'   \item{terms}{Names of tested coefficients}
#'   \item{n}{Sample size}
#'   \item{call}{Matched call}
#'   \item{method}{Description string}
#'
#' @details
#' For linear regression, the per-observation score is
#' \eqn{u_i = x_i (y_i - x_i^\top \hat\beta_{\text{unw}})} at the unweighted OLS estimate.
#' The test statistic is based on \eqn{R_i = (1 - q_i) u_i}, where
#' \eqn{q_i = w_i / E_s(w_i \mid x_i)}. The Hotelling F statistic is
#' \eqn{F = \frac{n-p}{p} \bar R^\top S^{-1} \bar R}, with df1 = p, df2 = n - p.
#'
#' @examples
#' \dontrun{
#' library(survey)
#' set.seed(123)
#' n <- 200
#' x <- rnorm(n)
#' y <- 1 + 2*x + rnorm(n)
#' w <- 1 + (y - min(y))
#' dat <- data.frame(y, x, w)
#' des <- svydesign(ids = ~1, weights = ~w, data = dat)
#' fit <- svyglm(y ~ x, design = des, family = gaussian())
#'
#' res <- estim_eq_test(fit, q_method = "linear")
#' print(res)
#' summary(res)
#' }
#'
#' @references
#' Pfeffermann, D., & Sverchkov, M. Y. (2003).
#' Fitting generalized linear models under informative sampling.
#' In R. L. Chambers & C. J. Skinner (Eds.), *Analysis of Survey Data* (Ch. 12). Wiley.
#'
#' @seealso \code{\link{diff_in_coef_test}}, \code{\link{lr_test}},
#'    \code{\link{wa_test}}, \code{\link{perm_test}}
#'
#' @export
estim_eq_test <- function(model, coef_subset = NULL, q_method = c("linear", "log"),
                          stabilize = TRUE, na.action = stats::na.omit) {
  if (!inherits(model, "svyglm")) stop("Model must be of class 'svyglm'.")
  fam <- model$family
  if (!(fam$family == "gaussian" && fam$link == "identity")) {
    stop("estim_eq_test currently supports only gaussian(identity) models.")
  }
  q_method <- match.arg(q_method)

  # Extract data
  w <- stats::weights(model$survey.design)
  X_full <- stats::model.matrix(model)
  y <- stats::model.response(stats::model.frame(model))

  # Handle missing
  dat <- data.frame(y = y, X_full, w = w)
  dat <- na.action(dat)
  y <- dat$y
  X_full <- as.matrix(dat[, setdiff(names(dat), c("y", "w"))])
  w <- dat$w

  # Subset coefficients
  X <- X_full
  terms <- colnames(X)
  if (!is.null(coef_subset)) {
    keep <- colnames(X) %in% coef_subset
    if (!any(keep)) stop("No matching coefficients found in model matrix.")
    X <- X[, keep, drop = FALSE]
    terms <- colnames(X)
  }

  n <- nrow(X)
  p <- ncol(X)

  # Unweighted OLS
  beta_hat <- solve(crossprod(X), crossprod(X, y))
  resid <- y - X %*% beta_hat
  U <- sweep(X, 1, resid, `*`)  # scores

  # Estimate E_s(w|x)
  W_df <- data.frame(w = w, X)
  if (q_method == "linear") {
    fit_w <- stats::lm(w ~ ., data = W_df)
    w_hat <- as.numeric(stats::predict(fit_w, newdata = W_df))
    w_hat <- pmax(w_hat, .Machine$double.eps)
  } else {
    fit_w <- stats::lm(log(w) ~ ., data = W_df)
    w_hat <- exp(as.numeric(stats::predict(fit_w, newdata = W_df)))
  }

  q <- w / w_hat
  if (stabilize) {
    lo <- stats::quantile(q, probs = 1e-3, na.rm = TRUE)
    hi <- stats::quantile(q, probs = 1 - 1e-3, na.rm = TRUE)
    q <- pmax(pmin(q, hi), lo)
  }

  # R_i = (1 - q_i) u_i
  R <- sweep(U, 1, (1 - q), `*`)
  Rbar <- colMeans(R)
  S <- stats::cov(R)

  # Hotelling F
  Fstat <- as.numeric(((n - p) / p) * (t(Rbar) %*% solve(S, Rbar)))
  pval <- stats::pf(Fstat, df1 = p, df2 = n - p, lower.tail = FALSE)

  structure(
    list(
      statistic = Fstat,
      p.value = pval,
      df1 = p,
      df2 = n - p,
      Rbar = as.numeric(Rbar),
      S = S,
      terms = terms,
      n = n,
      method = "Pfeffermann-Sverchkov Estimating Equations Test (linear case)",
      call = match.call()
    ),
    class = "estim_eq_test"
  )
}

#' @rdname estim_eq_test
#' @method print estim_eq_test
#' @param x An object of class estim_eq_test
#' @param ... Additional arguments passed to methods
#' @export
print.estim_eq_test <- function(x, ...) {
  cat("\n", x$method, "\n", sep = "")
  cat("F =", formatC(x$statistic, digits = 4, format = "f"),
      " df1 =", x$df1, " df2 =", x$df2,
      " p-value =", formatC(x$p.value, digits = 4, format = "f"), "\n")
  invisible(x)
}

#' @rdname estim_eq_test
#' @method summary estim_eq_test
#' @param object An object of class estim_eq_test
#' @param ... Additional arguments passed to methods
#' @export
summary.estim_eq_test <- function(object, ...) {
  cat("\nEstimating Equations Test (linear case)\n")
  cat("Call:\n")
  print(object$call)
  cat("\nMethod:\n ", object$method, "\n", sep = "")
  cat("\nTested terms:\n ", paste(object$terms, collapse = ", "), "\n", sep = "")
  cat("\nHotelling F:\n ",
      formatC(object$statistic, digits = 6, format = "f"),
      " (df1 =", object$df1, ", df2 =", object$df2, ")",
      ", p-value =", formatC(object$p.value, digits = 6, format = "f"), "\n", sep = "")
  invisible(object)
}

#' @rdname estim_eq_test
#' @method tidy estim_eq_test
#' @param x An object of class estim_eq_test
#' @param ... Additional arguments passed to methods
tidy.estim_eq_test <- function(x, ...) {
  tibble::tibble(
    term      = x$terms,
    Rbar      = x$Rbar,
    df1       = x$df1,
    df2       = x$df2,
    statistic = x$statistic,
    p.value   = x$p.value,
    method    = x$method,
    n         = x$n
  )
}

#' @rdname estim_eq_test
#' @method glance estim_eq_test
#' @param x An object of class estim_eq_test
#' @param ... Additional arguments passed to methods
glance.estim_eq_test <- function(x, ...) {
  tibble::tibble(
    statistic = x$statistic,
    p.value   = x$p.value,
    df1       = x$df1,
    df2       = x$df2,
    n         = x$n,
    method    = x$method
  )
}
