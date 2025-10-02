#' Likelihood-Ratio Test for Informative Survey Weights
#'
#' Implements the Breidt–Herndon likelihood-ratio test for assessing whether
#' survey weights are informative in linear regression models. The test compares
#' maximized log-likelihoods under equal weights (null) and survey weights
#' (alternative), with an asymptotic distribution given by a weighted chi-squared
#' mixture.
#'
#' @param model An object of class \code{svyglm}.
#' @param coef_subset Optional character vector of coefficient names to include
#'   in the test. Defaults to all coefficients.
#' @param na.action Function to handle missing data before testing.
#' @param likelihood Character string specifying the likelihood form:
#'   \code{"pseudo"} (default) for raw weighted likelihood,
#'   or \code{"scaled"} to normalize weights by their mean.
#'
#' @return An object of class \code{"lr_test"} containing:
#'   \item{statistic}{Likelihood-ratio test statistic (non-negative)}
#'   \item{p.value}{P-value for the test (Satterthwaite approximation)}
#'   \item{df}{Approximate degrees of freedom}
#'   \item{eigvals}{Eigenvalues of the Gamma matrix}
#'   \item{logLik_null}{Maximized log-likelihood under equal weights}
#'   \item{logLik_alt}{Maximized log-likelihood under survey weights}
#'   \item{method}{Name of the test performed}
#'   \item{call}{Function call}
#'
#' @details
#' The null hypothesis is that survey weights are not informative (equal weights
#' suffice). The alternative allows weights to affect the likelihood. The
#' asymptotic null distribution is a weighted chi-squared mixture; here we
#' approximate the p-value using a Satterthwaite moment-matching approach.
#'
#' @references
#' Breidt, F. J., & Opsomer, J. D. (1997).
#'   Testing for informativeness in analytic inference from complex surveys.
#'   *Survey Methodology*, 23(1), 1–11.
#'
#' Herndon, J. (2022).
#'   Testing and adjusting for informative sampling in survey data.
#'   *Journal of Survey Statistics and Methodology*, 10(3), 455–480.
#'
#' @seealso
#' \code{\link{diff_in_coef}}, \code{\link{wa_test}}, \code{\link{svytestCE}}
#'
#' @export
lr_test <- function(model,
                    coef_subset = NULL,
                    na.action = na.omit,
                    likelihood = c("pseudo","scaled")) {
  if (!inherits(model, "svyglm")) stop("Model must be of class 'svyglm'.")
  likelihood <- match.arg(likelihood)

  # Extract design matrix, response, weights
  wts <- weights(model$survey.design)
  X <- model.matrix(model)
  y <- model.response(model.frame(model))

  # Handle missing data
  dat <- data.frame(y = y, X, wts = wts)
  dat <- na.action(dat)
  y <- dat$y
  X <- as.matrix(dat[, setdiff(names(dat), c("y","wts"))])
  wts <- dat$wts

  # Optionally subset coefficients
  if (!is.null(coef_subset)) {
    keep_cols <- colnames(X) %in% coef_subset
    if (!any(keep_cols)) stop("No matching coefficients found in model.")
    X <- X[, keep_cols, drop = FALSE]
  }

  # Log-likelihood
  logLike <- function(params, y, X, wts, likelihood) {
    beta <- params[1:ncol(X)]
    sigma2 <- exp(params[ncol(X)+1])
    mu <- X %*% beta
    if (likelihood == "pseudo") {
      ll <- -0.5 * sum(wts * (log(2*pi*sigma2) + (y - mu)^2 / sigma2))
    } else { # scaled
      wts_scaled <- wts / mean(wts)
      ll <- -0.5 * sum(wts_scaled * (log(2*pi*sigma2) + (y - mu)^2 / sigma2))
    }
    return(-ll) # optim minimizes
  }

  init <- c(rep(0, ncol(X)), log(var(y)))

  # Null: equal weights
  wts_null <- rep(mean(wts), length(wts))
  fit_null <- optim(init, logLike, y=y, X=X, wts=wts_null,
                    likelihood=likelihood, method="BFGS", hessian=TRUE)
  ll_null <- -fit_null$value

  # Alt: actual weights
  fit_alt <- optim(init, logLike, y=y, X=X, wts=wts,
                   likelihood=likelihood, method="BFGS", hessian=TRUE)
  ll_alt <- -fit_alt$value

  # LR statistic
  LR <- 2 * (ll_alt - ll_null)
  if (LR < 0) {
    LR <- 0
    pval <- 1
    eigvals <- NA
    df <- NA
  } else {
    # Approximate null distribution via Satterthwaite
    beta_hat <- fit_alt$par[1:ncol(X)]
    sigma2_hat <- exp(fit_alt$par[ncol(X)+1])
    mu_hat <- X %*% beta_hat

    score_beta <- sweep(X, 1, (wts * (y - mu_hat) / sigma2_hat), `*`)
    score_sigma <- 0.5 * wts * ((y - mu_hat)^2 / sigma2_hat - 1)
    score_i <- cbind(score_beta, score_sigma)

    J_hat <- -fit_alt$hessian / length(y)
    K_hat <- crossprod(scale(score_i, scale=FALSE)) / length(y)
    Gamma <- solve(J_hat) %*% K_hat %*% solve(J_hat)
    eigvals <- eigen(Gamma, only.values=TRUE)$values

    df <- 2 * (sum(eigvals)^2) / sum(eigvals^2)
    scale <- sum(eigvals) / df
    pval <- 1 - pchisq(LR/scale, df=df)
  }

  structure(
    list(
      statistic = LR,
      p.value   = pval,
      df        = df,
      eigvals   = eigvals,
      logLik_null = ll_null,
      logLik_alt  = ll_alt,
      method    = paste("Breidt–Herndon Likelihood-Ratio Test (", likelihood, " likelihood)", sep=""),
      call      = match.call()
    ),
    class = "lr_test"
  )
}

#' @export
print.lr_test <- function(x, ...) {
  cat("\n", x$method, "\n", sep = "")
  cat("LR =", formatC(x$statistic, digits=4, format="f"),
      if (!is.na(x$df)) paste(" df ≈", formatC(x$df, digits=2, format="f")) else "",
      " p-value =", formatC(x$p.value, digits=4, format="f"), "\n")
  invisible(x)
}

#' @export
summary.lr_test <- function(object, ...) {
  cat("\nLikelihood-Ratio Test for Informative Weights\n")
  cat("Call:\n")
  print(object$call)
  cat("\nMethod:\n ", object$method, "\n", sep = "")
  cat("\nTest Statistic:\n")
  cat(" LR =", formatC(object$statistic, digits = 4, format = "f"))
  if (!is.na(object$df)) {
    cat(" on ~", formatC(object$df, digits = 2, format = "f"), "df (approx)")
  }
  cat(", p-value =", formatC(object$p.value, digits = 4, format = "f"), "\n")
  cat("\nLog-likelihoods:\n")
  cat(" Null =", formatC(object$logLik_null, digits = 4, format = "f"),
      " Alt =", formatC(object$logLik_alt, digits = 4, format = "f"), "\n")
  invisible(object)
}

#' @export
tidy.lr_test <- function(x, ...) {
  tibble::tibble(
    term     = "likelihood_ratio",
    statistic = x$statistic,
    df        = x$df,
    p.value   = x$p.value,
    method    = x$method
  )
}

#' @export
glance.lr_test <- function(x, ...) {
  tibble::tibble(
    statistic = x$statistic,
    df        = x$df,
    p.value   = x$p.value,
    logLik_null = x$logLik_null,
    logLik_alt  = x$logLik_alt,
    method    = x$method
  )
}
