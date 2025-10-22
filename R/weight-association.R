#' Weight-Association Tests for Survey Weights
#'
#' Implements several weight-association tests that examine whether survey
#' weights are informative about the response variable after conditioning
#' on covariates. Variants include DuMouchel-Duncan (DD), Pfeffermann-Sverchkov
#' (PS1 and PS2, with optional quadratic terms or user-supplied auxiliary designs),
#' and Wu-Fuller (WF).
#'
#' @param model An object of class \code{svyglm}.
#' @param type Character string specifying the test type:
#'   \code{"DD"}, \code{"PS1"}, \code{"PS1q"}, \code{"PS2"}, \code{"PS2q"}, \code{"WF"}.
#' @param coef_subset Optional character vector of coefficient names to include
#'   in the test. Defaults to all coefficients.
#' @param aux_design Optional matrix or function to generate auxiliary regressors
#'   for PS1/PS2 tests. If a function, it should take \code{X} and \code{y}
#'   and return a matrix of extra columns to include.
#' @param na.action Function to handle missing data before testing.
#'
#' @return An object of class \code{"wa_test"} containing:
#'   \item{statistic}{F-test statistic}
#'   \item{parameter}{Degrees of freedom (numerator, denominator)}
#'   \item{p.value}{P-value for the test}
#'   \item{method}{Name of the test performed}
#'   \item{call}{Function call}
#'
#' @details
#' Let \eqn{y} denote the response, \eqn{X} the design matrix of covariates,
#' and \eqn{w} the survey weights. The null hypothesis in all cases is that
#' the weights are \emph{non-informative} given \eqn{X}, i.e. they do not
#' provide additional information about \eqn{y} beyond the covariates.
#'
#' The following test variants are implemented:
#'
#' \itemize{
#'
#'   \item \strong{DuMouchel–Duncan (DD)}:
#'     After fitting the unweighted regression
#'     \deqn{\hat\beta = (X^\top X)^{-1} X^\top y,}
#'     compute residuals \eqn{e = y - X\hat\beta}.
#'     The DD test regresses \eqn{e} on the weights \eqn{w}:
#'     \deqn{e = \gamma_0 + \gamma_1 w + u.}
#'     A significant \eqn{\gamma_1} indicates association between weights
#'     and residuals, hence informativeness.
#'
#'   \item \strong{Pfeffermann–Sverchkov PS1}:
#'     Augments the outcome regression with functions of the weights as
#'     auxiliary regressors:
#'     \deqn{y = X\beta + f(w)\theta + \varepsilon.}
#'     Under the null, \eqn{\theta = 0}. Quadratic terms
#'     (\eqn{w^2}) can be included (\code{"PS1q"}), or the user may supply
#'     a custom auxiliary design matrix \eqn{f(w)}.
#'
#'   \item \strong{Pfeffermann–Sverchkov PS2}:
#'     First regress the weights on the covariates,
#'     \deqn{w = X\alpha + \eta,}
#'     and obtain fitted values \eqn{\hat w}.
#'     Then augment the outcome regression with \eqn{\hat w} (and optionally
#'     \eqn{\hat w^2} for \code{"PS2q"}):
#'     \deqn{y = X\beta + g(\hat w)\theta + \varepsilon.}
#'     Again, \eqn{\theta = 0} under the null.
#'
#'   \item \strong{Wu–Fuller (WF)}:
#'     Compares weighted and unweighted regression fits. Let
#'     \eqn{\hat\beta_W} and \eqn{\hat\beta_U} denote the weighted and
#'     unweighted estimators. The test statistic is based on
#'     \deqn{T = (\hat\beta_W - \hat\beta_U)^\top
#'                \widehat{\mathrm{Var}}^{-1}(\hat\beta_W - \hat\beta_U)
#'     }
#'     and follows an approximate \eqn{F} distribution. A large value
#'     indicates that weights materially affect the regression.
#'
#' }
#'
#' In all cases, the reported statistic is an \eqn{F}-test with numerator
#' degrees of freedom equal to the number of auxiliary regressors added,
#' and denominator degrees of freedom equal to the residual degrees of
#' freedom from the augmented regression.
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
#' # Run weight-association diagnostic test; reports F-stat, df's, and p-value
#' results <- wa_test(fit, type = "DD")
#' print(results)
#'
#' @references
#' DuMouchel, W. H., & Duncan, G. J. (1983).
#'   Using sample survey weights in multiple regression analyses of stratified samples.
#'   *Journal of the American Statistical Association*, 78(383), 535-543.
#'
#' Pfeffermann, D., & Sverchkov, M. (1999).
#'   Parametric and semi-parametric estimation of regression models fitted to survey data.
#'   *Sankhya: The Indian Journal of Statistics, Series B*, 61(1), 166-186.
#'
#' Pfeffermann, D., & Sverchkov, M. (2003).
#'   Fitting generalized linear models under informative sampling.
#'   In R. L. Chambers & C. J. Skinner (Eds.), *Analysis of Survey Data*
#'   (pp. 175-196). Wiley.
#'
#' Wu, Y., & Fuller, W. A. (2005).
#'   Preliminary testing procedures for regression with survey samples.
#'   In *Proceedings of the Joint Statistical Meetings, Survey Research Methods Section*
#'   (pp. 3683-3688). American Statistical Association.
#'
#' @seealso
#' \code{\link{diff_in_coef_test}} for the Hausman-Pfeffermann difference-in-coefficients test,
#' and \code{\link{svytestCE}} for the example dataset included in this package.
#'
#' @importFrom survey svyglm
#'
#' @export
wa_test <- function(model, type = c("DD", "PS1", "PS1q", "PS2", "PS2q", "WF"),
                    coef_subset = NULL, aux_design = NULL, na.action = stats::na.omit) {

  # Argument check for svyglm object
  if (!inherits(model, "svyglm")) {
    stop("`model` must be an object of class 'svyglm' (from the survey package).")
  }

  # Argument check for type
  valid_types <- c("DD", "PS1", "PS1q", "PS2", "PS2q", "WF")
  if (length(type) != 1L || !is.character(type)) {
    stop("`type` must be a single character string.")
  }
  if (!type %in% valid_types) {
    stop("Invalid `type`: must be one of ", paste(valid_types, collapse = ", "), ".")
  }

  # Argument check for coefficient subset, custom auxiliary design, and na.action function arguments
  if (!is.null(coef_subset) && !is.character(coef_subset)) {
    stop("`coef_subset` must be a character vector of coefficient names.")
  }
  if (!is.null(aux_design) && !(is.function(aux_design) || is.matrix(aux_design))) {
    stop("`aux_design` must be NULL, a function, or a numeric matrix.")
  }
  if (!is.function(na.action)) {
    stop("`na.action` must be a function (e.g., stats::na.omit).")
  }

  # Extract design matrices
  wts <- tryCatch(stats::weights(model$survey.design),
                  error = function(e) stop("Could not extract weights from `model$survey.design`."))

  X <- tryCatch(stats::model.matrix(model),
                error = function(e) stop("Could not extract model matrix from `model`."))

  y <- tryCatch(stats::model.response(stats::model.frame(model)),
                error = function(e) stop("Could not extract response from `model`."))

  if (length(y) != nrow(X) || length(y) != length(wts)) {
    stop("Mismatch in dimensions of response, design matrix, and weights.")
  }

  # Handle missing data and assign the objects
  dat <- data.frame(y = y, X, wts = wts)
  dat <- na.action(dat)
  y <- dat$y
  X <- as.matrix(dat[, setdiff(names(dat), c("y", "wts"))])
  wts <- dat$wts

  # Subset coefficients if specified
  if (!is.null(coef_subset)) {
    keep_cols <- colnames(X) %in% coef_subset
    if (!any(keep_cols)) {
      stop("None of the names in `coef_subset` matched the model coefficients: ",
           paste(colnames(X), collapse = ", "))
    }
    X <- X[, keep_cols, drop = FALSE]
  }

  # Check auxiliary design matrix if provided
  if (!is.null(aux_design) && type %in% c("PS1", "PS1q", "PS2", "PS2q")) {
    A <- if (is.function(aux_design)) {
      aux_design(model, X, y, wts)
    } else {
      aux_design
    }
    if (!is.matrix(A) || !is.numeric(A)) {
      stop("`aux_design` must be a numeric matrix (or a function returning one).")
    }
    if (nrow(A) != nrow(X)) {
      stop("`aux_design` must have the same number of rows as the model matrix.")
    }

    qrA <- qr(A, LAPACK = TRUE)
    if (qrA$rank < ncol(A)) {
      stop("The supplied auxiliary design matrix is rank-deficient (perfect collinearity). ",
           "Please remove redundant columns or supply a full-rank matrix.")
    }
  }

  # Dispatch to test type functions
  out <- switch(type,
                "DD"   = wa_DD(y, X, wts),
                "PS1"  = wa_PS1(y, X, wts, quadratic = FALSE, aux_design = aux_design),
                "PS1q" = wa_PS1(y, X, wts, quadratic = TRUE, aux_design = aux_design),
                "PS2"  = wa_PS2(y, X, wts, quadratic = FALSE, aux_design = aux_design),
                "PS2q" = wa_PS2(y, X, wts, quadratic = TRUE, aux_design = aux_design),
                "WF"   = wa_WF(y, X, wts))

  # Return structured result
  structure(list(
      statistic = out$statistic,
      parameter = out$df,
      p.value   = out$p.value,
      method    = out$method,
      call      = match.call()
    ), class = "wa_test")
}



#' @rdname wa_test
#' @method print wa_test
#' @param x An object of class wa_test
#' @param ... Additional arguments passed to methods
#' @export
print.wa_test <- function(x, ...) {
  cat("\n", x$method, "\n", sep = "")
  cat("F =", formatC(x$statistic, digits = 4, format = "f"),
      " df1 =", x$parameter[1],
      " df2 =", x$parameter[2],
      " p-value =", formatC(x$p.value, digits = 4, format = "f"), "\n")
  invisible(x)
}


#' @rdname wa_test
#' @method summary wa_test
#' @param object An object of class wa_test
#' @param ... Additional arguments passed to methods
#' @export
summary.wa_test <- function(object, ...) {
  cat("\nWeight-Association Test\n")
  cat("Call:\n")
  print(object$call)
  cat("\nMethod:\n ", object$method, "\n", sep = "")
  cat("\nTest Statistic:\n")
  cat(" F =", formatC(object$statistic, digits = 4, format = "f"),
      " on", object$parameter[1], "and", object$parameter[2], "df",
      ", p-value =", formatC(object$p.value, digits = 4, format = "f"), "\n")
  invisible(object)
}


#' DuMouchel-Duncan WA test
#' @keywords internal
wa_DD <- function(y, X, wts) {
  W <- diag(wts, nrow = length(wts))
  X_tilde <- W %*% X
  X_comb <- cbind(X, X_tilde)

  betas_full <- qr.solve(X_comb, y)
  RSS_full <- sum((y - X_comb %*% betas_full)^2)

  X_reduced <- X
  betas_reduced <- qr.solve(X_reduced, y)
  RSS_reduced <- sum((y - X_reduced %*% betas_reduced)^2)

  df1 <- ncol(X_tilde)
  df2 <- length(y) - ncol(X_comb)
  F_stat <- ((RSS_reduced - RSS_full) / df1) / (RSS_full / df2)
  list(statistic = F_stat, df = c(df1, df2),
       p.value = 1 - stats::pf(F_stat, df1, df2),
       method = "DuMouchel-Duncan Weight-Association Test")
}

#' Wu-Fuller test
#' @keywords internal
wa_WF <- function(y, X, wts) {
  X_main <- X[, -1, drop = FALSE]

  # Auxiliary: regress weights on 1, X, X^2 -> construct q = w / w_hat
  X_aux <- cbind(1, X_main, X_main^2)
  eta <- qr.solve(X_aux, wts)
  w_hat <- X_aux %*% eta
  q <- as.numeric(wts / w_hat)

  # Full: y ~ X + (diag(q) %*% X_main)
  X_tilde <- diag(q) %*% X_main
  X_full <- cbind(1, X_main, X_tilde)
  b_full <- qr.solve(X_full, y)
  RSS_full <- sum((y - X_full %*% b_full)^2)

  # Reduced: y ~ X
  b_red <- qr.solve(X, y)
  RSS_red <- sum((y - X %*% b_red)^2)

  df1 <- ncol(X_tilde)
  df2 <- length(y) - ncol(X_full)
  F_stat <- ((RSS_red - RSS_full) / df1) / (RSS_full / df2)
  list(statistic = F_stat, df = c(df1, df2),
       p.value = 1 - stats::pf(F_stat, df1, df2),
       method = "Wu-Fuller Weight-Association Test")
}

#' Pfeffermann-Sverchkov Test 1
#' @keywords internal
wa_PS1 <- function(y, X, wts, quadratic = FALSE, aux_design = NULL) {
  betas_u <- qr.solve(X, y)
  residuals <- y - X %*% betas_u
  res_diag <- diag(as.vector(residuals))
  X_main <- X[, -1, drop = FALSE]

  if (!is.null(aux_design)) {
    extra <- if (is.function(aux_design)) aux_design(X_main, y) else aux_design
    method <- "Pfeffermann-Sverchkov WA Test 1 (Custom Aux)"
  } else if (quadratic) {
    extra <- cbind(X_main^2)
    method <- "Pfeffermann-Sverchkov WA Test 1 (Quadratic)"
  } else {
    extra <- NULL
    method <- "Pfeffermann-Sverchkov WA Test 1"
  }

  # Full model
  X_design <- cbind(1, X_main, extra, residuals, residuals^2, res_diag %*% X_main)
  .check_full_rank(X_design, "full PS1 design")
  betas <- qr.solve(X_design, wts)
  W_hat <- X_design %*% betas
  RSS <- sum((wts - W_hat)^2)

  # Reduced model
  X_reduced <- cbind(1, X_main)
  .check_full_rank(X_reduced, "reduced PS1 design")
  betas_reduced <- qr.solve(X_reduced, wts)
  W_hat_reduced <- X_reduced %*% betas_reduced
  TSS <- sum((wts - W_hat_reduced)^2)

  df1 <- ncol(X_design) - (ncol(X_main) + 1)
  df2 <- nrow(X_design) - ncol(X_design)
  F_stat <- ((TSS - RSS) / df1) / (RSS / df2)
  list(statistic = F_stat, df = c(df1, df2),
       p.value = 1 - stats::pf(F_stat, df1, df2),
       method = method)
}


#' Pfeffermann-Sverchkov Test 2
#' @keywords internal
wa_PS2 <- function(y, X, wts, quadratic = FALSE, aux_design = NULL) {
  X_main <- X[, -1, drop = FALSE]

  if (!is.null(aux_design)) {
    extra <- if (is.function(aux_design)) aux_design(X_main, y) else aux_design
    XY_full    <- cbind(1, X_main, extra, y, y^2)
    XY_reduced <- cbind(1, X_main, extra)
    method <- "Pfeffermann-Sverchkov WA Test 2 (Custom Aux)"
    added_cols <- 2
  } else if (quadratic) {
    XY_full    <- cbind(1, X_main, X_main^2, y, y^2)
    XY_reduced <- cbind(1, X_main, X_main^2)
    method <- "Pfeffermann-Sverchkov WA Test 2 (Quadratic)"
    added_cols <- 2
  } else {
    XY_full    <- cbind(1, X_main, y, y^2)
    XY_reduced <- cbind(1, X_main)
    method <- "Pfeffermann-Sverchkov WA Test 2"
    added_cols <- 2
  }

  .check_full_rank(XY_full, "full PS2 design")
  betas_full <- qr.solve(XY_full, wts)
  RSS_full <- sum((wts - XY_full %*% betas_full)^2)

  .check_full_rank(XY_reduced, "reduced PS2 design")
  betas_reduced <- qr.solve(XY_reduced, wts)
  RSS_reduced <- sum((wts - XY_reduced %*% betas_reduced)^2)

  df1 <- added_cols
  df2 <- length(wts) - ncol(XY_full)
  F_stat <- ((RSS_reduced - RSS_full) / df1) / (RSS_full / df2)
  list(statistic = F_stat, df = c(df1, df2),
       p.value = 1 - stats::pf(F_stat, df1, df2),
       method = method)
}


#' @rdname wa_test
#' @method tidy wa_test
#' @param x An object of class wa_test
#' @param ... Additional arguments passed to methods
tidy.wa_test <- function(x, ...) {
  tibble::tibble(
    term     = "weights_association",
    statistic = x$statistic,
    df1       = x$parameter[1],
    df2       = x$parameter[2],
    p.value   = x$p.value,
    method    = x$method
  )
}

#' @rdname wa_test
#' @method glance wa_test
#' @param x An object of class wa_test
#' @param ... Additional arguments passed to methods
glance.wa_test <- function(x, ...) {
  tibble::tibble(
    statistic = x$statistic,
    df1       = x$parameter[1],
    df2       = x$parameter[2],
    p.value   = x$p.value,
    method    = x$method
  )
}

#' Internal: check full rank before solving
#' @keywords internal
.check_full_rank <- function(M, name = "design matrix") {
  qrM <- qr(M, LAPACK = TRUE)
  if (qrM$rank < ncol(M)) {
    stop(sprintf(
      "The %s is rank-deficient (perfect collinearity detected). ",
      name
    ),
    "Please remove redundant columns or supply a full-rank auxiliary design.",
    call. = FALSE)
  }
  invisible(TRUE)
}
