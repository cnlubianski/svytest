#' svytest: Survey Weight Regression Diagnostics
#'
#' The **svytest** package provides diagnostic tools for assessing the role of
#' survey weights in regression modeling. It implements formal tests such as
#' the Hausman-Pfeffermann Difference-in-Coefficients test, along with helper
#' datasets and utilities for reproducible survey analysis.
#'
#' @section Functions:
#' \itemize{
#'   \item \code{\link{diff_in_coef_test}}: Hausman-Pfeffermann Difference-in-Coefficients test.
#'   \item (future functions will be listed here as you add them).
#' }
#'
#' @section Datasets:
#' \itemize{
#'   \item \code{\link{svytestCE}}: Curated subset of the Consumer Expenditure dataset,
#'         useful for demonstrating survey-weighted regression diagnostics.
#' }
#'
#' @details
#' The package builds on the general Hausman specification test
#' (Hausman, 1978) and its adaptation to survey weights by Pfeffermann (1993).
#' These methods allow researchers to formally test whether survey weights
#' materially affect regression estimates.
#'
#' @references
#' Hausman, J. A. (1978). Specification Tests in Econometrics.
#'   *Econometrica*, 46(6), 1251-1271. \doi{10.2307/1913827}
#'
#' Pfeffermann, D. (1993). The Role of Sampling Weights When Modeling Survey Data.
#'   *International Statistical Review*, 61(2), 317-337. \doi{10.2307/1403631}
#'
#' @seealso
#' \code{\link{diff_in_coef_test}}, \code{\link{svytestCE}},
#' \code{\link{wa_test}}, \code{\link{perm_test}}
#'
#' @docType package
#' @name svytest
#' @keywords internal
"_PACKAGE"
