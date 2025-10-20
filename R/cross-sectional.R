#' Run All Diagnostic Tests for Informative Weights
#'
#' This function runs all implemented diagnostic tests:
#'   - \code{\link{wa_test}}() with types "DD","PS1","PS1q","PS2","PS2q","WF"
#'   - \code{\link{diff_in_coef_test}}()
#'   - \code{\link{estim_eq_test}}()
#'   - \code{\link{perm_test}}() with stats "pred_mean" and "coef_mahal"
#'
#' @param model A fitted \code{svyglm} object.
#' @param alpha Critical value for rejection (default 0.05).
#'
#' @return A list with:
#'   \item{results}{Data frame of test names, statistics, p-values, reject indicator}
#'   \item{recommendation}{Character string with suggested action}
#'   \item{raw}{List of raw test outputs, including permutation test objects}
#'
#' @examples
#' if (requireNamespace("survey", quietly = TRUE)) {
#'   # Load in survey package (required) and load in example data
#'   library(survey)
#'   data("svytestCE", package = "svytest")
#'
#'   # Create a survey design and fit a weighted regression model
#'   des <- svydesign(ids = ~1, weights = ~FINLWT21, data = svytestCE)
#'   fit <- svyglm(TOTEXPCQ ~ ROOMSQ + BATHRMQ + BEDROOMQ + FAM_SIZE + AGE, design = des)
#'
#'   # Run all diagnostic tests and return a list of statistics, including a recommendation
#'   results <- run_all_diagnostic_tests(fit)
#'   print(results)
#' }
#'
#' @seealso \code{\link{wa_test}}, \code{\link{diff_in_coef_test}},
#'   \code{\link{estim_eq_test}}, \code{\link{perm_test}}
#'
#' @export
run_all_diagnostic_tests <- function(model, alpha = 0.05, B = 1000) {
  # Run all tests
  tests <- list(
    DD   = wa_test(model, type = "DD"),
    PS1  = wa_test(model, type = "PS1"),
    PS1q = wa_test(model, type = "PS1q"),
    PS2  = wa_test(model, type = "PS2"),
    PS2q = wa_test(model, type = "PS2q"),
    WF   = wa_test(model, type = "WF"),
    HP   = diff_in_coef_test(model, var_equal = TRUE),
    PS3  = estim_eq_test(model),
    perm_mean  = perm_test(model, stat = "pred_mean",  B = B, engine = "R"),
    perm_mahal = perm_test(model, stat = "coef_mahal", B = B, engine = "R")
  )

  # Extract comparable results into a tidy data frame
  res <- lapply(names(tests), function(nm) {
    obj <- tests[[nm]]
    tibble::tibble(
      test = nm,
      statistic = if (!is.null(obj$statistic)) unname(obj$statistic) else NA_real_,
      p.value   = if (!is.null(obj$p.value)) obj$p.value else NA_real_,
      reject    = if (!is.null(obj$p.value)) obj$p.value <= alpha else NA
    )
  }) |> dplyr::bind_rows()

  # Recommendation logic
  if (any(res$reject, na.rm = TRUE)) {
    rec <- paste0("At least one test rejects H0 at significance level = ", alpha,
                  ". Recommendation: use survey weights in regression.")
  } else {
    rec <- paste0("No test rejects H0 at significance level = ", alpha,
                  ". Recommendation: unweighted regression may be acceptable.")
  }

  structure(
    list(
      results = res,
      recommendation = rec,
      raw = tests
    ),
    class = "run_all_diagnostic_tests"
  )
}

#' @rdname run_all_diagnostic_tests
#' @method print run_all_diagnostic_tests
#' @param x An object of class run_all_diagnostic_tests
#' @param ... Additional arguments passed to methods
#' @export
print.run_all_diagnostic_tests <- function(x, ...) {
  cat("\nDiagnostic Tests for Informative Weights\n")
  print(x$results)
  cat("\nRecommendation:\n", x$recommendation, "\n")
  invisible(x)
}
