#' Run All Diagnostic Tests for Informative Weights
#'
#' This function runs all implemented diagnostic tests:
#'   - \code{\link{wa_test}}() with types "DD", "PS1", "PS1q", "PS2", "PS2q", "WF"
#'   - \code{\link{diff_in_coef_test}}()
#'   - \code{\link{estim_eq_test}}()
#'   - \code{\link{perm_test}}() with stats "pred_mean" and "coef_mahal"
#'
#' @param model A fitted \code{svyglm} object.
#' @param alpha Critical value for rejection (default 0.05).
#' @param B Number of permutations for permutation tests (default 1000).
#'
#' @return A list with:
#'   \item{results}{Data frame of test names, statistics, p-values, reject indicator}
#'   \item{recommendation}{Character string with suggested action}
#'   \item{raw}{List of raw test outputs, including permutation test objects}
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
#' # Run all diagnostic tests and return a list of statistics, including a recommendation
#' results <- run_all_diagnostic_tests(fit)
#' print(results)
#'
#' @seealso \code{\link{wa_test}}, \code{\link{diff_in_coef_test}},
#'   \code{\link{estim_eq_test}}, \code{\link{perm_test}}
#'
#' @importFrom survey svyglm
#'
#' @export
run_all_diagnostic_tests <- function(model, alpha = 0.05, B = 1000) {
  # Argument checks
  if (!inherits(model, "svyglm")) {
    stop("`model` must be an object of class 'svyglm' (from the survey package).")
  }

  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single numeric value strictly between 0 and 1.")
  }

  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 1 || B != as.integer(B)) {
    stop("`B` must be a positive integer (number of permutations).")
  }
  B <- as.integer(B)

  # Function to ensure that all tests run safely
  safe_run <- function(expr) {
    tryCatch(expr, error = function(e) {
      warning("Test failed: ", conditionMessage(e))
      NULL
    })
  }

  # Run all tests
  tests <- list(
    DD   = safe_run(wa_test(model, type = "DD")),
    PS1  = safe_run(wa_test(model, type = "PS1")),
    PS1q = safe_run(wa_test(model, type = "PS1q")),
    PS2  = safe_run(wa_test(model, type = "PS2")),
    PS2q = safe_run(wa_test(model, type = "PS2q")),
    WF   = safe_run(wa_test(model, type = "WF")),
    HP   = safe_run(diff_in_coef_test(model, var_equal = TRUE)),
    PS3  = safe_run(estim_eq_test(model, q_method = "linear")),
    perm_mean  = safe_run(perm_test(model, stat = "pred_mean",  B = B, engine = "R")),
    perm_mahal = safe_run(perm_test(model, stat = "coef_mahal", B = B, engine = "R"))
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

  structure(list(
      results = res,
      recommendation = rec,
      raw = tests
    ), class = "run_all_diagnostic_tests")
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
