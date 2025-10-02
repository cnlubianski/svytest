#' Check Model Compatibility for DC Test
#'
#' Ensures the model is suitable for the Difference-in-Coefficients test.
#'
#' @param model An object of class \code{lm} or \code{svyglm}.
#' @return Throws an error if the model is incompatible.
#' @export
check_model_compatibility <- function(model) {
  if (!inherits(model, c("lm", "svyglm"))) {
    stop("Model must be of class 'lm' or 'svyglm'.")
  }
  if (inherits(model, "lm") && is.null(model$model)) {
    stop("Model frame not found. Refit with model = TRUE.")
  }
  if (inherits(model, "svyglm") && is.null(model$survey.design)) {
    stop("svyglm object missing survey.design.")
  }
  return(TRUE)
}
