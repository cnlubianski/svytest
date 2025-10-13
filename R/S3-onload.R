.onLoad <- function(libname, pkgname) {
  # Register glance/tidy methods only if broom is available
  if (requireNamespace("broom", quietly = TRUE)) {
    s3_register("broom::glance", "diff_in_coef_test")
    s3_register("broom::tidy",   "diff_in_coef_test")
    s3_register("broom::glance", "estim_eq_test")
    s3_register("broom::tidy",   "estim_eq_test")
    s3_register("broom::glance", "wa_test")
    s3_register("broom::tidy",   "wa_test")
    s3_register("broom::glance", "perm_test")
    s3_register("broom::tidy",   "perm_test")
    s3_register("broom::glance", "pred_power_test")
    s3_register("broom::tidy",   "pred_power_test")
  }
}

s3_register <- function(generic, class, method = NULL) {
  stopifnot(is.character(generic), length(generic) == 1)
  pieces <- strsplit(generic, "::")[[1]]
  stopifnot(length(pieces) == 2)
  package <- pieces[[1]]
  generic <- pieces[[2]]

  if (is.null(method)) {
    method <- paste0(generic, ".", class)
  }

  # Only register if the generic's package is available
  if (requireNamespace(package, quietly = TRUE)) {
    fun <- get(method, envir = parent.frame())
    registerS3method(generic, class, fun, envir = asNamespace(package))
  }
}
