test_that("wa_test rejects invalid arguments", {
  skip_if_not_installed("survey")
  library(survey)

  data(api, package = "survey")
  dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat)
  fit <- svyglm(api00 ~ ell + meals, design = dstrat)

  # Wrong model class
  expect_error(wa_test(lm(api00 ~ ell, data = apistrat)),
               "`model` must be an object of class 'svyglm'")

  # Invalid type
  expect_error(wa_test(fit, type = "BADTYPE"), "Invalid `type`")

  # Invalid coef_subset
  expect_error(wa_test(fit, type = "DD", coef_subset = 123),
               "`coef_subset` must be a character vector")

  # Invalid aux_design
  expect_error(wa_test(fit, type = "DD", aux_design = "not_a_matrix"),
               "`aux_design` must be NULL, a function, or a numeric matrix")

  # Invalid na.action
  expect_error(wa_test(fit, type = "DD", na.action = "omit"), "`na.action` must be a function")
})

test_that("wa_test returns expected structure for each type", {
  skip_if_not_installed("survey")
  library(survey)

  data(api, package = "survey")
  dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat)
  fit <- svyglm(api00 ~ ell + meals, design = dstrat)

  for (t in c("DD", "PS1", "PS1q", "PS2", "PS2q", "WF")) {
    res <- wa_test(fit, type = t)
    expect_s3_class(res, "wa_test")
    expect_true(is.numeric(res$statistic))
    expect_true(is.numeric(res$p.value))
    expect_true(res$p.value >= 0 && res$p.value <= 1)
    expect_true(is.character(res$method))
  }
})

test_that("wa_test works with coefficient subset", {
  skip_if_not_installed("survey")
  library(survey)

  data(api, package = "survey")
  dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat)
  fit <- svyglm(api00 ~ ell + meals + mobility, design = dstrat)

  res <- wa_test(fit, type = "DD", coef_subset = c("(Intercept)", "ell"))
  expect_s3_class(res, "wa_test")
  expect_true(is.numeric(res$statistic))
})

# More testing for auxiliary design matrices
