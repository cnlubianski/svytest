test_that("estim_eq_test rejects invalid arguments", {
  skip_if_not_installed("survey")
  library(survey)

  data(api, package = "survey")
  dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat)
  fit <- svyglm(api00 ~ ell + meals, design = dstrat)

  # Wrong model class
  expect_error(estim_eq_test(lm(api00 ~ ell, data = apistrat)),
               "`model` must be an object of class 'svyglm'")

  # Wrong family
  fit_binom <- svyglm(I(api00 > 700) ~ ell, design = dstrat, family = quasibinomial())
  expect_error(estim_eq_test(fit_binom), "supports only gaussian\\(identity\\) models")

  # Invalid q_method
  expect_error(estim_eq_test(fit, q_method = "bad"), "`q_method` must be one of")

  # Invalid coef_subset
  expect_error(estim_eq_test(fit, q_method = "linear", coef_subset = 123),
               "`coef_subset` must be a character vector")

  # Invalid stabilize
  expect_error(estim_eq_test(fit, q_method = "linear", stabilize = "yes"),
               "`stabilize` must be a single logical value")

  # Invalid na.action
  expect_error(estim_eq_test(fit, q_method = "linear", na.action = "omit"),
               "`na.action` must be a function")
})

test_that("estim_eq_test returns expected structure with linear q_method", {
  skip_if_not_installed("survey")
  library(survey)

  data(api, package = "survey")
  dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat)
  fit <- svyglm(api00 ~ ell + meals, design = dstrat)

  res <- estim_eq_test(fit, q_method = "linear")
  expect_s3_class(res, "estim_eq_test")
  expect_true(is.numeric(res$statistic))
  expect_true(is.numeric(res$p.value))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
  expect_equal(res$df1, length(res$terms))
  expect_equal(res$df2, res$n - length(res$terms))
})

test_that("estim_eq_test works with log q_method and stabilize = FALSE", {
  skip_if_not_installed("survey")
  library(survey)

  data(api, package = "survey")
  dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat)
  fit <- svyglm(api00 ~ ell + meals, design = dstrat)

  res <- estim_eq_test(fit, q_method = "log", stabilize = FALSE)
  expect_s3_class(res, "estim_eq_test")
  expect_true(is.numeric(res$statistic))
  expect_true(is.numeric(res$p.value))
})

test_that("estim_eq_test works with coefficient subset", {
  skip_if_not_installed("survey")
  library(survey)

  data(api, package = "survey")
  dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat)
  fit <- svyglm(api00 ~ ell + meals + mobility, design = dstrat)

  res <- estim_eq_test(fit, q_method = "linear", coef_subset = c("ell", "meals"))
  expect_s3_class(res, "estim_eq_test")
  expect_true(all(res$terms %in% c("ell", "meals")))
})
