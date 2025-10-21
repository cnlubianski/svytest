test_that("perm_test rejects invalid arguments", {
  skip_if_not_installed("survey")
  library(survey)

  data(api, package = "survey")
  dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat)
  fit <- svyglm(api00 ~ ell + meals, design = dstrat)

  # Wrong model class
  expect_error(perm_test(lm(api00 ~ ell, data = apistrat)),
               "Model must be of class 'svyglm'")

  # Invalid stat
  expect_error(perm_test(fit, stat = "foobar", B = 10),
               "Invalid `stat`")

  # Invalid engine
  expect_error(perm_test(fit, engine = "python", B = 10),
               "Invalid `engine`")

  # Invalid B
  expect_error(perm_test(fit, B = 0), "must be a positive integer")
  expect_error(perm_test(fit, B = 1.5), "must be a positive integer")

  # Invalid coef_subset
  expect_error(perm_test(fit, coef_subset = 123, B = 10),
               "`coef_subset` must be a character vector")

  # Invalid block length
  expect_error(perm_test(fit, block = rep(1, 5), B = 10),
               "`block` length must equal")
})

test_that("perm_test runs with R engine and returns expected structure", {
  skip_if_not_installed("survey")
  library(survey)

  data(api, package = "survey")
  dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat)
  fit <- svyglm(api00 ~ ell + meals, design = dstrat)

  res <- perm_test(fit, stat = "pred_mean", B = 20, engine = "R")

  expect_s3_class(res, "perm_test")
  expect_true(is.numeric(res$stat_obs))
  expect_true(is.numeric(res$stat_null))
  expect_length(res$perm_stats, 20)
  expect_true(is.numeric(res$p.value))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("perm_test works with block permutation", {
  skip_if_not_installed("survey")
  library(survey)

  data(api, package = "survey")
  dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat)
  fit <- svyglm(api00 ~ ell + meals, design = dstrat)

  block <- apistrat$stype
  res <- perm_test(fit, stat = "pred_mean", B = 1000, block = apistrat$cnum, engine = "R")

  expect_s3_class(res, "perm_test")
  expect_length(res$perm_stats, 1000)
})

test_that("perm_test works with custom_fun", {
  skip_if_not_installed("survey")
  library(survey)

  data(api, package = "survey")
  dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat)
  fit <- svyglm(api00 ~ ell + meals, design = dstrat)

  my_fun <- function(model, X, y, wts) stats::weighted.mean(y, wts)  # trivial custom stat
  res <- perm_test(fit, custom_fun = my_fun, B = 50)

  expect_s3_class(res, "perm_test")
  expect_length(res$perm_stats, 50)
})

test_that("perm_test runs with C++ engine, and faster than R", {
  skip_if_not_installed("survey")
  skip_if_not_installed("Rcpp")
  library(survey)

  data(api, package = "survey")
  dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat)
  fit <- svyglm(api00 ~ ell + meals, design = dstrat)

  res <- perm_test(fit, stat = "pred_mean", B = 1000, engine = "C++")
  expect_s3_class(res, "perm_test")
  expect_length(res$perm_stats, 1000)

  time_r <- system.time({
    res_r <- perm_test(fit, stat = "pred_mean", B = 1000, engine = "R")
  })
  time_cpp <- system.time({
    res_cpp <- perm_test(fit, stat = "pred_mean", B = 1000, engine = "C++")
  })
  expect_lt(time_cpp["elapsed"], time_r["elapsed"])
})

test_that("perm_test C++ and R engines give similar results", {
  skip_if_not_installed("survey")
  skip_if_not_installed("Rcpp")
  library(survey)

  data(api, package = "survey")
  dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat)
  fit <- svyglm(api00 ~ ell + meals, design = dstrat)

  set.seed(123)
  res_r <- perm_test(fit, stat = "pred_mean", B = 500, engine = "R")
  set.seed(123)
  res_cpp <- perm_test(fit, stat = "pred_mean", B = 500, engine = "C++")

  expect_equal(res_r$stat_obs, res_cpp$stat_obs)
  expect_equal(res_r$perm_stats, res_cpp$perm_stats)
  expect_equal(res_r$p.value, res_cpp$p.value)
})
