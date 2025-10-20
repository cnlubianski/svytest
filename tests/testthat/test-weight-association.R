# Tests using the provided example dataset
test_that("wa_test runs successfully with svytestCE dataset", {
  data("svytestCE", package = "svytest")

  # Create a survey design and fit a weighted regression model
  des <- svydesign(ids = ~1, weights = ~FINLWT21, data = svytestCE)
  fit <- svyglm(TOTEXPCQ ~ ROOMSQ + BATHRMQ + BEDROOMQ + FAM_SIZE + AGE, design = des)

  # Test all methods with the real dataset
  test_types <- c("DD", "PS1", "PS1q", "PS2", "PS2q", "WF")

  for (type in test_types) {
    result <- wa_test(fit, type = type)

    # Check structure
    expect_s3_class(result, "wa_test")
    expect_true(!is.null(result$statistic))
    expect_true(length(result$parameter) == 2)
    expect_true(is.numeric(result$p.value))
    expect_true(0 <= result$p.value && result$p.value <= 1)
    expect_true(!is.null(result$method))

    # Test print method doesn't error
    expect_output(print(result), "F =")

    # Test summary method doesn't error
    expect_output(summary(result), "Weight-Association Test")
  }
})

# Tests with simulated data - known informative weights
test_that("wa_test correctly identifies informative weights", {
  set.seed(1337)
  n <- 200

  # Create data where weights are STRONGLY associated with y
  X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
  weights <- exp(0.5 + rnorm(n, 0, 0.3))
  eps <- rnorm(n)

  # Y directly depends on weights
  y <- 1 + X1 + X2 + X3 + 0.7*weights + eps

  data <- data.frame(y = y, X1 = X1, X2 = X2, X3 = X3, weights = weights)
  des <- svydesign(ids = ~1, weights = ~weights, data = data)
  fit <- svyglm(y ~ X1, design = des)

  # All tests should reject the null (weights are informative)
  dd_test <- wa_test(fit, type = "DD")
  ps1_test <- wa_test(fit, type = "PS1")
  ps2_test <- wa_test(fit, type = "PS2")
  wf_test <- wa_test(fit, type = "WF")

  # With strong informativeness, p-values should be very small
  expect_lt(dd_test$p.value, 0.01)
  expect_lt(ps1_test$p.value, 0.01)
  expect_lt(ps2_test$p.value, 0.01)
  expect_lt(wf_test$p.value, 0.01)

  # The test statistics should be large
  expect_gt(dd_test$statistic, 3)
  expect_gt(ps1_test$statistic, 3)
  expect_gt(ps2_test$statistic, 3)
  expect_gt(wf_test$statistic, 3)
})

# Test with simulated data - known non-informative weights
test_that("wa_test correctly handles non-informative weights", {
  set.seed(789)
  n <- 200

  # Create data where weights are NOT associated with y
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  weights <- exp(0.5 + rnorm(n, 0, 0.3))  # Log-normal weights
  eps <- rnorm(n)

  # Y does NOT depend on weights
  y <- 1 + X1 + X2 + eps

  data <- data.frame(y=y, X1=X1, X2=X2, weights=weights)
  des <- svydesign(ids=~1, weights=~weights, data=data)
  fit <- svyglm(y ~ X1 + X2, design=des)

  # Tests should typically not reject the null (weights are not informative)
  # But p-values are random, so we can't always expect p > 0.05
  # Instead, test multiple times and ensure the distribution is reasonable
  p_values <- numeric(20)

  for (i in 1:20) {
    set.seed(1000 + i)
    X1 <- rnorm(n)
    X2 <- rnorm(n)
    weights <- exp(0.5 + rnorm(n, 0, 0.3))
    y <- 1 + X1 + X2 + rnorm(n)

    data <- data.frame(y=y, X1=X1, X2=X2, weights=weights)
    des <- svydesign(ids=~1, weights=~weights, data=data)
    fit <- svyglm(y ~ X1 + X2, design=des)

    p_values[i] <- wa_test(fit, type = "DD")$p.value
  }

  # Under the null, p-values should be roughly uniform
  # At least some p-values should be large
  expect_gte(mean(p_values > 0.25), 0.15)
})

# Test with subset of coefficients
test_that("wa_test handles coefficient subset correctly", {
  data("svytestCE", package = "svytest")
  des <- svydesign(ids = ~1, weights = ~FINLWT21, data = svytestCE)
  fit <- svyglm(TOTEXPCQ ~ ROOMSQ + BATHRMQ + BEDROOMQ + FAM_SIZE + AGE, design = des)

  # Test with all coefficients
  result_all <- wa_test(fit, type="DD")

  # Test with subset of coefficients
  result_subset <- wa_test(fit, type="DD", coef_subset=c("ROOMSQ", "BATHRMQ"))

  # Results should differ as we're testing different effects
  expect_false(identical(result_all$statistic, result_subset$statistic))
  expect_false(identical(result_all$parameter, result_subset$parameter))
})

# Test auxiliary design functionality
test_that("wa_test handles auxiliary design correctly", {
  set.seed(101)
  n <- 100
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  weights <- abs(1 + rnorm(n))
  y <- X1 + X2 + rnorm(n)

  data <- data.frame(y=y, X1=X1, X2=X2, weights=weights)
  des <- svydesign(ids=~1, weights=~weights, data=data)
  fit <- svyglm(y ~ X1 + X2, design=des)

  # Create a custom auxiliary design matrix
  aux_matrix <- cbind(sin(X1), cos(X2), X1*X2)
  colnames(aux_matrix) <- c("sin_X1", "cos_X2", "X1X2")

  # Test PS1 with auxiliary design
  result_aux <- wa_test(fit, type="PS1", aux_design=aux_matrix)
  expect_s3_class(result_aux, "wa_test")
  expect_true(grepl("Custom Aux", result_aux$method))

  # Test with auxiliary design function
  aux_fn <- function(X, y) {
    cbind(y^2, exp(y/10), X^2)
  }

  result_fn <- wa_test(fit, type="PS2", aux_design=aux_fn)
  expect_s3_class(result_fn, "wa_test")
  expect_true(grepl("Custom Aux", result_fn$method))
})

# Test error handling
test_that("wa_test produces appropriate errors", {
  set.seed(202)
  n <- 50
  X <- rnorm(n)
  weights <- abs(1 + rnorm(n))
  y <- X + rnorm(n)

  data <- data.frame(y=y, X=X, weights=weights)
  des <- svydesign(ids=~1, weights=~weights, data=data)
  fit <- svyglm(y ~ X, design=des)

  # Test with invalid model
  lm_fit <- lm(y ~ X, data=data)
  expect_error(wa_test(lm_fit), "Model must be of class 'svyglm'")

  # Test with invalid type
  expect_error(wa_test(fit, type="invalid"), "should be one of")

  # Test with non-existent coefficients
  expect_error(wa_test(fit, coef_subset="NonExistent"), "No matching coefficients found")
})

# Test with missing data
test_that("wa_test handles missing data correctly", {
  data("svytestCE", package = "svytest")

  # Introduce some missing values
  svytestCE_missing <- svytestCE
  idx <- sample(1:nrow(svytestCE), 20)
  svytestCE_missing$ROOMSQ[idx[1:5]] <- NA
  svytestCE_missing$BATHRMQ[idx[6:10]] <- NA
  svytestCE_missing$AGE[idx[11:15]] <- NA

  # Create design with missing data
  des <- svydesign(ids = ~1, weights = ~FINLWT21, data = svytestCE_missing)
  fit <- svyglm(TOTEXPCQ ~ ROOMSQ + BATHRMQ + BEDROOMQ + FAM_SIZE + AGE, design = des)

  # Test with default na.action
  result <- wa_test(fit, type="DD")
  expect_s3_class(result, "wa_test")

  # Custom na.action that removes only rows with NA in y or weights
  custom_na <- function(data) {
    data[!is.na(data$y) & !is.na(data$wts), ]
  }

  result_custom <- wa_test(fit, type="DD", na.action=custom_na)
  expect_s3_class(result_custom, "wa_test")
})

# Test tidy and glance methods
test_that("tidy and glance methods work correctly", {
  skip_if_not_installed("tibble")

  data("svytestCE", package = "svytest")
  des <- svydesign(ids = ~1, weights = ~FINLWT21, data = svytestCE)
  fit <- svyglm(TOTEXPCQ ~ ROOMSQ + BATHRMQ + BEDROOMQ, design = des)

  result <- wa_test(fit, type="DD")

  # Test tidy output
  tidy_output <- tidy(result)
  expect_s3_class(tidy_output, "tbl_df")
  expect_equal(nrow(tidy_output), 1)
  expect_true(all(c("term", "statistic", "p.value") %in% names(tidy_output)))

  # Test glance output
  glance_output <- glance(result)
  expect_s3_class(glance_output, "tbl_df")
  expect_equal(nrow(glance_output), 1)
  expect_true(all(c("statistic", "p.value", "method") %in% names(glance_output)))
})
