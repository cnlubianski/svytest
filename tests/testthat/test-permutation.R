# Load in functions
devtools::load_all()

test_that("perm_test works as expected", {
  # Load in survey package (required) and load in example data
  library(survey)
  data("svytestCE", package = "svytest")

  # Create a survey design and fit a weighted regression model
  formula <- TOTEXPCQ ~ ROOMSQ + BATHRMQ + BEDROOMQ + FAM_SIZE + AGE
  des <- svydesign(ids = ~1, weights = ~FINLWT21, data = svytestCE)
  fit <- svyglm(formula = formula, design = des)

  # Run the permutation test for predictive mean
  test_result_mean <- perm_test(fit, stat = "pred_mean", B = 500, engine = "R")

  # Check that the result is of class 'perm_test'
  expect_s3_class(test_result_mean, "perm_test")

})

test_that("perm_test handles errors correctly", {
  # Load in survey package (required) and load in example data
  library(survey)
  data("svytestCE", package = "svytest")

  # Create a survey design and fit a weighted regression model
  formula <- TOTEXPCQ ~ ROOMSQ + BATHRMQ + BEDROOMQ + FAM_SIZE + AGE n
  des <- svydesign(ids = ~1, weights = ~FINLWT21, data = svytestCE)
  fit <- svyglm(formula = formula, design = des)

  # Check that an error is thrown for unsupported model types
  fit_logistic <- svyglm(I(FAM_SIZE > 3) ~ ROOMSQ + BATHRMQ + BEDROOMQ + AGE,
                         design = des, family = binomial())
  expect_error(perm_test(fit_logistic, stat = "pred_mean", B = 500, engine = "R"),
               "perm_test currently supports only gaussian\\(identity\\) models.")

  # Check that an error is thrown for invalid statistic choice
  expect_error(perm_test(fit, stat = "invalid_stat", B = 500, engine = "R"),
               "Invalid 'stat' argument. Choose 'pred_mean' or 'coef_mahal'.")

})
