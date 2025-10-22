test_that("diff_in_coef_test works as expected", {
  # Load in survey package (required) and load in example data
  library(survey)
  data("svytestCE", package = "svytest")

  # Create a survey design and fit a weighted regression model
  formula <- TOTEXPCQ ~ ROOMSQ + BATHRMQ + BEDROOMQ + FAM_SIZE + AGE
  des <- svydesign(ids = ~1, weights = ~FINLWT21, data = svytestCE)
  fit <- svyglm(formula = formula, design = des)
  lm_fit <- lm(formula, data = svytestCE)

  # Run the difference in coefficients test
  test_result <- diff_in_coef_test(fit, var_equal = TRUE)

  # Check that the result is a list with expected components
  expect_type(test_result, "list")
  expect_true(all(c("statistic", "p.value") %in% names(test_result)))

  # Check that statistic is numeric and p-value is between 0 and 1
  expect_type(test_result$statistic, "double")
  expect_gte(test_result$p.value, 0)
  expect_lte(test_result$p.value, 1)

  # Check specific values (these values may change if the underlying data changes)
  expect_equal(round(test_result$statistic, 2), 15.33)
  expect_equal(round(test_result$p.value, 4), 0.0179)

  # The weighted coefficients should match exactly what svyglm produces
  expect_equal(coef(fit), test_result$betas_weighted)

  # The unweighted coefficients should match exactly what lm produces
  expect_equal(coef(lm_fit), test_result$betas_unweighted)
})
