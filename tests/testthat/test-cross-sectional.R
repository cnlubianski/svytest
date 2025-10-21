test_that("run_all_diagnostic_tests rejects invalid arguments", {
  skip_if_not_installed("survey")
  library(survey)

  data(api, package = "survey")
  dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat)
  fit <- svyglm(api00 ~ ell + meals, design = dstrat)

  # Wrong model class
  expect_error(run_all_diagnostic_tests(lm(api00 ~ ell, data = apistrat)),
               "`model` must be an object of class 'svyglm'")

  # Invalid alpha
  expect_error(run_all_diagnostic_tests(fit, alpha = -0.1),
               "`alpha` must be a single numeric value strictly between 0 and 1")
  expect_error(run_all_diagnostic_tests(fit, alpha = 2),
               "`alpha` must be a single numeric value strictly between 0 and 1")

  # Invalid B
  expect_error(run_all_diagnostic_tests(fit, B = 0),
               "`B` must be a positive integer")
  expect_error(run_all_diagnostic_tests(fit, B = 1.5),
               "`B` must be a positive integer")
})


test_that("run_all_diagnostic_tests returns expected structure", {
  skip_if_not_installed("survey")
  library(survey)

  data(api, package = "survey")
  dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat)
  fit <- svyglm(api00 ~ ell + meals, design = dstrat)

  res <- run_all_diagnostic_tests(fit, alpha = 0.05, B = 500)
  expect_s3_class(res, "run_all_diagnostic_tests")
  expect_true("results" %in% names(res))
  expect_true("recommendation" %in% names(res))
  expect_true("raw" %in% names(res))
  expect_true(all(c("test", "statistic", "p.value", "reject") %in% names(res$results)))
})

# A small segment of the large simulation to check monotone increasing rejection rates
test_that("each diagnostic's rejection rate increases with informativeness", {
  skip_on_cran()
  skip_if_not_installed("survey")
  library(survey)
  library(dplyr)
  library(sampling)
  library(future.apply) # for parallelization
  library(tidyr)

  set.seed(1337)

  # Simulation controls
  N <- 3000
  n <- 200
  sim_rep <- 100
  alpha <- c(0, 0.2, 0.4, 0.6)
  critical_value <- 0.05
  B_perm <- 500
  sigma <- 0.2
  delta <- 1
  cases <- tidyr::expand_grid(N, n, sigma, delta, alpha)
  test_names <- c("DD", "PS1", "PS1q", "PS2", "PS2q", "WF", "HP", "PS3", "perm_mean", "perm_mahal")

  # Read simulation functions in directly
  # Function to generate data from Wang et al.'s (2023) Study 1
  generate_data_study1 <- function(N, sigma, alpha, delta) {
    X <- runif(N, 0, 1)
    u <- runif(N, 0, 1)
    epsilon <- rnorm(N, 0, sd = sigma)
    Y <- 1 + X + epsilon
    w <- alpha * Y + 0.3 * X + delta * u
    data.frame(y = Y, x = X, w = w)
  }

  # Function to retrieve sample from population using Brewer sampling given n sample size
  generate_sample_brewer <- function(data, w, n) {
    pik <- sampling::inclusionprobabilities(w, n)
    choosen <- sampling::UPbrewer(pik)
    samp <- data[1:nrow(data) * choosen, ] |>
      mutate(w = 1 / pik[1:nrow(data) * choosen])
    samp
  }

  # Function to run one simulation replication
  run_one_rep <- function(N, n, sigma, alpha, delta) {
    # Generate population and sample
    pop <- generate_data_study1(N, sigma, alpha, delta)
    samp <- generate_sample_brewer(pop, w = pop$w, n = n)

    # Construct svyglm model object for diagnostic functions
    design <- svydesign(id = ~1, weights = ~w, data = samp)
    svyglm_model <- svyglm(y ~ x, design = design)

    # Call the package functions and return the p-values
    out <- run_all_diagnostic_tests(svyglm_model, alpha = critical_value, B = B_perm)
    return(out[['results']]$p.value)
  }

  # Function to perform the simulation
  simulate_study1 <- function(cases, B = sim_rep) {
    # Initiate a storage object
    results <- list()

    # Iterate through each case
    for (case in 1:nrow(cases)) {
      param <- cases[case, ]
      message("Running case ", case, " of ", nrow(cases))

      # Parallel replications
      reps <- future_lapply(seq_len(B), function(b) {
        run_one_rep(param$N, param$n, param$sigma, param$alpha, param$delta)
      }, future.seed = TRUE)

      case_df <- do.call(rbind, reps) |>
        as.data.frame() |>
        setNames(test_names) |>
        tibble::rownames_to_column("replication") |>
        dplyr::mutate(replication = as.integer(replication))

      results[[case]] <- case_df
    }

    results <- bind_rows(results)
    return(results)
  }

  res <- simulate_study1(cases, B = sim_rep)

  # Summarize rejection rates
  reject_table <- res %>%
    mutate(across(c(DD, PS1, PS1q, PS2, PS2q, WF, HP, PS3, perm_mean, perm_mahal),
                  ~ ifelse(. <= critical_value, 1, 0)),
           case = rep(seq_len(nrow(cases)), each = sim_rep)) %>%
    group_by(case) %>%
    summarize(across(-replication, mean)) %>%
    left_join(cases %>% mutate(case = row_number()), by = "case") %>%
    select(n, sigma, delta, alpha, DD, HP, PS1, PS1q, PS2, PS2q, PS3, perm_mean, perm_mahal)

  # Check monotone increasing rejection rates for each test
  for (test in test_names) {
    expect_true(all(diff(reject_table[[test]]) >= 0),
                info = paste("Column", nm, "is not monotone increasing"))
  }
})
