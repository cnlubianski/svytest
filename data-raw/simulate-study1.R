# ------------------------------------------------------------------------------
# Title: Parametric Diagnostic Tests Simulation
# Description: This script runs a simulation for the parametric diagnostic tests
# implemented throughout the package exported functions. This is used to
# motivate the use of the tests and compare it with the parametric alternatives.
# The simulation design is based off of one implemented by Wang et al. (2023).
# ------------------------------------------------------------------------------

# Preamble ---------------------------------------------------------------------
# Load in libraries
library(dplyr)
library(tidyr)
library(sampling)
library(survey)
library(future.apply) # for parallelization

# Set working directory and define script/storage paths
working_dir <- getwd()
if (!file.exists("svytest.Rproj")) stop("Working directory must be project root.")
sim_script_dir <- file.path(working_dir, "data-raw")
sim_storage_dir <- file.path(working_dir, "inst/extdata")

# Define simulation attributes
set.seed(1337)
sim_rep <- 1000 # How many times a specific case is run
critical_value <- 0.05

# Define cases and construct case grid
N <- 3000
n <- c(100, 200)
sigma <- c(0.1, 0.2)
alpha <- c(0, 0.2, 0.4, 0.6)
delta <- c(1.5, 1)
cases <- tidyr::expand_grid(N, n, sigma, delta, alpha)

# General simulation functions -------------------------------------------------
# Retrieve simulation utility functions
source(file.path(sim_script_dir, "simulation-utilities.R"))

# Function to run one simulation replication
run_one_rep <- function(N, n, sigma, alpha, delta) {
  # Generate population and sample
  pop <- generate_data_study1(N, sigma, alpha, delta)
  samp <- generate_sample_brewer(pop, w = pop$w, n = n)

  # Construct svyglm model object for diagnostic functions
  design <- svydesign(id = ~1, weights = ~w, data = samp)
  svyglm_model <- svyglm(y ~ x, design = design)

  # Call the package functions and return the p-values
  out <- list(
    DD   = wa_test(model = svyglm_model, type = "DD")$p.value,
    PS1  = wa_test(model = svyglm_model, type = "PS1")$p.value,
    PS1q = wa_test(model = svyglm_model, type = "PS1q")$p.value,
    PS2  = wa_test(model = svyglm_model, type = "PS2")$p.value,
    PS2q = wa_test(model = svyglm_model, type = "PS2q")$p.value,
    WF   = wa_test(model = svyglm_model, type = "WF")$p.value,
    HP   = diff_in_coef_test(model = svyglm_model, var_equal = TRUE)$p.value,
    PS3  = estim_eq_test(model = svyglm_model, q_method = "linear")$p.value
  )
  return(out)
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

    case_df <- tibble(row = seq_len(B), data = reps) |> unnest_wider(data)
    results[[case]] <- case_df
  }

  results <- bind_rows(results)
  return(results)
}

# Run simulation, wrangle results into table, then save table ------------------
res <- simulate_study1(cases, B = sim_rep)

# Summarize rejection rates
reject_table <- res %>%
  mutate(across(c(DD, PS1, PS1q, PS2, PS2q, WF, HP, PS3),
                ~ ifelse(. <= critical_value, 1, 0)),
         case = rep(seq_len(nrow(cases)), each = sim_rep)) %>%
  group_by(case) %>%
  summarize(across(-row, mean)) %>%
  left_join(cases %>% mutate(case = row_number()), by = "case") %>%
  select(n, sigma, delta, alpha, DD, HP, PS1, PS1q, PS2, PS2q, PS3)

# Save the results
saveRDS(reject_table, file.path(sim_storage_dir, "st1_reject_table.rds"))
