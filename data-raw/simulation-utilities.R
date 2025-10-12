# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------

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
