library(dplyr)
library(tidyr)
library(survey)
library(rpms)
library(sampling)

set.seed(51483464)

N <- 3000
n <- c(100, 200)
sigma <- c(0.1, 0.2)
alpha <- c(0.0, 0.2, 0.4, 0.6)
delta <- c(1.5, 1)
cases <- expand_grid(N, n, sigma, delta, alpha)
case <- 32

generate_sample_brewer = function(data, w, n, rescale = FALSE) {
  pik = inclusionprobabilities(w, n)
  choosen = UPbrewer(pik)
  samp = data[1:nrow(data) * choosen,] %>%
    mutate(w = 1 / pik[1:nrow(data) * choosen])
  if (rescale == TRUE) mutate(samp, w = w / sum(w))
  return(samp)
}

# Quantitative x variable
generate_data_quant = function(N, sigma, alpha, delta) {
  X <- runif(N, 0, 1)
  u <- runif(N, 0, 1)
  epsilon <- rnorm(N, 0, sd = sigma)

  Y <- 1 + X + epsilon
  w <- alpha * Y + 0.3 * X + delta * u
  data <- data.frame(y = Y, x = X, w)
  return(data)
}

# Categorical x variable
generate_data_cat = function(N, sigma, alpha, delta) {
  X <- factor(sample(c("A", "B", "C"), N, replace = TRUE))
  u <- runif(N, 0, 1)
  epsilon <- rnorm(N, 0, sd = sigma)

  Y <- 1 + as.numeric(X) + epsilon
  w <- alpha * Y + 0.3 * as.numeric(X) + delta * u
  data <- data.frame(y = Y, x = X, w)
  return(data)
}

# oth cat and quant x variables
generate_data_both = function(N, sigma, alpha, delta) {
  X1 <- runif(N, 0, 1)
  X2 <- factor(sample(c("A", "B", "C"), N, replace = TRUE))
  u <- runif(N, 0, 1)
  epsilon <- rnorm(N, 0, sd = sigma)

  Y <- 1 + X1 * 0.5 + as.numeric(X2) * 0.5 + epsilon
  w <- alpha * Y + 0.3 * (X1 * 0.5 + as.numeric(X2) * 0.5) + delta * u
  data <- data.frame(y = Y, x1 = X1, x2 = X2, w)
  return(data)
}

pop <- generate_data_quant(N = cases$N[case],
                           sigma = cases$sigma[case],
                           alpha = cases$alpha[case],
                           delta = cases$delta[case])
pop <- generate_data_cat(N = cases$N[case],
                         sigma = cases$sigma[case],
                         alpha = cases$alpha[case],
                         delta = cases$delta[case])
pop <- generate_data_both(N = cases$N[case],
                          sigma = cases$sigma[case],
                          alpha = cases$alpha[case],
                          delta = cases$delta[case])
samp <- generate_sample_brewer(pop, w = pop$w, n = cases$n[case])

# Difference in Coefficients Test ----------------------------------------------

## Hausman-Pfeffermann DC Test =================================================
HP_DC_test <- function(data, y, x, wts) {
  # Unweighted Regression
  X <- cbind(1, x)
  betas_u <- solve(t(X) %*% X) %*% t(X) %*% y

  # Weighted Regression
  W <- diag(x = wts, nrow = length(wts), ncol = length(wts))
  betas_w <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y

  # Calculate Test Statistic

  # Calculate sigma^2 estimate from OLS under null
  y_hat <- X %*% betas_u
  residuals <- y - y_hat
  SSE <- sum(residuals^2)
  sigma_sq_hat <- SSE / (length(y) - length(betas_u) - 1)

  # A (note that in paper, H = diag(w) which is our W)
  A <- (solve(t(X) %*% W %*% X) %*% t(X) %*% W) - (solve(t(X) %*% X) %*% t(X))

  # Calculate variance
  V_hat <- sigma_sq_hat * A %*% t(A)

  # Chi Squared Test Statistic
  Chi_statistic <- t(betas_u - betas_w) %*% solve(V_hat) %*% (betas_u - betas_w)
  p_value <- pchisq(q = Chi_statistic, df = length(betas_u), lower.tail = FALSE)
  return(as.numeric(p_value))
}

HP_DC_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)

samp <- samp |>
  rename(x1 = x) |>
  mutate(x2 = x1 + rnorm(nrow(samp)))

HP_DC_test(data = samp, y = samp$y, x = samp$x1, wts = samp$w)

HP_DC_test(data = samp, y = samp$y, x = cbind(samp$x1, samp$x2), wts = samp$w)



HP_DC_test <- function(formula, data, wts) {
  # Extract response and predictor variables
  y <- model.response(model.frame(formula, data))
  X <- model.matrix(formula, data)
  weights <- data[[deparse(substitute(wts))]]

  # Unweighted Regression
  betas_u <- solve(t(X) %*% X) %*% t(X) %*% y

  # Weighted Regression
  W <- diag(weights)
  betas_w <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y

  # Calculate sigma^2 estimate from OLS under null
  y_hat <- X %*% betas_u
  residuals <- y - y_hat
  SSE <- sum(residuals^2)
  sigma_sq_hat <- SSE / (length(y) - length(betas_u) - 1)

  # A (note that in paper, H = diag(w) which is our W)
  A <- (solve(t(X) %*% W %*% X) %*% t(X) %*% W) - (solve(t(X) %*% X) %*% t(X))

  # Calculate variance
  V_hat <- sigma_sq_hat * A %*% t(A)

  # Chi Squared Test Statistic
  Chi_statistic <- t(betas_u - betas_w) %*% solve(V_hat) %*% (betas_u - betas_w)
  p_value <- pchisq(q = Chi_statistic, df = length(betas_u), lower.tail = FALSE) |> as.numeric()

  # Define class for the object
  class(p_value) <- "HP_DC_test"
  return(p_value)
}

# Print method for custom_lm object
print.HP_DC_test <- function(obj) {
  cat("P-value:\n")
  print(unclass(obj))
}

# Does work!
result <- HP_DC_test(y ~ x, data = samp, wts = w)
print(result)


