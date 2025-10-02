### Set-up
if (!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
if (!require(rpms)) install.packages("rpms"); library(rpms)
if (!require(sampling)) install.packages("sampling"); library(sampling)
if (!require(survey)) install.packages("survey"); library(survey)

generate_data_study1 = function(N, sigma, alpha, delta) {
  X <- runif(N, 0, 1)
  u <- runif(N, 0, 1)
  epsilon <- rnorm(N, 0, sd = sigma)

  Y <- 1 + X + epsilon
  w <- alpha * Y + 0.3 * X + delta * u
  data = data.frame(y = Y, x = X, w)
  return(data)
}

generate_data_study2 = function(N, sigma, alpha) {
  X <- runif(N, 0, 1)
  u <- runif(N, 0, 1)
  epsilon <- rnorm(N, 0, sd = sigma)

  Y <- 1 + X + epsilon
  w <- alpha * (Y - 1.5 * alpha)^2 + 0.3 * X - 0.3 * X^2 + u
  data = data.frame(y = Y, x = X, w)
  return(data)
}

generate_data_study3 = function(N, alpha, psi) {
  X <- rnorm(N, 0, sd = sqrt(0.5))
  epsilon <- rnorm(N, 0, sd = sqrt(0.5))
  z <- rnorm(N, 0, sd = sqrt(0.5))
  beta = 2 - alpha

  eta = function(x) {
    ifelse(x < 0.2, 0.025,
           ifelse(0.2 <= x & x <= 1.2, 0.475 * (x - 0.2) + 0.025, 0.5))
  }
  Y <- 0.5 + X + epsilon
  w <- alpha * eta(X) + beta * eta(psi * epsilon + (1 - psi) * z)
  data = data.frame(y = Y, x = X, w)
  return(data)
}


generate_sample_brewer = function(data, w, n, rescale = FALSE) {
  pik = inclusionprobabilities(w, n)
  choosen = UPbrewer(pik)
  samp = data[1:nrow(data) * choosen,] %>%
    mutate(w = 1 / pik[1:nrow(data) * choosen])
  if (rescale == TRUE) mutate(samp, w = w / sum(w))
  return(samp)
}


generate_sample_poisson = function(data, w, n, rescale = FALSE) {
  choosen = as.numeric(runif(length(w), 0, (1 / n) * sum(w)) < w)
  samp = cbind(data, choosen) %>%
    filter(choosen == 1) %>%
    select(-choosen) %>%
    mutate(w = 1 / w) # Redefine from pi to weights w
  if (rescale == TRUE) mutate(samp, w = w / sum(w))
  return(samp)
}

generate_sample_poisson_not_presumed = function(data, w, n, rescale = FALSE) {
  choosen = as.numeric(runif(length(w), 0, (1 / n) * sum(w)) < w)
  samp = cbind(data, choosen) %>%
    filter(choosen == 1) %>%
    select(-choosen)
  if (rescale == TRUE) mutate(samp, w = w / sum(w))
  return(samp)
}

generate_sample_brewer_not_presumed = function(data, w, n, rescale = FALSE) {
  pik = inclusionprobabilities(w, n)
  choosen = UPbrewer(pik)
  samp = data[1:nrow(data) * choosen,]
  if (rescale == TRUE) mutate(samp, w = w / sum(w))
  return(samp)
}


## Study 1

set.seed(51483464)
B = 1000 # 10000 for large study # 1000 for normal

N <- 3000
n <- c(100, 200)
sigma <- c(0.1, 0.2)
alpha <- c(0.0, 0.2, 0.4, 0.6)
delta <- c(1.5, 1)
cases <- expand_grid(N, n, sigma, delta, alpha)

columns = c("case", "iteration", "DD", "PN", "HP", "PS1", "PS1q", "PS2", "PS2q",
            "PS3", "WF", "LR")
results = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(results) = columns

for (case in 1:nrow(cases)) {
  DD = PN = HP = PS1 = PS1q = PS2 = PS2q = PS3 = WF = LR = rep(NA, B)
  case_storage = data.frame(iteration = seq_len(B), DD, PN, HP, PS1, PS1q,
                            PS2, PS2q, PS3, WF, LR)
  for (b in 1:B) {
    pop = generate_data_study1(N = cases$N[case],
                               sigma = cases$sigma[case],
                               alpha = cases$alpha[case],
                               delta = cases$delta[case])
    samp = generate_sample_brewer(pop, w = pop$w, n = cases$n[case])

    case_storage$DD[b] = DD_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$PN[b] = PN_test(data = samp, y = samp$y, x = samp$x, wts = samp$w, est_split = 0.5)
    case_storage$HP[b] = HP_DC_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$PS1[b] = PS1_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$PS1q[b] = PS1q_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$PS2[b] = PS2_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$PS2q[b] = PS2q_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$PS3[b] = PS3_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$WF[b] = WF_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$LR[b] = LR_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
  }

  results = rbind(results, cbind(case, case_storage))
  print(case)
}
write.csv(results, "st1_results.csv")
#write.csv(results, "st1_results_large.csv")

results = read.csv("st1_results.csv")
#results = read.csv("st1_results_large.csv")

reject = results %>%
  mutate(HP = case_when(HP <= 0.05 ~ 1, TRUE ~ 0),
         DD = case_when(DD <= 0.05 ~ 1, TRUE ~ 0),
         PS1 = case_when(PS1 <= 0.05 ~ 1, TRUE ~ 0),
         PS1q = case_when(PS1q <= 0.05 ~ 1, TRUE ~ 0),
         PS2 = case_when(PS2 <= 0.05 ~ 1, TRUE ~ 0),
         PS2q = case_when(PS2q <= 0.05 ~ 1, TRUE ~ 0),
         PS3 = case_when(PS3 <= 0.05 ~ 1, TRUE ~ 0),
         WF = case_when(WF <= 0.05 ~ 1, TRUE ~ 0),
         LR = case_when(LR <= 0.05 ~ 1, TRUE ~ 0),
         PN = case_when(PN <= 0.05 ~ 1, TRUE ~ 0)) %>%
  select(-iteration) %>%
  group_by(case) %>%
  summarize(across(everything(), mean)) %>%
  mutate(HP = format(round(HP * 100, 1), nsmall = 1),
         DD = format(round(DD * 100, 1), nsmall = 1),
         PS1 = format(round(PS1 * 100, 1), nsmall = 1),
         PS1q = format(round(PS1q * 100, 1), nsmall = 1),
         PS2 = format(round(PS2 * 100, 1), nsmall = 1),
         PS2q = format(round(PS2q * 100, 1), nsmall = 1),
         PS3 = format(round(PS3 * 100, 1), nsmall = 1),
         WF = format(round(WF * 100, 1), nsmall = 1),
         LR = format(round(LR * 100, 1), nsmall = 1),
         PN = format(round(PN * 100, 1), nsmall = 1))

reject_table = cbind(cases, reject) %>% select(-c(N, case, X))
reject_table

write.csv(reject_table, "recent_sim1.csv")
#write.csv(reject_table, "recent_sim1_large.csv")


# Study 2

set.seed(51483464)
B = 10000 # 1000

N <- 3000
n <- c(100, 200)
sigma <- c(0.1)
alpha <- c(0, 0.5, 1.0, 1.5)
cases <- expand_grid(N, n, sigma, alpha)

columns = c("case", "iteration", "DD", "PN", "HP", "PS1", "PS1q", "PS2", "PS2q",
            "PS3", "WF", "LR")
results = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(results) = columns

for (case in 1:nrow(cases)) {
  DD = PN = HP = PS1 = PS1q = PS2 = PS2q = PS3 = WF = LR = rep(NA, B)
  case_storage = data.frame(iteration = seq_len(B), DD, PN, HP, PS1, PS1q,
                            PS2, PS2q, PS3, WF, LR)
  for (b in 1:B) {
    pop = generate_data_study2(N = cases$N[case],
                               sigma = cases$sigma[case],
                               alpha = cases$alpha[case])
    samp = generate_sample_brewer(pop, w = pop$w, n = cases$n[case])

    case_storage$DD[b] = DD_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$PN[b] = PN_test(data = samp, y = samp$y, x = samp$x, wts = samp$w, est_split = 0.5)
    case_storage$HP[b] = HP_DC_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$PS1[b] = PS1_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$PS1q[b] = PS1q_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$PS2[b] = PS2_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$PS2q[b] = PS2q_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$PS3[b] = PS3_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$WF[b] = WF_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$LR[b] = LR_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
  }

  results = rbind(results, cbind(case, case_storage))
  print(case)
}
# write.csv(results, "st2_results.csv")
write.csv(results, "st2_results_large.csv")

#results = read.csv("st2_results_large.csv")
results = read.csv("st2_results.csv")

reject = results %>%
  mutate(HP = case_when(HP <= 0.05 ~ 1, TRUE ~ 0),
         DD = case_when(DD <= 0.05 ~ 1, TRUE ~ 0),
         PS1 = case_when(PS1 <= 0.05 ~ 1, TRUE ~ 0),
         PS1q = case_when(PS1q <= 0.05 ~ 1, TRUE ~ 0),
         PS2 = case_when(PS2 <= 0.05 ~ 1, TRUE ~ 0),
         PS2q = case_when(PS2q <= 0.05 ~ 1, TRUE ~ 0),
         PS3 = case_when(PS3 <= 0.05 ~ 1, TRUE ~ 0),
         WF = case_when(WF <= 0.05 ~ 1, TRUE ~ 0),
         LR = case_when(LR <= 0.05 ~ 1, TRUE ~ 0),
         PN = case_when(PN <= 0.05 ~ 1, TRUE ~ 0)) %>%
  select(-iteration) %>%
  group_by(case) %>%
  summarize(across(everything(), mean)) %>%
  mutate(HP = format(round(HP * 100, 1), nsmall = 1),
         DD = format(round(DD * 100, 1), nsmall = 1),
         PS1 = format(round(PS1 * 100, 1), nsmall = 1),
         PS1q = format(round(PS1q * 100, 1), nsmall = 1),
         PS2 = format(round(PS2 * 100, 1), nsmall = 1),
         PS2q = format(round(PS2q * 100, 1), nsmall = 1),
         PS3 = format(round(PS3 * 100, 1), nsmall = 1),
         WF = format(round(WF * 100, 1), nsmall = 1),
         LR = format(round(LR * 100, 1), nsmall = 1),
         PN = format(round(PN * 100, 1), nsmall = 1))

reject_table = cbind(cases, reject) %>% select(-c(N, case))
reject_table

#write.csv(reject_table, "recent_sim2.csv")
write.csv(reject_table, "recent_sim2_large.csv")


# Study 3

set.seed(51483464)
B = 10000 # 1000

N <- 3000
n <- c(100, 200)
psi <- c(0.0, 0.1, 0.2, 0.3)
alpha <- c(1.00, 0.75, 0.50, 0.25)
cases <- expand_grid(N, n, alpha, psi)

columns = c("case", "iteration", "DD", "PN", "HP", "PS1", "PS1q", "PS2", "PS2q",
            "PS3", "PS3q", "WF", "LR")
results = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(results) = columns

for (case in 1:nrow(cases)) {
  DD = PN = HP = PS1 = PS1q = PS2 = PS2q = PS3 = WF = LR = rep(NA, B)
  case_storage = data.frame(iteration = seq_len(B), DD, PN, HP, PS1, PS1q,
                            PS2, PS2q, PS3, WF, LR)
  for (b in 1:B) {
    pop = generate_data_study3(N = cases$N[case],
                               alpha = cases$alpha[case],
                               psi = cases$psi[case])
    samp = generate_sample_poisson(pop, w = pop$w, n = cases$n[case])

    case_storage$DD[b] = DD_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$PN[b] = PN_test(data = samp, y = samp$y, x = samp$x, wts = samp$w, est_split = 0.5)
    case_storage$HP[b] = HP_DC_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$PS1[b] = PS1_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$PS1q[b] = PS1q_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$PS2[b] = PS2_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$PS2q[b] = PS2q_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$PS3[b] = PS3_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$WF[b] = WF_WA_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
    case_storage$LR[b] = LR_test(data = samp, y = samp$y, x = samp$x, wts = samp$w)
  }

  results = rbind(results, cbind(case, case_storage))
  print(case)
}
#write.csv(results, "st3_results.csv")
write.csv(results, "st3_results_large.csv")

results = read.csv("st3_results_large.csv")
# results = read.csv("st3_results.csv")

reject = results %>%
  mutate(HP = case_when(HP <= 0.05 ~ 1, TRUE ~ 0),
         DD = case_when(DD <= 0.05 ~ 1, TRUE ~ 0),
         PS1 = case_when(PS1 <= 0.05 ~ 1, TRUE ~ 0),
         PS1q = case_when(PS1q <= 0.05 ~ 1, TRUE ~ 0),
         PS2 = case_when(PS2 <= 0.05 ~ 1, TRUE ~ 0),
         PS2q = case_when(PS2q <= 0.05 ~ 1, TRUE ~ 0),
         PS3 = case_when(PS3 <= 0.05 ~ 1, TRUE ~ 0),
         WF = case_when(WF <= 0.05 ~ 1, TRUE ~ 0),
         LR = case_when(LR <= 0.05 ~ 1, TRUE ~ 0),
         PN = case_when(PN <= 0.05 ~ 1, TRUE ~ 0)) %>%
  select(-iteration) %>%
  group_by(case) %>%
  summarize(across(everything(), mean)) %>%
  mutate(HP = format(round(HP * 100, 1), nsmall = 1),
         DD = format(round(DD * 100, 1), nsmall = 1),
         PS1 = format(round(PS1 * 100, 1), nsmall = 1),
         PS1q = format(round(PS1q * 100, 1), nsmall = 1),
         PS2 = format(round(PS2 * 100, 1), nsmall = 1),
         PS2q = format(round(PS2q * 100, 1), nsmall = 1),
         PS3 = format(round(PS3 * 100, 1), nsmall = 1),
         WF = format(round(WF * 100, 1), nsmall = 1),
         LR = format(round(LR * 100, 1), nsmall = 1),
         PN = format(round(PN * 100, 1), nsmall = 1))

reject_table = cbind(cases, reject) %>% select(-c(N, case))
reject_table

#write.csv(reject_table, "recent_sim3.csv")
write.csv(reject_table, "recent_sim3_large.csv")
