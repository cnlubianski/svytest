# testing

N <- 3000
n <- c(100, 200)
sigma <- c(0.1, 0.2)
alpha <- c(0.0, 0.2, 0.4, 0.6)
delta <- c(1.5, 1)
cases <- expand_grid(N, n, sigma, delta, alpha)
case <- 32

pop <- generate_data_quant(N = cases$N[case],
                           sigma = cases$sigma[case],
                           alpha = cases$alpha[case],
                           delta = cases$delta[case])
samp <- generate_sample_brewer(pop, w = pop$w, n = cases$n[case])

















lm_object <- lm(y ~ x, data = samp, w = w)
weights(lm_object)

design <- svydesign(id = ~1, weights = samp$w, data = samp)
svyglm_object <- svyglm(y ~ x, design = design)






hi <- diff_in_coef(model = svyglm_object, var_equal = TRUE, robust_type = "HC1")
print(hi)
hi$p.value


bye <- wa_test(model = svyglm_object, type = "DD")
print(bye)
summary(bye)


die <- lr_test(model = svyglm_object, likelihood = "scaled")
print(die)
summary(die)

data("svytestCE")


svytestCE$AGE = as.factor(svytestCE$AGE)
des <- svydesign(ids = ~1, weights = ~FINLWT21, data = svytestCE)
fit_svy <- svyglm(TOTEXPCQ ~ ROOMSQ + BATHRMQ + BEDROOMQ + FAM_SIZE + AGE, design = des)
res <- diff_in_coef(model = fit_svy, var_equal = FALSE, robust = "HC0")
print(res)

lm_object <- lm(TOTEXPCQ ~ ROOMSQ + BATHRMQ + BEDROOMQ + FAM_SIZE + AGE,
                data = svytestCE)
summary(lm_object)

glance(res)


hop <- lr_test(model = fit_svy, likelihood = "scaled")
print(hop)


samp1 <- samp
samp1$w <- samp1$y + 1

design1 <- svydesign(id = ~1, weights = samp$w, data = samp1)
svyglm_object1 <- svyglm(y ~ x, design = design1)

die <- lr_test(model = svyglm_object1, likelihood = "scaled")
print(die)
summary(die)



svytestCE$FAM_SIZE <- as.factor(svytestCE$FAM_SIZE)

des <- svydesign(ids = ~1, weights = ~FINLWT21, data = svytestCE)
fit_svy <- svyglm(TOTEXPCQ ~ ROOMSQ + BATHRMQ + BEDROOMQ + FAM_SIZE + AGE, design = des)
res <- diff_in_coef(model = fit_svy, var_equal = FALSE)
print(res)

lm_object <- lm(TOTEXPCQ ~ ROOMSQ + BATHRMQ + BEDROOMQ + FAM_SIZE + AGE,
                data = svytestCE)
summary(lm_object)








library(survey)

#--- Simulation setup ---
set.seed(123)

n <- 200
x <- rnorm(n)
# True model: y depends on x
y <- 2 + 1.5 * x + rnorm(n, sd = 1)

# Make weights informative: larger weights for larger y
wts <- 1 + (y - min(y))   # strictly positive, correlated with y

# Put into a data frame
dat <- data.frame(y = y, x = x, w = wts)

#--- Survey design and model ---
des <- svydesign(ids = ~1, weights = ~w, data = dat)
fit_svy <- svyglm(y ~ x, design = des)

#--- Run your LR test ---
res <- lr_test(fit_svy)  # default likelihood = "pseudo"
print(res)
summary(res)
