README
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

svytest provides a unified suite of diagnostic tests for assessing the
informativeness of survey weights in regression models. It implements
both classical parametric tests (Hausman–Pfeffermann, DuMouchel–Duncan,
Pfeffermann–Sverchkov, Wu–Fuller, Estimating Equations) and novel
non‑parametric permutation tests.

# Why use `svytest`?

Survey weights are essential for unbiased descriptive statistics, but
their role in regression modeling is debated. Using weights
unnecessarily can inflate variance, while ignoring informative weights
can bias estimates.

svytest helps you decide whether weights matter in your regression by
offering:

- Difference‑in‑Coefficients tests (Hausman–Pfeffermann)
- Weight‑Association tests (DD, PS1/PS2, WF)
- Estimating Equations test (Pfeffermann–Sverchkov)
- Permutation tests (distribution‑free, robust alternatives)
- Simulation utilities to benchmark test performance

# Installation

You can install the released version of svytest from CRAN with:
`install.packages("svytest")` or install the development version from
GitHub with: `remotes::ionstall_github("cnlubianski/svytest")`.

# Usage

Load the package and your survey data.

``` r
library(svytest)
library(survey)

data("svytestCE", package = "svytest")
des <- svydesign(ids = ~1, weights = ~FINLWT21, data = svytestCE)
fit <- svyglm(TOTEXPCQ ~ ROOMSQ + BATHRMQ + BEDROOMQ + FAM_SIZE + AGE,
              design = des)
```

## Example 1: Difference-in-Coefficients Test

``` r
res_dc <- diff_in_coef_test(fit, var_equal = TRUE)
summary(res_dc)
#> 
#> Difference-in-Coefficients Test
#> Call:
#> diff_in_coef_test(model = fit, var_equal = TRUE)
#> 
#> Method:
#>  Hausman-Pfeffermann Difference-in-Coefficients Test
#> 
#> Test Statistic:
#>  Chi-sq = 15.3278  on 6 df , p-value = 0.0179
```

This compares weighted vs. unweighted regression coefficients. A small
p‑value suggests weights are informative.

## Example 2: Weight-Association Test (DuMouchel–Duncan)

``` r
res_dd <- wa_test(fit, type = "DD")
summary(res_dd)
#> 
#> Weight-Association Test
#> Call:
#> wa_test(model = fit, type = "DD")
#> 
#> Method:
#>  DuMouchel-Duncan Weight-Association Test
#> 
#> Test Statistic:
#>  F = 2.5557  on 6 and 21914 df , p-value = 0.0178
```

This regresses residuals on weights. A significant slope indicates
informativeness.

## Example 3: Estimating Equations Test

``` r
res_perm <- perm_test(fit, stat = "pred_mean", B = 500, engine = "R")
summary(res_perm)
```

# Getting Started

See the vignette for a comprehensive review of all tests, the
statistical underpinnings, and references: `vignette("svytest")`.
