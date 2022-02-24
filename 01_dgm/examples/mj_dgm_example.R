library(tidyverse)
library(mvtnorm)
#library(simhelpers)


# Data Generating Model------------------------------------------------------

# n = sample size total
# mu = mean vector for covariate generation
# sig = var cov matrix for covariate generation
# ps_formula = formula regressing covs on the treatment indicator (D)
# b = regression coefficients for ps_formula
# y_formula = formula regressing D and covs on the treatment indicator
# a = regression coefficients for y_formula
# delta = treatment effect


generate_data <- function(n,
                          mu,
                          sig,
                          ps_formula,
                          b,
                          y_formula,
                          a,
                          delta, ...) {

  # Covariates --------------------------------------------------------------
  X_dat <- data.frame(rmvnorm(n = n, mu, sig))
  Y <- rnorm(n) # a placeholder to create the model matrices
  X_dat$Y <- Y

  X <- model.matrix(ps_formula, X_dat)

  # Data generation ---------------------------------------------------------
  pr <- 1/ (1 + exp(-(X %*% b)))
  D <- rbinom(n, 1, pr)

  X_all <- as.data.frame(cbind(X, D, Y))
  outcome_dat <- model.matrix(y_formula, X_all)
  outcome_coef <- matrix(c(a, delta))

  Y <- outcome_dat %*% outcome_coef + rnorm(n)  # add tau and what not here for multi

  X <- as_tibble(X)
  dat <- tibble(Y = Y, D = D) %>%
    bind_cols(X)

  return(dat)
}

# Example -----------------------------------------------------------------

sig <- matrix(c(1, .2, .4, .5, .6,
                .2, 1, .5, .4, .7,
                .4, .5, 1, .3, .6,
                .5, .4, .3, 1, .2,
                .6, .7, .6, .2, 1), ncol = 5)

mu <- c(0,0,0,0,0)

b <- matrix(c(1.2, 0.8, -0.25, 0.6, -0.4, -0.8))
a <- matrix(c(.5, 0.3, -0.36, -0.73, -0.2, 0.71))
delta <- -0.4

ps_formula <- as.formula("Y ~ X1 + X2 + X3 + X4 + X5")
y_formula <- as.formula("Y ~ X1 + X2 + X3 + X4 + X5 + D")

set.seed(1222018)
example_dat <- generate_data(n = 1000, mu = mu, sig = sig, ps_formula, b, y_formula, a, delta)

table(example_dat$D)
summary(example_dat$Y)
# the intercept for outcome model should make sure Y
# is balanced enough so that we won't run into estimation issues


