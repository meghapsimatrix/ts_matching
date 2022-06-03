library(tidyverse)
library(MatchIt)
library(lme4)
library(cobalt)
library(broom)
library(matchMulti)
library(nlme)

# https://cran.r-project.org/web/packages/matchMulti/matchMulti.pdf

# load data ---------------------------------------------------------------

source("01_dgm/01_dgm_function.R")

set.seed(20220602)

k <- 100 #schools
j <-20  # teachers 
i <- 10 # students

icc3 <- 0.05  # icc for school 
icc2 <- 0.20 # icc for teachers
R2 <- .40


ps_coef <- matrix(c(.07, 0.8, -0.25, 0.6, -0.4, 1, 1))  #what values for pi 6 and 7 in the notes?
pr_star <- .5

outcome_coef <- matrix(c(1, 0.3, .5, .4, -0.2, 1, 1))

delta <- -0.4

example_dat <- generate_data(k = k,
                             j = j, 
                             i = i, 
                             icc3 = icc3, 
                             icc2 = icc2, 
                             R2 = R2,
                             ps_coef = ps_coef, 
                             pr_star = pr_star, 
                             outcome_coef = outcome_coef, 
                             delta = delta)

glimpse(example_dat)


# create cluster level data - cluster and site covariates


cluster_level_dat <- 
  example_dat %>%
  group_by(teacher_id) %>%
  summarize(school_id = mean(school_id),
            Z_k = mean(Z_k),
            W_jk = mean(W_jk),
            D = mean(D),
            Y_ijk = mean(Y_ijk)) %>%
  ungroup()


# quintiles based on site-level scores
# what should we do with these?
quint <- with(example_dat, quantile(Z_k, seq(0, 1, 0.2)))

example_dat$quintile <- cut(example_dat$Z_k, quint,
                            labels = c("A","B","C","D","E"),
                            include.lowest = TRUE)


glimpse(example_dat)


# Estimate propensity scores ----------------------------------------------

# random intercept only?
# what do we do with these propensity scores?

# unit level 
# ps in logit 
unit_ps_model <- glmer(D ~ X_ijk + W_jk + Z_k + 
                         (1 | school_id / teacher_id),
                       family = "binomial", 
                       data = example_dat)

example_dat$ps_unit <- predict(unit_ps_model, type = "link")
example_dat$ps_unit_pr <- predict(unit_ps_model, type = "response")

example_dat %>%
  mutate(D = as.character(D)) %>%
  ggplot(aes(x = ps_unit, fill = D)) +
  geom_density(alpha = 0.5) + 
  theme_minimal()


# check the predicted probabilites
check <- example_dat %>%
  select(ends_with("id"), D, ps_unit_pr) %>%
  mutate(ps_unit_pr = round(ps_unit_pr),
         match = D == ps_unit_pr)

table(check$match)
# perfect predictions for unit level bc 
# the dgm doesn't have unit level random error? for ps model

# cluster level 
cluster_ps_model <- glmer(D ~ W_jk + Z_k + 
                            (1 | school_id),
                         family = "binomial",
                         data = cluster_level_dat)

cluster_level_dat$ps_cluster <- predict(cluster_ps_model, type = "link")


cluster_level_dat %>%
  mutate(D = as.character(D)) %>%
  ggplot(aes(x = ps_cluster, fill = D)) +
  geom_density(alpha = 0.5) + 
  theme_minimal()

# Method 1 ----------------------------------------------------------------


# the following isn't working ps_unit is almost perfect prediction
# how do you do multilevel matchit ?
m_out_1 <- matchit(D ~ ps_unit, # is this right?
                   caliper = .25,
                   distance = example_dat$ps_unit,
                   data = example_dat)

# the following works match across all sites
m_out_1 <- matchit(D ~ X_ijk + W_jk + Z_k, 
                   caliper = .25,
                   distance = example_dat$ps_unit,
                   data = example_dat)


# exact match on site & group

m_out_2 <- matchit(D ~ X_ijk + W_jk + Z_k, 
                   caliper = .25,
                   exact = ~ school_id + quintile,
                   data = example_dat)


treatment_sites <- 
  example_dat %>%
  filter(D == 1)

# i think all sites are treatment sites
# within schools some teacher are trt and some control
# so we need to mess with data generation so some schools don't have a lot of treatment units or control 
sids <- sort(unique(treatment_sites$school_id))
length(sids)

