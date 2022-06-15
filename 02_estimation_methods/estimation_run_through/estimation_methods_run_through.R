library(tidyverse)
library(MatchIt)
library(lme4)
library(cobalt)
library(broom)
library(matchMulti)
library(nlme)
library(optmatch)

# https://cran.r-project.org/web/packages/matchMulti/matchMulti.pdf

# load data ---------------------------------------------------------------

source("01_dgm/01_dgm_function.R")

set.seed(20220602)

k <- 50 #schools
j <- 8  # teachers 
i <- 20 # students

icc3 <- 0.05  # icc for school 
icc2 <- 0.20 # icc for teachers
R2 <- .40

# 50 schools 6 teachers per school


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

# do we do hlm for ps model or not?


# unit --------------------------------------------------------------------

# Not multilevel
unit_ps_model_nm <- glm(D ~ X_ijk + W_jk + Z_k, 
                         family = "binomial", 
                         data = example_dat)

example_dat$ps_unit_nm <- predict(unit_ps_model_nm, type = "link")


# Multilevel 
unit_ps_model <- glmer(D ~ X_ijk + W_jk + Z_k + 
                         (1 | school_id),
                       family = "binomial", 
                       data = example_dat)

example_dat$ps_unit <- predict(unit_ps_model, type = "link")
example_dat$ps_unit_pr <- predict(unit_ps_model, type = "response")

# check common support
example_dat %>%
  mutate(D = as.character(D)) %>%
  ggplot(aes(x = ps_unit, fill = D)) +
  geom_density(alpha = 0.5) + 
  theme_minimal()



# cluster -----------------------------------------------------------------

# not multilevel
cluster_ps_model_nm <- glm(D ~ W_jk + Z_k,
                          family = "binomial",
                          data = cluster_level_dat)

cluster_level_dat$ps_cluster_nm <- predict(cluster_ps_model_nm, type = "link")


# multilevel
cluster_ps_model <- glmer(D ~ W_jk + Z_k + 
                            (1 | school_id),
                         family = "binomial",
                         data = cluster_level_dat)

cluster_level_dat$ps_cluster <- predict(cluster_ps_model, type = "link")

#check common support
cluster_level_dat %>%
  mutate(D = as.character(D)) %>%
  ggplot(aes(x = ps_cluster, fill = D)) +
  geom_density(alpha = 0.5) + 
  theme_minimal()

# Method 1 ----------------------------------------------------------------

#optimal method in matchit doesn't like exact or caliper

match_them <- function(dat, 
                       #equation_ps = "D ~ ps_unit",
                       ps_method = "nearest",
                       exact = NULL, 
                       ps = dat$ps_unit){
  

  m_out <- matchit(D ~ ps_unit, 
                   #caliper = .25,
                   method = ps_method,
                   #exact = exact,
                   distance = ps,
                   data = dat)
  
  match_dat <- match.data(m_out)
  
  return(match_dat)
  
}


# switch distance to appropriate ps (multi or not multilevel)
system.time(m_out_1 <- match_them(dat = example_dat, 
                                  ps = example_dat$ps_unit))

# exact match on site by looping over site

system.time(site_match_data <- 
  example_dat %>%
  group_by(school_id) %>%
  group_modify(~ match_them(.x)))
  


example_dat %>%
  group_by(school_id, D) %>%
  summarize(n_teach = n_distinct(teacher_id)) %>%
  View()
  

# are schools allowed to have 
# only treatment teachers or only control teachers


treatment_sites <- 
  example_dat %>%
  filter(D == 1)

# i think all sites are treatment sites
# within schools some teacher are trt and some control
# so we need to mess with data generation so some schools don't have a lot of treatment units or control 
sids <- sort(unique(treatment_sites$school_id))
length(sids)

