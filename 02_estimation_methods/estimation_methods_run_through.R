library(tidyverse)
library(MatchIt)
library(lme4)
library(cobalt)
library(broom)
library(matchMulti)
library(nlme)

# https://cran.r-project.org/web/packages/matchMulti/matchMulti.pdf


# load data ---------------------------------------------------------------

load("data/example_data.RData")

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

quint <- with(example_dat, quantile(Z_k, seq(0, 1, 0.2)))

example_dat$quintile <- cut(example_dat$Z_k, quint,
                            labels = c("A","B","C","D","E"),
                            include.lowest = TRUE)


glimpse(example_dat)


# Estimate propensity scores ----------------------------------------------

# random intercept only?
# what do we do with these propensity scores?

# unit level 
unit_ps_model <- glmer(D ~ X_ijk + W_jk + Z_k + (1 | teacher_id) + (1 | school_id),
                       family = "binomial", 
                       data = example_dat)

example_dat$ps_unit <- predict(unit_ps_model, type = "response")


# cluster level 
cluster_ps_model <- glmer(D ~  W_jk + Z_k + (1 | school_id),
                         family = "binomial",
                         data = cluster_level_dat)

cluster_level_dat$ps_cluster <- predict(cluster_ps_model, type = "response")



# Method 1 ----------------------------------------------------------------

# no distinction units only 
# by default nearest neighbor, no replacement

# just with matchit (not multi-level)
m_out_1 <- matchit(D ~ X_ijk + W_jk + Z_k,
                   caliper = .25,
                   data = example_dat)

m_out_1$match.matrix


balance_stats <- bal.tab(m_out_1)
balance_stats$Balance

m_data_1 <- match.data(m_out_1)
table(m_data_1$weights)  # if we get to a point were we have to use weights - lmer won't work

out_mod_1 <- lmer(Y_ijk ~ D  + X_ijk + W_jk + Z_k + (1 | teacher_id) + (1 | school_id),
                  data = m_data_1)

summary(out_mod_1)



# do I use the ps estimated using glmer to match?


# matchMulti --------------------------------------------------------------

match_simple <- matchMulti(example_dat, 
                           treatment = "D",
                           school.id = "teacher_id",
                           match.students = FALSE,
                           student.vars = "X_ijk",
                           verbose = TRUE)

# two level - level 2 is teachers where the treatment is 

#str(match_simple)

match_data <- as.data.frame(match_simple$matched)
head(match_data)

mod <- lme(Y_ijk ~ D, 
           random = ~ 1 | pair.id/ school_id,
           data = match_data)

summary(mod)
