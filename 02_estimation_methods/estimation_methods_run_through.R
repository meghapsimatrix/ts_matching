library(tidyverse)
library(MatchIt)
library(lme4)
library(cobalt)
library(broom)

load("data/example_data.RData")

glimpse(example_dat)


# Cluster level data ------------------------------------------------------

cluster_level <- 
  example_dat %>%
  group_by(teacher_id) %>%
  summarize(school_id = mean(school_id),
            Z_k = mean(Z_k),
            W_jk = mean(W_jk),
            D = mean(D),
            Y_ijk = mean(Y_ijk)) %>%
  ungroup()



# Method 1 ----------------------------------------------------------------

# no distinction units only 

# by default nearest neighbor, no replacement

# just with matchit
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
