library(tidyverse)
library(cobalt)


bal_dat <- fun.balance(dat_full = dat,
                       dat_match = m_1,
                       tx.var = "D",
                       vars = c("Z_k","r_k","W_jk","u_jk","X_jk","X_ijk","U_ijk"))

covs <- m_1 %>%
  select(Z_k, r_k, W_jk, u_jk, X_jk, X_ijk, U_ijk)

ps_eq <- f.build("D", covs)
  


bal_stats <-
  bal.tab(D ~ Z_k + r_k + W_jk + u_jk + X_jk + X_ijk + U_ijk,
          data = m_1,
          method = "weighting",
          weights = "weights",
          s.d.denom = "treated")

bal_dat <- 
  bal_stats$Balance %>%
  rownames_to_column("var") %>%
  select(var, smd = Diff.Adj)


system.time(bal_m1 <- calc_balance(m_1))
system.time(bal_m1_c <- calc_bal(m_1))
bal_m10 <- calc_balance(m_10)


dups_m10 <- m_10 %>%
  filter(duplicated())