library(simhelpers)
library(tidyverse)

#------------------------------------------------------
# Calculate performance measures
# (For some simulations, it may make more sense
# to do this as part of the simulation driver.)
#------------------------------------------------------

# want to calcualte rmse, bias, balance...
# use simhelpers to get mcse

calc_performance <- function(results) {
  
  abs_criteria <- 
    results %>%
    group_by(method) %>%
    do(calc_absolute(., estimates = estimate, true_param = true_effect,
                     perfm_criteria = c("bias", "rmse"))) %>%
    ungroup()
  
  mean_smd <- 
    results %>%
    group_by(method) %>%
    summarize(across(U_ijk:Z_k, ~  mean(.x, na.rm = TRUE))) %>%
    ungroup()
  
  performance_measures <- left_join(abs_criteria, mean_smd, by = "method") %>%
    mutate(method = as.numeric(method)) %>%
    arrange(method)

  
  return(performance_measures)
}
