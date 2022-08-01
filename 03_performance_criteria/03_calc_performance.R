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
  
  
  ci_cov <- 
    results %>%
    group_by(method) %>%
    do(calc_coverage(., lower_bound = ci_low, upper_bound = ci_high, true_param = true_effect)) %>%
    ungroup() %>%
    select(-K)
  
  mean_smd <- 
    results %>%
    group_by(method) %>%
    summarize(across(U_ijk:Z_k, ~  mean(.x, na.rm = TRUE))) %>%
    ungroup()
  
  mean_prop <- 
    results %>%
    group_by(method) %>%
    summarize(prop_t_m = mean(prop_t, na.rm = TRUE)) %>%
    ungroup
  
  performance_measures <- left_join(abs_criteria, ci_cov, by = "method") %>%
    left_join(mean_smd, by = "method") %>%
    left_join(mean_prop, by = "method") %>%
    mutate(method = as.numeric(method)) %>%
    arrange(method)


  
  return(performance_measures)
}
