library(dplyr)
library(purrr)
library(mvtnorm)
library(tidyr)
library(stringr)
library(tibble)
library(simhelpers)


#-----------------------------------------------------------
# Source the functions
#-----------------------------------------------------------

source("01_dgm/01_dgm_function.R")
source("02_estimation_methods/02_estimation_functions.R") 
source("03_performance_criteria/03_calc_performance.R")


#-----------------------------------------------------------
# Simulation Driver - should return a data.frame or tibble
#-----------------------------------------------------------

# need to fill this in 

run_sim <- function(iterations, model_params, design_params, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  results <-
    rerun(iterations, {
      

     # generate student level data ---------------------------------------------
      dat <- generate_data(k = k,
                           j = j, 
                           i = i, 
                           icc3 = icc3, 
                           icc2 = icc2, 
                           R2 = R2,
                           ps_coef = ps_coef, 
                           pr_star = pr_star, 
                           outcome_coef = outcome_coef, 
                           delta = delta)
      

     # aggregate to teacher level ----------------------------------------------
     cluster_level_dat <- 
        dat %>%
        group_by(teacher_id) %>%
        summarize(school_id = mean(school_id),
                  Z_k = mean(Z_k),
                  W_jk = mean(W_jk),
                  D = mean(D),
                  Y_ijk = mean(Y_ijk)) %>%
        ungroup()
     

     # match -------------------------------------------------------------------

     m_1 <- match_them(dat = dat, ps = dat$ps_unit)
     # m_2 <-
     # m_3 <- 
     

     # outcome models ----------------------------------------------------------

     # something that runs the hlm outcome model and saves estimates etc. 
     # and the balance stats for each method
     # then save as long format with method 1 - results then method 2 - results ...
     
     # o_1 <- 
     # o_2 <-
       
    
    }) %>%
    bind_rows()
  
  calc_performance(results, model_params)
}

#-------------------------------------
# Experimental Design
#-------------------------------------

# include design matrix, exclude to_test

set.seed(20220510) # change this seed value!

# now express the simulation parameters as vectors/lists

design_factors <- list(
  k = c(20, 50), # schools
  j = c(4, 8), # teachers
  icc3 = c(0.10, 0.20), #icc school level
  icc2 = c(0.10, 0.20) #icc teacher level
)

# combine into a design set
params <-
  cross_df(design_factors) %>%
  mutate(
    iterations = 2400, # change this to how many ever iterations
    seed = round(runif(1) * 2^30) + 1:n()
  )


#--------------------------------------------------------
# run simulations in parallel - future + furrr workflow
#--------------------------------------------------------

library(future)
library(furrr)

plan(multisession) # choose an appropriate plan from the future package
evaluate_by_row(params, run_sim)

# OR
plan(multisession)
system.time(
  results <-
    params %>%
    mutate(res = future_pmap(., .f = run_sim)) %>%
    unnest(cols = res)
)


#--------------------------------------------------------
# Save results and details
#--------------------------------------------------------

session_info <- sessionInfo()
run_date <- date()

save(params, results, session_info, run_date, file = "simulation_results.Rdata")