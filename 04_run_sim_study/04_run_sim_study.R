library(dplyr)
library(purrr)
library(mvtnorm)
library(tidyr)
library(stringr)
library(tibble)
library(simhelpers)
library(broom.mixed)
library(lme4)
library(lmerTest)
library(janitor)
library(parameters)


#-----------------------------------------------------------
# Source the functions
#-----------------------------------------------------------

source("01_dgm/01_dgm_function.R")
source("02_estimation_methods/02_estimation_functions.R") 
source("03_performance_criteria/03_calc_performance.R")


#-----------------------------------------------------------
# Simulation Driver - should return a data.frame or tibble
#-----------------------------------------------------------

# check if default parameters make sense

run_sim <- function(iterations, 
                    k,
                    j, 
                    i = 20, 
                    icc3, 
                    icc2, 
                    R2 = .40,
                    ps_coef = matrix(c(-1.25, 1, .25, 1.5, 1.5, .20, .20)), 
                    pr_star = .5, 
                    outcome_coef = matrix(c(1, 0.3, .5, .4, -0.2, 1, 1)), 
                    delta = -0.4,
                    seed = NULL) {
  
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
                  X_jk = mean(X_ijk),
                  Z_k = mean(Z_k),
                  W_jk = mean(W_jk),
                  Z_q5 = as.factor(as.character(mean(as.numeric(Z_q5)))),
                  W_q5 = as.factor(as.character(mean(as.numeric(W_q5)))),
                  D = mean(D),
                  Y_ijk = mean(Y_ijk)) %>%
        ungroup()
     
     # restrict to sites that have t and c clusters
     dat_m <- dat %>%
       group_by(school_id) %>%
       mutate(tmpm = mean(D)) %>%
       ungroup() %>%
       filter(tmpm > 0 & tmpm < 1) %>%
       select(-tmpm)

     # match -------------------------------------------------------------------

     m_1 <- match_them(dat = dat, 
                       equation =  D ~ X_ijk + W_jk + Z_k,
                       caliper = .25)
     
     m_2 <- match_them(dat = cluster_level_dat, 
                       equation =  D ~ X_jk + W_jk + Z_k,
                       caliper = .25,
                       student_level = TRUE,
                       student_dat = dat)
     
     m_3 <- multi_match(dat = dat,
                        trt = "D",
                        l1_cov = c("X_ijk"),
                        l2_cov = c("W_q5", "Z_q5"),
                        l2_id = "teacher_id",
                        caliper = .25)
     
     # do we need to make school_id a factor or char?
     m_4 <- match_them(dat = dat, 
                       equation =  D ~ X_ijk + W_jk + Z_k,
                       caliper = 1,
                       exact = "school_id")
     
     m_5 <- match_them(dat = cluster_level_dat, 
                       equation =  D ~ X_jk + W_jk + Z_k,
                       caliper = 1,
                       exact = "school_id",
                       student_level = TRUE,
                       student_dat = dat)
     
     m_6 <- dat_m %>%
              group_by(school_id) %>%
              do(multi_match(., 
                             trt = "D",
                             l1_cov = c("X_ijk"),
                             l2_cov = NULL,
                             l2_id = "teacher_id",
                             l3_id = "school_id",
                             caliper = 4,
                             add_id = "school"))
     
     m_7 <- match_them(dat = dat, 
                       equation =  D ~ X_ijk + W_jk + Z_k,
                       caliper = 1,
                       exact = "Z_q5")  # quintile based on Z_k
     
     
     m_8 <- match_them(dat = cluster_level_dat, 
                       equation =  D ~ X_jk + W_jk + Z_k,
                       caliper = 1,
                       exact = "Z_q5",  # quintile based on Z_k
                       student_level = TRUE,
                       student_dat = dat)  
     
    
     m_9 <- dat_m %>%
              group_by(Z_q5) %>%
              do(multi_match(., 
                             trt = "D",
                             l1_cov = c("X_ijk"),
                             l2_cov = c("W_q5", "Z_q5"),
                             l2_id = "teacher_id",
                             l3_id = "Z_q5",
                             caliper = 4,
                             add_id = "school"))
     

    m_10 <- match_hybrid(dat = dat, 
                         site_id = "school_id", 
                         teacher_id = "student_id", 
                         tx_var = "D", 
                         exact_vars = NULL, 
                         teacher_vars = c("X_ijk","W_jk"), 
                         site_vars = "Z_k",
                         group = "Z_q5", # changed the name here bc Jordan is using 2 different quintiles
                         crc = caliper,
                         replacement = FALSE,
                         ratio = 1,
                         seed = 1234,
                         by_student = TRUE,
                         student_dat = dat) # may need to change seed
    
    
    
    m_11 <- match_hybrid(dat = cluster_level_dat, 
                         site_id = "school_id", 
                         teacher_id = "teacher_id", 
                         tx_var = "D", 
                         exact_vars = NULL, 
                         teacher_vars = c("X_jk","W_jk"), 
                         site_vars = "Z_k",
                         group = "Z_q5",
                         crc = caliper,
                         replacement = FALSE,
                         ratio = 1,
                         seed = 1234,
                         by_student = FALSE,
                         student_dat = dat) 
    
    m_12_hold <- 
      dat_m %>%
      group_by(school_id) %>%
      do(multi_match(., 
                     trt = "D",
                     l1_cov = c("X_ijk"),
                     l2_cov = NULL,
                     l2_id = "teacher_id",
                     l3_id = "Z_q5",
                     caliper = 4,
                     add_id = "school")) %>%
      as.data.frame()
    
    
    # Identify Clusters NOT in a Within-Site Pair
    wsp <- m_12_hold %>% 
      ungroup() %>%
      select(teacher_id) %>% 
      mutate(inwsp = 1)
    
    
    dat_m2 <- dat %>% 
      left_join(wsp, by = "teacher_id") %>% 
      filter(is.na(inwsp) == TRUE)
    
    # restrict matching to groups that have T & C clusters
    
    dat_m2 <- dat_m2 %>% 
      group_by(Z_q5) %>% 
      mutate(tmpm = mean(D)) %>% 
      ungroup() %>%
      filter(tmpm > 0 & tmpm < 1) %>% 
      select(-tmpm)     
    
    
    m_12_hold_2 <- dat_m2 %>%
      mutate(Z_q5 = as.character(Z_q5)) %>%
      group_by(Z_q5) %>%
      do(multi_match(., 
                     trt = "D",
                     l1_cov = c("X_ijk"),
                     l2_cov = NULL,
                     l2_id = "teacher_id",
                     l3_id = "Z_q5",
                     caliper = 4,
                     add_id = "pair")) %>%
      as.data.frame() %>%
      mutate(Z_q5 = as.factor(Z_q5))
    
    
     m_12 <- bind_rows(m_12_hold, m_12_hold_2)
    
    

     # outcome models ----------------------------------------------------------

     # something that runs the hlm outcome model and saves estimates etc. 
     # and the balance stats for each method
     # then save as long format with method 1 - results then method 2 - results ...
    
     matched_sets <- tibble(matched_dat = list(m_1, m_2, m_3,
                                               m_4, m_5, m_6,
                                               m_7, m_8, m_9,
                                               m_10, m_11, m_12),
                            method = c("1", "2", "3", "4", "5", "6",
                                       "7", "8", "9", "10", "11", "12"),
                            n_t_all = sum(dat$D))
     
    pmap_dfr(matched_sets, estimate_effect) %>%
      mutate(true_effect = delta)
    
       
    
    }) %>%
    bind_rows()
  
  calc_performance(results)
}

#-------------------------------------
# Experimental Design
#-------------------------------------

# include design matrix, exclude to_test

set.seed(20220805) # change this seed value!

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
    iterations = 2, # change this to how many ever iterations
    seed = round(runif(1) * 2^30) + 1:n()
  )

# this is just to test
params <- params[1, ]


# run sim in serial -------------------------------------------------------


library(purrr)

system.time(
  results <-
    params %>%
    mutate(res = pmap(., .f = run_sim)) %>%
    unnest(cols = res)
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
