library(tidyverse)

# part of this is in base R bc that runs faster in sim studies
# part in dplyr tidyverse bc my brain works faster in that 
# will reconcile soon

# is there degree of heterogeneity in play in the dgm somewhere?

generate_data <- function(k, # schools
                          j, # teachers
                          i = 10, # students some default value
                          icc3, #school level icc
                          icc2, # teacher level icc
                          R2 = 0.40, # r-sq
                          ps_coef, # coefficients for ps model add default
                          pr_star, # overall proportion of teachers in trt group add default
                          outcome_coef, # coefficients for outcome model  add default
                          delta # treatment effect add default here
                          ){ 

  
  N <- i * j * k # total number of observations - 
  
  # covariates --------------------------------------------------------------
  
  # questions - should we do multivariate normal - 
  # the covariates aren't correlated with each other right now
  
  Z_k <- rnorm(k)  # mean 0 sd 1 by default
  W_jk <- rnorm(k * j)
  X_ijk <- rnorm(N)
  
  U_ijk <- rnorm(N)   # unobserved student-level covariate
  
  # Residuals ---------------------------------------------------------------
  
  r_k <- rnorm(k, 0, sqrt(icc3))
  u_jk <- rnorm(k * j, 0, sqrt(icc2))

  # data --------------------------------------------------------------------
  total_teachers <- k * j
  
  dat <- data.frame(student_id = 1:N,
                    teacher_id = rep(1:total_teachers, each = i),
                    school_id = rep(1:k, each = N/k))

  
  # covariates 
  dat$Z_k <- Z_k[dat$school_id]
  dat$W_jk <- W_jk[dat$teacher_id]
  dat$X_ijk <- X_ijk
  dat$U_ijk <- U_ijk
  
  # error terms
  dat$r_k <- r_k[dat$school_id]
  dat$u_jk <- u_jk[dat$teacher_id]

  # map(dat, ~ sum(is.na(.))) # check na 

  # design matrix -----------------------------------------------------------
  design_mat <- 
    dat %>%
    mutate(intercept = 1) %>%
    group_by(teacher_id) %>%
    mutate(X_jk = mean(X_ijk),
           U_jk = mean(U_ijk)) %>%
    ungroup()
  
  X_ps <- design_mat %>%
    select(intercept, X_jk, U_jk, W_jk, Z_k, r_k, u_jk) %>%
    distinct() %>%
    as.matrix()
  

 # Treatment indicator -----------------------------------------------------
  
  pr <- 1/ (1 + exp(-(X_ps %*% ps_coef)))  
  D <- ifelse(pr > pr_star, 1, 0) 
  #D <- rbinom(k * j, 1, pr)
  
  dat$D <- D[dat$teacher_id]
  

 # potential outcomes ------------------------------------------------------

 X_outcome <- 
    design_mat %>%
    select(intercept, X_ijk, U_ijk, W_jk, Z_k, r_k, u_jk) %>%
    as.matrix()
    
  Y_0_ijk <- rnorm(n = N,
                   mean = X_outcome %*% outcome_coef,
                   sd = sqrt(1 - icc2 - icc3 - R2))
  
  Y_1_ijk <- Y_0_ijk + delta
  

  # Observed outcome --------------------------------------------------------
  dat <- 
    dat %>%
    mutate(Y_ijk = D * Y_1_ijk + (1 - D) * Y_0_ijk)
  
  
  return(dat)
  
  
  
}


