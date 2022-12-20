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
  total_teachers <- k * j # total number of clusters
  
  # covariates --------------------------------------------------------------
  
  # questions - should we do multivariate normal - 
  # the covariates aren't correlated with each other right now
  
  # generate site-level data
  
  Z_k <- rnorm(k)  # mean 0 sd 1 by default
  G_k <- rnorm(k)
  H_k <- rnorm(k)
  r_k <- rnorm(k, 0, sqrt(icc3))
  
  
  
  # generate cluster & unit level data for each site
  
  hold <- NULL
  
  for(k. in 1:k) {
    
    # cluster level 
    school_id <- rep(k., i * j)
    W_jk <- rep(rnorm(j, -0.2 * Z_k[k.], 1), i) # observed cluster-level covariate
    C_jk <- rep(rnorm(j, -0.4 * Z_k[k.], 1), i) # do we need to change the #'s here
    L_jk <- rep(rnorm(j, -0.5 * Z_k[k.], 1), i)
    u_jk <- rep(rnorm(j, 0, sqrt(icc2)), i) # cluster-level residual
    
    # student level
    X_ijk <- rnorm(i * j, -0.45 * Z_k[k.] + -0.20 * r_k[k.], 1) # observed unit-level covariate
    A_ijk <- rnorm(i * j, -0.25 * Z_k[k.] + -0.40 * r_k[k.], 1) # observed unit-level covariate
    B_ijk <- rnorm(i * j, -0.35 * Z_k[k.] + -0.50 * r_k[k.], 1) # observed unit-level covariate
    U_ijk <- rnorm(i * j, -0.25 * Z_k[k.], 1)  # unobserved unit-level covariate 
    
    
    tmp <- cbind(school_id, 
                 rep(Z_k[k.], i * j), 
                 rep(G_k[k.], i * j), 
                 rep(H_k[k.], i * j), 
                 rep(r_k[k.], i * j), 
                 W_jk,C_jk, L_jk, u_jk, 
                 X_ijk, A_ijk, B_ijk, U_ijk)
    tmp <- tmp[order(W_jk, u_jk),]
    hold <- rbind(hold, tmp)
    
    
    
  }
  
  dat <- as.data.frame(hold)
  names(dat) <- c("school_id", "Z_k", "G_k", "H_k", "r_k", 
                  "W_jk", "C_jk", "L_jk", "u_jk",
                  "X_ijk", "A_ijk", "B_ijk", "U_ijk")
  dat$teacher_id <- rep(1:(total_teachers), each = i)
  dat$student_id <- 1:N
  dat$school_id <- rep(1:k, each = N/k)


  # design matrix -----------------------------------------------------------
  design_mat <- 
    dat %>%
    mutate(intercept = 1) %>%
    group_by(teacher_id) %>%
    mutate(X_jk = mean(X_ijk),
           A_jk = mean(A_ijk),
           B_jk = mean(B_ijk),
           U_jk = mean(U_ijk),
           X_jk_2 = mean(X_ijk^2),
           A_jk_2 = mean(A_ijk^2),
           X_A = mean(X_ijk * A_ijk),
           W_X = mean(W_jk * X_ijk),
           C_W = C_jk * W_jk,
           L_jk_2 = L_jk ^ 2,
           Z_X = mean(Z_k * X_ijk)) %>%
    ungroup()
  
  X_ps <- design_mat %>%
    select(intercept, 
           X_jk, A_jk, B_jk, U_jk,
           X_jk_2, A_jk_2, X_A,
           W_jk, C_jk, L_jk,
           L_jk_2, W_X, C_W,
           Z_k, G_k, H_k,
           Z_X,
           r_k, u_jk) %>%
    distinct() %>%
    as.matrix()
  

 # Treatment indicator -----------------------------------------------------
  
  pr <- 1/ (1 + exp(-(X_ps %*% ps_coef)))  
  D <- rbinom(k * j, 1, pr)
  
  dat$D <- D[dat$teacher_id]
  

 # potential outcomes ------------------------------------------------------

  X_outcome <- 
    design_mat %>%
    select(intercept, X_ijk, U_ijk, W_jk, Z_k, r_k, u_jk) %>%
    mutate(X_ijk_2 = X_ijk ^ 2, 
           X_ijk_3 = X_ijk ^ 3,
           W_jk_2 = W_jk ^ 2,
           W_jk_3 = W_jk ^3,
           X_W = X_ijk * W_jk,
           X_A = X_ijk * A_ijk) %>%
    select(intercept, X_ijk, X_ijk_2, X_ijk_3, U_ijk, 
           W_jk, W_jk_2, W_jk_3, Z_k, 
           X_W, X_A,
           r_k, u_jk) %>%
    as.matrix()
    
  Y_0_ijk <- rnorm(n = N,
                   mean = X_outcome %*% outcome_coef,
                   sd = sqrt(1 - icc2 - icc3 - R2))
  
  # standardizing Y_0_ijk
  Y_0_ijk <- (Y_0_ijk - mean(Y_0_ijk)) / sd(Y_0_ijk)
  
  Y_1_ijk <- Y_0_ijk + delta
  

  # Observed outcome --------------------------------------------------------
  dat <- 
    dat %>%
    mutate(Y_ijk = D * Y_1_ijk + (1 - D) * Y_0_ijk)
  
  
  dat <- 
    dat %>%
    group_by(teacher_id) %>%
    mutate(X_jk = mean(X_ijk)) %>%
    ungroup()
  
  
  dat$W_q5 <- cut(dat$W_jk, 5,
                  labels = c(1, 2, 3, 4, 5),
                  include.lowest = TRUE) 
  
  dat$Z_q5 <- cut(dat$Z_k, 5,
                  labels = c(1, 2, 3, 4, 5),
                  include.lowest = TRUE)
  
  dat <- as.data.frame(dat)
  
  
  return(dat)
  
  
  
}


