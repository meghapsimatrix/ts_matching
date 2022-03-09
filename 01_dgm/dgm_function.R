

generate_data <- function(k, # schools
                          j, # teachers
                          i, # students
                          icc3, #school level icc
                          icc2, # teacher level icc
                          ps_coef, # coefficients for ps model
                          pr_star, # overall proportion of teachers in trt group
                          outcome_coef # coefficients for outcome model 
                          ){ 

  
  N <- i * j * k # total numbr of observations - 
  
  # covariates --------------------------------------------------------------
  
  Z_k <-rnorm(k)  # mean 0 sd 1 by default
  W_jk <- rnorm(k *j)
  X_ijk <- rnorm(N)
  
  U_ijk <- rnorm(N)   # unobserved student-level covariate
  
  # Residuals ---------------------------------------------------------------
  
  u_k <- rnorm(k, 0, icc3)
  v_jk <- rnorm(k * j, 0, icc2)
  e_ijk <- rnorm(N, 1)

  # data --------------------------------------------------------------------
  
  dat <- data.frame(student_id = 1:N,
                    teacher_id = rep(rep(1:j, times = k), each = i),
                    school_id = rep(1:k, each = N/k))
  
  dat$teacher_id_spec <- as.numeric(paste0(dat$school_id, dat$teacher_id))
  
  
  # covariates 
  dat$Z_k <- Z_k[dat$school_id]
  dat$W_jk <- W_jk[dat$teacher_id_spec]
  dat$X_ijk <- X_ijk
  dat$U_ijk <- U_ijk
  
  # error terms
  dat$u_k <- u_k[dat$school_id]
  dat$v_jk <- v_jk[dat$teacher_id_spec]
  dat$e_ijk <- e_ijk
  

  # design matrix -----------------------------------------------------------
  design_mat_ps <- 
    dat %>%
    mutate(intercept = 1) %>%
    group_by(teacher_id_spec) %>%
    mutate(X_jk = mean(X_ijk),
           U_jk = mean(U_ijk)) %>%
    ungroup() %>%
    select(intercept, X_jk, U_jk, W_jk, Z_k, u_k, v_jk) %>%
    distinct()
    
  X_ps <- as.matrix(design_mat_ps)
  

 # Treatment indicator -----------------------------------------------------
  
  pr <- 1/ (1 + exp(-(X_ps %*% b)))  
  D <- ifelse(pr > pr_star, 1, 0) 
  
  dat$D <- D[dat$teacher_id_spec]
  

 # potential outcomes ------------------------------------------------------

  Y_0_ijk <- rnorm(outcome_dat %*% a )
  
  
  
}


pr_star <- .5
b <- matrix(c(1.2, 0.8, -0.25, 0.6, -0.4, 1, 1))  #what values for pi 6 and 7 in the notes?
                        
