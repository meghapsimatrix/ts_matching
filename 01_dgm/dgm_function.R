library(tidyverse)

# part of this is in base R bc that runs faster in sim studies
# part in dplyr tidyverse bc my brain works faster in that 
# will reconcile soon

# is there degree of heterogeneity in play in the dgm somewhere?

generate_data <- function(k, # schools
                          j, # teachers
                          i, # students
                          icc3, #school level icc
                          icc2, # teacher level icc
                          ps_coef, # coefficients for ps model
                          pr_star, # overall proportion of teachers in trt group
                          outcome_coef, # coefficients for outcome model 
                          delta # treatment effect
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
  
  r_k <- rnorm(k, 0, icc3)
  u_jk <- rnorm(k * j, 0, icc2)

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
                   sd = sqrt(1 - icc2 - icc3))
  
  Y_1_ijk <- Y_0_ijk + delta
  

  # Observed outcome --------------------------------------------------------
  dat <- 
    dat %>%
    mutate(Y_ijk = D * Y_1_ijk + (1 - D) * Y_0_ijk)
  
  
  return(dat)
  
  
  
}




# run function ------------------------------------------------------------

# memory exhaustion if I try to make all of these large 
# for checking the distributions - need to increased the numbers here
k <- 100 #schools
j <-20  # teachers 
i <- 10 # students

icc3 <- 0.05  # icc for school 
icc2 <- 0.20 # icc for teachers


ps_coef <- matrix(c(1.2, 0.8, -0.25, 0.6, -0.4, 1, 1))  #what values for pi 6 and 7 in the notes?
pr_star <- .5

outcome_coef <- matrix(c(1, 0.3, .5, .4, -0.2, 1, 1))

delta <- -0.4
                       

example_dat <- generate_data(k = k,
                             j = j, 
                             i = i, 
                             icc3 = icc3, 
                             icc2 = icc2, 
                             ps_coef = ps_coef, 
                             pr_star = pr_star, 
                             outcome_coef = outcome_coef, 
                             delta = delta)

save(example_dat, file = "data/example_data.RData")

# check if any is missing 
map(example_dat, ~ sum(is.na(.)))

# proportion of teachers in treatment group off a bit
# proportion still off
prop.table(table(example_dat$D))



# check covs --------------------------------------------------------------

# school level 

school_cov_dat <- 
  example_dat %>%
  select(school_id, Z_k) %>%
  distinct(.)

ggplot(school_cov_dat, aes(x = Z_k)) +
  geom_density() + 
  theme_minimal()

mean(school_cov_dat$Z_k)
sd(school_cov_dat$Z_k)

# teacher level 

teacher_cov_dat <- 
  example_dat %>%
  select(teacher_id, W_jk) %>%
  distinct(.)

ggplot(teacher_cov_dat, aes(x = W_jk)) +
  geom_density() + 
  theme_minimal()

mean(teacher_cov_dat$W_jk)
sd(teacher_cov_dat$W_jk)

# student level 

ggplot(example_dat, aes(x = X_ijk)) +
  geom_density() + 
  theme_minimal()

mean(example_dat$X_ijk)
sd(example_dat$X_ijk)

ggplot(example_dat, aes(x = U_ijk)) +
  geom_density() + 
  theme_minimal()

mean(example_dat$U_ijk)
sd(example_dat$U_ijk)

