library(tidyverse)
library(nlme)
library(lme4)
library(metafor)


set.seed(20220224)


# Sample size -------------------------------------------------------------

k <- 300 #schools
j <- 1000 # teachers 
i <- 100000 # students



N <- i * j * k # total number of obs



# ICC ---------------------------------------------------------------------

icc3 <- 0.1  # icc for school 
icc2 <- 0.2 # icc for teachers
icc1 <- 1


# covariates --------------------------------------------------------------

Z_k <-rnorm(k)  # mean 0 sd 1 by default
W_jk <- rnorm(k *j)





beta_0 <- 0.3

e_ij <- rnorm(N, 1)
v_ij <- rnorm(k * j, 0, icc2)
u_j <- rnorm(k, 0, icc3)

dat <- data.frame(student_id = 1:N,
                  teacher_id = rep(rep(1:j, times = k), each = i),
                  school_id = rep(1:k, each = N/k))

dat$teacher_id_spec <- as.numeric(paste0(dat$school_id, dat$teacher_id))
dat$v_ij_long <- v_ij[dat$teacher_id_spec]
dat$u_j_long <- u_j[dat$school_id]
dat$e_ijk <- rnorm(N)

dat$y <- beta_0 + dat$v_ij_long + dat$u_j_long + dat$e_ijk
# we can do X %*% beta



# check -------------------------------------------------------------------

mod <- lmer(y ~ 1 + (1 | teacher_id_spec) + (1|school_id),
            data = dat)

summary(mod)








