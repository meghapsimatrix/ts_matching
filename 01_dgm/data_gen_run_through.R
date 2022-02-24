library(tidyverse)
library(nlme)
library(lme4)
library(metafor)


set.seed(20220224)

i <- 100000
j <- 1000
k <- 300


N <- i * j * k

tau <- 0.1
omega <- 0.2

beta_0 <- 0.3


v_ij <- rnorm(k * j, 0, tau)
u_j <- rnorm(k, 0, omega)

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








