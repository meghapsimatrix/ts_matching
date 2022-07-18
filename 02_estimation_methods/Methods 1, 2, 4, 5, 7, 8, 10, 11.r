##############################################
# Read in packages and set working directory #
##############################################

# Read in packages #
options(java.parameters = "- Xmx1024m")
library(gdata)
library(gtools)
library(MatchIt)
library(lme4)
library(scclust)
library(distances)
library(tidyverse)
library(magrittr)
library(rstanarm)
library(Hmisc)
library(broom.mixed)

# Set working directory #
setwd("C:/Users/mcitkowicz/American Institutes for Research in the Behavioral Sciences/TFA AC Evaluation - !Methods Paper/simulation/code")

#################
# Generate data #
#################

# Read in function #
source("ts_matching-main 2/01_dgm/01_dgm_function.R")

# Generate data #
dat <- generate_data(k = 50, # schools
					 j = 8, # teachers
					 i = 20, # students some default value
					 icc3 = .10, #school level icc
					 icc2 = .10, # teacher level icc
					 R2 = 0.40, # r-sq
					 ps_coef = c(.07, 0.8, -0.25, 0.6, 0.001, 1, 1), # coefficients for ps model add default
					 pr_star = .55, # overall proportion of teachers in trt group add default
					 outcome_coef = c(.79,.04,.05,-.12,.44,.005,.03), # coefficients for outcome model  add default
					 delta = .2 # treatment effect add default here
					 )
table(dat$D)

# Define quintiles #
quint <- with(dat, quantile(Z_k, seq(0, 1, 0.2)))
dat$quintile <- cut(dat$Z_k, quint, labels = c("A","B","C","D","E"), include.lowest = TRUE)

# Subset data to the cluster level #
vars <- c("teacher_id","school_id","Z_k","W_jk","X_jk","D","quintile")
dat_agg <- dat[!duplicated(dat[,vars]),vars]

table(dat_agg$D)

# Save data #
write.csv(dat,file="test/student data.csv",row.names=FALSE)
write.csv(dat_agg,file="test/teacher data.csv",row.names=FALSE)

####################
# Method 1 (units) #
####################

caliper <- .25

res_M1 <- matchit(D ~ X_ijk + W_jk + Z_k, data = dat, method = "nearest", replace = FALSE, exact = NULL, caliper = caliper)
dat_M1 <- match.data(res_M1)

table(dat_M1$D)
summary(res_M1,standardize = TRUE)

#######################
# Method 2 (clusters) #
#######################

caliper <- .25

res_M2 <- matchit(D ~ X_jk + W_jk + Z_k, data = dat_agg, method = "nearest", replace = FALSE, exact = NULL, caliper = caliper)
dat_M2 <- match.data(res_M2)

table(dat_M2$D)
summary(res_M2,standardize = TRUE)

####################
# Method 4 (units) #
####################

caliper <- 1

res_M4 <- matchit(D ~ X_ijk + W_jk + Z_k, data = dat, method = "nearest", replace = FALSE, exact = "school_id", caliper = caliper)
dat_M4 <- match.data(res_M4)

table(dat_M4$D)
summary(res_M4,standardize = TRUE)

#######################
# Method 5 (clusters) #
#######################

caliper <- 1

res_M5 <- matchit(D ~ X_jk + W_jk + Z_k, data = dat_agg, method = "nearest", replace = FALSE, exact = "school_id", caliper = caliper)
dat_M5 <- match.data(res_M5)

table(dat_M5$D)
summary(res_M5,standardize = TRUE)

####################
# Method 7 (units) #
####################

caliper <- 1

res_M7 <- matchit(D ~ X_ijk + W_jk + Z_k, data = dat, method = "nearest", replace = FALSE, exact = "quintile", caliper = caliper)
dat_M7 <- match.data(res_M7)

table(dat_M7$D)
summary(res_M7,standardize = TRUE)

#######################
# Method 8 (clusters) #
#######################

caliper <- 1

res_M8 <- matchit(D ~ X_jk + W_jk + Z_k, data = dat_agg, method = "nearest", replace = FALSE, exact = "quintile", caliper = caliper)
dat_M8 <- match.data(res_M8)

table(dat_M8$D)
summary(res_M8,standardize = TRUE)

#####################
# Method 10 (units) #
#####################

source("Hybrid matching function.r")

caliper <- 1

res_M10 <- fun.match(dat=dat, 
                               SiteID="school_id", 
                               TeacherID="student_id", 
                               tx.var="D", 
                               exact.vars=NULL, 
                               teacher.vars=c("X_ijk","W_jk"), 
                               site.vars="Z_k",
                               group="quintile",
                               crc=caliper,
                               replacement=FALSE,
							   ratio=1,
							   seed=1234)
dat_M10 <- res_M10$`Matched Data`

table(dat_M10$D)

########################
# Method 11 (clusters) #
########################

source("Hybrid matching function.r")

caliper <- 1

res_M11 <- fun.match(dat=dat_agg, 
                               SiteID="school_id", 
                               TeacherID="teacher_id", 
                               tx.var="D", 
                               exact.vars=NULL, 
                               teacher.vars=c("X_jk","W_jk"), 
                               site.vars="Z_k",
                               group="quintile",
                               crc=caliper,
                               replacement=FALSE,
							   ratio=1,
							   seed=1234)
dat_M11 <- res_M11$`Matched Data`

table(dat_M11$D)

