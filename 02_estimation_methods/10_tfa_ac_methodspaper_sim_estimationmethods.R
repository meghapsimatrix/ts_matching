####################################################
## MATCHING METHODS FOR MULTISITE CLUSTER DESIGNS ##
## Simulation Code #################################
####################################################
## jrickles 06.16.22 ###############################

## NOTE: Adding to code developed by Megha Joshi
## estimation_methods_run_through.R

library(tidyverse)

library(MatchIt)
library(lme4)
library(cobalt)
library(broom)
library(matchMulti)
library(nlme)
library(optmatch)

# https://cran.r-project.org/web/packages/matchMulti/matchMulti.pdf

# load data ---------------------------------------------------------------

source("01_dgm/01_dgm_function.R")

set.seed(20220615)

k <- 50 #schools
j <- 8  # teachers 
i <- 20 # students

icc3 <- 0.20  # icc for school 
icc2 <- 0.20 # icc for teachers
R2 <- .40


ps_coef <- matrix(c(.07, 0.8, -0.25, 0.6, 0.001, 1, 1))  #what values for pi 6 and 7 in the notes?
pr_star <- .6

outcome_coef <- matrix(c(1, 0.3, .5, .4, -0.2, 1, 1))

delta <- -0.4

example_dat <- generate_data(k = k,
                             j = j, 
                             i = i, 
                             icc3 = icc3, 
                             icc2 = icc2, 
                             R2 = R2,
                             ps_coef = ps_coef, 
                             pr_star = pr_star, 
                             outcome_coef = outcome_coef, 
                             delta = delta)

glimpse(example_dat)

example_dat %>%
  group_by(school_id, D) %>%
  summarize(n_teach = n_distinct(teacher_id)) %>%
  View()



# create cluster level data - cluster and site covariates


cluster_level_dat <- 
  example_dat %>%
  group_by(teacher_id) %>%
  summarize(school_id = mean(school_id),
            Z_k = mean(Z_k),
            W_jk = mean(W_jk),
            D = mean(D),
            Y_ijk = mean(Y_ijk)) %>%
  ungroup()


# quintiles based on site-level scores
# what should we do with these?
quint <- with(example_dat, quantile(Z_k, seq(0, 1, 0.2)))

example_dat$quintile <- cut(example_dat$Z_k, quint,
                            labels = c(1,2,3,4,5),
                            include.lowest = TRUE)


glimpse(example_dat)
summary(example_dat)


# Estimate propensity scores ----------------------------------------------

# do we do hlm for ps model or not?


# unit --------------------------------------------------------------------

# Not multilevel
unit_ps_model_nm <- glm(D ~ X_ijk + W_jk + Z_k, 
                         family = "binomial", 
                         data = example_dat)

example_dat$ps_unit_nm <- predict(unit_ps_model_nm, type = "link")


# Multilevel 
unit_ps_model <- glmer(D ~ X_ijk + W_jk + Z_k + 
                         (1 | school_id),
                       family = "binomial", 
                       data = example_dat)

example_dat$ps_unit <- predict(unit_ps_model, type = "link")
example_dat$ps_unit_pr <- predict(unit_ps_model, type = "response")

# check common support
example_dat %>%
  mutate(D = as.character(D)) %>%
  ggplot(aes(x = ps_unit_nm, fill = D)) +
  geom_density(alpha = 0.5) + 
  theme_minimal()



# cluster -----------------------------------------------------------------

# not multilevel
cluster_ps_model_nm <- glm(D ~ W_jk + Z_k,
                          family = "binomial",
                          data = cluster_level_dat)

cluster_level_dat$ps_cluster_nm <- predict(cluster_ps_model_nm, type = "link")


# multilevel
cluster_ps_model <- glmer(D ~ W_jk + Z_k + 
                            (1 | school_id),
                         family = "binomial",
                         data = cluster_level_dat)

cluster_level_dat$ps_cluster <- predict(cluster_ps_model, type = "link")

#check common support
cluster_level_dat %>%
  mutate(D = as.character(D)) %>%
  ggplot(aes(x = ps_cluster_nm, fill = D)) +
  geom_density(alpha = 0.5) + 
  theme_minimal()



# Method 3 ----------------------------------------------------------------
# Clusters only (multimatch) & No site-level distinction

# prep data for matching
summary(example_dat)

tcn <- length(unique(example_dat$teacher_id[example_dat$D==1])) # get number of clusters in treatment group

std.cov <- c("X_ijk") # define unit-level covariates
tch.cov <- c("W_jk", "Z_k") # define cluster-level covariates

# create quintiles from cluster covariate for refined covariate balance
example_dat$W_q5 <- cut(example_dat$W_jk, 5,
                            labels = c(1,2,3,4,5),
                            include.lowest = TRUE) 
example_dat$Z_q5 <- cut(example_dat$Z_k, 5,
                            labels = c(1,2,3,4,5),
                            include.lowest = TRUE)

# Matching function
fun.matchmulti3 <- function(df, l1.cov, l2.cov, trt, l2.id) {
# Notes
# df = data frame
# l1.cov = vector of unit-level covariate names
# l2.cov = vector of cluster-level covariate names (categorical versions for refined covariate balance)
# trt = name of treatment variable
# l2.id = name of cluster-level identifier

 # set caliper for cluster-level pairing
 cluster.caliper <- buildCaliper(data = df, treatment = trt, ps.vars = l1.cov, group.id = l2.id, caliper = 0.25)

 # execute matching
 matchout <- matchMulti(df, treatment = trt, school.id = l2.id, school.fb = list(l2.cov), student.vars = l1.cov,
   match.students = FALSE, verbose=FALSE, school.caliper = cluster.caliper)

 # save matched data
 mdata <- as.data.frame(matchout$matched)

 return(mdata)

} 

test <- fun.matchmulti3(df = example_dat, l1.cov = std.cov, l2.cov = c("W_q5", "Z_q5"), trt = "D", l2.id = "teacher_id")
length(unique(test$teacher_id[test$D==1]))
length(unique(example_dat$teacher_id[example_dat$D==1]))


# Method 6 ----------------------------------------------------------------
# Clusters only (multimatch) & Within-site matching
# NOTE: to execute multimatch within sites, had to make two changes to the matching approach:
#       1. Removed refined covariate balance option for cluster-level covariates
#       2. Had to increase the caliper to 1.00 instead of 0.25

# Matching function
fun.matchmulti6 <- function(df, l1.cov, l2.cov, trt, l2.id, l3.id) {
# Notes
# df = data frame
# l1.cov = vector of unit-level covariate names
# l2.cov = vector of cluster-level covariate names (categorical versions for refined covariate balance)
  # NOTE: Not feasible to do l2 refined covariate balance with number of clusters within a site.
  # This part of the matching is removed for within-site matching.
# trt = name of treatment variable
# l2.id = name of cluster-level identifier
# l3.id = name of site-level identifier
  
 # restrict matching to sites that have T & C clusters
 df$tmp <- df[,trt]
 df$l3id <- df[,l3.id]
 dfm <- df %>% group_by(l3id) %>% mutate(tmpm = mean(tmp)) %>% ungroup() %>%
   filter(tmpm > 0 & tmpm < 1) %>% select(-tmp, -tmpm)

 # loop over each site to execute matching separately within each site
 sid <- unique(dfm$l3id) # index of sites to include in matching
 hold <- NULL # placeholder to store matched samples
 
 for(k in 1:length(sid)) {
   sid.k <- sid[k]
   df.k <- dfm[dfm$l3id == sid.k, ]
   df.k <- as.data.frame(df.k)

   # set caliper for cluster-level pairing
   cluster.caliper <- buildCaliper(data = df.k, treatment = trt, ps.vars = l1.cov, group.id = l2.id, caliper = 1.00)
   
   # match clusters within site k
   matchout <- matchMulti(df.k, treatment = trt, school.id = l2.id, student.vars = l1.cov,
    match.students = FALSE, verbose=FALSE, school.caliper = cluster.caliper)
   
   # save matched data
   mdata.k <- as.data.frame(matchout$matched)
   mdata.k$pair.id <- (mdata.k$l3id*100) + mdata.k$pair.id
   mdata.k <- subset(mdata.k, select = -c(l3id))
   
   hold <- rbind(hold, mdata.k)

 }
 
 return(hold)

} 


test <- fun.matchmulti6(df = example_dat, 
                        l1.cov = std.cov, 
                        l2.cov = c("W_q5"), 
                        trt = "D", 
                        l2.id = "teacher_id", 
                        l3.id = "school_id")

length(unique(test$teacher_id[test$D==1]))



# Method 9 ----------------------------------------------------------------
# Clusters only (multimatch) & Within site-defined group matching
# NOTE: to execute multimatch within sites, had to make two changes to the matching approach:
#       1. Removed refined covariate balance option for cluster-level covariates
#       2. Had to increase the caliper to 1.00 instead of 0.25

# Matching function
fun.matchmulti9 <- function(df, l1.cov, l2.cov, trt, l2.id, l3.id) {
# Notes
# df = data frame
# l1.cov = vector of unit-level covariate names
# l2.cov = vector of cluster-level covariate names (categorical versions for refined covariate balance)
  # NOTE: Not feasible to do l2 refined covariate balance with number of clusters within a site.
  # This part of the matching is removed for within-site matching.
# trt = name of treatment variable
# l2.id = name of cluster-level identifier
# l3.id = name of categorical variable that defines site groupings
  
 # restrict matching to sites that have T & C clusters
 df$tmp <- df[,trt]
 df$l3id <- df[,l3.id]
 dfm <- df %>% group_by(l3id) %>% mutate(tmpm = mean(tmp)) %>% ungroup() %>%
   filter(tmpm > 0 & tmpm < 1) %>% select(-tmp, -tmpm)

 # loop over each site to execute matching separately within each site
 sid <- unique(dfm$l3id) # index of groups to include in matching
 hold <- NULL # placeholder to store matched samples
 
 for(k in 1:length(sid)) {
   sid.k <- sid[k]
   df.k <- dfm[dfm$l3id == sid.k, ]
   df.k <- as.data.frame(df.k)

   # set caliper for cluster-level pairing
   cluster.caliper <- buildCaliper(data = df.k, treatment = trt, ps.vars = l1.cov, group.id = l2.id, caliper = 1.00)
   
   # match clusters within site k
   matchout <- matchMulti(df.k, treatment = trt, school.id = l2.id, student.vars = l1.cov, school.fb = list(l2.cov),
    match.students = FALSE, verbose=FALSE, school.caliper = cluster.caliper)
   
   # save matched data
   mdata.k <- as.data.frame(matchout$matched)
   mdata.k$pair.id <- (mdata.k$l3id*100) + mdata.k$pair.id
   mdata.k <- subset(mdata.k, select = -c(l3id))
   
   hold <- rbind(hold, mdata.k)

 }
 
 return(hold)

} 


test <- fun.matchmulti9(df = example_dat, l1.cov = std.cov, l2.cov = c("W_q5", "Z_q5"), trt = "D", l2.id = "teacher_id", l3.id = "quintile")
length(unique(test$teacher_id[test$D==1]))


# Method 12 ----------------------------------------------------------------
# Clusters only (multimatch) & hybrid matching
# NOTE: to execute multimatch within sites, had to make two changes to the matching approach:
#       1. Removed refined covariate balance option for cluster-level covariates
#       2. Had to increase the caliper to 1.00 instead of 0.25

# Matching function
fun.matchmulti12 <- function(df, l1.cov, l2.cov, trt, l2.id, l3.id, group) {
# Notes
# df = data frame
# l1.cov = vector of unit-level covariate names
# l2.cov = vector of cluster-level covariate names (categorical versions for refined covariate balance)
  # NOTE: Not feasible to do l2 refined covariate balance with number of clusters within a site.
  # This part of the matching is removed for within-site matching.
# trt = name of treatment variable
# l2.id = name of cluster-level identifier
# l3.id = name of site-level identifier
# group = name of categorical variable that defines site groupings  

# STEP 1. WITHIN-SITE MATCHING
  
 # restrict matching to sites that have T & C clusters
 df$tmp <- df[,trt]
 df$l2id <- df[,l2.id]
 df$l3id <- df[,l3.id]
 dfm <- df %>% group_by(l3id) %>% mutate(tmpm = mean(tmp)) %>% ungroup() %>%
   filter(tmpm > 0 & tmpm < 1) %>% select(-tmp, -tmpm)

 # loop over each site to execute matching separately within each site
 sid <- unique(dfm$l3id) # index of sites to include in matching
 hold1 <- NULL # placeholder to store matched samples
 
 for(k in 1:length(sid)) {
   sid.k <- sid[k]
   df.k <- dfm[dfm$l3id == sid.k, ]
   df.k <- as.data.frame(df.k)

   # set caliper for cluster-level pairing
   cluster.caliper <- buildCaliper(data = df.k, treatment = trt, ps.vars = l1.cov, group.id = l2.id, caliper = 1.00)
   
   # match clusters within site k
   matchout <- matchMulti(df.k, treatment = trt, school.id = l2.id, student.vars = l1.cov,
    match.students = FALSE, verbose=FALSE, school.caliper = cluster.caliper)
   
   # save matched data
   mdata.k <- as.data.frame(matchout$matched)
   mdata.k$pair.id <- (mdata.k$l3id*100) + mdata.k$pair.id
   mdata.k <- subset(mdata.k, select = -c(l3id))
   
   hold1 <- rbind(hold1, mdata.k)
 }


# STEP 2. WITHIN-GROUP MATCHING FOR CLUSTERS WITHOUT A WITHIN-SITE MATCH  
 
 # Identify Clusters NOT in a Within-Site Pair
 wsp <- hold1 %>% select(l2id) %>% mutate(inwsp = 1)
 dfm2 <- df %>% left_join(wsp, by = "l2id") %>% filter(is.na(inwsp) == TRUE)
    
 # restrict matching to groups that have T & C clusters
 dfm2$grp <- dfm2[,group]
 dfm2$tmp <- dfm2[,trt]
 dfm2 <- dfm2 %>% group_by(grp) %>% mutate(tmpm = mean(tmp)) %>% ungroup() %>%
   filter(tmpm > 0 & tmpm < 1) %>% select(-tmp, -tmpm)
 
 # loop over each group to execute matching separately within each site
 sid <- unique(dfm2$grp) # index of groups to include in matching
 hold <- NULL # placeholder to store matched samples
 
 for(k in 1:length(sid)) {
   sid.k <- sid[k]
   df.k <- dfm2[dfm2$grp == sid.k, ]
   df.k <- as.data.frame(df.k)

   # set caliper for cluster-level pairing
   cluster.caliper <- buildCaliper(data = df.k, treatment = trt, ps.vars = l1.cov, group.id = l2.id, caliper = 1.00)
   
   # match clusters within group k
   matchout <- matchMulti(df.k, treatment = trt, school.id = l2.id, student.vars = l1.cov,
    match.students = FALSE, verbose=FALSE, school.caliper = cluster.caliper)
   
   # save matched data
   mdata.k <- as.data.frame(matchout$matched)
   mdata.k$pair.id <- (as.numeric(as.character(mdata.k$grp))*100000) + mdata.k$pair.id
   mdata.k <- subset(mdata.k, select = -c(inwsp, l3id, grp))
   
   hold2 <- rbind(hold, mdata.k)

 }
 
 hold <- rbind(hold1, hold2)
 hold$match <- ifelse(hold$pair.id < 100000, "ws", "bs")
 hold <- subset(hold, select = -c(l2id))
 return(hold)

} 

test <- fun.matchmulti12(df = example_dat, l1.cov = std.cov, l2.cov = c("W_q5"), trt = "D", l2.id = "teacher_id", l3.id = "school_id", group = "quintile")
length(unique(test$teacher_id[test$D==1]))





#match.simple <- matchMulti(example_dat, treatment = 'D', school.id = 'teacher_id',
#  match.students = FALSE, student.vars = std.cov, verbose=TRUE, keep.target = round(tcn*0.8), tol = 0.25)

#balanceMulti(match.simple, student.cov = std.cov, school.cov = tch.cov)









summary(example_dat)
summary(match.data)

sum(example_dat$D)
sum(match.data$D)

length(unique(example_dat$teacher_id[example_dat$D==1]))
length(unique(match.data$teacher_id[match.data$D==1]))

# Method 1 ----------------------------------------------------------------

#optimal method in matchit doesn't like exact or caliper

match_them <- function(dat, 
                       ps_method = "nearest",
                       exact = NULL, 
                       ps = dat$ps_unit){
  

  m_out <- matchit(D ~ ps,  # the rhs doesn't matter I think 
                   caliper = .25,
                   method = ps_method,
                   exact = exact,
                   distance = ps,
                   data = dat)
  
  match_dat <- match.data(m_out)
  
  return(match_dat)
  
}



# unit level data ---------------------------------------------------------

example_dat <- example_dat %>%
  mutate(ps = ps_unit)

# ignoring sites and groups -----------------------------------------------

# switch distance to appropriate ps (multi or not multilevel)
system.time(m_out_1 <- match_them(dat = example_dat, 
                                  ps = example_dat$ps_unit))


# exact on site -----------------------------------------------------------

# this breaks down with this data - but was working with another seedless data i generated :D 
# exact match on site by looping over site
system.time(site_match_data <- 
  example_dat %>%
  group_by(school_id) %>%
  group_modify(~ match_them(.x)))

# exact match on site using exact
# this takes a while 
system.time(m_out_site <- match_them(dat = example_dat, 
                                     exact = ~ school_id,
                                     ps = example_dat$ps_unit))



# exact on group ----------------------------------------------------------

# exact match on site by looping over site
system.time(group_match_data <- 
              example_dat %>%
              group_by(quintile) %>%
              group_modify(~ match_them(.x)))
  
# exact match on site using exact
# this takes a while 
system.time(m_out_group <- match_them(dat = example_dat, 
                                     exact = ~ quintile,
                                     ps = example_dat$ps_unit))




# cluster level data ------------------------------------------------------

cluster_level_dat <- cluster_level_dat %>%
  mutate(ps = ps_cluster)

# is this the same as above but we are using the cluster level data?

# switch distance to appropriate ps (multi or not multilevel)
system.time(m_out_2 <- match_them(dat = cluster_level_dat, 
                                  ps = cluster_level_dat$ps_cluster))






# looping over treatment sites
# not complete yet - need to figure out what is happening
# all sites have treatment teachers? 
treatment_sites <- 
  example_dat %>%
  filter(D == 1)

# i think all sites are treatment sites
# within schools some teacher are trt and some control
# so we need to mess with data generation so some schools don't have a lot of treatment units or control 
sids <- sort(unique(treatment_sites$school_id))
length(sids)
