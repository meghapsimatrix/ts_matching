library(tidyverse)
library(MatchIt)
library(lme4)
library(cobalt)
library(broom)
library(matchMulti)
library(nlme)
library(optmatch)


# simple matching ---------------------------------------------------------

match_them <- function(dat, 
                       equation,
                       ps_method = "nearest",
                       caliper =  .25,
                       exact = NULL,
                       replace = FALSE){
  
  
  m_out <- matchit(equation,  
                   method = ps_method,
                   caliper = caliper,
                   exact = exact,
                   replace = replace,
                   data = dat)
  
  match_dat <- match.data(m_out)
  
  return(match_dat)
  
}


# method 3 ----------------------------------------------------------------

multi_match <- function(dat, 
                        l1_cov, 
                        l2_cov = NULL, 
                        trt, 
                        l2_id,
                        add_id = FALSE) {
  # Notes
  # dat = data frame
  # l1_cov = vector of unit-level covariate names
  # l2_cov = vector of cluster-level covariate names (categorical versions for refined covariate balance)
  # trt = name of treatment variable
  # l2_id = name of cluster-level identifier

  
  # set caliper for cluster-level pairing
  cluster_caliper <- buildCaliper(data = dat, 
                                  treatment = trt, 
                                  ps.vars = l1_cov, 
                                  group.id = l2_id, 
                                  caliper = 0.25)
  
  # execute matching
  matchout <- matchMulti(df, 
                         treatment = trt, 
                         school.id = l2_id, 
                         school.fb = list(l2_cov), 
                         student.vars = l1_cov,
                         match.students = FALSE, 
                         verbose = FALSE, 
                         school.caliper = cluster_caliper)
  
  # save matched data
  mdata <- as.data.frame(matchout$matched)
  
  if(add_id == TRUE){
    mdata$pair_id <- (mdata$l3id * 100) + mdata$pair_id
    mdata <- mdata %>% 
      select(-l3id)
    
  }
  
  return(mdata)
  
} 



# method 6 ----------------------------------------------------------------

mm_6 <- function(dat, 
                 l1_cov, 
                 l2_cov, 
                 trt, 
                 l2_id, 
                 l3_id) {
  # Notes
  # dat = data frame
  # l1_cov = vector of unit-level covariate names
  # l2_cov = vector of cluster-level covariate names (categorical versions for refined covariate balance)
  # NOTE: Not feasible to do l2 refined covariate balance with number of clusters within a site.
  # This part of the matching is removed for within-site matching.
  # trt = name of treatment variable
  # l2_id = name of cluster-level identifier
  # l3_id = name of site-level identifier
  
  # restrict matching to sites that have T & C clusters
  
  # this code is in the method 9 code too
  df$tmp <- df[, trt]
  df$l3id <- df[, l3_id]
  
  dat_m <- 
    dat %>% 
    group_by(l3id) %>% 
    mutate(tmpm = mean(tmp)) %>% 
    ungroup() %>%
    filter(tmpm > 0 & tmpm < 1) %>% 
    select(-tmp, -tmpm)
  
  hold <- 
    dat_m %>%
    group_by(l3id) %>%
    group_modify(~ (mm_them(.x)))  # add in other arguments
  
  return(hold)
  

} 



# method_9 ----------------------------------------------------------------

# Clusters only (multimatch) & Within site-defined group matching
# NOTE: to execute multimatch within sites, had to make two changes to the matching approach:
#       1. Removed refined covariate balance option for cluster-level covariates
#       2. Had to increase the caliper to 1.00 instead of 0.25

# Matching function
mm_9 <- function(df, 
                 l1_cov, 
                 l2_cov, 
                 trt, 
                 l2_id, 
                 l3_id) {
  # Notes
  # df = data frame
  # l1_cov = vector of unit-level covariate names
  # l2_cov = vector of cluster-level covariate names (categorical versions for refined covariate balance)
  # NOTE: Not feasible to do l2 refined covariate balance with number of clusters within a site.
  # This part of the matching is removed for within-site matching.
  # trt = name of treatment variable
  # l2_id = name of cluster-level identifier
  # l3_id = name of categorical variable that defines site groupings
  
  # restrict matching to sites that have T & C clusters
  df$tmp <- df[,trt]
  df$l3id <- df[,l3.id]
  
  dfm <- df %>% 
    group_by(l3id) %>% 
    mutate(tmpm = mean(tmp)) %>% 
    ungroup() %>%
    filter(tmpm > 0 & tmpm < 1) %>%
    select(-tmp, -tmpm)
  
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
