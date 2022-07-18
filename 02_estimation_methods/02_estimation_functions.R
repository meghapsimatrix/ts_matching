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


match_hybrid <-  function(dat, 
                          SiteID, 
                          TeacherID, 
                          tx.var, 
                          exact.vars, 
                          teacher.vars, 
                          site.vars, 
                          group, 
                          crc, 
                          replacement, 
                          ratio, 
                          seed){
  
  # set seed #
  set.seed(seed)
  
  # set data #
  if(is.null(site.vars)){
    df <- as.data.frame(dat[,c(SiteID, TeacherID, tx.var, exact.vars, teacher.vars, group)])
  }else{
    df <- as.data.frame(dat[,c(SiteID, TeacherID, tx.var, exact.vars, teacher.vars, site.vars, group)])
  }
  
  # estimate propensity score with two-level RE model #
  if(is.null(site.vars)){
    pfmla <- as.formula(paste(tx.var," ~ 1 + ",paste(teacher.vars,collapse="+")," + (1|",SiteID,")"))
  }else{
    pfmla <- as.formula(paste(tx.var," ~ 1 + ",paste(teacher.vars,collapse="+"),"+",paste(site.vars,collapse="+")," + (1|",SiteID,")"))
  }
  p3 <- glmer(pfmla, data=df, family=binomial(link="logit"))
  #p3 <- rstanarm::stan_glmer(pfmla, data=df, family=binomial(link="logit"))
  df$p3 <- log(fitted(p3)/(1-fitted(p3)))
  #df$p3 <- colMeans(rstanarm::posterior_linpred(p3))
  
  # save school-level EB coefficient estimates for between-site matching #
  
  #psrx <- ranef(p3)[SiteID]
  psrx <- broom.mixed::tidy(p3,effects='ran_vals')$estimate
  #psrx <- as.data.frame(cbind(c(1:summary(p3)$ngrps),c(rep(1,summary(p3)$ngrps)),as.data.frame(ranef(p3)[SiteID])))
  psrx <- as.data.frame(cbind(c(1:length(psrx)),c(rep(1,length(psrx))),as.data.frame(broom.mixed::tidy(p3,effects='ran_vals')$estimate)))
  names(psrx) <- c(SiteID,"ALL","U0")
  #psfx <- as.data.frame(cbind(1,t(as.matrix(fixef(p3)))))
  psfx <- as.data.frame(cbind(1,t(as.matrix(broom.mixed::tidy(p3)$estimate))))
  #names(psfx)[names(psfx) == "V1"] <- "ALL"
  names(psfx) <- c('ALL',broom.mixed::tidy(p3)$term)
  pspar3 <- merge(psfx,psrx, by="ALL", all.x=TRUE)
  pspar3$ps <- 3
  
  df$pscore <- df$p3
  psparx <- pspar3
  
  # conduct two-Stage matching method #
  m.data2 <- NULL # create placeholder for matched data
  m.data3 <- NULL # create placeholder for matched data
  
  # set treatment schools #
  sids <- sort(unique(df[df[,tx.var]==1,SiteID]))
  
  output <- list()
  output[['glm output']] <- p3
  output[['data']] <- df
  
  # loop over treatment schools #
  for(j in 1:length(sids)) {
    
    # set data #
    sid <- sids[j]
    t <- df[df[,SiteID]==sid & df[,tx.var]==1,] # use treatment units in school j
    c1 <- df[df[,SiteID]==sid & df[,tx.var]==0,] # use control units in school j
    if(any(df[,SiteID]!=sid & df[,group] %in% t[,group] & df[,tx.var]==0)) {
      c2 <- df[df[,SiteID]!=sid & df[,group] %in% t[,group] & df[,tx.var]==0,]
    }else{
      c2 <- df[df[,SiteID]!=sid & df[,group] %in% t[,group],][1,]
      c2[,2:dim(df)[2]] <- NA
    }		
    tmp1 <- rbind(t,c1)
    
    pspar.j <- as.vector(psparx[psparx[,SiteID]==j,]) # use ps model parameter estimates for school j
    
    # match 1: within-school #
    m.out1 <- tryCatch(matchit(tmp1[,tx.var] ~ tmp1$pscore, data = tmp1, method = "nearest", replace = replacement, exact = exact.vars, caliper=crc, distance=tmp1$pscore,ratio=ratio),error=function(e)NULL)
    output[['Within School Matches']][[j]] <- m.out1
    
    # create data frame of matched units from match 1 # 
    m.data3w <- tryCatch(match.data(m.out1, weights="PSW3W", distance="PS3W"),error=function(e)NULL)
    
    # save within-site matches for later analysis #
    m.data2.j <- tryCatch(match.data(m.out1, weights="PSW2", distance="PS2"),error=function(e)NULL)
    m.data2 <- rbind(m.data2,m.data2.j)
    
    # select unmatched treatment units in school j #
    t.m1 <- tryCatch(match.data(m.out1, weights="PSW3W", distance="PS3W", group="treat"),error=function(e)NULL)
    t.m1 <- t.m1[ ,c(TeacherID,exact.vars,"PSW3W")]
    t.x <- t; t.x$PSW3W <- NA # placeholder for schools with no within-school matches
    t2 <- tryCatch(merge(t, t.m1, by=c(TeacherID,exact.vars), all.x=TRUE),error=function(e)t.x)
    t2$PSW3W <- ifelse(is.na(t2$PSW3W),0,1)
    t2 <- t2[t2$PSW3W==0,]
    t2$pscore.j <- t2$pscore
    
    # set data for between-school match #
    # re-calculate control unit propensity score based on school j pscore model parameter estimates #
    c2$PSW3W <- 0
    
    if(is.null(site.vars)){
      tfmla <- NULL
      for(i in teacher.vars){
        tfmla.i <- paste0(pspar.j[,i],"*c2$",i)
        tfmla <- paste(tfmla,tfmla.i,sep=" + ")
      }
      
      newfmla <- paste(pspar.j[,"(Intercept)"],tfmla,"+",pspar.j[,"U0"])
    }else{
      tfmla <- NULL
      for(i in teacher.vars){
        tfmla.i <- paste0(pspar.j[,i],"*c2$",i)
        tfmla <- paste(tfmla,tfmla.i,sep=" + ")
      }
      
      sfmla <- NULL
      for(i in site.vars){
        sfmla.i <- paste0(pspar.j[,i],"*c2$",i)
        sfmla <- paste(sfmla,sfmla.i,sep=" + ")
      }
      newfmla <- paste(pspar.j[,"(Intercept)"],tfmla,sfmla,"+",pspar.j[,"U0"])
    }
    
    c2$pscore.j <- eval(parse(text=newfmla))
    
    tmp2 <- rbind(t2,c2)
    
    # match 2: between-school: only if control units fall within treatment unit range #
    m.out2 <- tryCatch(matchit(tmp2[,tx.var] ~ tmp2$pscore.j, data = tmp2, method = "nearest", replace = replacement, exact = exact.vars, caliper=crc, distance=tmp2$pscore.j,ratio=ratio),error=function(e)NULL)
    output[['Bewteen School Matches']][[j]] <- m.out2
    
    # create data frame of matched units from match 2 #
    m.data3b <- tryCatch(match.data(m.out2, weights="PSW3B", distance="PS3B"),error=function(e)NULL)
    if (!is.null(m.data3b)) m.data3b$PS3W <- NA else m.data3b <- NULL
    
    # combine matched files #
    if (!is.null(m.data3w)) m.data3w$PSW3B  <- 0 else m.data3w <- NULL
    if (!is.null(m.data3w)) m.data3w$PS3B <- NA else m.data3w <- NULL
    if (!is.null(m.data3w)) m.data3w$pscore.j <- m.data3w$pscore else m.data3w <- NULL
    
    if (!is.null(m.data3w)) m.data3w <- m.data3w[,order(names(m.data3w))] else m.data3w <- NULL
    if (!is.null(m.data3b)) m.data3b <- m.data3b[,order(names(m.data3b))] else m.data3b <- NULL
    m.datax <- rbind(m.data3w, m.data3b)
    
    # add treatment school id #
    if (!is.null(m.datax)) m.datax$Site_match <- sid else m.datax <- NULL
    
    # append to master file #
    if (!is.null(m.data3)) m.data3 <- m.data3[,order(names(m.data3))] else m.data3 <- NULL
    if (!is.null(m.datax)) m.datax <- m.datax[,order(names(m.datax))] else m.datax <- NULL
    m.data3 <- rbind(m.data3, m.datax)
    
  } # close 2 stage matching loop
  output[['Matched Data']] <- m.data3
  
  return(output)
}
