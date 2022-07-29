library(tidyverse)
library(MatchIt)
library(lme4)
library(cobalt)
library(broom)
library(matchMulti)
library(nlme)
library(optmatch)
library(lmerTest)
library(broom.mixed)
library(janitor)
library(parameters)


# simple matching ---------------------------------------------------------

match_them <- function(dat, 
                       equation,
                       ps_method = "nearest",
                       caliper =  .25,
                       exact = NULL,
                       replace = FALSE,
                       student_level = FALSE,
                       student_dat = NULL){
  
  
  m_out <- matchit(equation,  
                   method = ps_method,
                   caliper = caliper,
                   exact = exact,
                   replace = replace,
                   data = dat)
  
  match_dat <- match.data(m_out)
  
  
  if(student_level == TRUE){
    
    
    match_dat <- semi_join(student_dat, match_dat, by = c("teacher_id")) %>%
                 left_join(match_dat %>%
                             select(teacher_id, distance, weights, subclass), by = c("teacher_id")) 
  
  }
    
  return(match_dat)
    
  
}


# multi matching ----------------------------------------------------------



multi_match <- function(dat, # data 
                        trt, # name of treatment var
                        l1_cov, # student level cov
                        l2_cov = NULL, # teacher level cov
                        l2_id,
                        l3_id = NULL,
                        add_id = "none",
                        caliper,
                        match_students = FALSE
                        ) {

  dat <- as.data.frame(dat)
  
  # set caliper for cluster-level pairing
  cluster_caliper <- buildCaliper(data = dat, 
                                  treatment = trt, 
                                  ps.vars = l1_cov, 
                                  group.id = l2_id, 
                                  caliper = caliper)
  
  # execute matching
  matchout <- tryCatch(matchMulti(dat, 
                         treatment = trt, 
                         school.id = l2_id, 
                         school.fb = list(l2_cov), 
                         student.vars = l1_cov,
                         match.students = match_students, 
                         verbose = FALSE, 
                         school.caliper = cluster_caliper),
                        error = function(e) NULL)
  
  if(!is.null(matchout)){
  # save matched data
  mdata <- as.data.frame(matchout$matched) %>%
    rename(pair_id = pair.id)
  
  if(add_id == "school"){
    
    mdata <- as.data.frame(mdata) 
    
    mdata$pair_id <- (mdata$school_id * 100) + mdata$pair_id
    
  } else if(add_id == "group"){
    
    mdata$pair_id <- (as.numeric(as.character(mdata$Z_q5)) * 100000) + mdata$pair_id
    
  } else if(add_id == "none"){
    
    mdata <- mdata
  }
  
    mdata$weights <- 1
  
  } else{
    
    mdata <- dat
    mdata$weights <- 0
    
  }
  
  
  
  return(mdata)
  
  
} 


# hybrid matching ---------------------------------------------------------

match_hybrid <-  function(dat, 
                          site_id, 
                          teacher_id, 
                          tx_var, 
                          exact_vars, 
                          teacher_vars, 
                          site_vars, 
                          group, 
                          crc, 
                          replacement, 
                          ratio, 
                          seed,
                          by_student,
                          student_dat){
  
  # set seed #
  set.seed(seed)
  
  # set data #
  if(is.null(site_vars)){
    df <- as.data.frame(dat[,c(site_id, 
                               teacher_id, 
                               tx_var, 
                               exact_vars, 
                               teacher_vars, 
                               group)])
  }else{
    df <- as.data.frame(dat[,c(site_id, 
                               teacher_id, 
                               tx_var, 
                               exact_vars, 
                               teacher_vars, 
                               site_vars, 
                               group)])
  }
  
  # estimate propensity score with two-level RE model #
  if(is.null(site_vars)){
    pfmla <- as.formula(paste(tx_var," ~ 1 + ",paste(teacher_vars,collapse="+")," + (1|",site_id,")"))
  }else{
    pfmla <- as.formula(paste(tx_var," ~ 1 + ",paste(teacher_vars,collapse="+"),"+",paste(site_vars,collapse="+")," + (1|",site_id,")"))
  }
  p3 <- glmer(pfmla, data=df, family=binomial(link="logit"))
  #p3 <- rstanarm::stan_glmer(pfmla, data=df, family=binomial(link="logit"))
  df$p3 <- log(fitted(p3)/(1-fitted(p3)))
  #df$p3 <- colMeans(rstanarm::posterior_linpred(p3))
  
  # save school-level EB coefficient estimates for between-site matching #
  
  #psrx <- ranef(p3)[site_id]
  psrx <- broom.mixed::tidy(p3,effects='ran_vals')$estimate
  #psrx <- as.data.frame(cbind(c(1:summary(p3)$ngrps),c(rep(1,summary(p3)$ngrps)),as.data.frame(ranef(p3)[site_id])))
  psrx <- as.data.frame(cbind(c(1:length(psrx)),c(rep(1,length(psrx))),as.data.frame(broom.mixed::tidy(p3,effects='ran_vals')$estimate)))
  names(psrx) <- c(site_id,"ALL","U0")
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
  sids <- sort(unique(df[df[,tx_var]==1,site_id]))
  
  output <- list()
  output[['glm output']] <- p3
  output[['data']] <- df
  
  # loop over treatment schools #
  for(j in 1:length(sids)) {
    
    # set data #
    sid <- sids[j]
    t <- df[df[,site_id]==sid & df[,tx_var]==1,] # use treatment units in school j
    c1 <- df[df[,site_id]==sid & df[,tx_var]==0,] # use control units in school j
    if(any(df[,site_id]!=sid & df[,group] %in% t[,group] & df[,tx_var]==0)) {
      c2 <- df[df[,site_id]!=sid & df[,group] %in% t[,group] & df[,tx_var]==0,]
    }else{
      c2 <- df[df[,site_id]!=sid & df[,group] %in% t[,group],][1,]
      c2[,2:dim(df)[2]] <- NA
    }		
    tmp1 <- rbind(t,c1)
    
    pspar.j <- as.vector(psparx[psparx[,site_id]==j,]) # use ps model parameter estimates for school j
    
    # match 1: within-school #
    m.out1 <- tryCatch(matchit(tmp1[,tx_var] ~ tmp1$pscore, data = tmp1, method = "nearest", replace = replacement, exact = exact_vars, caliper=crc, distance=tmp1$pscore,ratio=ratio),error=function(e)NULL)
    output[['Within School Matches']][[j]] <- m.out1
    
    # create data frame of matched units from match 1 # 
    m.data3w <- tryCatch(match.data(m.out1, weights="PSW3W", distance="PS3W"),error=function(e)NULL)
    
    # save within-site matches for later analysis #
    m.data2.j <- tryCatch(match.data(m.out1, weights="PSW2", distance="PS2"),error=function(e)NULL)
    m.data2 <- rbind(m.data2,m.data2.j)
    
    # select unmatched treatment units in school j #
    t.m1 <- tryCatch(match.data(m.out1, weights="PSW3W", distance="PS3W", group="treat"),error=function(e)NULL)
    t.m1 <- t.m1[ ,c(teacher_id,exact_vars,"PSW3W")]
    t.x <- t; t.x$PSW3W <- NA # placeholder for schools with no within-school matches
    t2 <- tryCatch(merge(t, t.m1, by=c(teacher_id,exact_vars), all.x=TRUE),error=function(e)t.x)
    t2$PSW3W <- ifelse(is.na(t2$PSW3W),0,1)
    t2 <- t2[t2$PSW3W==0,]
    t2$pscore.j <- t2$pscore
    
    # set data for between-school match #
    # re-calculate control unit propensity score based on school j pscore model parameter estimates #
    c2$PSW3W <- 0
    
    if(is.null(site_vars)){
      tfmla <- NULL
      for(i in teacher_vars){
        tfmla.i <- paste0(pspar.j[,i],"*c2$",i)
        tfmla <- paste(tfmla,tfmla.i,sep=" + ")
      }
      
      newfmla <- paste(pspar.j[,"(Intercept)"],tfmla,"+",pspar.j[,"U0"])
    }else{
      tfmla <- NULL
      for(i in teacher_vars){
        tfmla.i <- paste0(pspar.j[,i],"*c2$",i)
        tfmla <- paste(tfmla,tfmla.i,sep=" + ")
      }
      
      sfmla <- NULL
      for(i in site_vars){
        sfmla.i <- paste0(pspar.j[,i],"*c2$",i)
        sfmla <- paste(sfmla,sfmla.i,sep=" + ")
      }
      newfmla <- paste(pspar.j[,"(Intercept)"],tfmla,sfmla,"+",pspar.j[,"U0"])
    }
    
    c2$pscore.j <- eval(parse(text=newfmla))
    
    tmp2 <- rbind(t2,c2)
    
    # match 2: between-school: only if control units fall within treatment unit range #
    m.out2 <- tryCatch(matchit(tmp2[,tx_var] ~ tmp2$pscore.j, data = tmp2, method = "nearest", replace = replacement, exact = exact_vars, caliper=crc, distance=tmp2$pscore.j,ratio=ratio),error=function(e)NULL)
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
  
  
  match_dat <- output$data
  
  if(by_student == TRUE){
    
  # Create weights for comparison students selected more than once #
    for(i in unique(match_dat$student_id)){
      match_dat[which(match_dat$student_id == i),"N"] <- nrow(match_dat[which(match_dat$student_id == i),])
    }
    match_dat$weights <- 1/match_dat$N
    
    # Drop the duplicate matched cases #
    match_dat <- match_dat[!duplicated(match_dat[c("student_id")]),]
    

    match_dat <- semi_join(student_dat, match_dat, by = c("student_id")) %>%
                 left_join(match_dat %>%
                      select(student_id, p3, pscore, weights), by = c("student_id"))
  
  
  
  
  } else if(by_student == FALSE){
    
    
    # Create weights for comparison teachers selected more than once #
    for(i in unique(match_dat$teacher_id)){
      match_dat[which(match_dat$teacher_id==i),"N"] <- nrow(match_dat[which(match_dat$teacher_id==i),])
    }
    match_dat$weights <- 1/match_dat$N
    
    # Drop the duplicate matched cases #
    match_dat <- match_dat[!duplicated(match_dat[c("teacher_id")]),]
    
    match_dat <- semi_join(student_dat, match_dat, by = c("teacher_id")) %>%
          left_join(match_dat %>%
                  select(teacher_id, p3, pscore, weights), by = c("teacher_id"))
    
  
  
  
  }
  

  
  return(match_dat)
  
}

# calculate balance -------------------------------------------------------


#Check balance Function #
calc_balance <- function(dat_match, 
                         tx_var = "D", 
                         vars = c("Z_k", "W_jk", "X_jk", "X_ijk", "U_ijk")){
  
  dat_match <- as.data.frame(dat_match)
  
  # write weighted SD function #
  weight_sd.fxn <- function(x,wt) sqrt(sum(wt*(x - weighted.mean(x, wt))^2)/(length(x)-1))
  
  # matched sample #
  table.trt_match <- rbind()
  table.ctrl_match <- rbind()
  table.sd_match <- rbind()
  
  for(i in vars){
    trt_match <- dat_match[dat_match[,tx_var]==1,]
    ctrl_match <- dat_match[dat_match[,tx_var]==0,]
    
    table.trt_match <- rbind(table.trt_match,c(i,weighted.mean(trt_match[,i],trt_match$weights)))
    table.ctrl_match <- rbind(table.ctrl_match,c(i,weighted.mean(ctrl_match[,i],ctrl_match$weights)))
    table.sd_match <- rbind(table.sd_match,c(i,weight_sd.fxn(dat_match[,i],dat_match$weights)))
  }
  table.match <- cbind(table.trt_match,table.ctrl_match[,2],table.sd_match[,2])
  colnames(table.match) <- c("Variable","Mean_trt","Mean_ctrl","SD_all")
  table.match <- data.frame(table.match)
  
  table.match$SMD <- (as.numeric(as.character(table.match$Mean_trt)) - as.numeric(as.character(table.match$Mean_ctrl)))/as.numeric(as.character(table.match$SD_all))
  

  # cutting out other stuff 
  res <- table.match %>%
    select(var = Variable,
           smd = SMD) %>%
    spread(var, smd)
  
  return(res)
  
}






# outcome_modeling --------------------------------------------------------


estimate_effect <- function(matched_dat,
                            method) {
  
  if(sum(matched_dat$weights == 0)){
    
    results <- data.frame(method = method,
                          U_ijk = NA,
                          W_jk = NA,
                          X_ijk = NA,
                          X_jk = NA,
                          Z_k = NA,
                          estimate = NA,
                          std_error = NA,
                          statistic = NA, 
                          df = NA,
                          p_value = NA,
                          ci_low = NA, 
                          ci_high = NA)
  } else{
 
  #Run model
  out_mod_1 <- lmer(Y_ijk ~ D  + X_ijk + X_jk + W_jk + Z_k + 
                      (1 | teacher_id) + (1 | school_id),
                    data = matched_dat)
  
  bal_res <- calc_balance(matched_dat)
  
  #Store results
  results <- tidy(out_mod_1) %>%
    filter(term == "D") %>% 
    mutate(method = method) %>%
    clean_names()

  ci <- ci(out_mod_1) %>%
    clean_names() %>%
    filter(parameter == "D") %>%
    select(term = parameter, ci_low, ci_high)
  
  results <- bind_cols(results, bal_res) %>%
    left_join(ci, by = "term") %>%
    select(method, U_ijk:Z_k, estimate:p_value, ci_low, ci_high)
  
  }
  
  return(results)
}


