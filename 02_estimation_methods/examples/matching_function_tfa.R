# Matching Function #
fun.match <- function(dat, 
                      SiteID, 
                      TeacherID, 
                      tx.var, 
                      exact.vars, 
                      teacher.vars, 
                      site.vars = NULL, 
                      score, # what is this?
                      crc, # caliper
                      replacement, # with or without replacement
                      ratio = 1,
                      seed = 1234){
  
  set.seed(seed)
  
  if(is.null(site.vars)){
    df <- as.data.frame(dat[,c(SiteID, TeacherID, tx.var, 
                               exact.vars, teacher.vars, score)])
  
    } else{
      
    df <- as.data.frame(dat[,c(SiteID, TeacherID, tx.var, exact.vars, 
                               teacher.vars, site.vars)]) # why is score not here?
    
  }
  
  if(length(unique(dat[,SiteID]))<12){
    
    df$site_cluster <- 1
    
  } else if(length(unique(dat[,SiteID]))>=12 & length(unique(dat[,SiteID]))<20){
    cluster_q <- quantile(df[,score], probs = seq(0, 1, 1/3))
    df$site_cluster <- ifelse(df[,score] <= cluster_q[2], 1,
                              ifelse(df[,score] > cluster_q[2] & df[,score] <= cluster_q[3],2,
                                     ifelse(df[,score] > cluster_q[3],3,NA)))
  }else{
    # define cluster sites based on subject score #
    cluster_q <- quantile(df[,score], probs = seq(0, 1, .2))
    df$site_cluster <- ifelse(df[,score] <= cluster_q[2],1,
                              ifelse(df[,score] > cluster_q[2] & df[,score] <= cluster_q[3],2,
                                     ifelse(df[,score] > cluster_q[3] & df[,score] <= cluster_q[4],3,
                                            ifelse(df[,score] > cluster_q[4] & df[,score] <= cluster_q[5],4,
                                                   ifelse(df[,score] > cluster_q[5],5,NA)))))
  }
  
  #cluster_q <- quantile(df[,score], probs = seq(0, 1, .25))
  #df$site_cluster <- ifelse(df[,score] <= cluster_q[2],1,
  #                          ifelse(df[,score] > cluster_q[2] & df[,score] <= cluster_q[3],2,
  #                                ifelse(df[,score] > cluster_q[3] & df[,score] <= cluster_q[4],3,
  #                                        ifelse(df[,score] > cluster_q[4] & df[,score] <= cluster_q[5],4,NA))))
  
  #cluster_scores <- df %>% dplyr::group_by(schoolid,sy) %>% dplyr::summarize_at(vars(score),unique) %>% data.frame
  #cluster_scores <- df %>% dplyr::group_by(schoolid) %>% 
  #  dplyr::summarise(tx.grps=case_when(any(tfa_flag_math==0)~'Ctr',
  #                                     TRUE~'Tx'),
  #                   score=unique(!!sym(score))) %>% data.frame
  
  #cluster_scores[,'site_cluster'] <- scclust::sc_clustering(distances::distances(cluster_scores$score),2,
  #                                                          type_labels=cluster_scores$tx.grps,
  #                                                          type_constraints = c('Ctr'=1))
  #cluster_scores <- df %>% dplyr::group_by(schoolid) %>% 
  #  dplyr::summarise(score=unique(!!sym(score))) %>% data.frame
  
  #cluster_scores[,'site_cluster'] <- scclust::sc_clustering(distances::distances(cluster_scores$score),2)
  #df <- dplyr::left_join(df,cluster_scores[,c('schoolid','site_cluster')],by=c('schoolid'))
  
  
  
  # multilevel formula for ps modeling 
  # treatment ~ teacher level variables and site level variables 
  # random intercept by site
  
  if(is.null(site.vars)){
    pfmla <- as.formula(paste(tx.var," ~ 1 + ",paste(teacher.vars,collapse="+")," + (1|",SiteID,")"))
  }else{
    # estimate propensity score with two-level RE model #
    pfmla <- as.formula(paste(tx.var," ~ 1 + ",paste(teacher.vars,collapse="+"),"+",paste(site.vars,collapse="+")," + (1|",SiteID,")"))
    pfmla <- as.formula(paste(c(paste(tx.var,"~1"),teacher.vars,site.vars,paste("(1|",SiteID,")")),collapse='+'))
  }
  #pfmla <- as.formula(paste0(c(tx.var,"~1",teacher.vars,site.vars,paste("(1|",SiteID,")")),collapse='+'))
  #p3 <- glmer(pfmla, data=df, family=binomial(link="logit"),control=glmerControl(optimizer='bobyqa',optCtrl = list(maxfun = 2e4)))
  
  # use bayesian glmer
  # and p3 - take out propensity scores and add it to the data
  
  p3 <- rstanarm::stan_glmer(pfmla, data=df, family=binomial(link="logit"))
  df$p3 <- log(fitted(p3)/(1-fitted(p3)))
  #df$p3 <- colMeans(rstanarm::posterior_linpred(p3))
  
  # save school-level EB coefficient estimates for between-site matching #
  
  #pspar3a <- coef(p3)[[SiteID]]
  #pspar3a$U0 <- pull(ranef(p3)[[SiteID]])
  #pspar3a$`(Intercept)` <- pspar3a$`(Intercept)`-pspar3a$U0
  #pspar3a[,SiteID] <- rownames(ranef(p3)[[SiteID]])
  
  
  # NOT SURE WHAT THE FOLLOWING IS DOING
  
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
  # site id - that are treatment 
  # unique and sorted
  
  sids <- sort(unique(df[df[,tx.var]==1,SiteID]))
  
  # list glm output is the ps p3
  # data used 
  
  output <- list()
  output[['glm output']] <- p3
  output[['data']] <- df
  
  
  # loop over treatment schools #
  
  # for each of the treatment schools
  
  for(j in 1:length(sids)) {
    #j <-4
    # set data #
    sid <- sids[j] # school id 
    t <- df[df[,SiteID]==sid & df[,tx.var]==1,] # use treatment units in school j
    c1 <- df[df[,SiteID]==sid & df[,tx.var]==0,] # use control units in school j
    
    
    if(any(df[,SiteID]!=sid & 
           df$site_cluster %in% t$site_cluster &
           df[,tx.var] == 0)) {
      
      c2 <- df[df[, SiteID]!=sid & 
                 df$site_cluster %in% t$site_cluster & 
                 df[,tx.var]==0,]
   
    }
    
    else{
      
      c2 <- df[df[,SiteID]!=sid & 
                 df$site_cluster %in% t$site_cluster,][1,]
      
      c2[,2:dim(df)[2]] <- NA
      
    }		
    
    tmp1 <- rbind(t, c1)  # treatment and control data? 
    
    pspar.j <- as.vector(psparx[psparx[,SiteID]==j,]) # use ps model parameter estimates for school j
    
    # match 1: within-school #
    m.out1 <- tryCatch(matchit(tmp1[,tx.var] ~ tmp1$pscore, # treatment ~ pscore
                               data = tmp1, 
                               method = "nearest", 
                               replace = replacement, 
                               exact = exact.vars, 
                               caliper = crc, 
                               distance = tmp1$pscore,
                               ratio = ratio),
                          error=function(e)NULL)
    
    output[['Within School Matches']][[j]] <- m.out1
    
    # create data frame of matched units from match 1 # 
    m.data3w <- tryCatch(match.data(m.out1, 
                                    weights = "PSW3W", 
                                    distance = "PS3W"),
                         error=function(e)NULL)
    
    # save within-site matches for later analysis #
    m.data2.j <- tryCatch(match.data(m.out1, 
                                     weights = "PSW2", 
                                     distance = "PS2"),
                          error=function(e)NULL)
    
    m.data2 <- rbind(m.data2,m.data2.j)
    
    # select unmatched treatment units in school j #
    t.m1 <- tryCatch(match.data(m.out1, 
                                weights="PSW3W", 
                                distance="PS3W", 
                                group="treat"),
                     error=function(e)NULL)
    
    t.m1 <- t.m1[ ,c(TeacherID,exact.vars,"PSW3W")]
    
    t.x <- t; t.x$PSW3W <- NA # placeholder for schools with no within-school matches
    
    t2 <- tryCatch(merge(t, t.m1, by=c(TeacherID,exact.vars), 
                         all.x = TRUE),
                   error=function(e)t.x)
    
    t2$PSW3W <- ifelse(is.na(t2$PSW3W),0,1)
    t2 <- t2[t2$PSW3W == 0,]
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
    
    #c2$pscore.j <- as.numeric(c(pspar.j[c(2:10,12)])) %*% t(as.matrix(cbind(1,c2[,c(teacher.vars,site.vars)],1)))
    
    tmp2 <- rbind(t2,c2)
    
    # match 2: between-school: only if control units fall within treatment unit range #
    m.out2 <- tryCatch(matchit(tmp2[,tx.var] ~ tmp2$pscore.j, 
                               data = tmp2, 
                               method = "nearest", 
                               replace = replacement,
                               exact = exact.vars, 
                               caliper=crc, 
                               distance=tmp2$pscore.j,
                               ratio=ratio),
                       error=function(e)NULL)
    output[['Bewteen School Matches']][[j]] <- m.out2
    
    # create data frame of matched units from match 2 #
    m.data3b <- tryCatch(match.data(m.out2, weights="PSW3B", distance="PS3B"),error=function(e)NULL)
    if (!is.null(m.data3b)) m.data3b$PS3W <- NA else m.data3b <- NULL
    
    # combine matched files #
    if (!is.null(m.data3w)) m.data3w$PSW3B  <- 0 else m.data3w <- NULL
    if (!is.null(m.data3w)) m.data3w$PS3B <- NA else m.data3w <- NULL
    if (!is.null(m.data3w)) m.data3w$pscore.j <- m.data3w$pscore else m.data3w <- NULL
    m.datax <- rbind(m.data3w, m.data3b)
    if (!is.null(m.datax)) m.datax$Site_match <- sid else m.datax <- NULL
    
    # append to master file #
    m.data3 <- rbind(m.data3, m.datax)
    
  } # close 2 stage matching loop
  output[['Matched Data']] <- m.data3
  
  return(output)
}




####################################### Math #######################################

# 1) Define the vector of teacher-level variables.
teacher.tchr_vars_math <- c("tlm_race_eth_min","tlm_sex_f","tlm_sped_y","tlm_ell_y",'tlm_z_math_py')

# 2) Define the vector of school-level variables.
site.tchr_vars_math <- NULL

# 3) Define variables on which to exact match.
exact.tchr_vars_math <- c("sy","tchr_exp_band_math","grade_band")

# 4) Read in the data.
dat_tch_math <- read.csv("./Data/CCSD teacher data (math).csv")

dat_tch_math %>%
  mutate(teacherid1_math = as.character(teacherid1_math)) %>%
  filter(is.na(sl_z_math_py)) 

# 6) Ensure there is no missingness on the matching variables.
for(i in teacher.tchr_vars_math) print(psych::describe(dat_tch_math[,i]))
for(i in site.tchr_vars_math) print(summary(dat_tch_math[,i]))
for(i in exact.tchr_vars_math) print(summary(dat_tch_math[,i]))

Hmisc::describe(dat_tch_math[,c("schoolid","teacherid1_math","tfa_flag_math",
                                exact.tchr_vars_math,
                                teacher.tchr_vars_math,
                                site.tchr_vars_math,
                                "sl_z_math_py")])

# create a summary table of the number of TFA and non-TFA teachers within each exact matching group
summary.exact.matches_math <- dat_tch_math %>% 
  group_by(sy,tchr_exp_band_math,grade_band,tfa_flag_math) %>% 
  summarise(N=n(),
            Mean=mean(tlm_z_math_py),
            Matches=paste0(tlm_z_math_py,collapse=', '))
write.csv(summary.exact.matches_math,file="./Data/Matched data/Exact match summary (math).csv",row.names=FALSE)

table(is.na(dat_tch_math$sl_z_math_py))

dat_tch_math <- dat_tch_math %>%
  filter(!is.na(sl_z_math_py))

# 7) Conduct matching.
# Conduct matching separately for each academic year (2018, 2019), caliper 4
dat_tch_math18 <- dat_tch_math %>% filter(sy==2018)




dat_math_output18 <- fun.match(dat=dat_tch_math18, 
                               SiteID="schoolid", 
                               TeacherID="teacherid1_math", 
                               tx.var="tfa_flag_math", 
                               exact.vars=exact.tchr_vars_math, 
                               teacher.vars=teacher.tchr_vars_math, 
                               site.vars=site.tchr_vars_math,
                               score="sl_z_math_py" , 
                               crc=4,
                               ratio=1,
                               replacement=TRUE)
dat_math18_matched <- dat_math_output18$`Matched Data`


dat_tch_math19 <- dat_tch_math %>% filter(sy==2019)

dat_math_output19 <- fun.match(dat=dat_tch_math19, 
                               SiteID="schoolid", 
                               TeacherID="teacherid1_math", 
                               tx.var="tfa_flag_math", 
                               exact.vars=exact.tchr_vars_math, 
                               teacher.vars=teacher.tchr_vars_math, 
                               #site.vars=site.tchr_vars_math,
                               score="sl_z_math_py" , 
                               crc=4,
                               ratio=1,
                               replacement=TRUE)
dat_math19_matched <- dat_math_output19$`Matched Data`


# 8) Save the matched data.
write.csv(bind_rows(dat_math18_matched,dat_math19_matched),file="Data/Matched data/CCSD matched teacher data (math).csv",row.names=FALSE)
write.csv(bind_rows(type.convert(dat_math_output18$data),type.convert(dat_math_output19$data)),"./Data/CCSD teacher data (math) pscore.csv",row.names = F)



####################################### ELA #######################################


# 1) Define the vector of teacher-level variables.
teacher.tchr_vars_ela <- c("tle_race_eth_min","tle_sex_f","tle_sped_y","tle_ell_y",'tle_z_ela_py')

# 2) Define the vector of school-level variables.
site.tchr_vars_ela <- NULL

# 3) Define variables on which to exact match.
exact.tchr_vars_ela <- c("sy","tchr_exp_band_ela","grade_band")

# 4) Read in the data.
dat_tch_ela <- read.csv("Data/CCSD teacher data (ela).csv")

table(is.na(dat_tch_ela$sl_z_ela_py))

dat_tch_ela <- dat_tch_ela %>%
  filter(!is.na(sl_z_ela_py), !is.na(schoolid)) 

# 6) Ensure there is no missingness on the matching variables.
for(i in teacher.tchr_vars_ela) print(summary(dat_tch_ela[,i]))
for(i in site.tchr_vars_ela) print(summary(dat_tch_ela[,i]))
for(i in exact.tchr_vars_ela) print(summary(dat_tch_ela[,i]))
describe(dat_tch_ela[,c("schoolid","teacherid1_ela","tfa_flag_ela",exact.tchr_vars_ela,teacher.tchr_vars_ela,site.tchr_vars_ela,"sl_z_ela_py")])

# create a summary table of the number of TFA and non-TFA teachers within each exact matching group
summary.exact.matches_ela <- dat_tch_ela %>% 
  group_by(sy,tchr_exp_band_ela,grade_band,tfa_flag_ela) %>% 
  summarise(N=n(),
            Mean=mean(tle_z_ela_py),
            Matches=paste0(tle_z_ela_py,collapse=', '))


write.csv(summary.exact.matches_ela,file="Data/Matched data/Exact match summary (ela).csv",row.names=FALSE)


# 7) Conduct matching.
# Conduct matching separately for each academic year (2018, 2019), caliper 4.0
dat_tch_ela18 <- dat_tch_ela %>% filter(sy==2018)

dat_tch_ela_output18 <- fun.match(dat=dat_tch_ela18, 
                                  SiteID="schoolid", 
                                  TeacherID="teacherid1_ela", 
                                  tx.var="tfa_flag_ela", 
                                  exact.vars=exact.tchr_vars_ela, 
                                  teacher.vars=teacher.tchr_vars_ela, 
                                  site.vars=site.tchr_vars_ela,
                                  score="sl_z_ela_py" , 
                                  crc=4,
                                  ratio=1,
                                  replacement=TRUE)
dat_tch_ela18_matched <- dat_tch_ela_output18$`Matched Data`
#10 MATCHES



dat_tch_ela19 <- dat_tch_ela %>% filter(sy==2019)

dat_tch_ela_output19 <- fun.match(dat=dat_tch_ela19, 
                                  SiteID="schoolid", 
                                  TeacherID="teacherid1_ela", 
                                  tx.var="tfa_flag_ela", 
                                  exact.vars=exact.tchr_vars_ela, 
                                  teacher.vars=teacher.tchr_vars_ela, 
                                  site.vars=site.tchr_vars_ela,
                                  score="sl_z_ela_py" , 
                                  crc=4,
                                  ratio=1,
                                  replacement=TRUE)
dat_tch_ela19_matched <- dat_tch_ela_output19$`Matched Data`


# 8) Save the matched data.
write.csv(bind_rows(dat_tch_ela18_matched,dat_tch_ela19_matched),file="Data/Matched data/CCSD matched teacher data (ela).csv",row.names=FALSE)
write.csv(bind_rows(type.convert(dat_tch_ela_output18$data),type.convert(dat_tch_ela_output19$data)),"./Data/CCSD teacher data (ela) pscore.csv",row.names = F)

