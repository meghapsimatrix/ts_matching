################################### Define matching function ###################################

# Matching Function #
fun.match <- function(dat, 
                      SiteID, 
                      TeacherID, tx.var, exact.vars, teacher.vars, site.vars, group, crc, replacement, ratio, seed){
  
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
