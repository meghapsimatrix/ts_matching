################################# Define check balance function ################################

# Check balance Function #
fun.balance <- function(dat_full, dat_match, tx.var, vars){
  
  # write weighted SD function #
  weight_sd.fxn <- function(x,wt) sqrt(sum(wt*(x - weighted.mean(x, wt))^2)/(length(x)-1))

  # full sample #
  table.trt_full <- rbind()
  table.ctrl_full <- rbind()
  table.sd_full <- rbind()
  
  for(i in vars){
    trt_full <- dat_full[dat_full[,tx.var]==1,]
    ctrl_full <- dat_full[dat_full[,tx.var]==0,]
    
    table.trt_full <- rbind(table.trt_full,c(i,mean(trt_full[,i],na.rm=TRUE)))
    table.ctrl_full <- rbind(table.ctrl_full,c(i,mean(ctrl_full[,i],na.rm=TRUE)))
    table.sd_full <- rbind(table.sd_full,c(i,sd(dat_full[,i],na.rm=TRUE)))
  }
  table.full <- cbind(table.trt_full,table.ctrl_full[,2],table.sd_full[,2])
  colnames(table.full) <- c("Variable","Mean_trt","Mean_ctrl","SD_all")
  table.full <- data.frame(table.full)
  
  table.full$SMD <- (as.numeric(as.character(table.full$Mean_trt)) - as.numeric(as.character(table.full$Mean_ctrl)))/as.numeric(as.character(table.full$SD_all))
  
  # matched sample #
  table.trt_match <- rbind()
  table.ctrl_match <- rbind()
  table.sd_match <- rbind()
  
  for(i in vars){
    trt_match <- dat_match[dat_match[,tx.var]==1,]
    ctrl_match <- dat_match[dat_match[,tx.var]==0,]
      
    table.trt_match <- rbind(table.trt_match,c(i,weighted.mean(trt_match[,i],trt_match$weights)))
    table.ctrl_match <- rbind(table.ctrl_match,c(i,weighted.mean(ctrl_match[,i],ctrl_match$weights)))
    table.sd_match <- rbind(table.sd_match,c(i,weight_sd.fxn(dat_match[,i],dat_match$weights)))
  }
  table.match <- cbind(table.trt_match,table.ctrl_match[,2],table.sd_match[,2])
  colnames(table.match) <- c("Variable","Mean_trt","Mean_ctrl","SD_all")
  table.match <- data.frame(table.match)
  
  table.match$SMD <- (as.numeric(as.character(table.match$Mean_trt)) - as.numeric(as.character(table.match$Mean_ctrl)))/as.numeric(as.character(table.match$SD_all))
  
  # combine tables #
  table.all <- cbind(table.full,table.match[,2:ncol(table.match)])
  colnames(table.all)[2:5] <- paste(colnames(table.all)[2:5],"Full",sep="_")
  colnames(table.all)[6:9] <- paste(colnames(table.all)[6:9],"Match",sep="_")
  
  return(table.all)
}
