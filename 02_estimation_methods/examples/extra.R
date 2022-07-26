

library(broom.mixed)



# no distinction units only 
# by default nearest neighbor, no replacement

# just with matchit (not multi-level)
m_out_1 <- matchit(D ~ X_ijk + W_jk + Z_k,
                   caliper = .25,
                   data = example_dat)

m_out_1$match.matrix


balance_stats <- bal.tab(m_out_1)
balance_stats$Balance

m_data_1 <- match.data(m_out_1)
table(m_data_1$weights)  # if we get to a point were we have to use weights - lmer won't work



#Function for running model and storing results

estimate_D <- function(MATCHED_SET,method) {
  #Create X_jk
  Xjk <- MATCHED_SET %>% group_by(teacher_id, school_id) %>% summarise_at(vars(X_ijk), list(X_jk=mean))
  MATCHED_SET <- merge(MATCHED_SET, Xjk, by=c("teacher_id","school_id"))
  
  #Run model
  out_mod_1 <- lmer(Y_ijk ~ D  + X_ijk + X_jk + W_jk + Z_k + 
                      (1 | teacher_id) + (1 | school_id),
                    data = MATCHED_SET)
  #Store results
  results <- tidy(out_mod_1)
  results <- results[2,]
  results <- results %>% mutate(method=method)
  
  return(results)
}


#Example of using the function
XX <- estimate_D(m_data_1,"Nearest Neighbor")








#Old code for saving dataset name

#name_data <- deparse(substitute(m_data_1))
#set2 <- character(1)
#set2[1] <- print(name_data)
#set2 <- as.data.frame(set2)
#YY <- data.frame(XX,set2)

