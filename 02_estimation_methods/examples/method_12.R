
match_within <- function(hold)


# Identify Clusters NOT in a Within-Site Pair
wsp <- hold %>% 
    ungroup() %>%
    select(teacher_id) %>% 
    mutate(inwsp = 1)


dat_m2 <- dat %>% 
  left_join(wsp, by = "teacher_id") %>% 
  filter(is.na(inwsp) == TRUE)

# restrict matching to groups that have T & C clusters

dat_m2 <- dat_m2 %>% 
  group_by(Z_q5) %>% 
  mutate(tmpm = mean(D)) %>% 
  ungroup() %>%
  filter(tmpm > 0 & tmpm < 1) %>% 
  select(-tmpm)

dat_m2 %>%
  group_by(Z_q5) %>%
  

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
