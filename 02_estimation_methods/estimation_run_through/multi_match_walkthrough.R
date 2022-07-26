m_3 <- multi_match(dat = dat,
                   trt = "D",
                   l1_cov = c("X_ijk"),
                   l2_cov = c("W_q5", "Z_q5"),
                   l2_id = "teacher_id",
                   caliper = .25)

l3_id <- "school_id"


dat$tmp <- dat[, trt]
dat$l3id <- dat[, l3_id]

dat_m <- dat %>%
  group_by(school_id) %>%
  mutate(tmpm = mean(D)) %>%
  ungroup() %>%
  filter(tmpm > 0 & tmpm < 1) %>%
  select(-tmpm)


m_6 <- dat_m %>%
  group_by(school_id) %>%
  do(multi_match(., 
                 trt = "D",
                 l1_cov = c("X_ijk"),
                 l2_cov = NULL,
                 l2_id = "teacher_id",
                 l3_id = "school_id",
                 caliper = 1,
                 add_id = FALSE))


m_9 <- dat_m %>%
  group_by(Z_q5) %>%
  do(multi_match(., 
                 trt = "D",
                 l1_cov = c("X_ijk"),
                 l2_cov = c("W_q5", "Z_q5"),
                 l2_id = "teacher_id",
                 l3_id = "Z_q5",
                 caliper = 1,
                 add_id = TRUE))


m_12_hold <- 
  dat_m %>%
  group_by(Z_q5) %>%
  do(multi_match(., 
                 trt = "D",
                 l1_cov = c("X_ijk"),
                 l2_cov = NULL,
                 l2_id = "teacher_id",
                 l3_id = "Z_q5",
                 caliper = 1,
                 add_id = TRUE))


# Identify Clusters NOT in a Within-Site Pair
wsp <- m_12_hold %>% 
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


m_12_hold_2 <- dat_m2 %>%
  mutate(Z_q5 = as.character(Z_q5)) %>%
  group_by(Z_q5) %>%
  do(multi_match(., 
                 trt = "D",
                 l1_cov = c("X_ijk"),
                 l2_cov = NULL,
                 l2_id = "teacher_id",
                 l3_id = "Z_q5",
                 caliper = 1,
                 add_id = TRUE))
