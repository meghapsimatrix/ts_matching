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
                 l2_cov = c("W_q5"),
                 l2_id = "teacher_id",
                 l3_id = "school_id",
                 caliper = 1,
                 add_id = FALSE))

     