source("test/dgm_params.R")

system.time(
  check <- run_sim(iterations = 2,
                   k = k,
                   j = j,
                   i = i,
                   icc3 = icc3,
                   icc2 = icc2,
                   R2 = R2,
                   ps_coef = ps_coef,
                   pr_star = pr_star,
                   outcome_coef = outcome_coef,
                   delta = delta,
                   seed = 20220726)
)

#163.818 


# checking perf criteria --------------------------------------------------

res <- pmap_dfr(matched_sets, estimate_effect) %>%
  mutate(true_effect = delta)


##check_res <- bind_rows(res, res)

perfm_res <- calc_performance(res)
