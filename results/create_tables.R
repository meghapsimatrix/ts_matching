load("results/results/results_clean.RData")

to_table <- results_clean %>%
  filter(icc2 == .1, icc3 == .1) %>%
  select(1:13) %>%
  select(-contains("icc"), -contains("mcse"), -width) %>%
  mutate_at(vars(bias, rmse, coverage), round, 3)

bias_res <- to_table %>%
  mutate(condition = paste0("k = ", k, ", ", "j = ", j)) %>%
  select(-k, -j) %>%
  gather(criteria, val, bias:coverage) %>%
  filter(criteria == "bias") %>%
  mutate(criteria = paste(condition, criteria)) %>%
  select(- condition) %>%
  spread(criteria, val)


write_csv(bias_res, "results/bias_res.csv")
