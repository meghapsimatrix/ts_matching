library(tidyverse)

load("results/results/results_clean.RData")

to_table <- results_clean %>%
  filter(icc2 == .1, icc3 == .1) %>%
  select(1:13) %>%
  select(-contains("icc"), -contains("mcse"), -width) %>%
  mutate_at(vars(bias, rmse, coverage), round, 2)

bias_res <- to_table %>%
  mutate(condition = paste0("k = ", k, ", ", "j = ", j)) %>%
  select(-k, -j) %>%
  gather(criteria, val, bias:coverage) %>%
  filter(criteria == "bias") %>%
  mutate(criteria = paste(condition, criteria)) %>%
  select(- condition) %>%
  spread(criteria, val)


write_csv(bias_res, "results/bias_res.csv")

to_table <- results_clean %>%
  filter(icc2 == .1, icc3 == .2) %>%
  select(1:13) %>%
  select(-contains("icc"), -contains("mcse"), -width) %>%
  mutate_at(vars(bias, rmse, coverage), round, 2)

bias_res <- to_table %>%
  mutate(condition = paste0("k = ", k, ", ", "j = ", j)) %>%
  select(-k, -j) %>%
  gather(criteria, val, bias:coverage) %>%
  filter(criteria == "bias") %>%
  mutate(criteria = paste(condition, criteria)) %>%
  select(- condition) %>%
  spread(criteria, val)


write_csv(bias_res, "results/bias_res_1_2.csv")


to_table <- results_clean %>%
  filter(icc2 == .2, icc3 == .1) %>%
  select(1:13) %>%
  select(-contains("icc"), -contains("mcse"), -width) %>%
  mutate_at(vars(bias, rmse, coverage), round, 2)

bias_res <- to_table %>%
  mutate(condition = paste0("k = ", k, ", ", "j = ", j)) %>%
  select(-k, -j) %>%
  gather(criteria, val, bias:coverage) %>%
  filter(criteria == "bias") %>%
  mutate(criteria = paste(condition, criteria)) %>%
  select(- condition) %>%
  spread(criteria, val)


write_csv(bias_res, "results/bias_res_2_1.csv")

to_table <- results_clean %>%
  filter(icc2 == .2, icc3 == .2) %>%
  select(1:13) %>%
  select(-contains("icc"), -contains("mcse"), -width) %>%
  mutate_at(vars(bias, rmse, coverage), round, 2)

bias_res <- to_table %>%
  mutate(condition = paste0("k = ", k, ", ", "j = ", j)) %>%
  select(-k, -j) %>%
  gather(criteria, val, bias:coverage) %>%
  filter(criteria == "bias") %>%
  mutate(criteria = paste(condition, criteria)) %>%
  select(- condition) %>%
  spread(criteria, val)


write_csv(bias_res, "results/bias_res.csv")



# for the x  --------------------------------------------------------------


to_table <- results_clean %>%
  filter(icc2 == .1, icc3 == .1) %>%
  select(1:5, bias, X_ijk) %>%
  select(-contains("icc")) %>%
  mutate_at(vars(bias, X_ijk), round, 2)

smd_res <- to_table %>%
  mutate(condition = paste0("k = ", k, ", ", "j = ", j)) %>%
  select(-k, -j) %>%
  gather(criteria, val, bias:X_ijk) %>%
  filter(criteria == "X_ijk") %>%
  mutate(criteria = paste(condition, criteria)) %>%
  select(- condition) %>%
  spread(criteria, val)


write_csv(smd_res, "results/smd_res_1_1.csv")

to_table <- results_clean %>%
  filter(icc2 == .1, icc3 == .2) %>%
  select(1:5, bias, X_ijk) %>%
  select(-contains("icc")) %>%
  mutate_at(vars(bias, X_ijk), round, 2)

smd_res <- to_table %>%
  mutate(condition = paste0("k = ", k, ", ", "j = ", j)) %>%
  select(-k, -j) %>%
  gather(criteria, val, bias:X_ijk) %>%
  filter(criteria == "X_ijk") %>%
  mutate(criteria = paste(condition, criteria)) %>%
  select(- condition) %>%
  spread(criteria, val)


write_csv(smd_res, "results/smd_res_1_2.csv")

to_table <- results_clean %>%
  filter(icc2 == .2, icc3 == .1) %>%
  select(1:5, bias, X_ijk) %>%
  select(-contains("icc")) %>%
  mutate_at(vars(bias, X_ijk), round, 2)

smd_res <- to_table %>%
  mutate(condition = paste0("k = ", k, ", ", "j = ", j)) %>%
  select(-k, -j) %>%
  gather(criteria, val, bias:X_ijk) %>%
  filter(criteria == "X_ijk") %>%
  mutate(criteria = paste(condition, criteria)) %>%
  select(- condition) %>%
  spread(criteria, val)


write_csv(smd_res, "results/smd_res_2_1.csv")



to_table <- results_clean %>%
  filter(icc2 == .2, icc3 == .2) %>%
  select(1:5, bias, X_ijk) %>%
  select(-contains("icc")) %>%
  mutate_at(vars(bias, X_ijk), round, 2)

smd_res <- to_table %>%
  mutate(condition = paste0("k = ", k, ", ", "j = ", j)) %>%
  select(-k, -j) %>%
  gather(criteria, val, bias:X_ijk) %>%
  filter(criteria == "X_ijk") %>%
  mutate(criteria = paste(condition, criteria)) %>%
  select(- condition) %>%
  spread(criteria, val)


write_csv(smd_res, "results/smd_res_2_2.csv")

