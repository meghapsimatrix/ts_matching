library(tidyverse)

load("results/simulation_results_mj_1.RData")
results_mj_1 <- results

load("results/simulation_results_mc1.RData")
results_mc_1 <- results

results_clean <- bind_rows(results_mj_1, results_mc_1) %>%
  group_by(method, k, j, icc3, icc2) %>%
  summarize_at(vars(bias:prop_t_stud_m), mean) %>%
  ungroup() %>%
  mutate(K = 400)


results_clean <- 
  results_clean %>%
  mutate(method = as.character(method)) %>%
  mutate(method = factor(method, levels = c("0.1", "0.2", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))) %>%
  mutate(k_j = paste("k = ", k, ",", "j = ", j),
         icc = paste("icc3 = ", icc3, ",", "icc2 = ", icc2))


ggplot(results_clean, aes(x = method, y = bias, color = method)) + 
  geom_point() +
  facet_grid(k_j ~ icc) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  guides(color = "none")

ggsave("results/prelim_graph_bias.png", device = "png", width = 12, height = 8)
