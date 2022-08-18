library(tidyverse)

# Load results ------------------------------------------------------------

files <- list.files("results/results", full.names = TRUE)
files

load_res <- function(file) {

  load(file)
  results

}

results <- map_dfr(files, load_res)


# clean results  ----------------------------------------------------------

K <- 200 * length(files)

results_clean <- 
  results %>%
  group_by(method, k, j, icc3, icc2) %>%
  summarize_at(vars(bias:prop_t_stud_m), mean) %>%
  ungroup() %>%
  mutate(K = K)


results_clean <- 
  results_clean %>%
  mutate(method = as.character(method)) %>%
  mutate(method = factor(method, levels = c("0.1", "0.2", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))) %>%
  mutate(k_j = paste("k = ", k, ",", "j = ", j),
         icc = paste("icc3 = ", icc3, ",", "icc2 = ", icc2))



# preliminary graphs ------------------------------------------------------

# bias --------------------------------------------------------------------

ggplot(results_clean, aes(x = method, y = bias, color = method)) + 
  geom_point() +
  facet_grid(k_j ~ icc) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  guides(color = "none")

ggsave("results/prelim_graph_bias.png", device = "png", width = 12, height = 8)


# rmse --------------------------------------------------------------------

ggplot(results_clean, aes(x = method, y = rmse, color = method)) + 
  geom_point() +
  facet_grid(k_j ~ icc) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  guides(color = "none")

ggsave("results/prelim_graph_mcse.png", device = "png", width = 12, height = 8)

# proportion  ---------------------------------------------------------------

ggplot(results_clean, aes(x = method, y = prop_t_m, color = method)) + 
  geom_point() +
  facet_grid(k_j ~ icc) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  guides(color = "none")

ggsave("results/prelim_graph_prop_t_m.png", device = "png", width = 12, height = 8)

# proportion  ---------------------------------------------------------------

ggplot(results_clean, aes(x = method, y = prop_t_stud_m, color = method)) + 
  geom_point() +
  facet_grid(k_j ~ icc) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  guides(color = "none")

ggsave("results/prelim_graph_prop_t_stud_m.png", device = "png", width = 12, height = 8)


# balance -----------------------------------------------------------------

# need to figure out how to graph this better

bal_res <- 
  results_clean %>%
  select(method, k_j, icc, U_ijk:Z_k) %>%
  gather(var, smd, U_ijk:Z_k)

bal_res %>%
  ggplot(aes(y = var, x = smd, color = method)) +
  geom_point() +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = c(-.1, .1), linetype = "dashed") +
  scale_y_discrete(limits = rev(levels(var))) +
  facet_grid(k_j ~ icc) +
  labs(x = "Standardized Mean Differences", y = "") +
  ggtitle("Covariate Balance") + 
  theme_bw()

