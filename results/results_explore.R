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
  filter(!(method %in% c("0.1", "0.2"))) %>%
  mutate(method = factor(method, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))) %>%
  mutate(k_j = paste0("k = ", k, ", ", "j = ", j),
         icc = paste("icc3 = ", icc3, ", ", "icc2 = ", icc2)) %>%
  mutate(matching_priority = fct_collapse(method,
                                          "No Distinction" = c("1", "2", "3"),
                                          "Within-site Matching" = c("4", "5", "6"),
                                          "Within-group Matching" = c("7", "8", "9"),
                                          "Within-site then -group Matching" = c("10", "11", "12")),
         matching_level = fct_collapse(method,
                                       "Units only" = c("1", "4", "7", "10"),
                                       "Clusters only" = c("2", "5", "8", "11"),
                                       "Units & Clusters" = c("3", "6", "9", "12"))) 



# preliminary graphs ------------------------------------------------------

# bias --------------------------------------------------------------------

# break off into sections
ggplot(results_clean, aes(x = method, 
                          y = bias, 
                          shape = matching_level, 
                          color = matching_priority)) + 
  geom_point() +
  annotate("rect",
           xmin = 0, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           fill = "blue", alpha = .2) +
  annotate("rect",
           xmin = 6.5, xmax = 9.5,
           ymin = -Inf, ymax = Inf,
           fill = "blue", alpha = .2) +
  facet_grid(k_j ~ icc) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() +
  ggtitle("Bias") +
  labs(shape = "", color = "", x = "Method", y = "Bias") +
  theme(legend.position = "bottom")

ggsave("results/graphs/prelim_graph_bias.png", device = "png", width = 12, height = 8)


# check_mcse 
# fairly small compared to bias
# do we want to include mcse in graph somehow? Like ci dot and line plot? - but it might now show bc the mcse is quite small
results_clean %>%
  group_by(method) %>%
  summarize(bias_mcse = max(bias_mcse)) %>%
  arrange(desc(bias_mcse))

# rmse --------------------------------------------------------------------

ggplot(results_clean, aes(x = method, 
                          y = rmse, 
                          shape = matching_level, 
                          color = matching_priority)) + 
  geom_point() +
  annotate("rect",
           xmin = 0, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           fill = "blue", alpha = .2) +
  annotate("rect",
           xmin = 6.5, xmax = 9.5,
           ymin = -Inf, ymax = Inf,
           fill = "blue", alpha = .2) +
  facet_grid(k_j ~ icc) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() +
  ggtitle("RMSE") +
  labs(shape = "", color = "", x = "Method", y = "RMSE") +
  theme(legend.position = "bottom")

ggsave("results/graphs/prelim_graph_rmse.png", device = "png", width = 12, height = 8)

# check_mcse 
# fairly small compared to rmse 
results_clean %>%
  group_by(method) %>%
  summarize(rmse_mcse = max(rmse_mcse)) %>%
  arrange(desc(rmse_mcse))

# proportion  ---------------------------------------------------------------

ggplot(results_clean, aes(x = method, 
                          y = prop_t_m, 
                          shape = matching_level, 
                          color = matching_priority)) + 
  geom_point() +
  annotate("rect",
           xmin = 0, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           fill = "blue", alpha = .2) +
  annotate("rect",
           xmin = 6.5, xmax = 9.5,
           ymin = -Inf, ymax = Inf,
           fill = "blue", alpha = .2) +
  facet_grid(k_j ~ icc) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() +
  ggtitle("Proportion Treatment Teachers Matched") +
  labs(shape = "", color = "", x = "Method", y = "Proportion of Treated Teachers Matched") +
  theme(legend.position = "bottom")


ggsave("results/graphs/prelim_graph_prop_t_m.png", device = "png", width = 12, height = 8)



# proportion  ---------------------------------------------------------------

ggplot(results_clean, aes(x = method, 
                          y = prop_t_stud_m, 
                          shape = matching_level, 
                          color = matching_priority)) + 
  geom_point() +
  annotate("rect",
           xmin = 0, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           fill = "blue", alpha = .2) +
  annotate("rect",
           xmin = 6.5, xmax = 9.5,
           ymin = -Inf, ymax = Inf,
           fill = "blue", alpha = .2) +
  facet_grid(k_j ~ icc) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() +
  ggtitle("Proportion Treatment Teachers Matched - Student Level") +
  labs(shape = "", color = "", x = "Method", y = "Proportion of Treated Students Matched") +
  theme(legend.position = "bottom")

ggsave("results/graphs/prelim_graph_prop_t_stud_m.png", device = "png", width = 12, height = 8)


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

