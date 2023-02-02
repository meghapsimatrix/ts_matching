library(tidyverse)

# Load results ------------------------------------------------------------

files <- list.files("results/take_4/results", full.names = TRUE)


load_res <- function(file) {

  load(file)
  results

}


results <- map_dfr(files, load_res)



# clean results  ----------------------------------------------------------

K <- 200 + 400 + 300

# check if seeds are different
results %>%
  group_by(seed) %>%
  count() %>%
  filter(n !=14)

results_clean <- 
  results %>%
  group_by(method, k, j, icc3, icc2) %>%
  summarize_at(vars(bias:prop_t_stud_m), mean) %>%
  ungroup() %>%
  mutate(K = K)


results_clean %>%
  group_by(method) %>%
  summarize(bias = mean(bias))

results_clean <- 
  results_clean %>%
  mutate(method = as.character(method)) %>%
  filter(!(method %in% c("0.1", "0.2"))) %>%
  mutate(method = factor(method, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))) %>%
  mutate(k_j = paste0("# Schools = ", k, ", ", "# Teachers = ", j),
         icc = paste0("icc3 = ", icc3, ", ", "icc2  = ", icc2)) %>%
  mutate(matching_priority = fct_collapse(method,
                                          "No Distinction" = c("1", "2", "3"),
                                          "Within-site Matching" = c("4", "5", "6"),
                                          "Within-group Matching" = c("7", "8", "9"),
                                          "Within-site then -group Matching" = c("10", "11", "12")),
         matching_level = fct_collapse(method,
                                       "Units only" = c("1", "4", "7", "10"),
                                       "Clusters only" = c("2", "5", "8", "11"),
                                       "Units & Clusters" = c("3", "6", "9", "12"))) 


save(results_clean, file = "results/results/results_clean.RData")


zoom_results <- results_clean %>% 
  filter(k == 20, icc3 == icc2)


# graph function ----------------------------------------------------------

make_plot <- function(dat, 
                      title,
                      y_label, 
                      yint = 0,
                      type_line = "dashed",
                      dot_size = 2,
                      text = "large"){
  
  p <- ggplot(dat, aes(x = method, 
                  y = outcome, 
                  shape = matching_level, 
                  color = matching_priority)) + 
    geom_point(size = dot_size) +
    annotate("rect",
             xmin = 0, xmax = 3.5,
             ymin = -Inf, ymax = Inf,
             fill = "blue", alpha = .2) +
    annotate("rect",
             xmin = 6.5, xmax = 9.5,
             ymin = -Inf, ymax = Inf,
             fill = "blue", alpha = .2) +
    facet_grid(icc ~ k_j) +
    geom_hline(yintercept = yint, linetype = type_line) +
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    ggtitle(title) +
    labs(shape = "", color = "", x = "Method", y = y_label) +
    theme(legend.position = "bottom")
  
  if(text == "large"){
    
    p <- p  +
      theme(strip.text.x = element_text(size = 12),
                   strip.text.y = element_text(size = 12),
                   axis.title.x = element_text(size = 12),
                   axis.title.y = element_text(size = 12),
                   legend.text = element_text(size = 10),
                   plot.title = element_text(size = 14, face = "bold"))
  }
  
  return(p)
}





# preliminary graphs ------------------------------------------------------

# bias --------------------------------------------------------------------

make_plot(dat = results_clean %>% mutate(outcome = bias), 
          title = "Bias",
          y_label = "Bias")

ggsave("results/graphs/full_graph_bias.png", device = "png", width = 12, height = 8)

make_plot(dat = zoom_results %>% mutate(outcome = bias), 
          title = "Bias",
          y_label = "Bias",
          dot_size = 3) 

ggsave("results/graphs/zoom_graph_bias.png", device = "png", width = 12, height = 8)


# check_mcse 
# fairly small compared to bias
# do we want to include mcse in graph somehow? Like ci dot and line plot? - but it might now show bc the mcse is quite small
results_clean %>%
  group_by(method) %>%
  summarize(bias_mcse = max(bias_mcse)) %>%
  arrange(desc(bias_mcse))

# rmse --------------------------------------------------------------------


make_plot(dat = results_clean %>% mutate(outcome = rmse), 
          title = "RMSE",
          y_label = "RMSE")

ggsave("results/graphs/full_graph_rmse.png", device = "png", width = 12, height = 8)

make_plot(dat = zoom_results %>% mutate(outcome = rmse), 
          title = "RMSE",
          y_label = "RMSE",
          dot_size = 3) 

ggsave("results/graphs/zoom_graph_rmse.png", device = "png", width = 12, height = 8)


     
# check_mcse 
# fairly small compared to rmse 
results_clean %>%
  group_by(method) %>%
  summarize(rmse_mcse = max(rmse_mcse)) %>%
  arrange(desc(rmse_mcse))

# proportion teachers  ---------------------------------------------------------------

make_plot(dat = results_clean %>% mutate(outcome = prop_t_m), 
          title = "Proportion of Treated Teachers Matched",
          y_label = "Proportion of Treated Teachers Matched",
          yint = 1)

ggsave("results/graphs/full_graph_prop_t_m.png", device = "png", width = 12, height = 8)

make_plot(dat = zoom_results %>% mutate(outcome = prop_t_m), 
          title = "Proportion of Treated Teachers Matched",
          y_label = "Proportion of Treated Teachers Matched",
          yint = 1,
          dot_size = 3)

ggsave("results/graphs/zoom_graph_prop_t_m.png", device = "png", width = 12, height = 8)



# proportion students  ---------------------------------------------------------------

make_plot(dat = results_clean %>% mutate(outcome = prop_t_stud_m), 
          title = "Proportion of Treated Students Matched",
          y_label = "Proportion of Treated Students Matched",
          yint = 1)

ggsave("results/graphs/full_graph_prop_t_stud_m.png", device = "png", width = 12, height = 8)

make_plot(dat = zoom_results %>% mutate(outcome = prop_t_stud_m), 
          title = "Proportion of Treated Students Matched",
          y_label = "Proportion of Treated Students Matched",
          yint = 1,
          dot_size = 3) 

ggsave("results/graphs/zoom_graph_prop_t_stud_m.png", device = "png", width = 12, height = 8)



# balance -----------------------------------------------------------------


make_plot(dat = results_clean %>% mutate(outcome = W_jk), 
          title = "SMD for W_jk",
          y_label = "Standardized Mean Difference",
          type_line = "solid") +
  geom_hline(yintercept = .25, linetype = "dashed")

ggsave("results/graphs/full_graph_smd_W_jk.png", device = "png", width = 12, height = 8)

make_plot(dat = zoom_results %>% mutate(outcome = W_jk), 
          title = "SMD for W_jk",
          y_label = "Standardized Mean Difference",
          type_line = "solid",
          dot_size = 3) +
  geom_hline(yintercept = .25, linetype = "dashed") 

ggsave("results/graphs/zoom_graph_smd_W_jk.png", device = "png", width = 12, height = 8)




make_plot(dat = results_clean %>% mutate(outcome = X_jk), 
          title = "SMD for X_jk",
          y_label = "Standardized Mean Difference",
          type_line = "solid") +
  geom_hline(yintercept = .25, linetype = "dashed")

ggsave("results/graphs/full_graph_smd_X_jk.png", device = "png", width = 12, height = 8)

make_plot(dat = zoom_results %>% mutate(outcome = X_jk), 
          title = "SMD for X_jk",
          y_label = "Standardized Mean Difference",
          type_line = "solid",
          dot_size = 3) +
  geom_hline(yintercept = .25, linetype = "dashed")

ggsave("results/graphs/zoom_graph_smd_X_jk.png", device = "png", width = 12, height = 8)



make_plot(dat = results_clean %>% mutate(outcome = X_ijk), 
          title = "SMD for X_ijk",
          y_label = "Standardized Mean Difference",
          type_line = "solid") +
  geom_hline(yintercept = .25, linetype = "dashed")

ggsave("results/graphs/full_graph_smd_X_ijk.png", device = "png", width = 12, height = 8)

make_plot(dat = zoom_results %>% mutate(outcome = X_ijk), 
          title = "SMD for X_ijk",
          y_label = "Standardized Mean Difference",
          type_line = "solid",
          dot_size = 3) +
  geom_hline(yintercept = .25, linetype = "dashed") 

ggsave("results/graphs/zoom_graph_smd_X_ijk.png", device = "png", width = 12, height = 8)


make_plot(dat = results_clean %>% mutate(outcome = Z_k), 
          title = "SMD for Z_k",
          y_label = "Standardized Mean Difference",
          type_line = "solid") +
  geom_hline(yintercept = .25, linetype = "dashed")

ggsave("results/graphs/full_graph_smd_Z_k.png", device = "png", width = 12, height = 8)

make_plot(dat = zoom_results %>% mutate(outcome = Z_k), 
          title = "SMD for Z_k",
          y_label = "Standardized Mean Difference",
          type_line = "solid",
          dot_size = 3) +
  geom_hline(yintercept = .25, linetype = "dashed") 

ggsave("results/graphs/zoom_graph_smd_Z_k.png", device = "png", width = 12, height = 8)



make_plot(dat = results_clean %>% mutate(outcome = U_ijk), 
          title = "SMD for U_ijk",
          y_label = "Standardized Mean Difference",
          type_line = "solid") +
  geom_hline(yintercept = .25, linetype = "dashed")

ggsave("results/graphs/full_graph_smd_U_ijk.png", device = "png", width = 12, height = 8)

make_plot(dat = zoom_results %>% mutate(outcome =  U_ijk), 
          title = "SMD for  U_ijk",
          y_label = "Standardized Mean Difference",
          type_line = "solid",
          dot_size = 3) +
  geom_hline(yintercept = .25, linetype = "dashed") 

ggsave("results/graphs/zoom_graph_smd_U_ijk.png", device = "png", width = 12, height = 8)



