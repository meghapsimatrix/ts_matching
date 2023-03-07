ggplot(results_clean, aes(x = bias, 
                          y = X_ijk, 
                          shape = matching_level, 
                          color = matching_priority)) + 
  geom_point() +
  facet_grid(icc ~ k_j) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() +
  theme(legend.position = "bottom")


ggsave("results/graphs/bias_smd.png", device = "png", width = 12, height = 9)
