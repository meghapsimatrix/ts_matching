
# estimate ps -------------------------------------------------------------

unit_ps_model <- glmer(D ~ X_ijk + W_jk + Z_k + 
                         (1 | school_id),
                       family = "binomial", 
                       data = dat)

dat$ps_unit <- predict(unit_ps_model, type = "link")
dat$ps_unit_pr <- predict(unit_ps_model, type = "response")


cluster_ps_model <- glmer(D ~ W_jk + Z_k + 
                            (1 | school_id),
                          family = "binomial",
                          data = cluster_level_dat)

cluster_level_dat$ps_cluster <- predict(cluster_ps_model, type = "link")
