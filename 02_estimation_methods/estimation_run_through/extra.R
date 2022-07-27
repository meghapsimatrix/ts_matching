# cobalt way

calc_bal <- function(dat){
  
  bal_stats <-
    bal.tab(D ~ Z_k + r_k + W_jk + u_jk + X_jk + X_ijk + U_ijk,
            data = dat,
            method = "weighting",
            weights = "weights",
            s.d.denom = "treated")
  
  bal_dat <- 
    bal_stats$Balance %>%
    rownames_to_column("var") %>%
    select(var, smd = Diff.Adj)
  
  return(bal_dat)
  
}





# no distinction units only 
# by default nearest neighbor, no replacement

# just with matchit (not multi-level)
m_out_1 <- matchit(D ~ X_ijk + W_jk + Z_k,
                   caliper = .25,
                   data = example_dat)

m_out_1$match.matrix


balance_stats <- bal.tab(m_out_1)
balance_stats$Balance

m_data_1 <- match.data(m_out_1)
table(m_data_1$weights)  # if we get to a point were we have to use weights - lmer won't work

out_mod_1 <- lmer(Y_ijk ~ D  + X_ijk + W_jk + Z_k + 
                    (1 | teacher_id) + (1 | school_id),
                  data = m_data_1)

summary(out_mod_1)



# do I use the ps estimated using glmer to match?


# matchMulti --------------------------------------------------------------

match_simple <- matchMulti(example_dat, 
                           treatment = "D",
                           school.id = "teacher_id",
                           match.students = FALSE,
                           student.vars = "X_ijk",
                           verbose = TRUE)

# two level - level 2 is teachers where the treatment is 

#str(match_simple)

match_data <- as.data.frame(match_simple$matched)
head(match_data)

mod <- lme(Y_ijk ~ D, 
           random = ~ 1 | pair.id/ school_id,
           data = match_data)

summary(mod)
