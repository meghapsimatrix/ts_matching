multi_match(dat = example_dat,
            trt = "D",
            l1_cov = c("X_ijk"),
            l2_cov = c("W_q5", "Z_q5"),
            l2_id = "teacher_id",
            caliper = .25)
