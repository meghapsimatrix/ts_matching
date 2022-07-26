
set.seed(20220615)

k <- 50 #schools
j <- 8  # teachers 
i <- 20 # students

icc3 <- 0.05  # icc for school 
icc2 <- 0.20 # icc for teachers
R2 <- .40


ps_coef <- matrix(c(.07, 0.8, -0.25, 0.6, 0.001, 1, 1))  #what values for pi 6 and 7 in the notes?
pr_star <- .5

outcome_coef <- matrix(c(1, 0.3, .5, .4, -0.2, 1, 1))

delta <- -0.4

