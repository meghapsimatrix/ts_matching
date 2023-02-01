method <- c(.1,.2,1,2,3,4,5,6,7,8,9,10,11,12)

# Bias
for(i in method) print(c(i,mean(results$bias[which(results$method==i)])))

# Proportion in treatment group (teacher)
for(i in method) print(c(i,mean(results$prop_t_m[which(results$method==i)])))

# Proportion in treatment group (student)
for(i in method) print(c(i,mean(results$prop_t_stud_m[which(results$method==i)])))

# X_ijk
for(i in method) print(c(i,mean(results$X_ijk[which(results$method==i)])))

# U_ijk
for(i in method) print(c(i,mean(results$U_ijk[which(results$method==i)])))

# X_jk
for(i in method) print(c(i,mean(results$X_jk[which(results$method==i)])))

# W_jk
for(i in method) print(c(i,mean(results$W_jk[which(results$method==i)])))

# V_jk
for(i in method) print(c(i,mean(results$V_jk[which(results$method==i)])))

# Z_k
for(i in method) print(c(i,mean(results$Z_k[which(results$method==i)])))

