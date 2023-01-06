G_k <- rnorm(k)
H_k <- rnorm(k)

L_jk <- rep(rnorm(j, -0.5 * Z_k[k.], 1), i)

A_ijk <- rnorm(i * j, -0.25 * Z_k[k.] + -0.40 * r_k[k.], 1) # observed unit-level covariate
B_ijk <- rnorm(i * j, -0.35 * Z_k[k.] + -0.50 * r_k[k.], 1) # observed unit-level covariate

rep(G_k[k.], i * j)
rep(H_k[k.], i * j)


A_jk_2 = mean(A_ijk^2)
X_A = mean(X_ijk * A_ijk)

#new l2 and w_jk - interaction