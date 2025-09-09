P <- 1  
T <- 300
N <- 10  
r <-3
d <- 1  

y <- as.matrix(Trade[, 2:(T + 1)])
X <- as.matrix(REER1[, 2:(T + 1)])  
H <- list()
for (t in (P+1):T) {
  H[[t]] <- cbind(y[, (t-1):(t-P)],X[, t])
}

e <- 1e-5
m <- 0
max_iter <- 1000  
alpha_0 <- matrix(c(0.5,0.5), (P + d), 1)
alpha_pre <- matrix(0, (P + d), 1)
alpha_next <- matrix(0, (P + d), 1)

Sxx_0 <- matrix(0, N, N)
Sxy_0 <- matrix(0, N, N)
for (t in (P+1):T) {
  Sxx_0 <- Sxx_0 + (1 / (T -P)) * H[[t]] %*% alpha_0 %*% t(alpha_0) %*% t(H[[t]])
  Sxy_0 <- Sxy_0 + (1 / (T -P)) * y[, t] %*% t(alpha_0) %*% t(H[[t]])
}

eigen_decomp_0 <- eigen(Sxy_0 %*% solve(Sxx_0) %*% t(Sxy_0))
eigen_values_0 <- eigen_decomp_0$values 
eigen_vectors_0 <- eigen_decomp_0$vectors
sorted_indices_0 <- order(eigen_values_0, decreasing = TRUE)[1:r]  
V_0 <- eigen_vectors_0[, sorted_indices_0]  

a_pre <- V_0
b_pre <- t(t(V_0) %*% Sxy_0 %*% solve(Sxx_0))

y_0 <- matrix(0, N, T)
F_norm <- function(matrix) {
  return(sqrt(sum(matrix^2)))
}

for (t in (P+1):T) {
  y_0[, t] <- y[,t]
}

sum1_0 <- matrix(0, (P + d), (P + d))
sum2_0 <- matrix(0, (P + d), 1)
for (t in (P+1):T) {
  sum1_0 <- sum1_0 + t(H[[t]]) %*% b_pre %*% t(a_pre) %*% a_pre %*% t(b_pre) %*% H[[t]]
  sum2_0 <- sum2_0 + t(H[[t]]) %*% b_pre %*% t(a_pre) %*% y_0[, t]
}

alpha_pre <- solve(sum1_0) %*% sum2_0

K1 <- kronecker(t(alpha_pre), a_pre %*% t(b_pre))
repeat {
  Sxx_pre <- matrix(0, N, N)
  Sxy_pre <- matrix(0, N, N)
  if (alpha_pre[1] > 1 || alpha_pre[P + 1] > 1) {
    alpha_pre <- alpha_0
  }
  for (t in (P+1):T) {
    Sxx_pre <- Sxx_pre + (1 / (T - P)) * H[[t]] %*% alpha_pre %*% t(alpha_pre) %*% t(H[[t]])
    Sxy_pre <- Sxy_pre + (1 / (T - P)) * y[, t] %*% t(alpha_pre) %*% t(H[[t]])
  }

  eigen_decomp_pre <- eigen(Sxy_pre %*% solve(Sxx_pre) %*% t(Sxy_pre))
  eigen_values_pre <- eigen_decomp_pre$values 
  eigen_vectors_pre <- eigen_decomp_pre$vectors
  sorted_indices_pre <- order(eigen_values_pre, decreasing = TRUE)[1:r] 
  V_pre <- eigen_vectors_pre[, sorted_indices_pre]  
  a_next <- V_pre
  b_next <- t(t(V_pre) %*% Sxy_pre %*% solve(Sxx_pre))
  
  y_pre <- matrix(0, N, T)
  for (t in (P+1):T) {
    y_pre[, t] <- a_pre%*% t(b_pre) %*% H[[t]] %*% alpha_pre}
  
  sum1_pre <- matrix(0, (P + d), (P + d))
  sum2_pre <- matrix(0, (P + d), 1)
  for (t in (P+1):T) {
    sum1_pre <- sum1_pre + t(H[[t]]) %*% b_next %*% t(a_next) %*% a_next %*% t(b_next) %*% H[[t]]
    sum2_pre <- sum2_pre + t(H[[t]]) %*% b_next %*% t(a_next) %*% y_pre[, t]
  }
  
  alpha_next <- solve(sum1_pre) %*% sum2_pre
  K2 <- kronecker(t(alpha_next), a_next %*% t(b_next))
  tolerance <- F_norm(K2 - K1)
 
  if (abs(tolerance) <= e ) {
    break
  }
  
  K1 <- K2
  alpha_pre <- alpha_next
  a_pre <- a_next
  b_pre <- b_next
  m <- m + 1
}
a_hat <- a_next
b_hat <- b_next
alpha_hat <- alpha_next
A_hat <- a_hat %*% t(b_hat) 
print(alpha_hat)

B <- 500  
alpha_boot <- matrix(0, nrow = (P + d), ncol = B) 

for (b in 1:B) {
  sampled_t <- sample((P + 1):T, T - P, replace = TRUE)

  y_boot <- y[, sampled_t]
  X_boot <- X[, sampled_t]
  
  H_boot <- list()
  for (i in 1:length(sampled_t)) {
    t_real <- sampled_t[i]
    H_boot[[i]] <- cbind(y[, t_real - 1], X[, t_real])
  }
  
  Sxx_boot <- matrix(0, N, N)
  Sxy_boot <- matrix(0, N, N)
  for (i in 1:length(sampled_t)) {
    Sxx_boot <- Sxx_boot + (1 / (T - P)) * H_boot[[i]] %*% alpha_hat %*% t(alpha_hat) %*% t(H_boot[[i]])
    Sxy_boot <- Sxy_boot + (1 / (T - P)) * y[, sampled_t[i]] %*% t(alpha_hat) %*% t(H_boot[[i]])
  }
  
  eigen_decomp_boot <- eigen(Sxy_boot %*% solve(Sxx_boot) %*% t(Sxy_boot))
  eigen_values_boot <- eigen_decomp_boot$values
  eigen_vectors_boot <- eigen_decomp_boot$vectors
  sorted_indices_boot <- order(eigen_values_boot, decreasing = TRUE)[1:r]
  V_boot <- eigen_vectors_boot[, sorted_indices_boot]
  
  a_boot <- V_boot
  b_boot <- t(t(V_boot) %*% Sxy_boot %*% solve(Sxx_boot))

  sum1_boot <- matrix(0, (P + d), (P + d))
  sum2_boot <- matrix(0, (P + d), 1)
  for (i in 1:length(sampled_t)) {
    sum1_boot <- sum1_boot + t(H_boot[[i]]) %*% b_boot %*% t(a_boot) %*% a_boot %*% t(b_boot) %*% H_boot[[i]]
    sum2_boot <- sum2_boot + t(H_boot[[i]]) %*% b_boot %*% t(a_boot) %*% y[, sampled_t[i]]
  }
  
  alpha_b <- solve(sum1_boot) %*% sum2_boot
  alpha_boot[, b] <- alpha_b[, 1]
}

alpha_se <- apply(alpha_boot, 1, sd)
alpha_tstat <- alpha_hat[,1] / alpha_se
alpha_pval <- 2 * (1 - pnorm(abs(alpha_tstat)))

alpha_result <- data.frame(
  Estimate = alpha_hat[,1],
  Std.Error = alpha_se,
  t.value = alpha_tstat,
  p.value = alpha_pval
)

print(alpha_result)

