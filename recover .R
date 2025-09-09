library(readxl)
library(Rsolnp)
library(mvtnorm)


P <- 1
T <- 300
N <- 10 
r <-3
d <- 1  


y <- as.matrix(Trade[, 2:(T + 1)])
X <- as.matrix(REER1[, 2:(T + 1)])  
Y <- list()
for (t in (P+1):T) {
  Y[[t]] <- y[, (t-1):(t-P)]
}



e <- 1e-5
m <- 0
max_iter <-1000
theta_0 <- matrix(c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5), (P+d), 1)
beta_0 <- matrix(theta_0[1:P],P,1)
rou_0 <- matrix(theta_0[(P+1):(P+d)],d,1)
theta_pre <- matrix(0,(P+d),1)
beta_pre <- matrix(0,P,1)
rou_pre <- matrix(0,d,1)
theta_next <- matrix(0,(P+d),1)
beta_next <- matrix(0,P,1)
rou_next <- matrix(0,d,1)


SXX_0 <- matrix(0,N,N)
SYX_0 <- matrix(0,N,N)
SZZ_0 <- matrix(0,N,N)
SYZ_0 <- matrix(0,N,N)
for (t in (P+1):T) {
  SXX_0 <- SXX_0 + (1/(T-P))*Y[[t]]%*%beta_0%*%t(beta_0)%*%t(Y[[t]])
  SYX_0 <- SYX_0 + (1/(T-P))*y[,t]%*%t(beta_0)%*%t(Y[[t]])
  SZZ_0 <- SZZ_0 + (1/(T-P))*X[,t]%*%rou_0%*%t(rou_0)%*%t(X[,t])
  SYZ_0 <- SYZ_0 + (1/(T-P))*y[,t]%*%t(rou_0)%*%t(X[,t])
}


eigen_decomp_V_0 <- eigen(SYX_0%*%(solve(SXX_0)+solve(SZZ_0))%*%t(SYX_0))
eigen_values_V_0 <- eigen_decomp_V_0$values
eigen_vectors_V_0 <- eigen_decomp_V_0$vectors
sorted_indices_V_0 <- order(eigen_values_V_0, decreasing = TRUE)[1:r]
V_0 <- eigen_vectors_V_0[, sorted_indices_V_0] 

a_pre <- V_0
b_pre <- t(t(V_0)%*%SYX_0%*%solve(SXX_0))

#U_0 a1_pre b1_pre
eigen_decomp_U_0 <- eigen(SYZ_0%*%(solve(SXX_0)+solve(SZZ_0))%*%t(SYZ_0))
eigen_values_U_0 <- eigen_decomp_U_0$values
eigen_vectors_U_0 <- eigen_decomp_U_0$vectors
sorted_indices_U_0 <- order(eigen_values_U_0, decreasing = TRUE)[1:r]
U_0 <- eigen_vectors_U_0[, sorted_indices_U_0]
a1_pre <- U_0
b1_pre <- t(t(U_0)%*%SYZ_0%*%solve(SZZ_0))

#Mt_pre
Mt_pre <- list()
for (t in (P+1):T) {
  Mt_pre[[t]] <- matrix(0,N,(P+d))
  Mt_pre[[t]] <- cbind(a_pre%*%t(b_pre)%*%Y[[t]],a1_pre%*%t(b1_pre)%*%X[,t])
}
#y_t^(m) y_0
y_0 <- matrix(0, N, T)
for (t in (P+1):T) {
  y_0[,t] <- y[,t]
}
F_norm <- function(matrix) {
  return(sqrt(sum(matrix^2)))
}


sum1_0 <- matrix(0,(P+d),(P+d))
sum2_0 <- matrix(0,(P+d),1)
for (t in(P+1):T) {
  sum1_0 <- sum1_0 + t(Mt_pre[[t]])%*%Mt_pre[[t]]
  sum2_0 <- sum2_0 + t(Mt_pre[[t]])%*%y_0[,t]
}
theta_pre <- solve(sum1_0) %*% sum2_0
beta_pre <- matrix(theta_pre[1:P],P,1)
rou_pre <- matrix(theta_pre[(P+1):(P+d)],d,1)

K1_beta <- kronecker(t(beta_pre), a_pre %*% t(b_pre))
K1_rou <- kronecker(t(rou_pre), a1_pre %*% t(b1_pre))

#迭代算法
repeat {

  SXX_pre <- matrix(0,N,N)
  SYX_pre <- matrix(0,N,N)
  SZZ_pre <- matrix(0,N,N)
  SYZ_pre <- matrix(0,N,N)
  if (theta_pre[1] > 1 || theta_pre[P+1] >1) {
    theta_pre <- theta_0
  }
  for (t in (P+1):T) {
    SXX_pre <- SXX_pre + (1/(T-P))*Y[[t]]%*%beta_pre%*%t(beta_pre)%*%t(Y[[t]])
    SYX_pre <- SYX_pre + (1/(T-P))*y[,t]%*%t(beta_pre)%*%t(Y[[t]])
    SZZ_pre <- SZZ_pre + (1/(T-P))*X[,t]%*%rou_pre%*%t(rou_pre)%*%t(X[,t])
    SYZ_pre <- SYZ_pre + (1/(T-P))*y[,t]%*%t(rou_pre)%*%t(X[,t])
  }
  

  eigen_decomp_V_pre <- eigen(SYX_pre%*%(solve(SXX_pre)+solve(SZZ_pre))%*%t(SYX_pre))
  eigen_values_V_pre <- eigen_decomp_V_pre$values
  eigen_vectors_V_pre<- eigen_decomp_V_pre$vectors
  sorted_indices_V_pre <- order(eigen_values_V_pre, decreasing = TRUE)[1:r]
  V_pre <- eigen_vectors_V_pre[, sorted_indices_V_pre] 
  a_next <- V_pre
  b_next <- t(t(V_pre)%*%SYX_pre%*%solve(SXX_pre))
  

  eigen_decomp_U_pre <- eigen(SYZ_pre%*%(solve(SXX_pre)+solve(SZZ_pre))%*%t(SYZ_pre))
  eigen_values_U_pre <- eigen_decomp_U_pre$values
  eigen_vectors_U_pre <- eigen_decomp_U_pre$vectors
  sorted_indices_U_pre <- order(eigen_values_U_pre, decreasing = TRUE)[1:r]
  U_pre <- eigen_vectors_U_pre[, sorted_indices_U_pre]
  a1_next <- U_pre
  b1_next <- t(t(U_pre)%*%SYZ_pre%*%solve(SZZ_pre))
  
  Mt_next <- list()
  for (t in (P+1):T) {
    Mt_next[[t]] <- matrix(0,N,(P+d))
    Mt_next[[t]] <- cbind(a_next%*%t(b_next)%*%Y[[t]],a1_next%*%t(b1_next)%*%X[,t])
  }
 
  y_pre <- matrix(0, N, T)
  for (t in (P+1):T) {
    y_pre[,t] <- a_pre%*%t(b_pre)%*%Y[[t]]%*%beta_pre+a1_pre%*%t(b1_pre)%*%X[,t]%*%rou_pre
  }
  
  sum1_pre <- matrix(0,(P+d),(P+d))
  sum2_pre <- matrix(0,(P+d),1)
  for (t in (P+1):T) {
    sum1_pre <- sum1_pre + t(Mt_next[[t]])%*%Mt_next[[t]]
    sum2_pre <- sum2_pre + t(Mt_next[[t]])%*%y_pre[,t]
  }
  
  theta_next <- solve(sum1_pre) %*% sum2_pre
  beta_next <- matrix(theta_next[1:P],P,1)
  rou_next <- matrix(theta_next[(P+1):(P+d)],d,1)
  
  K2_beta <- kronecker(t(beta_next), a_next %*% t(b_next))
  K2_rou <- kronecker(t(rou_next), a1_next %*% t(b1_next))
  
  tolerance <- F_norm(K2_beta - K1_beta)+F_norm(K2_rou - K1_rou)
  
  
  K1_beta <- K2_beta
  K1_rou <- K2_rou
  theta_pre <- theta_next
  beta_pre <- beta_next
  rou_pre <- rou_next
  
  
  a_pre <- a_next
  b_pre <- b_next
  a1_pre <- a1_next
  b1_pre <- b1_next
  
  m <- m+1
  if (abs(tolerance) <= e ){
    break
  }
}

a_hat <- a_next
b_hat <- b_next
a1_hat <- a1_next
b1_hat <- b1_next
theta_hat <- theta_next
beta_hat <- matrix(theta_hat[1:P],P,1)
rou_hat <- matrix(theta_hat[(P+1):(P+d)],d,1)
A_hat <- a_hat %*% t(b_hat) 
A1_hat <- a1_hat %*% t(b1_hat) 


I_N <- diag(N)
ut <- matrix(0, N, T)
for (t in (P+1):T) {
  ut[, t] <- y[, t] - A_hat %*% Y[[t]] %*% beta_hat - A1_hat %*% X[, t] %*% rou_hat
}

Sigma_hat <- cov(t(ut[, (P+1):T]))

beta0_init <- 0.1
a_tilde_init <- a_hat
b_tilde_init <- b_hat
W_init <- a_tilde_init %*% t(b_tilde_init)
W0_init <- W_init
diag(W0_init) <- 0

objective <- function(par) {
  beta0 <- par[1]
  a_tilde <- matrix(par[2:(1+N*r)], N, r)
  b_tilde <- matrix(par[(2+N*r):length(par)], N, r)
  
  W <- a_tilde %*% t(b_tilde)
  W0 <- W
  diag(W0) <- 0
  A_model <- solve(I_N - beta0 * W0) %*% W
  Sigma_model <- solve(I_N - beta0 * W0) %*% t(solve(I_N - beta0 * W0))
  loss <- sum((Sigma_hat - Sigma_model)^2) + sum((A_hat - A_model)^2)
  return(loss)
}

par_init <- c(beta0_init, as.vector(a_tilde_init), as.vector(b_tilde_init))

opt_result <- solnp(
  pars = par_init,
  fun = objective,
  control = list(trace = 0)
)


beta0_hat <- opt_result$pars[1]

a_tilde_hat <- matrix(opt_result$pars[2:(1+N*r)], N, r)
b_tilde_hat <- matrix(opt_result$pars[(2+N*r):length(opt_result$pars)], N, r)
W_hat <- a_tilde_hat %*% t(b_tilde_hat)
W_hat <- W_hat / norm(W_hat, type = "F")  
W1_hat <- (I_N-beta0_hat*W_hat)%*%A1_hat 
W1_hat <- W1_hat/  norm(W1_hat, type = "F")  

n_time <- T - P              
k_par  <- length(opt_result$pars)  

log_liks <- numeric(n_time)
for (t in (P+1):T) {
  log_liks[t - P] <- dmvnorm(ut[, t], mean = rep(0, N), sigma = Sigma_hat, log = TRUE)
}
logLik_total <- sum(log_liks)

AIC_value <- 2 * k_par - 2 * logLik_total
BIC_value <- log(n_time) * k_par - 2 * logLik_total

#输出结果
cat("------- 模型信息准则 -------\n")
cat("样本大小 n =", n_time, "\n")
cat("参数个数 k =", k_par, "\n")
cat("对数似然  =", round(logLik_total, 4), "\n")
cat("AIC       =", round(AIC_value, 4), "\n")
cat("BIC       =", round(BIC_value, 4), "\n")

B <- 500  
beta0_boot <- numeric(B)


for (b in 1:B) {
  t_idx <- sample((P+1):T, replace = TRUE)
  
  y_boot <- y
  X_boot <- X
  Y_boot <- list()
  for (i in seq_along(t_idx)) {
    t_star <- t_idx[i]
    y_boot[, P + i] <- y[, t_star]
    X_boot[, P + i] <- X[, t_star]
    Y_boot[[P + i]] <- Y[[t_star]]
  }
  
  ut_boot <- matrix(0, N, T)
  for (t in (P+1):T) {
    ut_boot[, t] <- y_boot[, t] - A_hat %*% Y_boot[[t]] %*% beta_hat - A1_hat %*% X_boot[, t] %*% rou_hat
  }

  Sigma_hat_boot <- cov(t(ut_boot[, (P+1):T]))
  objective_boot <- function(par) {
    beta0 <- par[1]
    a_tilde <- matrix(par[2:(1+N*r)], N, r)
    b_tilde <- matrix(par[(2+N*r):length(par)], N, r)
    
    W <- a_tilde %*% t(b_tilde)
    diag(W) <- 0
    W0 <- W
    
    A_model <- solve(I_N - beta0 * W0) %*% W
    Sigma_model <- solve(I_N - beta0 * W0) %*% t(solve(I_N - beta0 * W0))
    
    loss <- sum((Sigma_hat_boot - Sigma_model)^2) + sum((A_hat - A_model)^2)
    return(loss)
  }
  
  par_init_boot <- c(beta0_hat, as.vector(a_tilde_hat), as.vector(b_tilde_hat))
  
  opt_result_boot <- try(solnp(pars = par_init_boot, fun = objective_boot, control = list(trace = 0)), silent = TRUE)
  
  if (inherits(opt_result_boot, "try-error")) {
    beta0_boot[b] <- NA
  } else {
    beta0_boot[b] <- opt_result_boot$pars[1]
  }
}
beta0_boot <- na.omit(beta0_boot)

# 计算标准误与 p 值
se_beta0 <- sd(beta0_boot)
z_beta0 <- beta0_hat / se_beta0
p_beta0 <- 2 * (1 - pnorm(abs(z_beta0)))

cat("标准误 (SE):", round(se_beta0, 5), "\n")
cat("p 值:", round(p_beta0, 4), "\n")

