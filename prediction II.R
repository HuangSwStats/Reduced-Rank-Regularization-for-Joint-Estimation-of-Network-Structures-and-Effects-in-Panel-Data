
P <- 1 
T <- 300
N <- 10 
r <- 3
d <- 1


y <- as.matrix(Trade[, 2:(T + 1)])
X <- as.matrix(REER1[, 2:(T + 1)])

Y <- list()
for (t in (P+1):T) {
  Y[[t]] <- y[, (t-1):(t-P)]
}

estimate_params <- function(y, X, P, d, r, t_start, t_end) {
  theta_0 <- matrix(c(0.5, 0.5), nrow = P + d, ncol = 1)
  beta_0 <- matrix(theta_0[1:P], nrow = P, ncol = 1)
  rou_0 <- matrix(theta_0[(P+1):(P + d)], nrow = d, ncol = 1)
  
  e <- 1e-5
  max_iter <- 1000
  m <- 0
  
  n_samples <- t_end - t_start + 1
  if(n_samples <= 0) return(NULL)
  
  SXX_0 <- matrix(0, N, N)
  SYX_0 <- matrix(0, N, N)
  SZZ_0 <- matrix(0, N, N)
  SYZ_0 <- matrix(0, N, N)
  
  for (t in t_start:t_end) {
    Y_t <- Y[[t]]
    X_t <- X[, t]
    SXX_0 <- SXX_0 + (Y_t %*% beta_0 %*% t(beta_0) %*% t(Y_t))/n_samples
    SYX_0 <- SYX_0 + (y[, t] %*% t(beta_0) %*% t(Y_t))/n_samples
    SZZ_0 <- SZZ_0 + (X_t %*% rou_0 %*% t(rou_0) %*% t(X_t))/n_samples
    SYZ_0 <- SYZ_0 + (y[, t] %*% t(rou_0) %*% t(X_t))/n_samples
  }
  
  eigenV <- eigen(SYX_0 %*% (solve(SXX_0) + solve(SZZ_0)) %*% t(SYX_0))
  V <- eigenV$vectors[, 1:r]
  a_pre <- V
  b_pre <- t(t(V) %*% SYX_0 %*% solve(SXX_0))
  
  eigenU <- eigen(SYZ_0 %*% (solve(SXX_0) + solve(SZZ_0)) %*% t(SYZ_0))
  U <- eigenU$vectors[, 1:r]
  a1_pre <- U
  b1_pre <- t(t(U) %*% SYZ_0 %*% solve(SZZ_0))
  
  for(iter in 1:max_iter){
  
    Mt <- list()
    for(t in t_start:t_end){
      Mt[[t]] <- cbind(a_pre %*% t(b_pre) %*% Y[[t]], 
                       a1_pre %*% t(b1_pre) %*% X[, t])
    }
    
    sum1 <- matrix(0, P+d, P+d)
    sum2 <- matrix(0, P+d, 1)
    for(t in t_start:t_end){
      sum1 <- sum1 + t(Mt[[t]]) %*% Mt[[t]]
      sum2 <- sum2 + t(Mt[[t]]) %*% y[, t]
    }
    
    theta_next <- solve(sum1) %*% sum2
    beta_next <- matrix(theta_next[1:P], P, 1)
    rou_next <- matrix(theta_next[(P+1):(P+d)], d, 1)
    
    # 收敛判断
    if(norm(theta_next - theta_0, "F") < e) break
    theta_0 <- theta_next
  }
  
  return(list(
    beta = beta_next,
    rou = rou_next,
    a = a_pre,
    b = b_pre,
    a1 = a1_pre,
    b1 = b1_pre
  ))
}

window_size <- 20  # 训练窗口大小
pred_steps <- seq(P+1+window_size, T, by = 1)  # 预测点
errors <- numeric()

for(t_pred in pred_steps){
  t_start <- t_pred - window_size
  t_end <- t_pred - 1
  params <- estimate_params(y, X, P, d, r, t_start, t_end)
  if(is.null(params)) next
  Y_input <- Y[[t_pred]]
  X_input <- X[, t_pred]
  

  pred <- params$a %*% t(params$b) %*% Y_input %*% params$beta + 
    params$a1 %*% t(params$b1) %*% X_input %*% params$rou
  print(pred[1,])
  
  error <- norm(y[, t_pred] - pred, "F")^2
  errors <- c(errors, error)
}

# 输出结果
if(length(errors) > 0){
  cat("总预测点数:", length(errors), "\n")
  cat("MSE:", mean(errors), "\n")
  cat("RMSE:", sqrt(mean(errors)), "\n")
} else {
  cat("没有有效的预测点")
}

