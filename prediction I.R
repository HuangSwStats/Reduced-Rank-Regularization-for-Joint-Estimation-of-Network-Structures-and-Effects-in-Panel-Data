P <- 1  
T <- 300
N <- 10
r <- 3
d <- 1
window_size <- 20 
pred_step <- 1    


y <- as.matrix(Trade[, 2:(T + 1)])
X <- as.matrix(REER1[, 2:(T + 1)])


estimate_params <- function(train_start, train_end) {
  H_train <- list()
  for (t in train_start:train_end) {
    H_train[[t]] <- cbind(y[, (t-1):(t-P)], X[, t])
  }
  
  e <- 1e-5
  max_iter <- 1000
  alpha_0 <- matrix(c(0.5, 0.5), (P + d), 1)
  
  Sxx_0 <- matrix(0, N, N)
  Sxy_0 <- matrix(0, N, N)
  for (t in train_start:train_end) {
    Sxx_0 <- Sxx_0 + (1/(train_end - train_start + 1)) * H_train[[t]] %*% alpha_0 %*% t(alpha_0) %*% t(H_train[[t]])
    Sxy_0 <- Sxy_0 + (1/(train_end - train_start + 1)) * y[, t] %*% t(alpha_0) %*% t(H_train[[t]])
  }
  

  eigen_decomp_0 <- eigen(Sxy_0 %*% solve(Sxx_0) %*% t(Sxy_0))
  V_0 <- eigen_decomp_0$vectors[, order(eigen_decomp_0$values, decreasing=TRUE)[1:r]]
  
  a_pre <- V_0
  b_pre <- t(t(V_0) %*% Sxy_0 %*% solve(Sxx_0))
  
  for(iter in 1:max_iter) {
    sum1 <- matrix(0, (P + d), (P + d))
    sum2 <- matrix(0, (P + d), 1)
    for(t in train_start:train_end) {
      H_t <- H_train[[t]]
      sum1 <- sum1 + t(H_t) %*% b_pre %*% t(a_pre) %*% a_pre %*% t(b_pre) %*% H_t
      sum2 <- sum2 + t(H_t) %*% b_pre %*% t(a_pre) %*% y[, t]
    }
    
    alpha_next <- solve(sum1) %*% sum2
    if(norm(alpha_next - alpha_0, "F") < e) break
    alpha_0 <- alpha_next
  }
  return(list(a=a_pre, b=b_pre, alpha=alpha_next))
}

total_error <- 0
pred_count <- 0

pred_points <- seq(P + window_size + 1, T, by = pred_step)

for(pred_t in pred_points) {
  train_start <- pred_t - window_size
  train_end <- pred_t - 1
  params <- estimate_params(train_start, train_end)
  H_pred <- cbind(y[, (pred_t-1):(pred_t-P)], X[, pred_t])

  pred <- params$a %*% t(params$b) %*% H_pred %*% params$alpha
  print(pred[1,])

  actual <- y[, pred_t]
  error <- norm(actual - pred, "F")^2
  total_error <- total_error + error
  pred_count <- pred_count + 1
  
  #cat(sprintf("在时间点 %d 预测，MSE: %.4f\n", pred_t, error))
}

if(pred_count > 0) {
  avg_mse <- total_error / pred_count
  cat(sprintf("\n总预测点数: %d\n平均MSE: %.4f\nRMSE: %.4f",
              pred_count, avg_mse, sqrt(avg_mse)))
} else {
  cat("没有可用的预测点")
}  

