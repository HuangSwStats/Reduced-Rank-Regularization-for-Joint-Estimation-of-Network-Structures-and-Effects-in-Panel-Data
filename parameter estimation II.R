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
  max_iter <-1000#
  theta_0 <- matrix(c(0.5,0.5), (P+d), 1)
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
  

  eigen_decomp_U_0 <- eigen(SYZ_0%*%(solve(SXX_0)+solve(SZZ_0))%*%t(SYZ_0))
  eigen_values_U_0 <- eigen_decomp_U_0$values
  eigen_vectors_U_0 <- eigen_decomp_U_0$vectors
  sorted_indices_U_0 <- order(eigen_values_U_0, decreasing = TRUE)[1:r]
  U_0 <- eigen_vectors_U_0[, sorted_indices_U_0]
  a1_pre <- U_0
  b1_pre <- t(t(U_0)%*%SYZ_0%*%solve(SZZ_0))
  
  Mt_pre <- list()
  for (t in (P+1):T) {
    Mt_pre[[t]] <- matrix(0,N,(P+d))
    Mt_pre[[t]] <- cbind(a_pre%*%t(b_pre)%*%Y[[t]],a1_pre%*%t(b1_pre)%*%X[,t])
  }

  y_0 <- matrix(0, N, T)
  for (t in (P+1):T) {
    y_0[,t] <- y[,t]
  }
  F_norm <- function(matrix) {
    return(sqrt(sum(matrix^2)))
  }
  
  #更新theta(更新beta rou)
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
    
    #V_pre a_next b_next
    eigen_decomp_V_pre <- eigen(SYX_pre%*%(solve(SXX_pre)+solve(SZZ_pre))%*%t(SYX_pre))
    eigen_values_V_pre <- eigen_decomp_V_pre$values
    eigen_vectors_V_pre<- eigen_decomp_V_pre$vectors
    sorted_indices_V_pre <- order(eigen_values_V_pre, decreasing = TRUE)[1:r]
    V_pre <- eigen_vectors_V_pre[, sorted_indices_V_pre] 
    a_next <- V_pre
    b_next <- t(t(V_pre)%*%SYX_pre%*%solve(SXX_pre))
    
    #U_pre a1_next b1_next
    eigen_decomp_U_pre <- eigen(SYZ_pre%*%(solve(SXX_pre)+solve(SZZ_pre))%*%t(SYZ_pre))
    eigen_values_U_pre <- eigen_decomp_U_pre$values
    eigen_vectors_U_pre <- eigen_decomp_U_pre$vectors
    sorted_indices_U_pre <- order(eigen_values_U_pre, decreasing = TRUE)[1:r]
    U_pre <- eigen_vectors_U_pre[, sorted_indices_U_pre]
    a1_next <- U_pre
    b1_next <- t(t(U_pre)%*%SYZ_pre%*%solve(SZZ_pre))
    #Mt_next
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
  A_hat <- a_hat %*% t(b_hat) #估计矩阵A
  A1_hat <- a1_hat %*% t(b1_hat) #估计矩阵A1
  print(theta_hat)

  
  B <- 500  
  theta_boot <- matrix(0, nrow = P + d, ncol = B)
  
  for (b in 1:B) {
    sampled_t <- sample((P+1):T, T-P, replace = TRUE)
    
    Yb <- Y[sampled_t]
    Xb <- X[, sampled_t]
    yb <- y[, sampled_t]
    
  
    SXX_b <- matrix(0, N, N)
    SYX_b <- matrix(0, N, N)
    SZZ_b <- matrix(0, N, N)
    SYZ_b <- matrix(0, N, N)
    
    for (i in 1:length(sampled_t)) {
      SXX_b <- SXX_b + (1 / (T - P)) * Yb[[i]] %*% beta_hat %*% t(beta_hat) %*% t(Yb[[i]])
      SYX_b <- SYX_b + (1 / (T - P)) * yb[, i] %*% t(beta_hat) %*% t(Yb[[i]])
      SZZ_b <- SZZ_b + (1 / (T - P)) * Xb[, i] %*% rou_hat %*% t(rou_hat) %*% t(Xb[, i])
      SYZ_b <- SYZ_b + (1 / (T - P)) * yb[, i] %*% t(rou_hat) %*% t(Xb[, i])
    }
   
    Vb <- eigen(SYX_b %*% (solve(SXX_b) + solve(SZZ_b)) %*% t(SYX_b))$vectors[, 1:r]
    Ub <- eigen(SYZ_b %*% (solve(SXX_b) + solve(SZZ_b)) %*% t(SYZ_b))$vectors[, 1:r]
    ab <- Vb; ab1 <- Ub
    bb <- t(t(Vb) %*% SYX_b %*% solve(SXX_b))
    bb1 <- t(t(Ub) %*% SYZ_b %*% solve(SZZ_b))
    

    Mb <- list()
    for (i in 1:length(sampled_t)) {
      Mb[[i]] <- cbind(ab %*% t(bb) %*% Yb[[i]], ab1 %*% t(bb1) %*% Xb[, i])
    }
    
    
    sum1b <- matrix(0, P + d, P + d)
    sum2b <- matrix(0, P + d, 1)
    for (i in 1:length(sampled_t)) {
      sum1b <- sum1b + t(Mb[[i]]) %*% Mb[[i]]
      sum2b <- sum2b + t(Mb[[i]]) %*% yb[, i]
    }
    thetab <- solve(sum1b) %*% sum2b
    theta_boot[, b] <- thetab[, 1]
  }
  

  theta_se <- apply(theta_boot, 1, sd)
  theta_tstat <- theta_hat[, 1] / theta_se
  theta_pval <- 2 * (1 - pnorm(abs(theta_tstat)))
  
 
  theta_result <- data.frame(
    Estimate = theta_hat[, 1],
    Std.Error = theta_se,
    t.value = theta_tstat,
    p.value = theta_pval
  )
  
  print(theta_result)
  
  
 