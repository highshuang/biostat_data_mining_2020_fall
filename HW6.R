##### PHS597 HW6
# assess the conference interval for the model Y=Xb + e; 
# where e is standard normal distribution


library(MASS)
library(metaSEM)
library(matlib)

set.seed(1)

# data 
n <- 100
p <- 4
sigma <- vec2symMat(c(1, 0.1, 0.1, 1, 0.1, 1))
X_np <- cbind(rep(1,n), mvrnorm(n = n, mu = c(0.1, 0.5, -1), Sigma = sigma))

# simulate error terms and calculate y with true b = c(0.1, 1.3, 0.8, -1)
b <- c(0.1, 1.3, 0.8, -1) #c(1, 2.3, -0.8, 1.4)
Y <- X_np %*% b + rnorm(n = n, mean = 0, sd = 1)


# Part I: implementation of Non-Parametric Bootstrap 
Boot_nonpar <- function(B=10000, Y=Y, X=X_np) {
  set.seed(597)
  # get the size of the orignal data 
  N <- nrow(X)
  
  # matrix to store beta estimate 
  beta <- matrix(rep(0, B*p), nrow = p)
 
  for (i in 1:B) {
    # resample a bootstrap sample 
    index <- sample(N, replace = T)
    sample_X <- X[index, ]
    sample_Y <- Y[index, ]
  
    beta[, i] <- ginv(t(sample_X) %*% sample_X) %*% t(sample_X) %*% sample_Y
  }

  # get 2.5, 50, 97.5 percentile of beta
  CI <- apply(beta, 1, quantile, probs = c(0.025, 0.5, 0.975))
  
  return(CI)
}

a <- Boot_nonpar(B=1000, Y=Y, X=X_np)
a[3,]-a[1,]


# Part II: implementation of Non-parametric Bootstrap 
Boot_par <- function(B=1000, Y=Y, X=X_np) {
  set.seed(597)
  # get the size of the orignal data 
  N <- nrow(X)
  p <- ncol(X)
  
  # matrix to store beta estimate 
  beta <- matrix(rep(0, B*p), nrow = p)  
  
  # calculate b_hat and sigma_hat from raw data 
  b_hat <- ginv(t(X) %*% X) %*% t(X) %*% Y
  sigma_hat <- sqrt(sum((Y- X %*% b_hat)^2) / (N-p))
  
  for (i in 1:B) {
    # simulate bootstrp data based on model y*
    sample_Y <- X %*% b_hat + rnorm(n = N, mean = 0, sd = sigma_hat)
    
    
    beta[, i] <- ginv(t(X) %*% X) %*% t(X) %*% sample_Y
    
  }
  
  
  # get 2.5, 50, 97.5 percentile of beta
  CI <- apply(beta, 1, quantile, probs = c(0.025, 0.5, 0.975))
  
  return(CI)

}

c <- Boot_par(B=1000, Y=Y, X=X_np)
c[3,] -c[1,]








