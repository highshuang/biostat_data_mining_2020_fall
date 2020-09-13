####### HW2 ##########
library(pls)
library(MASS)
library(metaSEM)
library(matlib)
library(glmnet)
library(robustHD)
library(clusterGeneration)

set.seed(1)
#### Question1: Implement cyclic gradient descent algorithm for LASSO

n <- 50
p <- 20

# generate random positive definite matrix for covariance matrix 
pdmat <- genPositiveDefMat(dim = p)$Sigma

# simulation data 
X <- mvrnorm(n = n, mu = rep(0, p), Sigma = pdmat)
b_true <- rnorm(p, mean = 0, sd = 2)
Y <- X %*% b_true + rnorm(n, mean = 0, sd = 0.01)

lasso <- function(X, Y, n_iter=100, lambda=5) {
  # initialization 
  n <- length(Y)
  
  # add intercept to design matrix
  X <- cbind(rep(1, n), X)
  
  p <- dim(X)[2]
  
  # initialize beta with all 0
  b <- matrix(rep(0, p), nrow = p)
  
  # calculate the numerator = 1/n * X_j^2 outside the loop
  numerator <- apply(X, MARGIN=2, FUN = function(a) sum(a^2)/n)
  
  for (iter in 1:n_iter) {
    for (j in 1:p) {
      
      # store sum[Y_i - sum(Xij' b_j')]
      R <- t(as.vector(X[, j])) %*% (Y - X[, -j] %*% b[-j]) / n 
      
      # update intercept 
      if (j == 1) {
        b[j] = R 
        next 
      }

      if (R > lambda) {
        # update b_positive
        b[j] <- (R - lambda)/numerator[j]
      } else if (R < -lambda) {
        # update b_negative
        b[j] <- (R+lambda)/numerator[j]
      } else {
        b[j] <- 0
      }
    }
  }
  
  # return 
  return(b)
  
}
# standardize test data 
standardize_data <- function(X) {
  # intialize 
  n <- length(X[,1])
  
  X_mean <- colMeans(X)
  
  X_centered <- t(t(X) - X_mean)
  
  X_sd <- apply(X_centered, 2,  function(x) sqrt(sum(x^2)/n))
  
  X_stand <- t(t(X_centered) / X_sd)
  
  return(X_stand)
}



# LASSO with standardized X 
b1.1 <- lasso(standardize_data(X), Y, lambda = 3, n_iter = 100)
print(b1.1)


# glmnet() with standardized X
fit1.1 <- glmnet(standardize_data(X), Y, lambda = 3, 
              standardize = FALSE)
print(coef(fit1.1))

# comparison between self defined LASSO and glmnet: (difference is 0.00033)
print(paste0("The difference between self defined LASSO and glmnet with standardized input is  ", 
             sum(abs(b1.1 - as.vector(coef(fit1.1))))))


# LASSO with original X 
b1.2 <- lasso(X, Y, lambda = 3, n_iter = 100)
print(b1.2)


# glmnet() with original X 
fit1.2 <- glmnet(X, Y, lambda = 3, standardize = FALSE)
print(coef(fit1.2))

# comparison between self defined LASSO and glmnet: (difference is 0.0026)
print(paste0("The difference between self defined LASSO and glmnet with standardized input is  ", 
             sum(abs(b1.2 - as.vector(coef(fit1.2))))))







#### Question 2: Partial least square

PLS_self <- function(X, Y) {
  # initialization 
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  Z <- matrix(rep(0, n*p), ncol = p)
  
  for (m in 1:p) {
    
    # step a:
    phi_m <- as.vector(t(X) %*% Y)  
    Z[, m] <- rowSums(X %*% diag(phi_m))
    
    # step b:
    for (j in 1:p) {
       # update Xj(m)
       X[, j] <- X[, j] - as.vector(Z[, m] %*% X[, j] / Z[, m] %*% Z[, m]) * Z[ ,m]
    }
    
  }
  return(Z)
  
}


# test data 
X <- scale(as.matrix(mtcars[, -1]))
Y<- scale(as.matrix(mtcars[, 1]))

res <- PLS_self(a, b)

fit1 <- lm(b~res[,c(1:3)])
summary(unlist(fit1))

fit <- plsr(b~a, scale = FALSE)

fit$loadings
summary(fit$coefficients)
m1 <- fit$scores[,1]
sum(standardize(res)[,1]^2)

sum(m1^2)

fit$scores[,1] %*% fit$scores[,2]
t(res) %*% res
t(res_no) %*% res_no

res_no <- apply(res, 2, normal_vec)

fit_no <- apply(fit$scores, 2, normal_vec)

t(fit_no) %*% fit_no

t(res) %*% X[, 10] 


