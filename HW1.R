########## PHS597 HW1: simulation ##########
####### Author: Shuang Gao ##########
library(MASS)
library(metaSEM)
library(matlib)

set.seed(1)

# simulate data 
n <- 100
p <- 4
sigma <- vec2symMat(c(1, 0.1, 0.1, 1, 0.1, 1))

# add 1s to design matrix 
X_np <- cbind(rep(1,n), mvrnorm(n = n, mu = rep(0, p-1), Sigma = sigma))

b_true <- rnorm(p, mean = 0, sd = 2)

Y <- X_np %*% b_true + rnorm(n, mean = 0, sd = 0.01)

# Question 1: 2-stage regression verification 

# function return Least square estiamtes and residuals 
lm_ls <- function(Y, X_np) {
  if (ncol(X_np)>1) {
    # use generialized inverse 
    b_hat <- Ginv(t(X_np) %*% X_np) %*% t(X_np) %*% Y 
    residual <- Y - X_np %*% b_hat 
    
    return(list("beta_hat" = b_hat, "residual" = residual))
  }
  # case when t(X_np) %*% X_np dim 1x1 and cause error in Ginv()
  
  b_hat <- 1/ (t(X_np) %*% X_np) * t(X_np) %*% Y 
  
  residual <- Y - X_np %*% b_hat 
  
  return(list("beta_hat" = b_hat, "residual" = residual))

}

# step 1: regress Y over X1,..., X_p-1 and store resiudal
reg1 <- lm_ls(Y, X_np[, -p])
res1 <- reg1$residual

# step 2: regress Xp over X1,..., Xp-1 and store residual
reg2 <- lm_ls(X_np[, p], X_np[, -p])
res2 <- reg2$residual 

# step 3: regress residual 1 over resiudal 2
reg3 <- lm_ls(res1, cbind(rep(1, n), res2))

# compare with the multivariate linear regression 
reg4 <- lm_ls(Y, X_np)

# comparison 
print(paste0("The result from multivaraite linear regression is ", reg4$beta_hat[p]))
print(paste0("The result from 2-stage linear regression is ", reg3$beta_hat[2]))
print(paste0("The difference between these two estiamtes is ", 
             all.equal.numeric(reg4$beta_hat[p], reg3$beta_hat[2]),
             ", which is small enough to be considered the same."))


### Question 2: successive orthogonalization verification 

# remove the intercept col in the design matrix 
X2 <- X_np[, -1]
p2 <- 3

# initialize Z matrix with dim n by p
Z <- matrix(c(rep(1, n), rep(0, p2*n)), ncol = p2+1)

for (j in 1:p2) {
  # regress xj on z0,..., zj-1
  result<- lm_ls(X2[, j], as.matrix(Z[, 1:j]))
  
  # store residuals 
  Z[, (j+1)] <- result$residual 
}

# regress y on Zp
res_Q2 <- lm_ls(Y, as.matrix(Z[, p2+1]))


# comparison for Q2
print(paste0("The result from multivaraite linear regression is ", reg4$beta_hat[p]))
print(paste0("The result from successive orthogalinization is ", res_Q2$beta_hat))
print(paste0("The difference between these two estiamtes is ", 
             all.equal.numeric(reg4$beta_hat[p], as.numeric(res_Q2$beta_hat)),
             ", which is small enough to be considered the same."))




















