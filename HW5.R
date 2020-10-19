#### HW5 #######
library(quadprog)
library(ggplot2)

set.seed(1)

# simulated data with separable case 
# generate data 
x1 <- c(rnorm(n = 30, mean = 2, sd = 0.5), 
        rnorm(n = 30, mean = 4, sd = 0.5))

x2 <- c(rnorm(n = 30, mean = 2, sd = 0.5), 
        rnorm(n = 30, mean = 4, sd = 0.5))

data <- data.frame(x1, x2)
X <- cbind(x1, x2)

y <- c(rep(-1, 30), rep(1, 30))
data <- cbind(x1, x2, y)

plot(data[,1], data[,2])

# SVM dual formulation: min 0.5*t(a)*D*a - sum(a) 
# constraints: t(y)*a = 0 and all a >= 0 

# calculate D matrix with Dij = yiyjt(xj)xi symmetric 
D <- matrix(0, 60, 60)

# fill in the symmetric matrix D 
for (i in 1:60) {
  for (j in 1:i) {
    D[i, j] <- sum(data[i, -3]*data[j, -3]) * data[i, 3] * data[j, 3]
    D[j, i] <- D[i, j]
  }
}

# solve with quadratic solver
# solve.QP(Dmat, dvec, Amat, bvec, meq=0, factorized=FALSE)
dvec <- rep(1, 60)

# add perturbance to work around the semi-positive definite matrix D, 
# which cannot be calculated through quadprog
eps <- 1e-7
Dmat <- D + diag(rep(eps, 60))
is.positive.definite(Dmat)

# set equality constrait: t(y)*a = 0
A.eq <- rbind(data[,3])
b.eq <- 0

A.ineq <-diag(60) 
b.ineq <- rep(0, 60)

sol <- solve.QP(Dmat = Dmat,
                dvec = dvec,
                Amat = t(rbind(A.eq, A.ineq)),
                bvec = c(b.eq, b.ineq),
                meq = 1) # this argument says the first "meq" rows of Amat are equalities

# calculate intercept 
alpha <- sol$solution

# get the data point where alpha > 0 and thus the support vectors
sv.index <- which(alpha > 1e-6)
alpha_pos <- alpha[alpha>1e-6]

# support vector samples
sv <- data[c(sv.index), ]

# calculate beta from alpha > 0
b <- t(sv[,-3]) %*% (alpha_pos * sv[, 3])

# vector length of b 
len_b <- sqrt(sum(b^2))

# calculate b0 for each sample 
b0_all <- sv[, 3] - sv[, -3] %*% b

# average all the b0 
b0 <- mean(b0_all)


# label the support vector 
group <- c(rep(0, 60))
group[c(8, 54)] <- 1

a <- data.frame(cbind(data, group))
ggplot(data = a, aes(x = x1, y = x2, color = as.factor(group))) +
  geom_point() +
  geom_abline(intercept = -b0/b[2,1], slope = -b[1,1]/b[2,1]) +
  geom_abline(intercept = -(b0-1)/b[2,1], slope = -b[1,1]/b[2,1], color = "lightblue") +
  geom_abline(intercept = -(b0+1)/b[2,1], slope = -b[1,1]/b[2,1], color = "lightblue") +
  labs(color = "Support Vector") +
  theme_minimal() + 
  labs(title = "Support Vector Classifier with Separable Case")



