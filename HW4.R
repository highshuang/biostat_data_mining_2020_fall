###### HW4 ########
library(mvtnorm)
library(pracma)
library(msos)
library(matrixcalc)
library(ggplot2)
library(matlib)

set.seed(1)

# generate  data for 3 classes 
K <- 3
n1 <- 100
n2 <- 100
n3 <- 100
N <- n1 + n2+ n3

# each class with common covariance matrix I*0.5^2
x1 <- c(rnorm(n = n1, mean = 2, sd = 0.5), 
        rnorm(n = n2, mean = 4, sd = 0.5),
        rnorm(n = n3, mean = 6.5, sd = 0.5))

x2 <- c(rnorm(n = n1, mean = 2, sd = 0.5), 
        rnorm(n = n2, mean = 4, sd = 0.5),
        rnorm(n = n3, mean = 6, sd = 0.5))

data <- data.frame(x1, x2)

# add class label 
data$y <- c(rep(1, n1), rep(2, n2), rep(3, n3))

# PART 1: Implemente of LDA
## estimate parameters 
p1 <- n1 / N
p2 <- n2 / N
p3 <- n3 / N

# mean for each class 
u1 <- colMeans(data[c(1:n1), c(1:2)]) 
u2 <- colMeans(data[c((n1+1):(n1+n2)), c(1:2)]) 
u3 <- colMeans(data[c((n1+n2+1):N), c(1:2)]) 

# create mean matrix 
M1 <- as.matrix(rep(1, n1), ncol = 1) %*% u1
M2 <- as.matrix(rep(1, n2), ncol = 1) %*% u2
M3 <- as.matrix(rep(1, n3), ncol = 1) %*% u3
M <- rbind(M1, M2, M3)

# common cov 
X <- as.matrix(data[, c(1:2)])
Sig <- t(X-M) %*% (X-M) / (N-K)

# decision boundary 
# class 1 vs class2 
a1 <- log(p1/p2) - 0.5*t(u1+u2) %*% ginv(Sig) %*% (u1-u2)
a2 <- (ginv(Sig) %*% (u1-u2))[1]
a3 <- (ginv(Sig) %*% (u1-u2))[2]

# class 2 vs class3 
b1 <- log(p2/p3) - 0.5*t(u2+u3) %*% ginv(Sig) %*% (u2-u3)
b2 <- (ginv(Sig) %*% (u2-u3))[1]
b3 <- (ginv(Sig) %*% (u2-u3))[2]

# class 1 vs class3 
c1 <- log(p1/p3) - 0.5*t(u1+u3) %*% ginv(Sig) %*% (u1-u3)
c2 <- (ginv(Sig) %*% (u1-u3))[1]
c3 <- (ginv(Sig) %*% (u1-u3))[2]


# plot 
ggplot(data, aes(x = x1, y = x2, color = as.factor(y))) +
  geom_point() +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  theme_minimal() +
  geom_abline(intercept = -a1/a3, slope = -a2/a3) +
  geom_abline(intercept = -b1/b3, slope = -b2/b3) +
  geom_abline(intercept = -c1/c3, slope = -c2/c3) +
  labs(title = "Implementation of LDA")

# PART 2: Implemente of QDA
# sigma_k estimation 
Sig1 <- t(X[c(1:n1),]-M[c(1:n1),]) %*% (X[c(1:n1),]-M[c(1:n1),]) / (n1-1)

Sig2 <- t(X[c((n1+1):(n1+n2)),]-M[c((n1+1):(n1+n2)),]) %*% (X[c((n1+1):(n1+n2)),]-M[c((n1+1):(n1+n2)),]) / (n2-1)

Sig3 <- t(X[c((n1+n2+1):N),]-M[c((n1+n2+1):N),]) %*% (X[c((n1+n2+1):N),]-M[c((n1+n2+1):N),]) / (n3-1)

# decision boundary 
# class 1 vs class2 
res11 <- ginv(Sig2) - ginv(Sig1)
res12 <- t(u2) %*% ginv(Sig2) - t(u1) %*% ginv(Sig1)
res13 <- 0.5*t(u2) %*% ginv(Sig2)%*% u2 - 0.5*t(u1) %*% ginv(Sig1)%*% u1 + 
  0.5*logdet(Sig2)-0.5*logdet(Sig1) + log(p1/p2)

# general equation for conic section 
a1 <- 0.5*res11[1,1]
a2 <- 0.5*res11[2,2]
a3 <- 0.5*(res11[1,2]+res11[2,1])
a4 <- -res12[1]
a5 <- -res12[2]
a6 <- res13

# get conic sections for each level 
d1 <- transform(expand.grid(x=seq(0,10,length=300),
                           y=seq(0,10,length=300)),z=a1*x^2+a3*x*y+a4*x+a5*y+a2*y^2+as.vector(a6))

# class 2 vs class3 
res21 <- ginv(Sig3) - ginv(Sig2)
res22 <- t(u3) %*% ginv(Sig3) - t(u2) %*% ginv(Sig2)
res23 <- 0.5*t(u3) %*% ginv(Sig3)%*% u3 - 0.5*t(u2) %*% ginv(Sig2)%*% u2 + 0.5*logdet(Sig3)-0.5*logdet(Sig2) + log(p2/p3)

# general equation for conic section 
a1 <- 0.5*res21[1,1]
a2 <- 0.5*res21[2,2]
a3 <- 0.5*(res21[1,2]+res21[2,1])
a4 <- -res22[1]
a5 <- -res22[2]
a6 <- res23

d2 <- transform(expand.grid(x=seq(0,10,length=300),
                            y=seq(0,10,length=300)),z=a1*x^2+a3*x*y+a4*x+a5*y+a2*y^2+as.vector(a6))

# class 1 vs class3 
res31 <- ginv(Sig3) - ginv(Sig1)
res32 <- t(u3) %*% ginv(Sig3) - t(u1) %*% ginv(Sig1)
res33 <- 0.5*t(u3) %*% ginv(Sig3)%*% u3 - 0.5*t(u1) %*% ginv(Sig1)%*% u1 + 0.5*logdet(Sig3)-0.5*logdet(Sig1) + log(p1/p3)

# general equation for conic section 
a1 <- 0.5*res31[1,1]
a2 <- 0.5*res31[2,2]
a3 <- 0.5*(res31[1,2]+res31[2,1])
a4 <- -res32[1]
a5 <- -res32[2]
a6 <- res33

d3 <- transform(expand.grid(x=seq(0,10,length=300),
                            y=seq(0,10,length=300)),z=a1*x^2+a3*x*y+a4*x+a5*y+a2*y^2+as.vector(a6))


# plot 
ggplot(data, aes(x = x1, y = x2, color = as.factor(y))) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#9999CC"))+
  labs(title = "Implementation of QDA") +
  stat_contour(aes(x,y,z=z,colour=factor(..level..==0)),data=d1,breaks = c(0), linetype = 1) +
  stat_contour(aes(x,y,z=z,colour=factor(..level..==0)),data=d2,breaks = c(0), linetype = 2) +
  stat_contour(aes(x,y,z=z,colour=factor(..level..==0)),data=d3,breaks = c(0), linetype = 6) +
  theme(legend.position = "none")



# Part 3: LDA by Sphering data 
set.seed(1)
# generate data for 2 classes 
m1 <- mvrnorm(n = 200, mu = c(1, 3), Sigma = matrix(c(1, -0.6, -0.6, 1)*0.5, ncol = 2))
m2 <- mvrnorm(n = 200, mu = c(4.1, 3), Sigma = matrix(c(1, -0.6, -0.6, 1)*0.5, ncol = 2))

data2 <- data.frame(rbind(m1, m2))

data2$y <- c(rep(1, 200), rep(2, 200))

# mean for each class 
u1 <- colMeans(m1)
u2 <- colMeans(m2)


# create mean matrix 
M1 <- as.matrix(rep(1, 200), ncol = 1) %*% u1
M2 <- as.matrix(rep(1, 200), ncol = 1) %*% u2
M <- rbind(M1, M2)

# common cov 
X <- as.matrix(data2[, -3])
Sig <- t(X-M) %*% (X-M) / (400-2)

# sphering the data with common sigma 
U <- eigen(Sig)$vectors
D_diag <- eigen(Sig)$values

# transform the original data 
x_new <- as.matrix(data2[, -3]) %*% U %*% diag(D_diag^(-0.5))
mu1_new <- diag(D_diag^(-0.5)) %*% t(U) %*% u1
mu2_new <- diag(D_diag^(-0.5)) %*% t(U) %*% u2

# decision boundary ax1+bx2+c=0
c <- log(50/50) - 0.5* t(mu1_new + mu2_new) %*% (mu1_new - mu2_new)
a <- (mu1_new - mu2_new)[1]
b <- (mu1_new - mu2_new)[2]

# the transformed decision boundary 
# slope 
-a/b

# intercept 
-c/b

# transform the decison boundary to the original space 
dec <- U %*% diag(D_diag^(-0.5)) %*% as.vector(c(a, b))

# plot 
  ggplot() +
    geom_point(aes(x = data2[,1], y = data2[,2], color = as.factor(data2[,3]))) +
    geom_point(aes(x = x_new[,1], y = x_new[,2], color = as.factor(data2[,3]))) +
    scale_color_manual(values = c("#00AFBB", "#E7B800")) +
    theme_minimal() +
    geom_abline(intercept = -c/dec[2,1], slope = -dec[1,1]/dec[2,1]) +
    geom_abline(intercept = -c/b, slope = -a/b) +
    labs(title = "LDA via sphering")

  
  
# Part 4: QDA with sphering 
  set.seed(1)
  # generate data for 2 classes 
  m1 <- mvrnorm(n = 200, mu = c(1, 3), Sigma = matrix(c(1, -0.6, -0.6, 1)*0.5, ncol = 2))
  m2 <- mvrnorm(n = 200, mu = c(4.1, 3), Sigma = matrix(c(0.4, 0, 0, 0.4), ncol = 2))
  
  data3 <- data.frame(rbind(m1, m2))
  
  data3$y <- c(rep(1, 200), rep(2, 200))
  
  # mean for each class 
  u1 <- colMeans(m1)
  u2 <- colMeans(m2)
  
  
  # create mean matrix 
  M1 <- as.matrix(rep(1, 200), ncol = 1) %*% u1
  M2 <- as.matrix(rep(1, 200), ncol = 1) %*% u2
  M <- rbind(M1, M2)
  
  # common cov 
  X <- as.matrix(data3[, -3])
  
# cov matrix for each class
Sig1 <- t(X[c(1:200), ]-M[c(1:200), ]) %*% (X[c(1:200),] -M[c(1:200), ]) / (200-1)
Sig2 <- t(X[c(201:400), ]-M[c(201:400), ]) %*% (X[c(201:400),] -M[c(201:400), ]) / (200-1)

U1 <- eigen(Sig1)$vectors
D1_diag <- eigen(Sig1)$values 

U2 <- eigen(Sig2)$vectors
D2_diag <- eigen(Sig2)$values 


X_new <- rbind(as.matrix(data2[c(1:200), -3]) %*% U1 %*% diag(D1_diag^(-0.5)), 
                as.matrix(data2[c(201:400), -3]) %*% U2 %*% diag(D2_diag^(-0.5)))

mu1_new <- diag(D1_diag^(-0.5)) %*% t(U1) %*% u1
mu2_new <- diag(D2_diag^(-0.5)) %*% t(U2) %*% u2

# decision boundary ax1+bx2+c=0
c <- -0.5*sum(log(D1_diag)) + 0.5*sum(log(D2_diag)) - 0.5*sum(mu1_new^2) + 0.5*sum(mu2_new^2) + log(50/50)
a <- (mu1_new - mu2_new)[1]
b <- (mu1_new - mu2_new)[2]


# transform back 
# class 1 U1 D1
t1 <- U1 %*% diag(D1_diag^(-0.5)) %*% as.vector(c(a,b))

t2 <- U2 %*% diag(D2_diag^(-0.5)) %*% as.vector(c(a,b))

data_trans <- data.frame(cbind(X_new, c(rep(1,200), rep(2,200))))

# plot
ggplot() +
  geom_point(aes(x = data3[,1], y = data3[,2], color = as.factor(data3[,3]))) +
  geom_point(aes(x = data_trans[,1], y = data_trans[,2], color = as.factor(data_trans[,3]))) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  #theme_minimal() +
  geom_abline(intercept = -c/t1[2,1], slope = -t1[1,1]/t1[2,1]) +
  geom_abline(intercept = -c/t2[2,1], slope = -t2[1,1]/t2[2,1]) 
  labs(title = "QDA via sphering")









  
  
  
  


