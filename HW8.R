##### PHS597 HW8
library(MASS)
library(splines)


set.seed(1)

# simulate data 
n <- 100
x <- rnorm(n, 0, 1)
y <- exp(x)/(1+exp(x)) + rnorm(n, 0, 1)
data <- data.frame(x,y)
# break points 0 -1 1 
knot <- c(-1, 0, 1)
plot(x, y)

### Part I: Cubic spline 
# consider it as f(x) spanned by 1, X, X^2, X^3, (X-)

# create the design matrix 
Z <- data.frame(cbind(rep(1, n), x, x^2, x^3))
Z$knot1 <- ifelse(Z$x>= knot[1], (Z$x-knot[1])^3, 0)
Z$knot2 <- ifelse(Z$x>= knot[2], (Z$x-knot[2])^3, 0)
Z$knot3 <- ifelse(Z$x>= knot[3], (Z$x-knot[3])^3, 0)

Z <- as.matrix(Z)

# least square estimates 
beta1 <- ginv(t(Z) %*% Z) %*% t(Z) %*% y

# interval[-inf, -1]: beta0+beta1X+beta2X^2+beta3X^3
x1 <- seq(-2.5, -1, by = 0.05)
n_x1 <- length(x1)
z1 <- cbind(rep(1, n_x1), x1, x1^2, x1^3)
y1 <- z1 %*% beta1[1:4]
data1 <- data.frame("x"= x1, "y"=y1)

# interval [-1, 0]:  (beta0+beta4)+(beta1-3*beta4)X+(beta2+3*beta4)X^2+(beta3+beta4)X^3
x2 <- seq(-1, 0, by = 0.05)
n_x2 <- length(x2)
z2 <- cbind(rep(1, n_x2), x2, x2^2, x2^3, (x2+1)^3)
y2 <- z2 %*% beta1[1:5]
data2 <- data.frame("x"= x2, "y"=y2)

# interval [0, 1]
x3 <- seq(0, 1, by = 0.05)
n_x3 <- length(x3)
z3 <- cbind(rep(1, n_x3), x3, x3^2, x3^3, (x3+1)^3, (x3^3))
y3 <- z3 %*% beta1[1:6]
data3 <- data.frame("x"= x3, "y"=y3)

# interval [1, +inf]
x4 <- seq(1, 2.5, by = 0.05)
n_x4 <- length(x4)
z4 <- cbind(rep(1, n_x4), x4, x4^2, x4^3, (x4+1)^3, (x4^3), (x4-1)^3)
y4 <- z4 %*% beta1
data4 <- data.frame("x"= x4, "y"=y4)


# the plot for cubic spline 
ggplot(data = data,aes(x, y)) +
  geom_point() +
  geom_line(data = data1, color = 'red') +
  geom_line(data = data2, color = 'blue') +
  geom_line(data = data3, color = 'green') +
  geom_line(data = data4, color = 'purple') +
  theme_minimal() +
  ggtitle("Cublic splines fitting")


### Part II: natural cublic spline 
# consider it as f(x) spanned by 1, X, (X+1)+^3-(X-1)+^3, (X)+^3-(X-1)+^3

# create the design matrix 
M <- data.frame(cbind(rep(1, n), x))
M$knot1 <- ifelse(M$x>= knot[1], (M$x-knot[1])^3, 0)
M$knot2 <- ifelse(M$x>= knot[2], (M$x-knot[2])^3, 0)
M$knot3 <- ifelse(M$x>= knot[3], (M$x-knot[3])^3, 0)
M$factor3 <- M$knot1 - M$knot3 - 2*M$knot2 + 2*M$knot3

M <- as.matrix(M[, -c(3:5)])

# least square estimates 
beta2 <- ginv(t(M) %*% M) %*% t(M) %*% y

# interval[-inf, -1]: beta1+beta2X
x1 <- seq(-2.5, -1, by = 0.05)
n_x1 <- length(x1)
z1 <- cbind(rep(1, n_x1), x1)
y1 <- z1 %*% beta2[1:2]
data1 <- data.frame("x"= x1, "y"=y1)


# interval [-1, 0]: 
x2 <- seq(-1, 0, by = 0.05)
n_x2 <- length(x2)
z2 <- cbind(rep(1, n_x2), x2, (x2+1)^3)
y2 <- z2 %*% beta2
data2 <- data.frame("x"= x2, "y"=y2)

# interval [0, 1]
x3 <- seq(0, 1, by = 0.05)
n_x3 <- length(x3)
z3 <- cbind(rep(1, n_x3), x3, (x3+1)^3-2*x3^3)
y3 <- z3 %*% beta2
data3 <- data.frame("x"= x3, "y"=y3)

# interval [1, +inf]
x4 <- seq(1, 2.5, by = 0.05)
n_x4 <- length(x4)
z4 <- cbind(rep(1, n_x4), x4, (x4+1)^3-(x4-1)^3-2*(x4)^3+2*(x4-1)^3)
y4 <- z4 %*% beta2
data4 <- data.frame("x"= x4, "y"=y4)



# the plot for natural cubic spline 
ggplot(data = data,aes(x, y)) +
  geom_point() +
  geom_line(data = data1, color = 'red') +
  geom_line(data = data2, color = 'blue') +
  geom_line(data = data3, color = 'green') +
  geom_line(data = data4, color = 'purple') +
  theme_minimal() +
  ggtitle("Natural Cublic splines fitting")

### Part III: smoothing splines 
# calculate integral Nk''Nj'' 
int_kj <- function(knot_k, knot_j, knot_K_1, knot_K) {
  
  # pre-define some constants 
  C_k1 <- 1/(knot_k-knot_K)
  C_k2 <- 1/(knot_K_1-knot_K)
  C_j1 <- 1/(knot_j-knot_K)
  C_j2 <- 1/(knot_K_1-knot_K)
    
  # interval 1: < knot_j, integral = 0
  # interval 2: [knot_j, knot_K_1]
  integrand1 <- function(x) {36*C_k1*(x-knot_k)*C_j1*(x-knot_j)}
  a1 <- integrate(integrand1, lower = knot_j, upper=knot_K_1)  
  
  # interval 3: [knot_k_1, knot_K]
  integrand2 <- function(x) {(6*C_k1*(x-knot_k)- 6*C_k2*(x-knot_K_1))*(6*C_j1*(x-knot_j)- 6*C_j2*(x-knot_K_1))}
  a2 <- integrate(integrand2, lower = knot_K_1, upper=knot_K)
  
  # interval 4: >=knot_K
  integrand3 <- function(x) {(6*C_k1*(x-knot_k)- 6*C_k2*(x-knot_K_1) + 6*(C_k2-C_k1)*(x-knot_k))*(6*C_j1*(x-knot_j)- 6*C_j2*(x-knot_K_1)+6*(C_j2-C_j1)*(x-knot_j))}
  a3 <- integrate(integrand3, lower = knot_K, upper=knot_K) # diverge if upper is Inf, thus set upper to 10000, relatively large respect to the simulated data range 
  
  # integral of piecewise function 
  return(a1$value+a2$value+a3$value)
}

# calculate Omega except for the first two basis N1 = 1 and N2 = x 
A <- matrix(rep(0, 98*98), ncol = 98)

# sorted x used as knots for smoothing spline
x_sort <- sort(x)

# fill in A 
for (k in 1:98) {
  for (j in k:98) {

      A[k, j] <- int_kj(knot_k=x_sort[k], knot_j=x_sort[j], knot_K_1=x_sort[(n-1)], knot_K=x_sort[n])
      A[j, k] <- A[i, k]
    
  }
}

# get Omega for integral 
Omega <- matrix(rep(0, n*n), ncol = n)
Omega[3:n, 3:n] <- A

# calcualte basis function N 
N <- matrix(rep(0, n*n), ncol = n)

# N1(x) = 1
N[, 1] <- rep(1, n)
N[, 2] <- x

cal_Nk <- function(x, knot_k, knot_K_1=x_sort[(n-1)], knot_K=x_sort[n]) {
  
  # pre-define some constants 
  C1 <- 1/(knot_k-knot_K)
  C2 <- 1/(knot_K_1-knot_K)
  
  res <- C1* ((max(0, x-knot_k))^3-(max(0, x-knot_K))^3) - C2* ((max(0, x-knot_K_1))^3-(max(0, x-knot_K))^3) 
  
  return(res)
}


# fill in basis matrix N 
for (i in 3:n) {
  for (j in 1:n) {
    N[j, i] <- cal_Nk(x = x_sort[j], knot_k = x_sort[(i-2)])
  }
}

# set a small lambda
lambda <- 0.0001
# calculate theta_hat 
theta_hat <- ginv(t(N)%*%N+lambda*Omega) %*% t(N) %*% y 


# generate data for plotting 
m <- seq(-2.5, 2.5, by = 0.05)
y_m <- rep(0, length(m))
length(m)
# calcualte basis function n
new_N <- matrix(rep(0, 101*100), ncol = 100)

# N1(x) = 1
new_N[, 1] <- rep(1, 101)
new_N[, 2] <- m


# fill in new N 
for (i in 3:length(theta_hat)) {
  for (j in 1:length(m)) {
    new_N[j, i] <- cal_Nk(x = m[j], knot_k = x_sort[(i-2)])
  }
}

y_m <- new_N %*% theta_hat
data_s <- data.frame("x"= m, "y"=y_m)


# the plot Smoothing splines fitting with small lamda show more curves
# larger lambda leads to more linear fit 
ggplot(data = data,aes(x, y)) +
  geom_point() +
  geom_line(data = data_s, color = 'red') +
  theme_minimal() +
  ggtitle(" Smoothing splines fitting with small lamda show more curves") 



###### Part IV: B-spline
ggplot(data, aes(x, y)) +
  geom_point() +
  stat_smooth(method = lm, formula = y ~ splines::bs(x, df = 3)) +
  theme_minimal() +
  ggtitle("B-spline regression fitting")



