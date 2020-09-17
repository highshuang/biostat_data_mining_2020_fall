#### HW3 #########
library(ggplot2)

set.seed(1)
# generate data 
x1 <- c(rnorm(n = 30, mean = 2, sd = 1), 
        rnorm(n = 30, mean = 4, sd = 1),
        rnorm(n = 30, mean = 6.5, sd = 1))

x2 <- c(rnorm(n = 30, mean = 2, sd = 1), 
        rnorm(n = 30, mean = 4, sd = 1),
        rnorm(n = 30, mean = 6, sd = 1))

data <- data.frame(x1, x2)

# create dummy variable for each category 
data$y1 <- c(rep(1, 30), rep(0, 60))
data$y2 <- c(rep(0, 30), rep(1, 30), rep(0, 30))
data$y3 <- c(rep(0, 60), rep(1, 30))

# calculate model parameter estimates from the training data 
X <- as.matrix(cbind(rep(1, 90), data[, c(1:2)]), ncol = 3)
Y <- as.matrix(data[, c(3:5)])
B_hat <- Ginv(t(X) %*% X) %*% t(X) %*% Y

# decision boundary 
# f1 = f2 line 
l1 <- B_hat[, 1] - B_hat[, 2]

# f2 = f3 line 
l2 <- B_hat[, 2] - B_hat[, 3]

# f1 = f3 line
l3 <- B_hat[, 1] - B_hat[, 3]

# prepare group label for categories 
data$group <- c(rep(1, 30), rep(2, 30), rep(3, 30))

ggplot(data, aes(x = x1, y = x2, color = as.factor(group))) +
  geom_point() +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  theme_minimal() +
  geom_abline(intercept = -l1[1]/l1[3], slope = -l1[2]/l1[3]) +
  geom_abline(intercept = -l2[1]/l2[3], slope = -l2[2]/l2[3]) +
  geom_abline(intercept = -l3[1]/l3[3], slope = -l3[2]/l3[3]) 

# graph: the ambiguity shown when classify cluster 2 against other clusters 
  
  


