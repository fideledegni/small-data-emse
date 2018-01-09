#####################################################
#  Compte rendu de TP Gaussian Process Regression   #
#              DEGNI Fidèle                         #
#####################################################

rm(list=ls()) # to clear the environment
#setwd('~/Dropbox/ICM-EMSE/0_3A/Data_Science/4_Small_Data/gaussian_process/TP_GPR')
setwd('C:/Users/Fidèle DEGNI/Dropbox/ICM-EMSE/0_3A/Data_Science/4_Small_Data/gaussian_process/TP_GPR')

#### loading some packages and functions ####
library("plot3D")
library("MASS")
source("kernFun.R")

#1 Kernels functions created in kernFun.R


#2 and 3
x <- seq(0, 1, 0.01) # regular grid
# covariance parameters
param <- c(1,1) # c(1,0.2), c(0.5,0.2), c(1,0.1), c(1,0.5), c(1,1)
k1 <- expKern(x, x, param) # computing the covariance matrix using an exp. kernel
k2 <- brownKern(x,x,param)
k3 <- sincKern(x,x,param)
image2D(k1, theta = 90, xlab = "x", ylab = "y", main = "Plotting the expKern covariance matrix") # plotting the covariance matrix
image2D(k2, theta = 90, xlab = "x", ylab = "y", main = "Plotting the brownKern covariance matrix")
image2D(k3, theta = 90, xlab = "x", ylab = "y", main = "Plotting the sincKern covariance matrix")
# Q: what can you observe from the covariance matrix?
# R: We observe that the covariance matrix is symetric, which is normal and expected!

#?mvrnorm # using the help from RStudio
## to complete  ##
## simulating some samples using the "mvrnorm" function
samples1 <- t(mvrnorm(n=10, mu= rep(0,nrow(k1)),Sigma = k1))
samples2 <- t(mvrnorm(n=10, mu= rep(0,nrow(k2)),Sigma = k2))
samples3 <- t(mvrnorm(n=10, mu= rep(0,nrow(k3)),Sigma = k3))
#?matplot # a function to plot the samples. The samples are indexed by columns
matplot(x,samples1[,seq(5)],type = "l", lwd = 2, xlab = "x", ylab = "y", 
        main = expression(paste("expKern sample paths with parameter = ", c(1,1))))
matplot(x,samples2[,seq(5)],type = "l", lwd = 2, xlab = "x", ylab = "y", 
        main = expression(paste("brownKern sample paths with parameter = ", c(1,1))))
matplot(x,samples3[,seq(5)],type = "l", lwd = 2, xlab = "x", ylab = "y", 
        main = expression(paste("sincKern sample paths with parameter = ", c(1,1))))
# Q: what can you observe from the samples?
# R: we observe that the sample paths are different depending on the kernel
#   and the parameters


#4
# We use the brownian kernel : brownKern
samples <- t(mvrnorm(n=200, mu= rep(0,nrow(k2)),Sigma = k2))
pos_vector1 <- samples[12,]
pos_vector2 <- samples[18,]
pos_vector3 <- samples[87,]

plot(pos_vector1, col = "red", pch = 20, xlab = "sample rank", ylab = "values", 
     main = "Cloud of points for close by input points\n Position 12 in red and position 18 in blue")
points(pos_vector2, col = "blue", pch = 20)

plot(pos_vector1, col = "red", pch = 20, xlab = "sample rank", ylab = "values", 
     main = "Cloud of points for far away input points\n Position 12 in red and position 87 in blue")
points(pos_vector3, col = "blue", pch = 20)




#5 Kernels funcstion created in kernFun.R
x1 <- seq(0, 4*pi, 0.01) # regular grid
x2 <- seq(-20, 20, 0.1)
# covariance parameters
param <- c(1,1) # c(1,0.2), c(0.5,0.2), c(1,0.1), c(1,0.5), c(1,1)
kk <- pi_periodicKern(x1, x1, param) # computing the covariance matrix using an exp. kernel
kkk <- sym_1Kern(x2, x2, param)

## simulating some samples using the "mvrnorm" function
sampleskk <- t(mvrnorm(n=10, mu= rep(0,nrow(kk)),Sigma = kk))
sampleskkk <- t(mvrnorm(n=10, mu= rep(0,nrow(kkk)),Sigma = kkk))
matplot(x1,sampleskk[,seq(5)],type = "l", lwd = 2, xlab = "x", ylab = "y", 
        main = expression(paste("pi_periodicKern sample paths with parameter = ", c(1,1))))
matplot(x2,sampleskkk[,seq(5)],type = "l", lwd = 2, xlab = "x", ylab = "y", 
        main = expression(paste("sym_1Kern sample paths with parameter = ", c(1,1))))

# Large number of samples with pi_periodicKern
sampleskk2 <- t(mvrnorm(n=200, mu= rep(0,nrow(kk)),Sigma = kk))
pos_vectorkk1 <- sampleskk2[12,]
pos_vectorkk2 <- sampleskk2[18,]
pos_vectorkk3 <- sampleskk2[87,]

plot(pos_vectorkk1, col = "red", pch = 20, xlab = "sample rank", ylab = "values", 
     main = "Cloud of points for close by input points\n Position 12 in red and position 18 in blue")
points(pos_vectorkk2, col = "blue", pch = 20)

plot(pos_vectorkk1, col = "red", pch = 20, xlab = "sample rank", ylab = "values", 
     main = "Cloud of points for far away input points\n Position 12 in red and position 87 in blue")
points(pos_vectorkk3, col = "blue", pch = 20)

# Large number of samples with sym_1Kern
sampleskkk2 <- t(mvrnorm(n=200, mu= rep(0,nrow(kkk)),Sigma = kkk))
pos_vectorkkk1 <- sampleskkk2[12,]
pos_vectorkkk2 <- sampleskkk2[18,]
pos_vectorkkk3 <- sampleskkk2[87,]

plot(pos_vectorkkk1, col = "red", pch = 20, xlab = "sample rank", ylab = "values", 
     main = "Cloud of points for close by input points\n Position 12 in red and position 18 in blue")
points(pos_vectorkkk2, col = "blue", pch = 20)

plot(pos_vectorkkk1, col = "red", pch = 20, xlab = "sample rank", ylab = "values", 
     main = "Cloud of points for far away input points\n Position 12 in red and position 87 in blue")
points(pos_vectorkkk3, col = "blue", pch = 20)



#6 Bonus



#7
x <- seq(0, 1, 0.01) # regular grid
f <- function(x){
  return(x + sin(4*pi*x))
}
n <- 20
X <- runif(n)
Y <- f(X)

#8
m <- function(x, X, Y, kern, param){
  return(kern(x,X, param)%*%solve(kern(X,X, param))%*%Y)
}

cc <- function(x, y, X, kern, param){
  return(kern(x,y, param) - kern(x,X, param)%*%solve(kern(X,X, param))%*%kern(X,y, param))
}

#9 and 10
#kern <- brownKern
kern <- brownKern
param <- c(1, 0.2)
mx <- m(x, X, Y, kern, param)
intx = x
for (i in 1:length(x)) {
  intx[i] = 1.96*sqrt(cc(x[i], x[i], X, kern, param))
}
plot(x, f(x), col = "blue", pch = 20, 
     main = "f(x) in blue, m(x) in red with 95% confidence intervals\n Experiment points in green")
points(x, mx, col = "red", pch = 20, lwd = 1)
points(X, f(X), col = "green", pch = 24)
points(x, mx+intx, col = "gray", type = "l")
points(x, mx-intx, col = "gray", type = "l")


#11 Generate samples from the conditional process.
covx = matrix(0, nrow = length(x), ncol = length(x))
for (i in 1:length(x)) {
  for (j in 1:length(x))
    covx[i,j] <- cc(x[i], x[j], X, kern, param)
}
samples2 <- t(mvrnorm(n=10, mu= mx, Sigma = covx))


#12 Bonus


