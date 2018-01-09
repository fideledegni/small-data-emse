rm(list=ls()) # to clear the environment

#### loading some packages and functions ####
library("plot3D")
library("MASS")
source("kernFun.R")

#### Example with the Exp. kernel  ####
x <- seq(0, 1, 0.01) # regular grid
param <- c(1,0.2) # covariance parameters
k1 <- expKern(x, x, param) # computing the covariance matrix using an exp. kernel
k2 <- brownKern(x,x,param)
k3 <- sincKern(x,x,param)

image2D(k1, theta = 90, xlab = "x", ylab = "y") # plotting the covariance matrix
image2D(k2, theta = 90, xlab = "x", ylab = "y") # plotting the covariance matrix
image2D(k3, theta = 90, xlab = "x", ylab = "y") # plotting the covariance matrix

# Q: what can you observe from the covariance matrix?
?mvrnorm # using the help from RStudio

## to complete  ##
## simulating some samples using the "mvrnorm" function
samples <- t(mvrnorm(n=100, mu= rep(0,nrow(k3)),Sigma = k3))
# ?matplot # a function to plot the samples. The samples are indexed by columns
# Q: what can you observe from the samples?
matplot(x,samples[,seq(7)],type = "l")
