#####################################################
#  Compte rendu de TP Sensibility analysis          #
#              DEGNI Fidèle                         #
#              RODRIGUES Leticia                    #
#####################################################

rm(list=ls()) #  cleaning up
#setwd('~/Dropbox/ICM-EMSE/0_3A/Data_Science/4_Small_Data/sensitivity_analysis/TP')
setwd('C:/Users/Fidèle DEGNI/Dropbox/ICM-EMSE/0_3A/Data_Science/4_Small_Data/sensitivity_analysis/TP')


n <- 1000
X1 <- runif(n, min = -pi, max = pi)
X2 <- runif(n, min = -pi, max = pi)
X3 <- runif(n, min = -pi, max = pi)

# Ishigami
f <- function(X) {
  return( sin(X[,1]) + 7*(sin(X[,2]))^2 + 0.1*(X[,3]^4)*sin(X[,1]) )
}

Y <- f(cbind(X1, X2, X3))
mu <- mean(Y)

# Splines
ss1 <- smooth.spline(X1, Y-mu)
ss2 <- smooth.spline(X2, Y-mu)
ss3 <- smooth.spline(X3, Y-mu)

op <-par(mfrow = c(1,3))
plot(X1, Y-mu)
lines(ss1, col = "blue", lwd = 3)
plot(X2, Y-mu)
lines(ss2, col = "blue", lwd = 3)
plot(X3, Y-mu)
lines(ss3, col = "blue", lwd = 3)
par(op)




library(sensitivity)
#library(DiceView)
#library(DiceOptim)

f2 <- function(X, b12, b11) {
  return( X[,1] - 2*X[,2] + b12*X[,1]*X[,2] + b11*X[,1]^2 )
}


mMooris <- morris(model = f2, b12 = 10, b11 = 1, factors = 3, r = 10,
            design = list(type = "oat", levels = 5, grid.jump = 3), binf = -0.5, bsup = 0.5)
plot(mMooris)

# With Ishigami function
mMooris2 <- morris(model = f, factors = 3, r = 10,
                  design = list(type = "oat", levels = 5, grid.jump = 3), binf = -pi, bsup = pi)
plot(mMooris2)

f3 <- function(X) {
  return( X[,1]*X[,2] )
}
# With product function
mMooris3 <- morris(model = f3, factors = 2, r = 100,
                   design = list(type = "oat", levels = 5, grid.jump = 3), binf = 0, bsup = 3)
plot(mMooris3)



# Sobol indices
fast1 <- fast99(model = f2, b12 = 10, b11 = 1, factors = 3, n = 1000, q = "qunif", q.arg = list(min = -0.5, max = 0.5))
plot(fast1)

fast2 <- fast99(model = f, factors = 3, n = 1000, q = "qunif", q.arg = list(min = -pi, max = pi))
plot(fast2)

fast3 <- fast99(model = f3, factors = 2, n = 1000, q = "qunif", q.arg = list(min = 0, max = 3))
plot(fast3)


