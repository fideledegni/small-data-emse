#####################################################
#  Compte rendu de TP Sensibility analysis          #
#              DEGNI Fidèle                         #
#              RODRIGUES Leticia                    #
#####################################################

rm(list=ls()) #  cleaning up
#setwd('~/Dropbox/ICM-EMSE/0_3A/Data_Science/4_Small_Data/sensitivity_analysis/TP')
setwd('C:/Users/Fidèle DEGNI/Dropbox/ICM-EMSE/0_3A/Data_Science/4_Small_Data/sensitivity_analysis/TP')
library(sensitivity)

n <- 100
X1 <- runif(n, min = 0, max = 1)
X2 <- runif(n, min = 0, max = 1)
X3 <- runif(n, min = 0, max = 1)


# Product function
f.product <- function(X) {
  return( X[,1] * X[,2] * X[,3] )
}

# Morris
morris.product <- morris(model = f.product, factors = 3, r = 1000,
                          design = list(type = "oat", levels = 5, grid.jump = 3), binf = 0, bsup = 1)
plot(morris.product, xlim = c(0, 1), ylim = c(0, 1), main = "Morris")

# Sobol indices
sobol.product <- fast99(model = f.product,  factors = 3, n = 1000, q = "qunif", q.arg = list(min = 0, max = 1))
plot(sobol.product)
title(main = "Indices de Sobol")

# Offet des variables sur la réponse
Y.product <- f.product(cbind(X1, X2, X3))
mu0.p <- mean(Y.product)

# Splines pour lissage
ss1.p <- smooth.spline(X1, Y.product-mu0.p)
ss2.p <- smooth.spline(X2, Y.product-mu0.p)
ss3.p <- smooth.spline(X3, Y.product-mu0.p)

op <-par(mfrow = c(1,3))
plot(X1, Y.product-mu0.p)
lines(ss1.p, col = "blue", lwd = 3)
plot(X2, Y.product-mu0.p)
lines(ss2.p, col = "blue", lwd = 3)
plot(X3, Y.product-mu0.p)
lines(ss3.p, col = "blue", lwd = 3)
par(op)
title(main = "Effet des variables")




####### Cas test Volcan #######
