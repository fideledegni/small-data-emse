############################################################
#  Compte rendu de TP Advanced numerical engineering       #
#                 DEGNI Fidèle                             #
#                 RODRIGUES Leticia                        #
############################################################

rm(list=ls()) #  cleaning up
#setwd('~/Dropbox/ICM-EMSE/0_3A/Data_Science/4_Small_Data/advanced_num_engineering/TP')
setwd('C:/Users/Fidèle DEGNI/Dropbox/ICM-EMSE/0_3A/Data_Science/4_Small_Data/advanced_num_engineering/TP')

# reminer of detailed information available in the "volcan_test_case.pdf" file
source("mainScript_DatascienceClass.R")

library(DiceOptim)
library(sensitivity)


# Hypercubes latins
generateLHS <- function(npts, dim) {
  X <- matrix(nrow = npts, ncol = dim)
  for (i in 1:dim) {
    X[, i] <- sample.int(npts) - 1
  }
  #Normalisation pour avoir les coordonnees entre 0 et 1
  X <- X/(npts-1)
  return(X)
}

# Plan d'expériences avec hypercubes latins
npts <- 50
dim <- 5
design.fact <- generateLHS(npts, dim)

response = compute_wls(design.fact)
plot(response, ylab = "Observation", main = "Obsevation des valeurs calculées par compute_wls")

# Modèle de krigeage
fitted.model1 <- km(~1, design=design.fact, response=response,
                    covtype="matern5_2", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))
plot(fitted.model1)



# Prediction function
kriging.mean <- function(Xnew, m) {
  predict(m, Xnew, "UK", se.compute = FALSE, checkNames = FALSE)$mean
}
f.volcano <- function(X) {
  kriging.mean(X, fitted.model1)
}

# Morris
morris.volcano <- morris(model = f.volcano, factors = 5, r = 1000,
                         design = list(type = "oat", levels = 5, grid.jump = 3), binf = 0, bsup = 1)
plot(morris.volcano, xlim = c(0, 3), ylim = c(0, 1), main = "Morris")

# Sobol indices
sobol.volcano <- fast99(model = f.volcano,  factors = 5, n = 1000, q = "qunif", q.arg = list(min = 0, max = 1))
plot(sobol.volcano)
title(main = "Indices de Sobol")


# Effet des variables sur la réponse

# Plan d'expériences de test pour faire des tracés
n <- 10000
Xs <- runif(n = n)
Ys <- runif(n = n)
Zs <- runif(n = n)
AP <- expand.grid(seq(0, 1, length=100), seq(0, 1, length=100))

design.fact2 <- cbind(Xs, Ys, Zs, AP)
names(design.fact2) <- c("Xs", "Zs", "Ys", "A", "P")


Y.volcano <- f.volcano(design.fact2)
mu0.v <- mean(Y.volcano)

# Splines pour lissage
ss1.v <- smooth.spline(Xs, Y.volcano-mu0.v)
ss2.v <- smooth.spline(Ys, Y.volcano-mu0.v)
ss3.v <- smooth.spline(Zs, Y.volcano-mu0.v)
ss4.v <- smooth.spline(AP[,1], Y.volcano-mu0.v)
ss5.v <- smooth.spline(AP[,2], Y.volcano-mu0.v)

op <-par(mfrow = c(1,5))
plot(Xs, Y.volcano-mu0.v)
lines(ss1.v, col = "blue", lwd = 3)
plot(Ys, Y.volcano-mu0.v)
lines(ss2.v, col = "blue", lwd = 3)
plot(Zs, Y.volcano-mu0.v)
lines(ss3.v, col = "blue", lwd = 3)
plot(AP[,1], Y.volcano-mu0.v, xlab = "a (rayon)")
lines(ss4.v, col = "blue", lwd = 3)
plot(AP[,2], Y.volcano-mu0.v, xlab = "p (pression)")
lines(ss5.v, col = "blue", lwd = 3)
par(op)
title(main = "Effet des variables")

print(fitted.model1)

# Les interactions
inter <- predict(fitted.model1, design.fact2, "UK", se.compute = FALSE, checkNames = FALSE)$mean -
  predict(ss1.v, Xs)$y -
  predict(ss2.v, Ys)$y -
  predict(ss3.v, Zs)$y -
  predict(ss4.v, AP[,1])$y -
  predict(ss5.v, AP[,2])$y


image(seq(0, 1, length=100), seq(0, 1, length=100), matrix(inter, ncol = 100), 
      xlab = "rayon", ylab = "pression", main = "Interaction")


contour(seq(0, 1, length=100), seq(0, 1, length=10), matrix(inter, ncol = 100), 
        xlab = "rayon", ylab = "pression", main = "Interaction")
