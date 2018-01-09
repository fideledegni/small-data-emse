###########################
# MAIN CODE for
# Inversion of a punctual displacements source from 3D data
# The data used are under-sampled with the quadtree method (irregular grid)
#
# Rodolphe Le Riche, Nicolas Durrande, Valerie Cayol, Victor Picheny
#
#
# Notes: 
# * use of global variables, starting with Glb_...
###########################

####### load utilities ##########################

rm(list=ls()) #  cleaning up
#setwd('~/Dropbox/ICM-EMSE/0_3A/Data_Science/4_Small_Data/test_case_volcano/VolcanoTestCase')
setwd('C:/Users/Fidèle DEGNI/Dropbox/ICM-EMSE/0_3A/Data_Science/4_Small_Data/test_case_volcano/VolcanoTestCase')


library(R.matlab)
source("./mogi_3D.R")
source("./wls_ulos.R")
source('kernels.R')

###### input for variables identification ###########
n2i <- list(xs=1, ys=2, zs=3, a=4, p=5) # name to index for variables
varnames <- c("xs","ys","zs","a","p")
nbvar <- 5

# optimum
xstar <- NA
xstar[n2i$xs] <- 367000 # X location of source in m UTM coordinates
xstar[n2i$ys] <- 7650300 # Y location of source in m UTM
xstar[n2i$zs] <- 0 # Elevation of source with respect to sea level in m
xstar[n2i$a] <- 500 # source radius in m
xstar[n2i$p] <- 20 # Source overpressure in MPa

# order of magnitude of the variables 
# (different from 0, to be used when scaling, setting bounds...)
xmag <- NA
xmag[n2i$xs] <- 367000 # X location of source in UTM coordinates
xmag[n2i$ys] <- 7650300 # Y location of source in UTM
xmag[n2i$zs] <- 1000 # Elevation of source with respect to sea level
xmag[n2i$a] <- 500 # source radius
xmag[n2i$p] <- 100 # Source overpressure in MPa

# bounds on variables
xmax <- NA
xmin <- NA
xmin[n2i$xs]<-364000
xmax[n2i$xs]<-368000
xmin[n2i$ys]<-7649000
xmax[n2i$ys]<-7651000
xmin[n2i$zs]<- -3000
xmax[n2i$zs]<-1000
xmin[n2i$a]<- 50
xmax[n2i$a]<- 1000
xmin[n2i$p]<- -500
xmax[n2i$p]<-500

Glb_var <<- list(n2i=n2i,nbvar=nbvar,xmag=xmag,xmax=xmax,xmin=xmin) # always useful stuff

####### load data ###########

data <- readMat('data_nonoise.mat') # TODO: do a read from cvs version (rodo)
Glb_xi <<- as.matrix(data$locdata[,1])
Glb_yi <<- as.matrix(data$locdata[,2])
Glb_zi <<- as.matrix(data$locdata[,3])
Glb_ulos <<- as.matrix(data$locdata[,4])

# calculate data Covariance matrix, store it in a Global variable
# covariance from exponential kernel, var = 5e-4m2, cor_length = 850 m
# and invert it
Xdata <- data$locdata[,1:2] # z's are not accounted for in Xdata

Glb_CXinv <<- solve(kExp(Xdata,Xdata,c(5e-4,850,850))) # calculated once for all, used in wls_ulos
# # has been compared to the matlab covariance matrix below, max difference 1e-6
rm(data)
rm(Xdata)

#######  useful functions #############################
# Scale from [0 1] to [xmin xmax]
unnorm_var <- function(Xnorm){
  if (is.null(dim(Xnorm))) Xnorm <- matrix(data = Xnorm, nrow=1) # numeric vector
  nbrep <- nrow(Xnorm)
  Xu <- matrix(rep(xmin,times=nbrep),byrow = T,ncol=nbvar) + 
    Xnorm * matrix(rep((xmax-xmin),times=nbrep),byrow = T,ncol=nbvar)
  colnames(Xu) <- varnames
  return(Xu)
}
# Scale from [xmin xmax] to [0 1] 
norm_var <- function(X){
  if (is.null(dim(X))) X <- matrix(data = X, nrow=1)
  nbrep <- nrow(X)
  Xu <- (X - matrix(rep(xmin,times=nbrep),byrow = T,ncol=nbvar)) / 
    matrix(rep((xmax-xmin),times=nbrep),byrow = T,ncol=nbvar)  
  colnames(Xu) <- varnames
  return(Xu)
}

# normalize the output so that it is centered with a unit std dev
# because wls ranges from 0 to 10^9, do a scaling in log(1+wls)
# wls normalization
normalizeWLS <- function(awls){
  lawls <- log(1+awls)
  return((lawls - 8.54)/3.2)
}
# objective function
compute_wls <- function(x) {
  if (is.null(dim(x))) {
    x <- matrix(data = x, nrow=1)
  }  else if (ncol(x)!=nbvar) {
    x <- t(x)
  }
  xx <- unnorm_var(Xnorm=x)
  y <- apply(xx, 1, wls_ulos)
  return(normalizeWLS(y))
}





#### Our code ######
library(DiceView)
library(DiceOptim)
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

#3. Krieging models
npts <- 70
dim <- 5
design.fact <- generateLHS(npts, dim)
response = compute_wls(design.fact)
plot(response, ylab = "Observation", main = "Obsevation des valeurs calculées par compute_wls")

# Unnorm and visualise reponses with respect to the different variables
unnormed_fact <- unnorm_var(design.fact)
plot(unnormed_fact[,1], response, xlab = "xs", ylab = "Observation", 
     main = "Observation en fonction de la longitude")
plot(unnormed_fact[,2], response, xlab = "ys", ylab = "Observation", 
     main = "Observation en fonction de la latitude")
plot(unnormed_fact[,3], response, xlab = "zs", ylab = "Observation", 
     main = "Observation en fonction de l'élevation")
plot(unnormed_fact[,4], response, xlab = "a", ylab = "Observation", 
     main = "Observation en fonction du rayon")
plot(unnormed_fact[,5], response, xlab = "p", ylab = "Observation", 
     main = "Observation en fonction de la pression")


### Constructiondes modèles ###

# tendance constante
fitted.model1 <- km(~1, design=design.fact, response=response,
                    covtype="matern5_2", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))
plot(fitted.model1)

a <- design.fact[,4] # rayon
p <- design.fact[,5] # pression

# tendance lineaire avec le rayon
fitted.model2 <- km(~1+a, design=design.fact, response=response,
                    covtype="matern5_2", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))
plot(fitted.model2)

# ajout de la tendance quadratique avec la pression
fitted.model3 <- km(~1+a+p^2, design=design.fact, response=response,
                    covtype="matern5_2", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))
plot(fitted.model3)

# ajout de la tendance quadratique avec la pression et suppression de la constante
fitted.model4 <- km(~a+p^2, design=design.fact, response=response,
                    covtype="matern5_2", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))
plot(fitted.model4)

# ajout de la tendance quadratique avec la pression, noyau gauss
fitted.model5 <- km(~1+a+p^2, design=design.fact, response=response,
                    covtype="gauss", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))
plot(fitted.model5)

# ajout de la tendance quadratique avec la pression, noyau exp
fitted.model6 <- km(~1+a+p^2, design=design.fact, response=response,
                    covtype="exp", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))
plot(fitted.model6)

# fitted.model2 est le modèle retenu pour la suite



#4. Identification of the volcano source
#4.1 EGO
### EGO, 5 steps ##################
library(rgenoud)
nsteps <- 25
lower <- rep(0,dim)
upper <- rep(1,dim)
oEGO <- EGO.nsteps(model=fitted.model1, fun=compute_wls, nsteps=nsteps,
                   lower=lower, upper=upper, control=list(pop.size=20, BFGSburnin=2))
print(oEGO$par)
print(oEGO$value)

# Unnorm variables
unnormed_var <- unnorm_var(oEGO$par)

plot(unnormed_var[,1], oEGO$value, xlab = "xs", ylab = "Observation", 
     main = "Nouveaux points en fonction de la longitude")
plot(unnormed_var[,2], oEGO$value, xlab = "ys", ylab = "Observation", 
     main = "Nouveaux points en fonction de la latitude")
plot(unnormed_var[,3], oEGO$value, xlab = "zs", ylab = "Observation", 
     main = "Nouveaux points en fonction de l'élevation")
plot(unnormed_var[,4], oEGO$value, xlab = "a", ylab = "Observation", 
     main = "Nouveaux points en fonction du rayon")
plot(unnormed_var[,5], oEGO$value, xlab = "p", ylab = "Observation", 
     main = "Nouveaux points en fonction de la pression")



plot(response, ylab = "Observation",
     main = "Convergence de la fonction objectif \n Nouveaux points en bleu")
points(oEGO$value, col = "blue")
pairs(oEGO$lastmodel@X, main = "Distribution des points dans le plan")
hist(oEGO$value)

min(oEGO$value)
ind <- which.min(oEGO$value)
argMin <- unnormed_var[ind,]
argMin

