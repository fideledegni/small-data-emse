#####################################################
#  Compte rendu de TP Plan d'experiences            #
#              DEGNI Fidèle                         #
#              RODRIGUES Leticia                    #
#####################################################

rm(list=ls()) # to clear the environment
#setwd('~/Dropbox/ICM-EMSE/0_3A/Data_Science/4_Small_Data/experiments_design/TP1')
setwd('C:/Users/Fidèle DEGNI/Dropbox/ICM-EMSE/0_3A/Data_Science/4_Small_Data/experiments_design/TP1')

#1. Construction de plans
#1.1 Tessellations centroidales de Voronoi
# Algorithme de MacQueen
generateCTV <- function(npts, dim, nite) {
  X <- matrix(data = runif(npts*dim), nrow = npts, ncol = dim)
  j <- rep(1, npts)
  for (i in 1:nite) {
    w <- runif(dim)
    norme <- rep(0, npts)
    for (k in 1:npts) {
      norme[k] <- norm(t(X[k,]-w))
    }
    min_ind <- which.min(norme)
    X[min_ind,] <- (j[min_ind]*X[min_ind,] + w)/(j[min_ind]+1)
    j[min_ind] <- j[min_ind] +1
  }
  return(X)
}
# Test
npts <- 20
dim <- 2
nite <- 100
X <- generateCTV(npts, dim, nite)
(X)
plot(X, xlab = "x1", ylab = "x2", col = "blue", pch = 20, lwd = 5, 
     main = "Test : generateCTV avec :\n npts = 20, dim = 2 et nite = 100")


#1.2 Hypercubes latins
generateLHS <- function(npts, dim) {
  X <- matrix(nrow = npts, ncol = dim)
  for (i in 1:dim) {
    X[, i] <- sample.int(npts) - 1
  }
  #Normalisation pour avoir les coordonnees entre 0 et 1
  X <- X/(npts-1)
  return(X)
}
npts <- 30
dim <- 2
X <- generateLHS(npts, dim)
(X)
plot(X, xlab = "x1", ylab = "x2", col = "blue", pch = 20, lwd = 5, 
     main = "Test : generateLHS avec :\n npts = 30, dim = 2")


#1.3 Criteres
evalMinDist <- function(X) {
  D <- diag(NA, nrow(X))
  for (i in 1:(nrow(X)-1)) {
    for (j in (i+1):nrow(X)) {
      norme = norm(t(X[i,]-X[j,]))
      D[i, j] <- norme
      D[j, i] <- norme
    }
  }
  d <- min(D, na.rm = TRUE)
  diag(D) <- 0
  return( list(minDist = d, allDist = D) )
}
evalMinDist(X)



#1.4 Hypercubes latins
# Recherche aleatoire
generateLHSOA <- function(npts, dim, nbtest) {
  LHSO <- generateLHS(npts, dim)
  critere <- evalMinDist(LHSO)$minDist
  for (i in 1:nbtest) {
    LHSO_temp <- generateLHS(npts, dim)
    critere_temp <- evalMinDist(LHSO_temp)$minDist
    if (critere_temp < critere) {
      LHSO <- LHSO_temp
      critere <- critere_temp
    }
  }
  return(LHSO)
}
nbtest <- 30
npts <- 30
dim <- 5
X <- generateLHSOA(npts, dim, nbtest)
evalMinDist(X)$minDist



# Algorithme d'echange avec recuit simule
# Rq : On accepte la modification seulement lorsque l'echange est benefique !
generateLHSOE <- function(npts, dim, nbtest) {
  LHSO <- generateLHS(npts, dim)
  critere <- evalMinDist(LHSO)$minDist
  for (i in 1:nbtest) {
    LHSO_temp <- LHSO
    ind <- which(LHSO_temp==max(LHSO_temp), arr.ind = TRUE) # les points critiques
    pc_l <- ind[1,1] # numero de ligne d'un point critique
    pc_c <- sample.int(dim, 1) # choix aleatoire d'une colonne du point critique
    other_l <- sample.int(npts, 1) # choix aleatoire d'une autre cellule
    s <- LHSO_temp[pc_l, pc_c] # savegarde pour la permutation
    LHSO_temp[pc_l, pc_c] <- LHSO_temp[other_l, pc_c]
    LHSO_temp[other_l, pc_c] <- s
    critere_temp <- evalMinDist(LHSO_temp)$minDist
    # Si l'echange est benefique on accepte la modification
    if (critere_temp < critere) {
      LHSO <- LHSO_temp
      critere <- critere_temp
    }
  }
  return(LHSO)
}
nbtest <- 30
npts <- 30
dim <- 5
X <- generateLHSOE(npts, dim, nbtest)
evalMinDist(X)$minDist




#2. Analyse
nite <- 100
CVT_2_10 <- generateCTV(10, 2, nite)
CVT_5_70 <- generateCTV(70, 5, nite)
CVT_10_150 <- generateCTV(150, 10, nite)
LHS_2_10 <- generateLHS(10, 2)
LHS_5_70 <- generateLHS(70, 5)
LHS_10_150 <- generateLHS(150, 10)
nbtest <- 30
LHSOA_2_10 <- generateLHSOA(10, 2, nbtest)
LHSOA_5_70 <- generateLHSOA(70, 5, nbtest)
LHSOA_10_150 <- generateLHSOA(150, 10, nbtest)
LHSOE_2_10 <- generateLHSOE(10, 2, nbtest)
LHSOE_5_70 <- generateLHSOE(70, 5, nbtest)
LHSOE_10_150 <- generateLHSOE(150, 10, nbtest)

source("discrepancy.R")

# dim = 2 et npts = 10
# Repartition sur les marginales de dim 1
op <-par(mfrow = c(4,1))
hist(CVT_2_10[,1])
hist(LHS_2_10[,1])
hist(LHSOA_2_10[,1])
hist(LHSOE_2_10[,1])
par(op)

op <-par(mfrow = c(4,1))
hist(CVT_2_10[,2])
hist(LHS_2_10[,2])
hist(LHSOA_2_10[,2])
hist(LHSOE_2_10[,2])
par(op)

# Repartition sur les marginales de dim 2
pairs(CVT_2_10, main = "CVT_2_10")
pairs(LHS_2_10, main = "LHS_2_10")
pairs(LHSOA_2_10, main = "LHSOA_2_10")
pairs(LHSOE_2_10, main = "LHSOE_2_10")

# Valeurs du critere maximin
evalMinDist(CVT_2_10)$minDist
evalMinDist(LHS_2_10)$minDist
evalMinDist(LHSOA_2_10)$minDist
evalMinDist(LHSOE_2_10)$minDist

# Valeurs de la discrepance
discrepancy(CVT_2_10)
discrepancy(LHS_2_10)
discrepancy(LHSOA_2_10)
discrepancy(LHSOE_2_10)


# dim = 5 et npts = 70
# Repartition sur les marginales de dim 1
op <-par(mfrow = c(4,1))
hist(CVT_5_70[,1])
hist(LHS_5_70[,1])
hist(LHSOA_5_70[,1])
hist(LHSOE_5_70[,1])
par(op)

op <-par(mfrow = c(4,1))
hist(CVT_5_70[,2])
hist(LHS_5_70[,2])
hist(LHSOA_5_70[,2])
hist(LHSOE_5_70[,2])
par(op)

op <-par(mfrow = c(4,1))
hist(CVT_5_70[,3])
hist(LHS_5_70[,3])
hist(LHSOA_5_70[,3])
hist(LHSOE_5_70[,3])
par(op)

op <-par(mfrow = c(4,1))
hist(CVT_5_70[,4])
hist(LHS_5_70[,4])
hist(LHSOA_5_70[,4])
hist(LHSOE_5_70[,4])
par(op)

# Repartition sur les marginales de dim 2
pairs(CVT_5_70, main = "CVT_5_70")
pairs(LHS_5_70, main = "LHS_5_70")
pairs(LHSOA_5_70, main = "LHSOA_5_70")
pairs(LHSOE_5_70, main = "LHSOE_5_70")

# Valeurs du critere maximin
evalMinDist(CVT_5_70)$minDist
evalMinDist(LHS_5_70)$minDist
evalMinDist(LHSOA_5_70)$minDist
evalMinDist(LHSOE_5_70)$minDist

# Valeurs de la discrepance
discrepancy(CVT_5_70)
discrepancy(LHS_5_70)
discrepancy(LHSOA_5_70)
discrepancy(LHSOE_5_70)



# dim = 10 et npts = 150
# Repartition sur les marginales de dim 1
op <-par(mfrow = c(4,1))
hist(CVT_10_150[,1])
hist(LHS_10_150[,1])
hist(LHSOA_10_150[,1])
hist(LHSOE_10_150[,1])
par(op)

op <-par(mfrow = c(4,1))
hist(CVT_10_150[,2])
hist(LHS_10_150[,2])
hist(LHSOA_10_150[,2])
hist(LHSOE_10_150[,2])
par(op)

# Repartition sur les marginales de dim 2
pairs(CVT_10_150, main = "CVT_10_150")
pairs(LHS_10_150, main = "LHS_10_150")
pairs(LHSOA_10_150, main = "LHSOA_10_150")
pairs(LHSOE_10_150, main = "LHSOE_10_150")

# Valeurs du critere maximin
evalMinDist(CVT_10_150)$minDist
evalMinDist(LHS_10_150)$minDist
evalMinDist(LHSOA_10_150)$minDist
evalMinDist(LHSOE_10_150)$minDist

# Valeurs de la discrepance
discrepancy(CVT_10_150)
discrepancy(LHS_10_150)
discrepancy(LHSOA_10_150)
discrepancy(LHSOE_10_150)



# Histogramme des interdistances pour dim = 2 et npts = 10
op <-par(mfrow = c(2,1))
hist(evalMinDist(CVT_2_10)$allDist)
hist(evalMinDist(LHS_2_10)$allDist)
op <-par(mfrow = c(2,1))
hist(evalMinDist(LHSOA_2_10)$allDist)
hist(evalMinDist(LHSOE_2_10)$allDist)
par(op)

# Histogramme des interdistances pour dim = 5 et npts = 70
op <-par(mfrow = c(2,1))
hist(evalMinDist(CVT_5_70)$allDist)
hist(evalMinDist(LHS_5_70)$allDist)
op <-par(mfrow = c(2,1))
hist(evalMinDist(LHSOA_5_70)$allDist)
hist(evalMinDist(LHSOE_5_70)$allDist)
par(op)

# Histogramme des interdistances pour dim = 10 et npts = 150
op <-par(mfrow = c(2,1))
hist(evalMinDist(CVT_10_150)$allDist)
hist(evalMinDist(LHS_10_150)$allDist)
op <-par(mfrow = c(2,1))
hist(evalMinDist(LHSOA_10_150)$allDist)
hist(evalMinDist(LHSOE_10_150)$allDist)
par(op)


