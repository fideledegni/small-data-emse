discrepancy <- function (design)
{
  X <- as.matrix(design)
  dimension <- dim(X)[2]
  n <- dim(X)[1]
  
  dL2 <- 0
  for (j in 1:n) {
    for (i in 1:n) {
      if (i != j) {
        t <- c()
        for (l in 1:dimension) t <- c(t, 1 - max(X[i, 
                                                   l], X[j, l]))
        t <- (prod(t))/(n^2)
      }
      else {
        t1 <- 1 - X[i, ]
        t1 <- prod(t1)
        t2 <- 1 - X[i, ]^2
        t2 <- prod(t2)
        t <- t1/(n^2) - ((2^(1 - dimension))/n) * t2
      }
      dL2 <- dL2 + t
    }
  }
  DisL2star <- sqrt(3^(-dimension) + dL2)
  return(DisL2star)
}