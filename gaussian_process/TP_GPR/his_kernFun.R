linKern <- function(x, y, param){
  # input:
  #  x,y: input vectors
  #  param: parameters (sigma)
  # output:
  #  kern: covariance matrix cov(x,y)
  sigma <- param[1]
  kern <- sigma^2*outer(x, y, '*')
  return(kern)
}

cosKern <- function(x, y, param){
  # input:
  #  x,y: input vectors
  #  param: parameters (sigma,theta)
  # output:
  #  kern: covariance matrix cov(x,y)
  sigma <- param[1]
  theta <- param[2]
  dist <- outer(x/theta, y/theta, '-')
  kern <- sigma^2*cos(dist)
  return(kern)
}

expKern <- function(x, y, param){
  # input:
  #  x,y: input vectors
  #  param: parameters (sigma,theta)
  # output:
  #  kern: covariance matrix cov(x,y)
  sigma <- param[1]
  theta <- param[2]
  dist <- outer(x/theta, y/theta, '-')
  kern <- sigma^2*exp(-abs(dist))
  return(kern)
}

seKern <- function(x,y,param){
  # input:
  #  x,y: input vectors
  #  param: parameters (sigma,theta)
  # output:
  #  kern: covariance matrix cov(x,y)
  sigma <- param[1]
  theta <- param[2]
  dist<- outer(x,y, '-')
  dist2<- dist^2
  kern<-sigma^2*exp((-(dist2))/(2*theta^2))
  return(kern)
}

mat5_2Kern <- function (x,y,param){
  # input:
  #  x,y: input vectors
  #  param: parameters (sigma,theta)
  # output:
  #  kern: covariance matrix cov(x,y)
  sigma <- param[1]
  theta <- param[2]
  dist<- outer(x,y, '-')
  dist2<- dist^2
  term<-(sqrt(5)*dist)/theta
  kern<- sigma^2*(1+ term + (term^2/3))*exp(-term)
  return(kern)
}

brownKern <- function(x,y,param){
  # input:
  #  x,y: input vectors
  #  param: parameters (sigma,theta)
  # output:
  #  kern: covariance matrix cov(x,y)
  sigma <- param[1]
  mini <- matrix(0,nrow = length(x),ncol = length(y))
  for(i in 1:length(x)){
    for(j in 1:length(y)){
      mini[i,][j] <- min(x[i],y[j])
    }
  }
  kern<- sigma^2 * mini
  return(kern)
}

sincKern <- function(x,y,param){
  # input:
  #  x,y: input vectors
  #  param: parameters (sigma,theta)
  # output:
  #  kern: covariance matrix cov(x,y)
  sigma <- param[1]
  theta <- param[2]
  dist<- outer(x,y, '-')
  kern<- (sigma^2*theta*sin(dist/theta))/dist
  diag(kern) <- 1
  return(kern)
}



