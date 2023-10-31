####Question 4 - GP####

#Part 1
X<-seq(0,10,.1)
Yfun<-function(x){
  return (x*(sin(x)+sin(3*x))+rnorm(length(x),0,2))
}
plot(X,Yfun(X),xlim=c(0,10),ylim=c(-15,15), type = 'l')

SquaredExpKernel <- function(x1,x2,sigmaF=1,l=0.3){
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- sigmaF^2*exp(-0.5*( (x1-x2[i])/l )^2 )#Hyperparameter, S4
  }
  return(K)
}

#Gets the posterior mean and variance
posteriorGP <- function(X, y, XStar, sigmaNoise, k, ...){
  #sigmaNoise is the sd
  K <- k(X, X, ...)
  n <- length(X)
  #Hyperparameter, S11
  L <- t(chol(K+((sigmaNoise**2) * diag(n))))#Upper to Lower Triangular matrix
  #Predicting mean
  alpha <- solve(t(L) ,solve(L,y))
  KStar <- k(X, XStar, ...)
  FStar <- t(KStar) %*% alpha
  #Predicting variance
  v <- solve(L, KStar)
  V_FStart <- k(XStar, XStar, ...) - (t(v) %*% v)
  return(list(pos_mean = FStar,
              pos_var = V_FStart))
}

sigmaNoise <- 2; sigmaF <- 1.5; l <- .6
post_GP <- posteriorGP(X, Yfun(X), X, sigmaNoise, SquaredExpKernel, sigmaF, l)

plot(X,post_GP$pos_mean, type = 'l', ylab = 'Posterior mean',ylim=c(-15,15))
lines(X, post_GP$pos_mean - 1.96*sqrt(diag(post_GP$pos_var)),
      col = "blue", lty = 2)
lines(X, post_GP$pos_mean + 1.96*sqrt(diag(post_GP$pos_var)),
      col = "blue", lty = 2)
lines(X,Yfun(X),col='red')
legend("topright", legend = c("Posterior Mean", "95% Confidence Bands", "Observation Points"),
       col = c("black", "red", "yellow"), lty = c(1, 2, 1),
       cex = 0.6,
       bty = 'n')



#Part 2
X<-seq(0,10,2)
Yfun<-function(x){
  return (x*(sin(x)+sin(3*x))+rnorm(length(x),0,.2))
}
plot(X,Yfun(X),xlim=c(0,10),ylim=c(-15,15))

sigmaNoise <- 2; sigmaF <- 1.5; l <- .6
XStar <- seq(0,10,.1)
post_GP <- posteriorGP(X, Yfun(X), XStar, sigmaNoise, SquaredExpKernel, sigmaF, l)

plot(XStar,post_GP$pos_mean, type = 'l', ylab = 'Posterior mean',ylim=c(-15,15))
lines(XStar, post_GP$pos_mean - 1.96*sqrt(diag(post_GP$pos_var)),
      col = "blue", lty = 2)
lines(XStar, post_GP$pos_mean + 1.96*sqrt(diag(post_GP$pos_var)),
      col = "blue", lty = 2)
points(XStar,Yfun(XStar),col='red')
legend("topright", legend = c("Posterior Mean", "95% Confidence Bands", "Observation Points"),
       col = c("black", "red", "red"), lty = c(1, 2, NA),
       cex = 0.6,
       bty = 'n')

X
Yfun(X)
post_GP$pos_mean
upper_band <- post_GP$pos_mean + 1.96*sqrt(diag(post_GP$pos_var))
lower_band <- post_GP$pos_mean - 1.96*sqrt(diag(post_GP$pos_var))
plot(XStar, upper_band-lower_band)
X <- c(X, XStar[which.max(upper_band-lower_band)])


####Jose's code####
X<-seq(0,10,.1)
Yfun<-function(x){
  return (x*(sin(x)+sin(3*x))+rnorm(length(x),0,2))
}
plot(X,Yfun(X),xlim=c(0,10),ylim=c(-15,15))

SEKernel <- function(x1,x2){
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- (sigmaF^2)*exp(-0.5*( (x1-x2[i])/l)^2 )
  }
  return(K)
}

posteriorGP <- function(X,y,k,sigmaNoise,xStar){
  n <- length(y)
  L <- t(chol(k(X,X)+((sigmaNoise^2)*diag(n))))
  a <- solve(t(L),solve(L,y))
  kStar <- k(X,xStar)
  mu <- t(kStar)%*%a
  v <- solve(L,kStar)
  var <- k(xStar,xStar)-(t(v)%*%v)
  logmar <- -0.5*(t(y)%*%a)-sum(log(diag(L)))-(n/2)*log(2*pi)
  return(list("mu"=mu,"var"=var,"logmar"=logmar))
}

sigmaF <- 1.5
l <- .5 # These hyperparameter values work well as the y probability bands cover around 95% of the y values in the dataset.
xData <- X
yData <- Yfun(X)
sigmaN <- 2
xGrid <- seq(0,10,.1)
res<-posteriorGP(X=xData,y=yData,k=SEKernel,sigmaNoise=sigmaN,xStar=xGrid)
plot(xData,yData,xlim=c(0,10),ylim=c(-15,15))
lines(xGrid, res$mu, col="blue", lwd = 2)
lines(xGrid, res$mu - 1.96*sqrt(diag(res$var)), col = "red")
lines(xGrid, res$mu + 1.96*sqrt(diag(res$var)), col = "red")
lines(xGrid, res$mu - 1.96*sqrt(diag(res$var)+sigmaN^2), col = "red", lwd=2)
lines(xGrid, res$mu + 1.96*sqrt(diag(res$var)+sigmaN^2), col = "red", lwd=2)


#Part 2
X<-seq(0,10,2)
Yfun<-function(x){
  return (x*(sin(x)+sin(3*x))+rnorm(length(x),0,.2))
}
plot(X,Yfun(X),xlim=c(0,10),ylim=c(-15,15))

for(i in 1:4){
  sigmaF <- 1.5
  l <- .5 # These hyperparameter values work well as the y probability bands cover around 95% of the y values in the dataset.
  xData <- X
  yData <- Yfun(X)
  sigmaN <- .2
  xGrid <- seq(0,10,.1)
  res<-posteriorGP(X=xData,y=yData,k=SEKernel,sigmaNoise=sigmaN,xStar=xGrid)
  plot(xData,yData,xlim=c(0,10),ylim=c(-15,15))
  lines(xGrid, res$mu, col="blue", lwd = 2)
  lines(xGrid, res$mu - 1.96*sqrt(diag(res$var)), col = "red")
  lines(xGrid, res$mu + 1.96*sqrt(diag(res$var)), col = "red")
  lines(xGrid, res$mu - 1.96*sqrt(diag(res$var)+sigmaN^2), col = "red", lwd=2)
  lines(xGrid, res$mu + 1.96*sqrt(diag(res$var)+sigmaN^2), col = "red", lwd=2)
  
  foo<-.1*which.max(res$mu + 1.96*sqrt(diag(res$var)+sigmaN^2) - (res$mu - 1.96*sqrt(diag(res$var)+sigmaN^2)))
  X<-c(X,foo)
}
