####GP #####
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

SquaredExpKernel <- function(x1,x2,sigmaF=1,l=0.3){
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- sigmaF^2*exp(-0.5*( (x1-x2[i])/l )^2 )#Hyperparameter, S4
  }
  return(K)
}

xTrain <- c(-1.0, -0.6, -0.2, 0.4, 0.8)
y <- c(0.768, -0.044, -0.940, 0.719, -0.664)

sigmaF <- 1#Prior hyperparameters
l <- 1#Prior hyperparameters
sigmaNoise <- 1#Assumed

XStar <- seq(-1, 1, length = 100)#New points for prediction
post_GP <- posteriorGP(xTrain, y, XStar, sigmaNoise, SquaredExpKernel, sigmaF, l)
plot(XStar,post_GP$pos_mean, type = 'l', ylab = 'Posterior mean',ylim=c(-2.5,2.5))
lines(XStar, post_GP$pos_mean - 1.96*sqrt(diag(post_GP$pos_var)),
      col = "blue", lty = 2)
lines(XStar, post_GP$pos_mean + 1.96*sqrt(diag(post_GP$pos_var)),
      col = "blue", lty = 2)
lines(XStar, post_GP$pos_mean - 1.96*sqrt(diag(post_GP$pos_var)+sigmaNoise^2), 
      col = "red", lty=2)
lines(XStar, post_GP$pos_mean + 1.96*sqrt(diag(post_GP$pos_var)+sigmaNoise^2), 
      col = "red", lty=2)
points(xTrain,y,col='black', pch = 19)
legend("topright", legend = c("Posterior Mean", "95% Confidence Bands", "Observation Points"),
       col = c("black", "blue", "black"), lty = c(1, 2, NA), pch = c(NA, NA, 19),
       cex = 0.6,
       bty = 'n')


plot(xData,yData,xlim=c(0,10),ylim=c(-15,15))
lines(xGrid, res$mu, col="blue", lwd = 2)
lines(xGrid, res$mu - 1.96*sqrt(diag(res$var)), col = "red")
lines(xGrid, res$mu + 1.96*sqrt(diag(res$var)), col = "red")
lines(xGrid, res$mu - 1.96*sqrt(diag(res$var)+sigmaN^2), col = "red", lwd=2)
lines(xGrid, res$mu + 1.96*sqrt(diag(res$var)+sigmaN^2), col = "red", lwd=2)