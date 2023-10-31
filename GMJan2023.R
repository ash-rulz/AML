library(kernlab)

# Matern32  kernel
k <- function(sigmaf = 1, ell = 1)  
{   
  rval <- function(x, y = NULL) 
  {	r = sqrt(crossprod(x-y))
  return(sigmaf^2*(1+sqrt(3)*r/ell)*exp(-sqrt(3)*r/ell))   
  }   
  class(rval) <- "kernel"   
  return(rval) 
} 

sigmaF <- 1
l <- 0.5

#Effect of ell
zGrid <- seq(0.01,1,by=0.01)
MaternFunc <- k(sigmaf = sigmaF, ell = l)
dist_0.5 <- sapply(zGrid, function(x) MaternFunc(0,x))

MaternFunc <- k(sigmaf = sigmaF, ell = 0.1)
dist_0.1 <- sapply(zGrid, function(x) MaternFunc(0,x))

MaternFunc <- k(sigmaf = sigmaF, ell = 1)
dist_1 <- sapply(zGrid, function(x) MaternFunc(0,x))

plot(zGrid, dist_0.5, type = 'l', col = 'red')
lines(zGrid, dist_0.1, type = 'l', col = 'blue')
lines(zGrid, dist_1, type = 'l', col = 'black')

#Effect of sigmaF
sigmaF <- 0.5
l <- 0.5

MaternFunc <- k(sigmaf = 0.1, ell = l)
dist_0.1 <- sapply(zGrid, function(x) MaternFunc(0,x))

MaternFunc <- k(sigmaf = 0.5, ell = l)
dist_0.5 <- sapply(zGrid, function(x) MaternFunc(0,x))

MaternFunc <- k(sigmaf = 0.9, ell = l)
dist_0.9 <- sapply(zGrid, function(x) MaternFunc(0,x))

plot(zGrid, dist_0.5, type = 'l', col = 'red', ylim = c(0,1))
lines(zGrid, dist_0.1, type = 'l', col = 'blue')
lines(zGrid, dist_0.9, type = 'l', col = 'black')

#2nd question
sigmaF <- 1
l <- 1
MaternFunc <- k(sigmaf = sigmaF, ell = l)
sigma2Noise = 0.05

GPfit <- gausspr(x = distance,
                 y = logratio,
                 scaled = TRUE,
                 type = 'regression',
                 kernel = MaternFunc,
                 var = sigma2Noise^2)

x <- distance
n <- length(x)
xs <- seq(min(x), max(x), 0.001)
Kxs <- kernelMatrix(kernel = MaternFunc, x = x, y = xs)
Kxx <- kernelMatrix(kernel = MaternFunc, x = x, y = x)
Kxs <- kernelMatrix(kernel = MaternFunc, x = x, y = xs)
Kss <- kernelMatrix(kernel = MaternFunc, x = xs, y = xs)
cov_mat <- Kss-t(Kxs)%*%
  solve(Kxx + sigma2Noise^2*diag(n), Kxs)

posMeanF <- predict(GPfit, xs)
plot(distance, logratio, pch = 19, ylim = c(-0.7, 0))
lines(xs, posMeanF, col="red", lwd = 2)

lines(xs, posMeanF - 1.96*sqrt(diag(cov_mat)),
      col = "blue", lty = 2)
lines(xs, posMeanF + 1.96*sqrt(diag(cov_mat)),
      col = "blue", lty = 2)

lines(xs, posMeanF - 1.96*sqrt((diag(cov_mat) + sigma2Noise^2)), col = "purple")
lines(xs, posMeanF + 1.96*sqrt((diag(cov_mat) + sigma2Noise^2)), col = "purple")

#Jose's implementation
load("lidar.RData") # loading the data
sigmaNoise = 0.05
x = distance
y = logratio

# Set up the kernel function
kernelFunc <- k(sigmaf = 1, ell = 5)

# Plot the data and the true 
plot(x, y, main = "", cex = 0.5)

GPfit <- gausspr(x, y, kernel = kernelFunc, var = sigmaNoise^2)
xs = seq(min(x),max(x), length.out = 100)
meanPred <- predict(GPfit, xs) # Predicting the training data. To plot the fit.
lines(xs, meanPred, col="blue", lwd = 2)

# Compute the covariance matrix Cov(f)
n <- length(x)
Kss <- kernelMatrix(kernel = kernelFunc, x = xs, y = xs)
Kxx <- kernelMatrix(kernel = kernelFunc, x = x, y = x)
Kxs <- kernelMatrix(kernel = kernelFunc, x = x, y = xs)
Covf = Kss-t(Kxs)%*%solve(Kxx + sigmaNoise^2*diag(n), Kxs)

# Probability intervals for f
lines(xs, meanPred - 1.96*sqrt(diag(Covf)), col = "red")
lines(xs, meanPred + 1.96*sqrt(diag(Covf)), col = "red")

# Prediction intervals for y 
lines(xs, meanPred - 1.96*sqrt((diag(Covf) + sigmaNoise^2)), col = "purple")
lines(xs, meanPred + 1.96*sqrt((diag(Covf) + sigmaNoise^2)), col = "purple")


legend("topright", inset = 0.02, legend = c("data","post mean","95% intervals for f", "95% predictive intervals for y"), 
       col = c("black", "blue", "red", "purple"), 
       pch = c('o',NA,NA,NA), lty = c(NA,1,1,1), lwd = 2, cex = 0.55)
