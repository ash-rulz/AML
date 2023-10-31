####GP####

#Part 1
SquaredExpKernel <- function(x1,x2,sigmaF=1,l=0.3){
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- sigmaF^2*exp(-0.5*( (x1-x2[i])/l )^2 )#Hyperparameter, S4
  }
  return(K)
}
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
  logmar <- -0.5*(t(y)%*%alpha)-sum(log(diag(L)))-(n/2)*log(2*pi)
  return(list(pos_mean = FStar,
              pos_var = V_FStart,
              log_mar = logmar))
}

xTrain <- c(-1.0, -0.6, -0.2, 0.4, 0.8)
y <- c(0.768, -0.044, -0.940, 0.719, -0.664)
sigmaNoise <- 0; sigmaF <- 1; l <- 0.3
XStar <- seq(-1, 1, 0.01)#New points for prediction
post_GP <- posteriorGP(xTrain, y, XStar, sigmaNoise, SquaredExpKernel, sigmaF, l)
plot(XStar,post_GP$pos_mean, type = 'l', ylab = 'Posterior mean',ylim=c(-2.5,2.5))
lines(XStar, post_GP$pos_mean - 1.96*sqrt(diag(post_GP$pos_var)),
      col = "blue", lty = 2)
lines(XStar, post_GP$pos_mean + 1.96*sqrt(diag(post_GP$pos_var)),
      col = "blue", lty = 2)
points(xTrain,y,col='black', pch = 19)
legend("topright", legend = c("Posterior Mean", "95% Confidence Bands", "Observation Points"),
       col = c("black", "blue", "black"), lty = c(1, 2, NA), pch = c(NA, NA, 19),
       cex = 0.6,
       bty = 'n')
which(XStar == 0.00)
plot(XStar, post_GP$pos_var[which(XStar == 0.00),])
abline(v = xTrain)
abline(h = 0)


#Part 2

data <- read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/TempTullinge.csv",
                 header=TRUE, sep=";")
time <- 1:nrow(data) #1... 2190
temp <- data[,2]
day <- (1:nrow(data))%%365 #1...365 1...365
day[day == 0] <- 365
#Consider only time = 1, 6, 11, . . ., 2186
indx <- seq(1, length(time), 5)
time_sub <- time[indx]
day_sub <- day[indx]
temp_sub <- data[indx,2]

sigmaf <- 20;ell <- 0.2
SquaredExpKernel <- function(x1,x2,sigmaF=20,l=0.2){
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- sigmaF^2*exp(-0.5*( (x1-x2[i])/l )^2 )#Hyperparameter, S4
  }
  return(K)
}

posteriorGP <- function(X, y, k, par){
  #sigmaNoise is the sd
  K <- k(X, X)
  n <- length(X)
  #Hyperparameter, S11
  L <- t(chol(K+((par**2) * diag(n))))#Upper to Lower Triangular matrix
  #Predicting mean
  alpha <- solve(t(L) ,solve(L,y))
  logmar <- -0.5*(t(y)%*%alpha)-sum(log(diag(L)))-(n/2)*log(2*pi)
  return(list(log_mar = logmar))
}

#Scale the x and y -> (val-mean)/sd
scaled_time_sub <- scale(time_sub)
#scaled_time <- scale(time)
scaled_temp_sub <- scale(temp_sub)

sigma2NoiseScaled <- seq(0.1, 10, 0.01)
sigmaF <- 20; l <- 0.2
log_marginal <- c()
for (i in sigma2NoiseScaled) {
  post_GP <- posteriorGP(scaled_time_sub,
                         scaled_temp_sub,
                         scaled_time_sub, #Pass the entire time for prediction
                         i, #Variance
                         SquaredExpKernel, sigmaF, l)  
  log_marginal <- c(log_marginal, as.vector(post_GP$log_mar))
}

plot(sigma2NoiseScaled, log_marginal)

sigma2NoiseScaled[which.max(log_marginal)]#0.19

foo<-optim(par = 0.1, fn = posteriorGP, 
           X=scaled_time_sub,
           y=scaled_temp_sub,
           k=SquaredExpKernel, 
           method="L-BFGS-B",
           lower = 0.01,
           upper = 10,
           control=list(fnscale=-1))
foo$value
foo$par

#Jose's code
# Question 3: GPs.

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

sigmaF <- 1
l <- 0.3
xData <- c(-1,-0.6,-0.2,0.4,0.8)
yData <- c(0.768,-0.044,-0.94,0.719,-0.664)
xGrid <- seq(-1,1,0.01)
res<-posteriorGP(X=xData,y=yData,k=SEKernel,sigmaNoise=0,xStar=xGrid)
plot(xData,yData,xlim=c(-1,1),ylim=c(-0.5,0.5))
xGrid[101]
lines(xGrid, res$var[101,], col = "green")
abline(h=0)
abline(v=-1)
abline(v=-0.6)
abline(v=-0.2)
abline(v=0.4)
abline(v=0.8)
abline(v=0)

#Part 2
tempData <- read.csv('https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/TempTullinge.csv', header=TRUE, sep=';')
temp <- tempData$temp
plot(temp, type="l")
time = 1:length(temp)
day = rep(1:365,6)

# Extract every 5:th observation
subset <- seq(1, length(temp), by = 5)
temp <- temp[subset]
time = time[subset]
plot(time,temp, type="l")

sigmaF <- 20
l <- 0.2
polyFit <- lm(scale(temp) ~  scale(time) + I(scale(time)^2))
sigmaNoiseFit = sd(polyFit$residuals)
res<-posteriorGP(X=scale(time),y=scale(temp),k=SEKernel,sigmaNoise=sigmaNoiseFit,xStar=scale(time))
lines(time, res$mu*sd(temp)+mean(temp), col="green", lwd = 2)
lines(time, res$mu*sd(temp)+mean(temp) - 1.96*sd(temp)*sqrt(diag(res$var)), col = "red")
lines(time, res$mu*sd(temp)+mean(temp) + 1.96*sd(temp)*sqrt(diag(res$var)), col = "red")

LM <- function(X,y,k,par){
  n <- length(y)
  L <- t(chol(k(X,X)+((par^2)*diag(n))))
  a <- solve(t(L),solve(L,y))
  logmar <- -0.5*(t(y)%*%a)-sum(log(diag(L)))-(n/2)*log(2*pi)
  return(logmar)
}

# Grid search.

besti<- 0.1
bestLM<-LM(X=scale(time),y=scale(temp),k=SEKernel,par=besti)
bestLM
besti
log_m <- c()
sigma_seq <- seq(0.2,10,0.1)
for(i in sigma_seq){
  aux<-LM(X=scale(time),y=scale(temp),k=SEKernel,par=i)
  log_m <- c(log_m, aux)
}
bestLM
besti
plot(sigma_seq, log_m)
sigma_seq[which.max(log_m)]

res<-posteriorGP(X=scale(time),y=scale(temp),k=SEKernel,sigmaNoise=besti,xStar=scale(time))
lines(time, res$mu*sd(temp)+mean(temp), col="green", lwd = 2)
lines(time, res$mu*sd(temp)+mean(temp) - 1.96*sd(temp)*sqrt(diag(res$var)), col = "blue")
lines(time, res$mu*sd(temp)+mean(temp) + 1.96*sd(temp)*sqrt(diag(res$var)), col = "blue")
