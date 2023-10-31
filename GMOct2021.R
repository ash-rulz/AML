true_dag = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")
plot(true_dag)

#Sample 1000 points - using approximate inference
bn_fit_true <- bn.fit(true_dag, data = asia)
sample_d <- c()
sample_s <- c()
for (i in 1:1000) {
  a <- sample(1:2, 1, prob = bn_fit_true$A$prob)
  s <- sample(1:2, 1, prob = bn_fit_true$S$prob)
  t <- sample(1:2, 1, prob = bn_fit_true$T$prob[,a])
  l <- sample(1:2, 1, prob = bn_fit_true$L$prob[,s])
  b <- sample(1:2, 1, prob = bn_fit_true$B$prob[,s])
  e <- sample(1:2, 1, prob = bn_fit_true$E$prob[,l,t])
  d <- sample(1:2, 1, prob = bn_fit_true$D$prob[,b,e])
  x <- sample(1:2, 1, prob = bn_fit_true$X$prob[,e])  
  sample_s <- c(sample_s, s)
  sample_d <- c(sample_d, d)
}
cond_s <- sample_s[which(sample_d == 2)]
sum(cond_s == 2)/length(cond_s)#P(S='no'|d='yes')
sum(cond_s == 1)/length(cond_s)#P(S='yes'|d='yes')

#Compare with true
junction = compile(as.grain(bn_fit_true))
evidence <- setEvidence(junction,
                        nodes = c("D"),
                        states = "yes")
posterior <- querygrain(evidence,
                        nodes = "S",
                        type = "marginal")
#It matches
#Jose's code
library(bnlearn)
library(gRain)

data("asia")
net<-model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")
net<-bn.fit(net,asia,method="bayes")

mydata<-matrix(NA,1000,8)
for(i in 1:1000){
  a<-sample(1:2,1,prob=net$A$prob)
  s<-sample(1:2,1,prob=net$S$prob)
  t<-sample(1:2,1,prob=net$T$prob[,a])
  
  l<-sample(1:2,1,prob=net$L$prob[,s])
  b<-sample(1:2,1,prob=net$B$prob[,s])
  
  e<-sample(1:2,1,prob=net$E$prob[,l,t])
  d<-sample(1:2,1,prob=net$D$prob[,b,e])
  
  x<-sample(1:2,1,prob=net$X$prob[,e])
  
  mydata[i,]<-c(a,s,t,l,b,d,e,x)
}

foo<-mydata[which(mydata[,6]==2),2]
table(foo)/length(foo)

net<-as.grain(net)
net<-compile(net)
net<-setEvidence(net,nodes=c("D"),states=c("yes"))
querygrain(net,c("S"))
