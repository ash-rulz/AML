library(HMM)
library(entropy)
set.seed(567)
States=1:10
Symbols=1:10
transProbs=matrix(c(.5,.5,0,0,0,0,0,0,0,0,
                    0,.5,.5,0,0,0,0,0,0,0,
                    0,0,.5,.5,0,0,0,0,0,0,
                    0,0,0,.5,.5,0,0,0,0,0,
                    0,0,0,0,.5,.5,0,0,0,0,
                    0,0,0,0,0,.5,.5,0,0,0,
                    0,0,0,0,0,0,.5,.5,0,0,
                    0,0,0,0,0,0,0,.5,.5,0,
                    0,0,0,0,0,0,0,0,.5,.5,
                    .5,0,0,0,0,0,0,0,0,.5), nrow=length(States), ncol=length(States), byrow = TRUE)
emissionProbs=matrix(c(.2,.2,.2,0,0,0,0,0,.2,.2,
                       .2,.2,.2,.2,0,0,0,0,0,.2,
                       .2,.2,.2,.2,.2,0,0,0,0,0,
                       0,.2,.2,.2,.2,.2,0,0,0,0,
                       0,0,.2,.2,.2,.2,.2,0,0,0,
                       0,0,0,.2,.2,.2,.2,.2,0,0,
                       0,0,0,0,.2,.2,.2,.2,.2,0,
                       0,0,0,0,0,.2,.2,.2,.2,.2,
                       .2,0,0,0,0,0,.2,.2,.2,.2,
                       .2,.2,0,0,0,0,0,.2,.2,.2), nrow=length(States), ncol=length(States), byrow = TRUE)


startProbs=c(.1,.1,.1,.1,.1,.1,.1,.1,.1,.1)
hmm=initHMM(States,Symbols,startProbs,transProbs,emissionProbs)
sim=simHMM(hmm,100)
logf=forward(hmm,sim$observation[1:100])
ef=exp(logf)
pt=prop.table(ef,2)
maxpt=apply(pt,2,which.max)
table(maxpt==sim$states)
post=posterior(hmm,sim$observation[1:100])
maxpost=apply(post,2,which.max)
table(maxpost==sim$states)

# Forward phase
a<-matrix(NA,nrow=100, ncol=length(States))
for(i in States){
  a[1,i]<-emissionProbs[sim$observation[1],i]*startProbs[i]
}

for(t in 2:100){
  for(i in States){
    a[t,i]<-emissionProbs[i,sim$observation[t]]*sum(a[t-1,]*transProbs[,i])
  }
}

for(t in 1:100){
  a[t,]<-a[t,]/sum(a[t,])
}

maxa=apply(a,1,which.max)
table(maxa==sim$states)


# Backward phase

b<-matrix(NA,nrow=100, ncol=length(States))
for(i in States){
  b[100,i]<-1
}

for(t in 99:1){
  for(i in States){
    b[t,i]<-sum(b[t+1,]*emissionProbs[,sim$observation[t+1]]*transProbs[i,])
  }
}

for(t in 1:100){
  for(i in States){
    b[t,i]<-b[t,i]*a[t,i]
  }
  b[t,]<-b[t,]/sum(b[t,])
}

maxb=apply(b,1,which.max)
table(maxb==sim$states)
