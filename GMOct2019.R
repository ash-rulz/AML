library(bnlearn)
library(gRain)
net = model2network("[C][D][M|C:D]")
plot(net)

# discrete Bayesian network from expert knowledge.
cptC = matrix(c(1/3, 1/3, 1/3), ncol = 3, dimnames = list(NULL, c("1", "2", "3")))
cptD = matrix(c(1/3, 1/3, 1/3), ncol = 3, dimnames = list(NULL, c("1", "2", "3")))
cptM = c(0,0.5,0.5,
         0,0,1,
         0,1,0,
         0,0,1,
         .5,0,.5,
         1,0,0,
         0,1,0,
         1,0,0,
         .5,.5,0)
dim(cptM) = c(3, 3, 3)
dimnames(cptM) = list("M" = c("1", "2", "3"),
                      "D" = c("1", "2", "3"),
                      "C" = c("1", "2", "3"))
cfit = custom.fit(net, dist = list(M = cptM, D = cptD, C = cptC))
junction = compile(as.grain(cfit))
#p(C∣D = d1,M = d2)
evidence <- setEvidence(junction,
                        nodes = c("D", "M"),
                        states = c('1', '2'))
posterior <- querygrain(evidence,
                        nodes = "C",
                        type = "marginal")
#Switch to door 3
#p(C∣D = d1,M = d3)
evidence <- setEvidence(junction,
                        nodes = c("D", "M"),
                        states = c('1', '3'))
posterior <- querygrain(evidence,
                        nodes = "C",
                        type = "marginal")
#Switch to door 2

#Approximate using cpdist
table(cpdist(cfit,"C",evidence=TRUE))
prop.table(table(cpdist(cfit,"C",(D=="1" & M=="2"))))
prop.table(table(cpdist(cfit,"C",(D=="1" & M=="3"))))

#Question 2
net = model2network("[A][B][C|A:B]")
plot(net)

cptA <- matrix(runif(2), ncol = 2, dimnames = list(NULL, c("0", "1")))
cptA <- prop.table(cptA)
cptB <- matrix(runif(2), ncol = 2, dimnames = list(NULL, c("0", "1")))
cptB <- prop.table(cptB)
cptC = c(1,0,0,1,
         0,1,1,0)
dim(cptC) = c(2, 2, 2)
dimnames(cptC) = list("C" = c("1", "0"),
                      "B" = c("1", "0"),
                      "A" = c("0", "1"))

cfit = custom.fit(net, dist = list(C = cptC, B = cptB, A = cptA))

data_samples <- rbn(cfit,1000)

graph <- hc(data_samples)
plot(graph)


#Jose's solution
library(bnlearn)
library(gRain)

xor<-model2network("[A][B][C|A:B]") # Structure
cptA = matrix(c(0.5, 0.5), ncol = 2, dimnames = list(NULL, c("0", "1"))) # Parameters
cptB = matrix(c(0.5, 0.5), ncol = 2, dimnames = list(NULL, c("0", "1")))
cptC = c(1,0,0,1,0,1,1,0)
dim(cptC) = c(2, 2, 2)
dimnames(cptC) = list("C" = c("0", "1"), "A" =  c("0", "1"), "B" = c("0", "1"))
xorfit<-custom.fit(xor,list(A=cptA,B=cptB,C=cptC))

for(i in 1:10){
  xordata<-rbn(xorfit,1000) # Sample
  xorhc<-hc(xordata) # Learn
  plot(xorhc)
}

# The HC fails because the data comes from a distribution that is not faithful to the graph.
# In other words, the graph has an edge A -> C but A is marginally independent of C in the
# distribution. The same for B and C. And, of course, the same for A and B since there is no 
# edge between them. So, the HC cannot find two dependent variables to link with an edge and,
# thus, it adds no edge to the initial empty graph.
