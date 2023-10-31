# discrete Bayesian network from expert knowledge.
net = model2network("[C][A|C][D|C][Y|A:C]")
plot(net)

results <- matrix(NA, nrow = 1000, ncol = 4)
for (i in 1:1000) {
  cond_prob_C_0 <- runif(1)
  cptC <-  matrix(runif(2), 
                ncol = 2, dimnames = list(NULL, c("0", "1")))
  cptC <- prop.table(cptC)
  
  cptA = matrix(runif(4), 
                ncol = 2, 
                nrow = 2,
                dimnames = list("A" = c("0", "1"),
                                "C" = c("0", "1")))
  cptA <- prop.table(cptA,2)
  
  cptD = matrix(runif(4), 
                ncol = 2, 
                nrow = 2,
                dimnames = list("D" = c("0", "1"),
                                "C" = c("0", "1")))
  cptD <- prop.table(cptD,2)
  
  cptY = runif(8)
  dim(cptY) = c(2, 2, 2)
  dimnames(cptY) = list("Y" = c("0", "1"), "C" =  c("0", "1"),
                        "A" = c("0", "1"))
  cptY[,,1] <- prop.table(cptY[,,1],2)
  cptY[,,2] <- prop.table(cptY[,,2],2)
  
  cfit = custom.fit(net, dist = list(A = cptA, D = cptD, C = cptC,
                                     Y = cptY))
  junction = compile(as.grain(cfit))
  evidence <- setEvidence(junction,
                          nodes = c("A", "C"),
                          states = c('1','1'))
  pos_c_1 <- querygrain(evidence,
                        nodes = "Y",
                        type = "marginal")$Y[2]
  evidence <- setEvidence(junction,
                          nodes = c("A", "C"),
                          states = c('1','0'))
  pos_c_2 <- querygrain(evidence,
                        nodes = "Y",
                        type = "marginal")$Y[2]
  evidence <- setEvidence(junction,
                          nodes = c("A", "C"),
                          states = c('0','1'))
  pos_c_3 <- querygrain(evidence,
                        nodes = "Y",
                        type = "marginal")$Y[2]
  evidence <- setEvidence(junction,
                          nodes = c("A", "C"),
                          states = c('0','0'))
  pos_c_4 <- querygrain(evidence,
                        nodes = "Y",
                        type = "marginal")$Y[2]
  non_dec_c <- (pos_c_1 >= pos_c_2) & (pos_c_3 >= pos_c_4)
  non_inc_c <- (pos_c_1 <= pos_c_2) & (pos_c_3 <= pos_c_4)
  mon_c <- ifelse(sum(non_dec_c, non_inc_c) > 0,
                  TRUE, FALSE)
  
  evidence <- setEvidence(junction,
                          nodes = c("A", "D"),
                          states = c('1','1'))
  pos_d_1 <- querygrain(evidence,
                        nodes = "Y",
                        type = "marginal")$Y[2]
  evidence <- setEvidence(junction,
                          nodes = c("A", "D"),
                          states = c('1','0'))
  pos_d_2 <- querygrain(evidence,
                        nodes = "Y",
                        type = "marginal")$Y[2]
  evidence <- setEvidence(junction,
                          nodes = c("A", "D"),
                          states = c('0','1'))
  pos_d_3 <- querygrain(evidence,
                        nodes = "Y",
                        type = "marginal")$Y[2]
  evidence <- setEvidence(junction,
                          nodes = c("A", "D"),
                          states = c('0','0'))
  pos_d_4 <- querygrain(evidence,
                        nodes = "Y",
                        type = "marginal")$Y[2]
  non_dec_d <- (pos_d_1 >= pos_d_2) & (pos_d_3 >= pos_d_4)
  non_inc_d <- (pos_d_1 <= pos_d_2) & (pos_d_3 <= pos_d_4)
  mon_d <- ifelse(sum(non_dec_d, non_inc_d) > 0,
                  TRUE, FALSE)
  results[i,] <- c(non_dec_c,non_inc_c,non_dec_d,non_inc_d)
}


#Jose's solution
library(bnlearn)
library(gRain)

trials <- 1000
results <- matrix(data = NA, nrow = trials, ncol = 4)

for(i in 1:trials){
  net <- model2network("[D|C][C][A|C][Y|A:C]")
  cptC <- runif(2)
  dim(cptC) <- c(2)
  dimnames(cptC) <- list(c("1", "0"))
  cptC <- cptC/sum(cptC)
  cptD <- runif(4) 
  dim(cptD) <- c(2,2)
  dimnames(cptD) <- list("D" = c("1", "0"), "C" =  c("1", "0"))
  cptD <- prop.table(cptD,2)
  cptA <- runif(4) 
  dim(cptA) <- c(2,2)
  dimnames(cptA) <- list("A" = c("1", "0"), "C" =  c("1", "0"))
  cptA <- prop.table(cptA,2)
  cptY <- runif(8) 
  dim(cptY) <- c(2,2,2)
  dimnames(cptY) <- list("Y" = c("1", "0"), "A" =  c("1", "0"), "C" =  c("1", "0"))
  cptY <- prop.table(cptY,2:3)
  netfit <- custom.fit(net,list(C=cptC, D=cptD, A=cptA, Y=cptY))
  netcom <- compile(as.grain(netfit))
  
  pYAC <- querygrain(setEvidence(netcom,nodes=c("A","C"),states=c("1","1")),c("Y"))
  pYAc <- querygrain(setEvidence(netcom,nodes=c("A","C"),states=c("1","0")),c("Y"))
  pYaC <- querygrain(setEvidence(netcom,nodes=c("A","C"),states=c("0","1")),c("Y"))
  pYac <- querygrain(setEvidence(netcom,nodes=c("A","C"),states=c("0","0")),c("Y"))
  
  nondecC <- (pYAc$Y[1] <= pYAC$Y[1] & pYac$Y[1] <= pYaC$Y[1])
  nonincC <- (pYAc$Y[1] >= pYAC$Y[1] & pYac$Y[1] >= pYaC$Y[1])
  
  pYAD <- querygrain(setEvidence(netcom,nodes=c("A","D"),states=c("1","1")),c("Y"))
  pYAd <- querygrain(setEvidence(netcom,nodes=c("A","D"),states=c("1","0")),c("Y"))
  pYaD <- querygrain(setEvidence(netcom,nodes=c("A","D"),states=c("0","1")),c("Y"))
  pYad <- querygrain(setEvidence(netcom,nodes=c("A","D"),states=c("0","0")),c("Y"))
  
  nondecD <- (pYAd$Y[1] <= pYAD$Y[1] & pYad$Y[1] <= pYaD$Y[1])
  nonincD <- (pYAd$Y[1] >= pYAD$Y[1] & pYad$Y[1] >= pYaD$Y[1])
  
  results[i,] <- c(nondecC,nonincC,nondecD,nonincD)
}

colSums(results[which(results[,1]==FALSE & results[,2]==FALSE),])
colSums(results[which(results[,1]==FALSE & results[,2]==TRUE),])
colSums(results[which(results[,1]==TRUE & results[,2]==FALSE),])

0 0 0 0
0 254 128 126
233   0 119 114