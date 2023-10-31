library(bnlearn)
library(gRain)
net = model2network(paste0("[Z1][C1][X1|Z1]",
                           "[Z2|Z1:C1][C2|C1][X2|Z2]",
                           "[Z3|Z2:C2][C3|C2][X3|Z3]",
                           "[Z4|Z3:C3][C4|C3][X4|Z4]"))
graphviz.plot(net)

init_Z_prob <- c(0.5,0.6)#Z1
init_C_prob <- c(0,0,0,1)#C1
C_trans <- c(0,0,0,1,
             1,0,0,0,
             0,1,0,0,
             0,0,1,0)
X <- c(0.6,0.4,0.7,0.3)# Z-X remains the same
Z_trans1 <- c(0.9,0.1,0,1)#If C > 1
Z_trans2 <- c(0.9,0.1,0.2,0.8)#If C = 1

cptZ1 = matrix(init_Z_prob, ncol = 2, byrow = FALSE)
dimnames(cptZ1) = list(NULL, c("H", "S"))
cptC1 = matrix(init_C_prob, ncol = 4, , byrow = FALSE)
dimnames(cptC1) = list(NULL, c("1", "2", "3", "4"))
cptX1 = matrix(X, ncol = 2, nrow = 2, byrow = FALSE)
dimnames(cptX1) = list(c("H", "S"), c("H", "S"))

cptZ2 <- matrix(c(Z_trans2, Z_trans1, Z_trans1, Z_trans1), byrow = FALSE)
dim(cptZ2) <- c(2,2,4)
dimnames(cptZ2) = list("Z2" = c("H", "S"), "Z1" =  c("H", "S"),
                      "C1" = c("1", "2", "3", "4"))
cptC2 = matrix(C_trans, byrow = FALSE)
dim(cptC2) <- c(4,4)
dimnames(cptC2) = list(c("1", "2", "3","4"), c("1", "2", "3", "4"))
cptX2 = matrix(X, ncol = 2, nrow = 2, byrow = FALSE)
dimnames(cptX2) = list(c("H", "S"), c("H", "S"))

cptZ3 <- matrix(c(Z_trans2, Z_trans1, Z_trans1, Z_trans1), byrow = FALSE)
dim(cptZ3) <- c(2,2,4)
dimnames(cptZ3) = list("Z3" = c("H", "S"), "Z2" =  c("H", "S"),
                       "C2" = c("1", "2", "3", "4"))
cptC3 = matrix(C_trans, byrow = FALSE)
dim(cptC3) <- c(4,4)
dimnames(cptC3) = list(c("1", "2", "3","4"), c("1", "2", "3", "4"))
cptX3 = matrix(X, ncol = 2, nrow = 2, byrow = FALSE)
dimnames(cptX3) = list(c("H", "S"), c("H", "S"))

cptZ4 <- matrix(c(Z_trans2, Z_trans1, Z_trans1, Z_trans1), byrow = FALSE)
dim(cptZ4) <- c(2,2,4)
dimnames(cptZ4) = list("Z4" = c("H", "S"), "Z3" =  c("H", "S"),
                       "C3" = c("1", "2", "3", "4"))
cptC4 = matrix(C_trans, byrow = FALSE)
dim(cptC4) <- c(4,4)
dimnames(cptC4) = list(c("1", "2", "3","4"), c("1", "2", "3", "4"))
cptX4 = matrix(X, ncol = 2, nrow = 2, byrow = FALSE)
dimnames(cptX4) = list(c("H", "S"), c("H", "S"))

cfit = custom.fit(net, dist = list(C1 = cptC1,Z1 = cptZ1, X1 = cptX1,
                                   C2 = cptC2, Z2 = cptZ2, X2 = cptX2,
                                   C3 = cptC3, Z3 = cptZ3, X3 = cptX3,
                                   C4 = cptC4, Z4 = cptZ4, X4 = cptX4))
junction = compile(as.grain(cfit))
evidence <- setEvidence(junction,
                        nodes = c(""),
                        states = c(""))
posterior <- querygrain(evidence,
                        nodes = c("Z1", "Z2", "Z3", "Z4",
                                  "C1", "C2", "C3", "C4"),
                        type = "marginal")
as.vector(ifelse(posterior$Z1[1] > posterior$Z1[2], 'H', 'S'))
as.vector(ifelse(posterior$Z2[1] > posterior$Z2[2], 'H', 'S'))
as.vector(ifelse(posterior$Z3[1] > posterior$Z3[2], 'H', 'S'))
as.vector(ifelse(posterior$Z4[1] > posterior$Z4[2], 'H', 'S'))
