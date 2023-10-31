library(bnlearn)
library(gRain)
data(chest_cpt)

init_prob <- c(.5,.5)
emission_prob <- c(.6,.4,
                   .3,.7)
transition_prob <- c(.9,.1,
                     .8,.2)

cptZ1 = matrix(init_prob)
dim(cptZ1) = c(1, 2)
dimnames(cptZ1) = list(NULL, c("H", "S"))
cptX1 = matrix(emission_prob, byrow = FALSE)
dim(cptX1) = c(2, 2)
dimnames(cptX1) = list(X1 = c("H", "S"), 
                       Z1 = c("H", "S"))

cptZ2 = matrix(transition_prob, byrow = FALSE)
dim(cptZ2) = c(2, 2)
dimnames(cptZ2) = list(Z2 = c("H", "S"), 
                       Z1 = c("H", "S"))
cptX2 = matrix(emission_prob, byrow = FALSE)
dim(cptX2) = c(2, 2)
dimnames(cptX2) = list(X2 = c("H", "S"), 
                       Z2 = c("H", "S"))

cptZ3 = matrix(transition_prob, byrow = FALSE)
dim(cptZ3) = c(2, 2)
dimnames(cptZ3) = list(Z3 = c("H", "S"), 
                       Z2 = c("H", "S"))
cptX3 = matrix(emission_prob, byrow = FALSE)
dim(cptX3) = c(2, 2)
dimnames(cptX3) = list(X3 = c("H", "S"), 
                       Z3 = c("H", "S"))

cptZ4 = matrix(transition_prob, byrow = FALSE)
dim(cptZ4) = c(2, 2)
dimnames(cptZ4) = list(Z4 = c("H", "S"), 
                       Z3 = c("H", "S"))
cptX4 = matrix(emission_prob, byrow = FALSE)
dim(cptX4) = c(2, 2)
dimnames(cptX4) = list(X4 = c("H", "S"), 
                       Z4 = c("H", "S"))

net <- model2network(paste0("[Z1][X1|Z1]",
                            "[Z2|Z1][X2|Z2]",
                            "[Z3|Z2][X3|Z3]",
                            "[Z4|Z3][X4|Z4]"
                            ))

cfit = custom.fit(net, dist = list(Z1 = cptZ1, X1 = cptX1, 
                                   Z2 = cptZ2, X2 = cptX2,
                                   Z3 = cptZ3, X3 = cptX3,
                                   Z4 = cptZ4, X4 = cptX4
                                   ))
junction = compile(as.grain(cfit))
evidence <- setEvidence(junction,c("Z1"),list(c(0.5, 0.5)))
posterior <- querygrain(evidence,
                        nodes = "Z4",
                        type = "marginal")
posterior

evidence <- setEvidence(junction,c("X2"), c('H'))
posterior <- querygrain(evidence,
                        nodes = "Z4",
                        type = "marginal")
posterior

evidence <- setEvidence(junction,c("X2", "X3"), c('H', 'H'))
posterior <- querygrain(evidence,
                        nodes = "Z4",
                        type = "marginal")
posterior
