#### Problem 1 ####

library(bnlearn)
data("asia")

set.seed(123)
graph1 <- bnlearn::hc(asia)
plot(graph1)

graph2 <- bnlearn::hc(asia, restart = 5)
plot(graph2)
#Restart: when it hits the local optimum, it restarts again to 
#rxplore new areas in the graph.


#Change the initial structure
custom_initial_structure <- bnlearn::model2network("[A][S][T|A][E][L][B][D][X]")
graph3 <- bnlearn::hc(asia, start = custom_initial_structure)
plot(graph3)
#G doesn't have A -> T link, as this link may be removed to achieve the
#best possible score

# Calculate the BIC score for the learned Bayesian network
score_1 <- bnlearn::score(graph1, data = asia, type = "bde")
score_2 <- bnlearn::score(graph2, data = asia, type = "bde", iss = 10)
score_3 <- bnlearn::score(graph3, data = asia, type = "bde")
#Question: BIC scores are the same

#Change the equivalent sample size (a.k.a imaginary sample size) in the BDeu score
library(graphviz)
fitted_bn <- bn.fit(asia, method = "hc", score = "bdeu", ess = 10)
graphviz.plot(fitted_bn)
# Check if the two BN structures are equivalent
are_equivalent <- all.equal(graph1, graph2)

#### Problem 2: Exact Inference####
library(bnlearn)#To create structure from BN
library(gRain)#For probabilistic inference of Graph

#Split the data
set.seed(123)
train_idx <- sample(1:nrow(asia), size = nrow(asia) * 0.8)
data_train <- asia[train_idx,]
data_test <- asia[-train_idx,]

#Learn the structure
learned_graph <- bnlearn::hc(data_train)
plot(learned_graph)

#Learn the parameters/conditional probability distributions
bn_fit <- bn.fit(learned_graph, data = data_train)

#Perform probabilistic inference like computing posterior probabilities,
#conditional probabilities, classifications ...
#Convert the bn.fit object to a grain
#bn_grain <- as.grain(bn_fit)
junction = compile(as.grain(bn_fit))

pred_exact_inf <- character(length = nrow(data_test))
for (i in 1:nrow(data_test)) {
  # Create an evidence object with the observed variables from the test data
  evidence <- setEvidence(junction, 
                          nodes = c("A", "T", "L", "B", "E", "X", "D"), 
                          states = t(data_test[i ,-2]))
  
  # Query the posterior probability of S(S = yes and S = no)
  posterior <- querygrain(evidence, 
                          nodes = "S",
                          type = "marginal")
  
  # Predict the correct class
  pred_exact_inf[i] <- ifelse(posterior[[1]][1] > posterior[[1]][2], 
                              "no", "yes")
}

#Confusion matrix
true_labels <- as.vector(data_test$S)
confusion_matrix <- table(Actual = true_labels, Predicted = pred_exact_inf)
sum(diag(confusion_matrix))/sum(confusion_matrix)

#Compare with the true DAG
true_dag = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")
plot(true_dag)

#Learn the parameters/conditional probability distributions for true dag
bn_fit_true <- bn.fit(true_dag, data = data_train)
junction = compile(as.grain(bn_fit_true))
pred_true_dag <- character(length = nrow(data_test))
for (i in 1:nrow(data_test)) {
  # Create an evidence object with the observed variables from the test data
  evidence <- setEvidence(junction, 
                          nodes = c("A", "T", "L", "B", "E", "X", "D"), 
                          states = t(data_test[i ,-2]))
  
  # Query the posterior probability of S(S = yes and S = no)
  posterior <- querygrain(evidence, 
                          nodes = "S",
                          type = "marginal")
  
  # Predict the correct class
  pred_true_dag[i] <- ifelse(posterior[[1]][1] > posterior[[1]][2], 
                             "no", "yes")
}

#Confusion matrix
true_labels <- as.vector(data_test$S)
confusion_matrix_true <- table(Actual = true_labels, Predicted = pred_true_dag)
sum(diag(confusion_matrix_true))/sum(confusion_matrix_true)

#### Problem 2: Approximate Inference - cpquery####
library(bnlearn)#To create structure from BN
library(gRain)#For probabilistic inference of Graph

#Split the data
set.seed(123)
train_idx <- sample(1:nrow(asia), size = nrow(asia) * 0.8)
data_train <- asia[train_idx,]
data_test <- asia[-train_idx,]

#Learn the structure
learned_graph <- bnlearn::hc(data_train)
plot(learned_graph)

#Learn the parameters/conditional probability distributions
bn_fit <- bn.fit(learned_graph, data = data_train)

#Perform probabilistic inference like computing posterior probabilities,
#conditional probabilities, classifications ...
pred_approx_inf <- character(length = nrow(data_test))
for (i in 1:nrow(data_test)) {
  # Create an evidence object with the observed variables from the test data
  evidence <- data_test[i ,-2]
  
  # Query the posterior probability of S(S = yes and S = no)
  posterior <- cpquery(bn_fit, 
                       event = (S == "yes"),
                       evidence = as.list(evidence),
                       method = "lw")
  
  # Predict the correct class
  pred_approx_inf[i] <- ifelse(posterior >= 0.5, "yes", "no")
}

#Confusion matrix
true_labels <- as.vector(data_test$S)
confusion_matrix <- table(Actual = true_labels, Predicted = pred_approx_inf)
sum(diag(confusion_matrix))/sum(confusion_matrix)

#### Problem 2: Approximate Inference - cpdist####
library(bnlearn)#To create structure from BN
library(gRain)#For probabilistic inference of Graph

#Split the data
set.seed(123)
train_idx <- sample(1:nrow(asia), size = nrow(asia) * 0.8)
data_train <- asia[train_idx,]
data_test <- asia[-train_idx,]

#Learn the structure
learned_graph <- bnlearn::hc(data_train)
plot(learned_graph)

#Learn the parameters/conditional probability distributions
bn_fit <- bn.fit(learned_graph, data = data_train)

#Perform probabilistic inference like computing posterior probabilities,
#conditional probabilities, classifications ...
pred_approx_inf <- character(length = nrow(data_test))
for (i in 1:nrow(data_test)) {
  # Create an evidence object with the observed variables from the test data
  evidence <- data_test[1 ,-2]
  
  # Query the posterior probability of S(S = yes and S = no)
  posterior_sim <- cpdist(bn_fit, 
                      nodes = "S",
                       evidence = as.list(evidence),
                       method = "lw")
  
  # Predict the correct class
  pred_approx_inf[i] <- ifelse(mean(posterior_sim == 'yes') >= 0.5, "yes", "no")
}

#Confusion matrix
true_labels <- as.vector(data_test$S)
confusion_matrix <- table(Actual = true_labels, Predicted = pred_approx_inf)
sum(diag(confusion_matrix))/sum(confusion_matrix)

#### Problem 3: Classify S given MB####
library(bnlearn)#To create structure from BN
library(gRain)#For probabilistic inference of Graph

#Split the data
set.seed(123456789)
train_idx <- sample(1:nrow(asia), size = nrow(asia) * 0.8)
data_train <- asia[train_idx,]
data_test <- asia[-train_idx,]

#Learn the structure using training
learned_graph <- bnlearn::hc(data_train)
plot(learned_graph)

# Identify the Markov blanket of variable S
mb_S <- mb(learned_graph, node = "S")
#Question: No 'D' in MB

#Learn the parameters/conditional probability distributions
bn_fit <- bn.fit(learned_graph, data = data_train)

#Perform probabilistic inference like computing posterior probabilities,
#conditional probabilities, classifications ...
#Convert the bn.fit object to a grain
#bn_grain <- as.grain(bn_fit)
junction = compile(as.grain(bn_fit))

pred_mb_s <- character(length = nrow(data_test))
for (i in 1:nrow(data_test)) {
  # Create an evidence object with the observed variables from the test data
  evidence <- setEvidence(junction, 
                          nodes = names(data_test[i, mb_S]), 
                          states = t(data_test[i ,mb_S]))
  
  # Query the posterior probability of S(S = yes and S = no)
  posterior <- querygrain(evidence, 
                          nodes = "S",
                          type = "marginal")
  
  # Predict the correct class
  pred_mb_s[i] <- ifelse(posterior[[1]][1] > posterior[[1]][2], "no", "yes")
}

#Confusion matrix
true_labels <- as.vector(data_test$S)
confusion_matrix <- table(Actual = true_labels, Predicted = pred_mb_s)
sum(diag(confusion_matrix))/sum(confusion_matrix)


#### Problem 4: Naive Bayes classifier as BN####
library(bnlearn)#To create structure from BN
library(gRain)#For probabilistic inference of Graph

#S is the parent in Naive Bayes, all observed data is independent when conditioned
#on the class label S
naive_bayes_graph <- model2network("[S][A|S][T|S][L|S][B|S][D|S][E|S][X|S]")
graphviz.plot(naive_bayes_graph)

#Split the data
set.seed(123456789)
train_idx <- sample(1:nrow(asia), size = nrow(asia) * 0.8)
data_train <- asia[train_idx,]
data_test <- asia[-train_idx,]

#Learn the parameters/conditional probability distributions
bn_fit <- bn.fit(naive_bayes_graph, data = data_train)

#Perform probabilistic inference like computing posterior probabilities,
#conditional probabilities, classifications ...
#Convert the bn.fit object to a grain
#bn_grain <- as.grain(bn_fit)
junction = compile(as.grain(bn_fit))

pred_bayes_inf <- character(length = nrow(data_test))
for (i in 1:nrow(data_test)) {
  # Create an evidence object with the observed variables from the test data
  evidence <- setEvidence(junction, 
                          nodes = c("A", "T", "L", "B", "E", "X", "D"), 
                          states = t(data_test[i ,-2]))
  
  # Query the posterior probability of S(S = yes and S = no)
  posterior <- querygrain(evidence, 
                          nodes = "S",
                          type = "marginal")
  
  # Predict the correct class
  pred_bayes_inf[i] <- ifelse(posterior[[1]][1] > posterior[[1]][2], 
                              "no", "yes")
}

#Confusion matrix
true_labels <- as.vector(data_test$S)
confusion_matrix <- table(Actual = true_labels, Predicted = pred_bayes_inf)
sum(diag(confusion_matrix))/sum(confusion_matrix)

