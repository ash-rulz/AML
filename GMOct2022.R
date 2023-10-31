library(bnlearn)
library(gRain)
data <- bnlearn::asia

#Estimate the BN parameter values for the True BN - True param
true_dag = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")
plot(true_dag)
#Learn the parameters/conditional probability distributions for true dag
bn_fit_true <- bn.fit(true_dag, data = data)

junction = compile(as.grain(bn_fit_true))
pred_True_D <- character(length = nrow(data))

for (i in 1:nrow(data)) {
  # Create an evidence object with the observed variables from the test data
  evidence <- setEvidence(junction,
                          nodes = c("A", "S", "T", "L", "B", "E", "X"),
                          states = t(data[i ,-8]))
  # Query the posterior probability of S(S = yes and S = no)
  posterior <- querygrain(evidence,
                          nodes = "D",
                          type = "marginal")
  # Predict the correct class
  pred_True_D[i] <- ifelse(posterior[[1]][1] > posterior[[1]][2],
                             "no", "yes")
}
#Confusion matrix
true_labels <- as.vector(data$D)
confusion_matrix_true <- table(Actual = true_labels, Predicted = pred_True_D)
acc_D_True <- sum(diag(confusion_matrix_true))/sum(confusion_matrix_true)

#Use the 10 complete cases to estimate the BN parameter values
data_train <- data[1:10,]
data_test <- data[11:nrow(data),]
data_test$B <- NULL
data_test$E <- NULL
#Learn the parameters/conditional probability distributions for true dag
bn_fit_10 <- bn.fit(true_dag, data = data_train, method = 'bayes')

#Impute the missing values of B and E for the remaining 4990 cases
#Compute the posterior distribution of B, E conditioned on the values for the rest of the variables
#Impute the value for B, E by sampling the posterior distribution
junction = compile(as.grain(bn_fit_10))
pred_B <- character(length = nrow(data_test))
pred_E <- character(length = nrow(data_test))

for (i in 1:nrow(data_test)) {
  # Create an evidence object with the observed variables from the test data
  evidence <- setEvidence(junction,
                          nodes = c("A", "S", "T", "L", "X", "D"),
                          states = t(data_test[i , ]))
  # Query the posterior probability of B
  posterior_B <- querygrain(evidence,
                          nodes = "B",
                          type = "marginal")
  # Query the posterior probability of E
  posterior_E <- querygrain(evidence,
                            nodes = "E",
                            type = "marginal")
  
  # Predict the correct class
  pred_B[i] <- ifelse(posterior_B[[1]][1] > posterior_B[[1]][2],
                             "no", "yes")
  pred_E[i] <- ifelse(posterior_E[[1]][1] > posterior_E[[1]][2],
                      "no", "yes")
}
data_test$B <- pred_B
data_test$E <- pred_E

data_train$B <- NULL
data_train$E <- NULL
data_train$B <- data$B[1:10]
data_train$E <- data$E[1:10]

imputed_data <- rbind(data_train, data_test)

#Use imputed dataset to estimate the BN parameter values.
#Learn the parameters/conditional probability distributions
bn_fit_imputed <- bn.fit(true_dag, data = imputed_data)


bn_fit_true$D
bn_fit_10$D
bn_fit_imputed$D

junction = compile(as.grain(bn_fit_10))
pred_10_D <- character(length = nrow(data))

for (i in 1:nrow(data)) {
  # Create an evidence object with the observed variables from the test data
  evidence <- setEvidence(junction,
                          nodes = c("A", "S", "T", "L", "B", "E", "X"),
                          states = t(data[i ,-8]))
  # Query the posterior probability of S(S = yes and S = no)
  posterior <- querygrain(evidence,
                          nodes = "D",
                          type = "marginal")
  # Predict the correct class
  pred_10_D[i] <- ifelse(posterior[[1]][1] > posterior[[1]][2],
                           "no", "yes")
}
#Confusion matrix
true_labels <- as.vector(data$D)
confusion_matrix_10 <- table(Actual = true_labels, Predicted = pred_10_D)
acc_D_10 <- sum(diag(confusion_matrix_10))/sum(confusion_matrix_10)


junction = compile(as.grain(bn_fit_imputed))
pred_imputed_D <- character(length = nrow(data))

for (i in 1:nrow(data)) {
  # Create an evidence object with the observed variables from the test data
  evidence <- setEvidence(junction,
                          nodes = c("A", "S", "T", "L", "B", "E", "X"),
                          states = t(data[i ,-8]))
  # Query the posterior probability of S(S = yes and S = no)
  posterior <- querygrain(evidence,
                          nodes = "D",
                          type = "marginal")
  # Predict the correct class
  pred_imputed_D[i] <- ifelse(posterior[[1]][1] > posterior[[1]][2],
                         "no", "yes")
}
#Confusion matrix
true_labels <- as.vector(data$D)
confusion_matrix_imputed <- table(Actual = true_labels, Predicted = pred_imputed_D)
acc_D_imputed <- sum(diag(confusion_matrix_imputed))/sum(confusion_matrix_imputed)

acc_D_True
acc_D_10
acc_D_imputed


#Jose's solution
# BNs

library(bnlearn)
library(gRain)

set.seed(123)
data("asia")
tr<-asia
hc<-model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]") # True Asia DAG
hc<-bn.fit(hc,tr[1:10,],method="bayes") # Learn the parameters
hc$D
hc2<-as.grain(hc) # Convert to gRain object
hc2<-compile(hc2) # Compile for inference

for(i in 11:5000){
  z<-NULL
  for(j in c("A","S","T","L","X","D")){
    if(tr[i,j]=="no"){
      z<-c(z,"no")
    }
    else{
      z<-c(z,"yes")
    }
  }
  
  hc3<-setEvidence(hc2,nodes=c("A","S","T","L","X","D"),states=z) # Enter the evidence
  b<-querygrain(hc3,c("B")) # Get posterior distribution
  tr[i,"B"]<-sample(c("no","yes"),size=1,prob=b$B)
  
  hc3<-setEvidence(hc2,nodes=c("A","S","T","L","X","D"),states=z) # Enter the evidence
  e<-querygrain(hc3,c("E")) # Get posterior distribution
  tr[i,"E"]<-sample(c("no","yes"),size=1,prob=e$E)
}

hc<-model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]") # True Asia DAG
hc<-bn.fit(hc,tr,method="bayes") # Learn the parameters
hc$D

hc<-model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]") # True Asia DAG
hc<-bn.fit(hc,asia,method="bayes") # "True" Asia parameters
hc$D

# The parameters learned from the completed dataset are closer to the "true" ones.