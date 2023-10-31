library(HMM)

####1: Build a HMM####
num_states <- 10#There are 10 possible States(unobserved)
num_symbols <- 10#There are 5 possible Symbols, but 10 states(observed)

#Uniform initial probability for each state
initial_probs <- rep(1/num_states, num_states)

#Transition probabilities(Unobserved states to unobserved states)
transition_probs <- matrix(0, num_states, num_states)
for (i in 1:num_states) {
  #It has 2 choices: move right or stay. Cannot move backward
  transition_probs[i,i] <- 1/2
  
  right <- ifelse(i+1 == 11, 1, i+1)#Considers circular reference
  transition_probs[i,right] <- 1/2
}
transition_probs

#Emission probabilities(Unobserved states to observed symbols)
emission_probs <- matrix(
  c(0.2,0.2,0.2,0,0,0,0,0,0.2,0.2,
    0.2,0.2,0.2,0.2,0,0,0,0,0,0.2,
    0.2,0.2,0.2,0.2,0.2,0,0,0,0,0,
    0,0.2,0.2,0.2,0.2,0.2,0,0,0,0,
    0,0,0.2,0.2,0.2,0.2,0.2,0,0,0,
    0,0,0,0.2,0.2,0.2,0.2,0.2,0,0,
    0,0,0,0,0.2,0.2,0.2,0.2,0.2,0,
    0,0,0,0,0,0.2,0.2,0.2,0.2,0.2,
    0.2,0,0,0,0,0,0.2,0.2,0.2,0.2,
    0.2,0.2,0,0,0,0,0,0.2,0.2,0.2), num_states, num_states)
emission_probs

#Create the HMM
hmm_model <- initHMM(States = 1:num_states, Symbols = 1:num_symbols,
                     startProbs = initial_probs,
                     transProbs = transition_probs,
                     emissionProbs = emission_probs)

####2: Simulate HMM for 100 steps####
#Simulate the HMM for 100 time steps
run_simulation <- function(num_steps, hmm_model){
  #set.seed(123)  
  simulated_data <- simHMM(hmm_model, length = num_steps)
  
  # All the observations are within [i − 2, i + 2] distance 
  return(simulated_data)
}
num_steps <- 100 #Time steps
simulated_data <- run_simulation(num_steps, hmm_model)

#Extract the sequence of states and observations
states_sequence <- simulated_data$states
observations_sequence <- simulated_data$observation

#Print the sequences (states and observations)
print("States Sequence:")
print(states_sequence)
print("Observations Sequence:")
print(observations_sequence)

####3: Compute filtered and smoothed probability distributions. Also most probable path####
get_probabilities <- function(hmm_model, simulated_data){
  
  #Extract the sequence of observations
  observations_sequence <- simulated_data$observation
  
  #Get the non-normalized forward probabilities(state*time steps)
  alpha <- exp(forward(hmm_model, observations_sequence))
  
  #Get the backward probabilities
  beta <- exp(backward(hmm_model, observations_sequence))
  
  #Compute filtered probability distribution
  filtered_probs <- apply(alpha, 2, prop.table)
  
  #Compute smoothed probability distribution
  smoothed_probs <- alpha * beta #Unnormalized
  smoothed_probs <- apply(smoothed_probs, 2, prop.table)#Normalized
  
  # Compute the most probable path
  most_probable_path <- viterbi(hmm_model, observations_sequence)
  
  return(list(filtered_probs = filtered_probs,
              smoothed_probs = smoothed_probs,
              most_probable_path = most_probable_path))
}
prob_list <- get_probabilities(hmm_model, simulated_data)
prob_list$filtered_probs[,c(98,99,100)]
prob_list$smoothed_probs[,c(98,99,100)]
prob_list$most_probable_path

####4: Compute the accuracy of filtered, smoothed and probable path####
compute_accuracy <- function(prob_distr, true_path, probable_path = NULL){
  if(is.null(probable_path)){
    path <- apply(prob_distr,2, which.max)
    return((sum(path == true_path)/num_steps))
  }else{
    return((sum(probable_path == true_path)/num_steps))
  }
}

states_sequence <- simulated_data$states

#Get filtered probability distribution accuracy
prob_distr <- prob_list$filtered_probs
filtered_accuracy <- compute_accuracy(prob_distr, states_sequence)
print(paste0('Filtered accuracy:', filtered_accuracy))

#Get smoothed probability distribution accuracy
prob_distr <- prob_list$smoothed_probs
smoothed_accuracy <- compute_accuracy(prob_distr, states_sequence)
print(paste0('Smoothed accuracy:', smoothed_accuracy))

#Get most probable path accuracy
most_probable_path <- prob_list$most_probable_path
probable_path_accuracy <- compute_accuracy(prob_distr, states_sequence,
                                           most_probable_path)
print(paste0('Probable path accuracy:', probable_path_accuracy))

#5: Repeat the previous exercise with different simulated samples
filtered_accuracies <- c()
smoothed_accuracies <- c()
probable_path_accuracies <- c()
for (i in seq(1,100,1)) {
  
  num_steps <- 100 #Time steps
  simulated_data <- run_simulation(num_steps, hmm_model)
  states_sequence <- simulated_data$states
  
  prob_list <- get_probabilities(hmm_model, simulated_data)
  
  #Get filtered probability distribution accuracy
  prob_distr <- prob_list$filtered_probs
  filtered_accuracies <- c(filtered_accuracies,
                           compute_accuracy(prob_distr, states_sequence))
  
  #Get smoothed probability distribution accuracy
  prob_distr <- prob_list$smoothed_probs
  smoothed_accuracies <- c(smoothed_accuracies,
                           compute_accuracy(prob_distr, states_sequence))
  
  #Get most probable path accuracy
  most_probable_path <- prob_list$most_probable_path
  probable_path_accuracies <- c(probable_path_accuracies,
                                compute_accuracy(prob_distr, states_sequence,
                                             most_probable_path))
}
plot(filtered_accuracies,type = 'l', col = 'red',
     ylim = c(0.2, 1))
lines(smoothed_accuracies,type = 'l', col = 'blue')
lines(probable_path_accuracies,type = 'l', col = 'green')
# Adding legends
legend("topright", 
       legend = c("Filtered", "Smoothed", "Probable Path"),
       col = c("red", "blue", "green"),
       lty = 1,
       cex = 0.5)

#6: Is it always true that the later in time (i.e., the more observations you have received)
# the better you know where the robot is ?
library(entropy)

num_steps <- 300 #Time steps
simulated_data <- run_simulation(num_steps, hmm_model)

prob_list <- get_probabilities(hmm_model, simulated_data)

#Get filtered probability distribution accuracy
filtered_probs <- prob_list$filtered_probs

v_entropy <- c()
for (i in seq(30,num_steps,5)) {
  v_entropy <- c(v_entropy, entropy.empirical(filtered_probs[,i]))
}
plot(v_entropy,type = 'l', xlab = 'Time steps', ylab = 'Entropy')

#7: Compute the probabilities of the hidden states for the time step 101 given
# previous observations
num_steps <- 100 #Time steps
simulated_data <- run_simulation(num_steps, hmm_model)

prob_list <- get_probabilities(hmm_model, simulated_data)

#Get filtered probability distribution accuracy
filtered_probs <- prob_list$filtered_probs

#Product of the filtered prob at state 100 and the transition probabilities
filtered_probs[,100] %*% transition_probs
