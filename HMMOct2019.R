num_states <- 10#There are 10 possible States(unobserved)
num_symbols <- 11#There are 6 possible Symbols, but 11 states(observed)

#We want it to start from 1 always
initial_probs <-  rep(1/num_states, num_states)

transition_probs <- matrix(0, num_states, num_states)
for (i in 1:num_states) {
  #It has 2 choices: move right or stay. Cannot move backward
  transition_probs[i,i] <- 1/2
  right <- ifelse(i+1 == 11, 1, i+1)#Considers circular reference
  transition_probs[i,right] <- 1/2
}
#Emission probabilities(Unobserved states to observed symbols)
emission_probs <- matrix(
  c(0.1,0.1,0.1,0,0,0,0,0,0.1,0.1,0.5,
    0.1,0.1,0.1,0.1,0,0,0,0,0,0.1,0.5,
    0.1,0.1,0.1,0.1,0.1,0,0,0,0,0,0.5,
    0,0.1,0.1,0.1,0.1,0.1,0,0,0,0,0.5,
    0,0,0.1,0.1,0.1,0.1,0.1,0,0,0,0.5,
    0,0,0,0.1,0.1,0.1,0.1,0.1,0,0,0.5,
    0,0,0,0,0.1,0.1,0.1,0.1,0.1,0,0.5,
    0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.5,
    0.1,0,0,0,0,0,0.1,0.1,0.1,0.1,0.5,
    0.1,0.1,0,0,0,0,0,0.1,0.1,0.1,0.5), 
  byrow = TRUE,
  num_states, num_symbols)

hmm_model <- initHMM(States = 1:num_states, Symbols = 1:num_symbols,
                     startProbs = initial_probs,
                     transProbs = transition_probs,
                     emissionProbs = emission_probs)

hmm_model

simulated_data <- simHMM(hmm_model, length = 10)
simulated_data$states
simulated_data$observation

observations_sequence <- c(1,11,11,11)
alpha <- exp(forward(hmm_model, observations_sequence))
beta <- exp(backward(hmm_model, observations_sequence))
smoothed_probs <- alpha * beta #Unnormalized
smoothed_probs <- apply(smoothed_probs, 2, prop.table)#Normalized
apply(smoothed_probs,2,which.max)
posterior(hmm_model,observations_sequence)

most_probable_path <- viterbi(hmm_model, observations_sequence)
