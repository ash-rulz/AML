num_states <- 2
num_symbols <- 2

#Uniform initial probability for each state
initial_probs <- rep(1/num_states, num_states)

#Transition probabilities(Unobserved states to unobserved states)
transition_probs <- matrix(
  c(.9,.1,
    .2,.8), num_states, num_states, byrow = TRUE)

emission_probs <- matrix(
  c(.7,.3,
    .4,.6), num_states, num_symbols, byrow = TRUE)

hmm_model <- initHMM(States = 1:num_states, Symbols = 1:num_symbols,
                     startProbs = initial_probs,
                     transProbs = transition_probs,
                     emissionProbs = emission_probs)
simulated_data <- simHMM(hmm_model, length = 10)
simulated_data$states
simulated_data$observation
get_probabilities <- function(hmm_model, simulated_data){
  observations_sequence <- simulated_data$observation
  alpha <- exp(forward(hmm_model, observations_sequence))
  beta <- exp(backward(hmm_model, observations_sequence))
  filtered_probs <- apply(alpha, 2, prop.table)
  smoothed_probs <- alpha * beta #Unnormalized
  smoothed_probs <- apply(smoothed_probs, 2, prop.table)#Normalized
  return(list(filtered_probs = filtered_probs,
              smoothed_probs = smoothed_probs))
}
prob_list <- get_probabilities(hmm_model, simulated_data)
prob_list$filtered_probs[,c(1:4)]
