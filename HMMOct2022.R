library(HMM)
####1: Build a HMM####
num_states <- 10#There are 10 possible States(unobserved)
num_symbols <- 10#There are 5 possible Symbols, but 10 states(observed)
#Uniform initial probability for each state
initial_probs <- rep(1/num_states, num_states)

#Transition probabilities(Unobserved states to unobserved states)
transition_probs <- matrix(
  c(0.5,0.5,0,0,0,0,0,0,0,0,
    0,0.5,0.5,0,0,0,0,0,0,0,
    0,0,0.5,0.5,0,0,0,0,0,0,
    0,0,0,0.5,0.5,0,0,0,0,0,
    0,0,0,0,0.5,0.5,0,0,0,0,
    0,0,0,0,0,0.5,0.5,0,0,0,
    0,0,0,0,0,0,0.5,0.5,0,0,
    0,0,0,0,0,0,0,0.5,0.5,0,
    0,0,0,0,0,0,0,0,0.5,0.5,
    0.5,0,0,0,0,0,0,0,0,0.5
    ), byrow = TRUE,
  num_states, num_states)

#Emission probabilities(Unobserved states to observed symbols)
emission_probs <- matrix(
  c(1/6,1/6,1/9,1/9,1/9,0,0,0,1/6,1/6,#1.1
    1/6,1/6,1/9,1/9,1/9,0,0,0,1/6,1/6,#1.2
    1/6,1/6,1/9,1/9,1/9,1/6,1/6,0,0,0,#2.1
    1/6,1/6,1/9,1/9,1/9,1/6,1/6,0,0,0,#2.2
    1/6,1/6,1/9,1/9,1/9,1/6,1/6,0,0,0,#2.3
    0,0,1/9,1/9,1/9,1/6,1/6,1/3,0,0,#3.1
    0,0,1/9,1/9,1/9,1/6,1/6,1/3,0,0,#3.2
    0,0,0,0,0,1/6,1/6,1/3,1/6,1/6,#4
    1/6,1/6,0,0,0,0,0,1/3,1/6,1/6,#5.1
    1/6,1/6,0,0,0,0,0,1/3,1/6,1/6#5.1
    ), byrow = TRUE,
  num_states, num_states)

#Create the HMM
hmm_model <- initHMM(States = 1:num_states, Symbols = 1:num_symbols,
                     startProbs = initial_probs,
                     transProbs = transition_probs,
                     emissionProbs = emission_probs)

simulated_data <- simHMM(hmm_model, length = 20)
# All the observations are within [i - 1, i + 1] distance

states_sequence <- simulated_data$states
observations_sequence <- simulated_data$observation

