library(HMM)
####1: Build a HMM####
num_states <- 4#There are 10 possible States(unobserved)
num_symbols <- 2#

initial_probs <- rep(1/num_states, num_states)

transition_probs <- matrix(c(
  .75,.25,0,0,
  0,0,.5,.5,
  .5,.5,0,0,
  0,0,.25,.75
), num_states, num_states, byrow = TRUE)

emission_probs <- matrix(
  c(.9,.1,
    .1,.9,
    .9,.1,
    .1,.9), num_states, num_symbols, byrow = TRUE)
hmm_model <- initHMM(States = 1:num_states, Symbols = 1:num_symbols,
                     startProbs = initial_probs,
                     transProbs = transition_probs,
                     emissionProbs = emission_probs)
simulated_data <- simHMM(hmm_model, length = 10)
simulated_data$observation
simulated_data$states
