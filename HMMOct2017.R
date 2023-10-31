num_states <- 100#There are 100 segments(unobserved)
num_symbols <- 2#There are 2 possible Symbols, in front or not in front(observed)

initial_probs <- rep(1/num_states, num_states)

transition_probs <- matrix(0, num_states, num_states)
for (i in 1:num_states) {
  #It has 2 choices: move right or stay. Cannot move backward
  transition_probs[i,i] <- 0.1
  right <- ifelse(i+1 == 101, 1, i+1)#Considers circular reference
  transition_probs[i,right] <- 0.9
}

emission_probs <- matrix(
  c(rep(0.1,num_states), rep(0.9,num_states)), 
  byrow = FALSE, num_states, num_symbols)
emission_probs[10,] <- c(0.9,0.1)
emission_probs[11,] <- c(0.9,0.1)
emission_probs[12,] <- c(0.9,0.1)

emission_probs[20,] <- c(0.9,0.1)
emission_probs[21,] <- c(0.9,0.1)
emission_probs[22,] <- c(0.9,0.1)

emission_probs[30,] <- c(0.9,0.1)
emission_probs[31,] <- c(0.9,0.1)
emission_probs[32,] <- c(0.9,0.1)

hmm_model <- initHMM(States = 1:num_states, Symbols = 1:num_symbols,
                     startProbs = initial_probs,
                     transProbs = transition_probs,
                     emissionProbs = emission_probs)
