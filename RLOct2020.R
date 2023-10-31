####1st attempt####
library(ggplot2)

arrows <- c("ˆ", ">", "v", "<")
action_deltas <- list(c(1,0), # up
                      c(0,1), # right
                      c(-1,0), # down
                      c(0,-1)) # left
vis_environment <- function(iterations=0, epsilon = 0.5, alpha = 0.1, gamma = 0.95, beta = 0){
  
  # Visualize an environment with rewards. 
  # Q-values for all actions are displayed on the edges of each tile.
  # The (greedy) policy for each state is also displayed.
  # 
  # Args:
  #   iterations, epsilon, alpha, gamma, beta (optional): for the figure title.
  #   reward_map (global variable): a HxW array containing the reward given at each state.
  #   q_table (global variable): a HxWx4 array containing Q-values for each state-action pair.
  #   H, W (global variables): environment dimensions.
  
  df <- expand.grid(x=1:H,y=1:W)
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == 0,q_table[x,y,1],NA),df$x,df$y)
  df$val1 <- as.vector(round(foo, 2))
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == 0,q_table[x,y,2],NA),df$x,df$y)
  df$val2 <- as.vector(round(foo, 2))
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == 0,q_table[x,y,3],NA),df$x,df$y)
  df$val3 <- as.vector(round(foo, 2))
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == 0,q_table[x,y,4],NA),df$x,df$y)
  df$val4 <- as.vector(round(foo, 2))
  foo <- mapply(function(x,y) 
    ifelse(reward_map[x,y] == 0,arrows[GreedyPolicy(x,y)],reward_map[x,y]),df$x,df$y)
  df$val5 <- as.vector(foo)
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == 0,max(q_table[x,y,]),
                                     ifelse(reward_map[x,y]<0,NA,reward_map[x,y])),df$x,df$y)
  df$val6 <- as.vector(foo)
  
  print(ggplot(df,aes(x = y,y = x)) +
          scale_fill_gradient(low = "white", high = "green", na.value = "red", name = "") +
          geom_tile(aes(fill=val6)) +
          geom_text(aes(label = val1),size = 4,nudge_y = .35,na.rm = TRUE) +
          geom_text(aes(label = val2),size = 4,nudge_x = .35,na.rm = TRUE) +
          geom_text(aes(label = val3),size = 4,nudge_y = -.35,na.rm = TRUE) +
          geom_text(aes(label = val4),size = 4,nudge_x = -.35,na.rm = TRUE) +
          geom_text(aes(label = val5),size = 10) +
          geom_tile(fill = 'transparent', colour = 'black') + 
          ggtitle(paste("Q-table after ",iterations," iterations\n",
                        "(epsilon = ",epsilon,", alpha = ",alpha,"gamma = ",gamma,", beta = ",beta,")")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_x_continuous(breaks = c(1:W),labels = c(1:W)) +
          scale_y_continuous(breaks = c(1:H),labels = c(1:H)))
  
}

transition_model <- function(x, y, action, beta){
  # Computes the new state after given action is taken. The agent will follow the action
  # with probability (1-beta) and slip to the right or left with probability beta/2 each.
  #
  # Args:
  # x, y: state coordinates.
  # action: which action the agent takes (in {1,2,3,4}).
  # beta: probability of the agent slipping to the side when trying to move.
  # H, W (global variables): environment dimensions.
  #
  # Returns:
  # The new state after the action has been taken.
  delta <- sample(-1:1, size = 1, prob = c(0.5*beta,1-beta,0.5*beta))
  final_action <- ((action + delta + 3) %% 4) + 1
  foo <- c(x,y) + unlist(action_deltas[final_action])
  foo <- pmax(c(1,1),pmin(foo,c(H,W)))
  return (foo)
}
GreedyPolicy <- function(x, y){
  # Get a greedy action for state (x,y) from q_table.
  #
  # Args:
  # x, y: state coordinates.
  # q_table (global variable): a HxWx4 array containing Q-values for each state-action pair.
  #
  # Returns:
  # An action, i.e. integer in {1,2,3,4}.
  #Get the max q-value from the q_table for the state (x,y)
  max_q = max(q_table[x,y,])
  #Get the indexes which has maximum values
  max_index = which(q_table[x,y,] == max_q)#There can be more than 1 index
  #Question: can the agent stay in the same location?
  #In case of ties, sample randomly from those ties
  if (length(max_index) > 1) {
    action <- sample(max_index, size = 1)
  }else{
    action <- max_index
  }
  return(action)
}
# greedy_action <- GreedyPolicy(1,1)
# greedy_action
EpsilonGreedyPolicy <- function(x, y, epsilon){
  # Get an epsilon-greedy action for state (x,y) from q_table.
  #
  # Args:
  # x, y: state coordinates.
  # epsilon: probability of acting randomly.
  #
  # Returns:
  # An action, i.e. integer in {1,2,3,4}.
  #Explore with epsilon probability or exploit
  explore <- ifelse(epsilon > runif(1), 1, 0)
  if(explore){
    action <- sample(1:4, size = 1)#Exploration
  }else{
    action <- GreedyPolicy(x,y)#Exploitation
  }
  return(action)
}
transition_model <- function(x, y, action, beta){
  # Computes the new state after given action is taken. The agent will follow the action
  # with probability (1-beta) and slip to the right or left with probability beta/2 each.
  #
  # Args:
  # x, y: state coordinates.
  # action: which action the agent takes (in {1,2,3,4}).
  # beta: probability of the agent slipping to the side when trying to move.
  # H, W (global variables): environment dimensions.
  #
  # Returns:
  # The new state after the action has been taken.
  delta <- sample(-1:1, size = 1, prob = c(0.5*beta,1-beta,0.5*beta))
  final_action <- ((action + delta + 3) %% 4) + 1
  #Question: why is action_deltas(y,x)
  foo <- c(x,y) + unlist(action_deltas[final_action])
  foo <- pmax(c(1,1),pmin(foo,c(H,W)))
  return (foo)
}
q_learning <- function(start_state, epsilon = 0.5, alpha = 0.1, gamma = 0.95,
                       beta = 0){
  # Perform one episode of Q-learning. The agent should move around in the
  # environment using the given transition model and update the Q-table.
  # The episode ends when the agent reaches a terminal state.
  #
  # Args:
  # start_state: array with two entries, describing the starting position of the agent.
  # epsilon (optional): probability of acting randomly.
  # alpha (optional): learning rate.
  # gamma (optional): discount factor.
  # beta (optional): slipping factor.
  # reward_map (global variable): a HxW array containing the reward given at each state.
  # q_table (global variable): a HxWx4 array containing Q-values for each state-action pair.
  #
  # Returns:
  # reward: reward received in the episode.
  # correction: sum of the temporal difference correction terms over the episode.
  # q_table (global variable): Recall that R passes arguments by value. So, q_table being
  # a global variable can be modified with the superassigment operator <<-.
  #RL_L1_QLearning, slide 14
  #Initialize S
  old_x <- start_state[1]
  old_y <- start_state[2]
  episode_correction <- 0
  reward <- 0
  repeat{
    #print(paste0('State - (', old_x, ',', old_y,')'))
    #Choose A using policy
    action <- EpsilonGreedyPolicy(old_x,old_y,epsilon)
    #Take A and get the next state(S')
    s_prime <- transition_model(old_x, old_y, action, beta)
    new_x <- s_prime[1]
    new_y <- s_prime[2]
    #Observe the reward for moving into new state
    current_reward <- reward_map[new_x, new_y]
    #Calculate the Q-value temporal difference
    old_q <- q_table[old_x, old_y, action]
    new_q <- max(q_table[new_x, new_y, ])
    temporal_diff <- current_reward + (gamma * new_q) - old_q
    #Update Q-table
    q_table[old_x, old_y, action] <<- old_q + alpha * temporal_diff
    #Update S <- S'
    old_x <- new_x
    old_y <- new_y
    #Compute the episode correction - sum of temporal difference
    #Question: is it absolute or w/o
    episode_correction <- episode_correction + temporal_diff
    #Update the reward for the episode
    reward <- reward + current_reward
    if(current_reward == 10 | current_reward == -10)
      # End episode.
      return (c(reward,episode_correction))
  }
}

sarsa <- function(start_state, epsilon = 0.5, alpha = 0.1, gamma = 0.95,
                       beta = 0){
  #Initialize S
  old_x <- start_state[1]
  old_y <- start_state[2]
  episode_correction <- 0
  reward <- 0
  #Choose A using policy
  action <- EpsilonGreedyPolicy(old_x,old_y,epsilon)
  repeat{
    print(paste0('State - (', old_x, ',', old_y,',', action,')'))
    
    #Take A and get the next state(S')
    s_prime <- transition_model(old_x, old_y, action, beta)
    new_x <- s_prime[1]
    new_y <- s_prime[2]
    
    #Observe the reward for moving into new state
    current_reward <- reward_map[new_x, new_y]
    
    #Calculate the Q-value temporal difference
    old_q <- q_table[old_x, old_y, action]
    new_action <- EpsilonGreedyPolicy(old_x,old_y,epsilon)
    new_q <- q_table[new_x, new_y, new_action]
    temporal_diff <- current_reward + (gamma * new_q) - old_q
    #Update Q-table
    q_table[old_x, old_y, action] <<- old_q + alpha * temporal_diff
    
    #Update S <- S'
    old_x <- new_x
    old_y <- new_y
    
    #Update the action to the new action
    action <- new_action
    
    #Compute the episode correction - sum of temporal difference
    #Question: is it absolute or w/o
    episode_correction <- episode_correction + temporal_diff
    #Update the reward for the episode
    reward <- reward + current_reward
    if(current_reward == 10 | current_reward == -10)
      # End episode.
      return (c(reward,episode_correction))
  }
}

MovingAverage <- function(x, n){
  cx <- c(0,cumsum(x))
  rsum <- (cx[(n+1):length(cx)] - cx[1:(length(cx) - n)]) / n
  return (rsum)
}

# Environment C (modified).
H <- 3
W <- 6
reward_map <- matrix(-1, nrow = H, ncol = W)
reward_map[1,2:5] <- -10
reward_map[1,6] <- 10
#Q-Learning
q_table <- array(0,dim = c(H,W,4))
vis_environment()

#Q-Learning
q_table <- array(0,dim = c(H,W,4))
reward_1 <- NULL
for(i in 1:5000){
  foo <- q_learning(gamma = 1, epsilon = 0.5,
                    beta = 0, alpha = 0.1,
                    start_state = c(1,1))
  reward_1 <- c(reward_1,foo[1])
}
q_table
vis_environment(i, gamma = 1, epsilon = 0.5,
                beta = 0, alpha = 0.1)

#New alg
q_table <- array(0,dim = c(H,W,4))
reward_2 <- NULL
for(i in 1:5000){
  foo <- sarsa(gamma = 1, epsilon = 0.5,
                    beta = 0, alpha = 0.1,
                    start_state = c(1,1))
  reward_2 <- c(reward_2,foo[1])
}
q_table
vis_environment(i, gamma = 1, epsilon = 0.5,
                beta = 0, alpha = 0.1)

#plot(MovingAverage(reward,100),type = "l")

####2nd attempt####
library(ggplot2)
arrows <- c("ˆ", ">", "v", "<")
action_deltas <- list(c(1,0), # up
                      c(0,1), # right
                      c(-1,0), # down
                      c(0,-1)) # left
vis_environment <- function(iterations=0, epsilon = 0.5, alpha = 0.1, gamma = 0.95, beta = 0){
  # Visualize an environment with rewards.
  # Q-values for all actions are displayed on the edges of each tile.
  # The (greedy) policy for each state is also displayed.
  #
  #Args:
  # iterations, epsilon, alpha, gamma, beta (optional): for the figure title.
  # reward_map (global variable): a HxW array containing the reward given at each state.
  # q_table (global variable): a HxWx4 array containing Q-values for each state-action pair.
  # H, W (global variables): environment dimensions.
  df <- expand.grid(x=1:H,y=1:W)
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == -1,q_table[x,y,1],NA),df$x,df$y)
  df$val1 <- as.vector(round(foo, 2))
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == -1,q_table[x,y,2],NA),df$x,df$y)
  df$val2 <- as.vector(round(foo, 2))
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == -1,q_table[x,y,3],NA),df$x,df$y)
  df$val3 <- as.vector(round(foo, 2))
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == -1,q_table[x,y,4],NA),df$x,df$y)
  df$val4 <- as.vector(round(foo, 2))
  foo <- mapply(function(x,y)
    ifelse(reward_map[x,y] == -1,arrows[GreedyPolicy(x,y)],reward_map[x,y]),df$x,df$y)
  df$val5 <- as.vector(foo)
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == -1,max(q_table[x,y,]),
                                     ifelse(reward_map[x,y]<0,NA,reward_map[x,y])),df$x,df$y)
  df$val6 <- as.vector(foo)
  print(ggplot(df,aes(x = y,y = x)) +
          scale_fill_gradient(low = "white", high = "green", na.value = "red", name = "") +
          geom_tile(aes(fill=val6)) +
          geom_text(aes(label = val1),size = 4,nudge_y = .35,na.rm = TRUE) +
          geom_text(aes(label = val2),size = 4,nudge_x = .35,na.rm = TRUE) +
          geom_text(aes(label = val3),size = 4,nudge_y = -.35,na.rm = TRUE) +
          geom_text(aes(label = val4),size = 4,nudge_x = -.35,na.rm = TRUE) +
          geom_text(aes(label = val5),size = 10) +
          geom_tile(fill = 'transparent', colour = 'black') +
          ggtitle(paste("Q-table after ",iterations," iterations\n",
                        "(epsilon = ",epsilon,", alpha = ",alpha,"gamma = ",gamma,", beta = ",beta,")")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_x_continuous(breaks = c(1:W),labels = c(1:W)) +
          scale_y_continuous(breaks = c(1:H),labels = c(1:H)))
}
  
transition_model <- function(x, y, action, beta){
  # Computes the new state after given action is taken. The agent will follow the action
  # with probability (1-beta) and slip to the right or left with probability beta/2 each.
  ##
  Args:
    # x, y: state coordinates.
    # action: which action the agent takes (in {1,2,3,4}).
    # beta: probability of the agent slipping to the side when trying to move.
    # H, W (global variables): environment dimensions.
    ##
    Returns:
    # The new state after the action has been taken.
    delta <- sample(-1:1, size = 1, prob = c(0.5*beta,1-beta,0.5*beta))
    final_action <- ((action + delta + 3) %% 4) + 1
    foo <- c(x,y) + unlist(action_deltas[final_action])
    foo <- pmax(c(1,1),pmin(foo,c(H,W)))
    return (foo)
}
GreedyPolicy <- function(x, y){
  #Get the max q-value from the q_table for the state (x,y)
  max_q = max(q_table[x,y,])
  #Get the indexes which has maximum values
  max_index = which(q_table[x,y,] == max_q)#There can be more than 1 index
  #Question: can the agent stay in the same location?
  #In case of ties, sample randomly from those ties
  if (length(max_index) > 1) {
    action <- sample(max_index, size = 1)
  }else{
    action <- max_index
  }
  return(action)
}
EpsilonGreedyPolicy <- function(x, y, epsilon){
  #Explore with epsilon probability or exploit
  explore <- ifelse(epsilon > runif(1), 1, 0)
  if(explore){
    action <- sample(1:4, size = 1)#Exploration
  }else{
    action <- GreedyPolicy(x,y)#Exploitation
  }
  return(action)
}

transition_model <- function(x, y, action, beta){
  delta <- sample(-1:1, size = 1, prob = c(0.5*beta,1-beta,0.5*beta))
  final_action <- ((action + delta + 3) %% 4) + 1
  #Question: why is action_deltas(y,x)
  foo <- c(x,y) + unlist(action_deltas[final_action])
  foo <- pmax(c(1,1),pmin(foo,c(H,W)))
  return (foo)
}
q_learning <- function(start_state, epsilon = 0.5, alpha = 0.1, gamma = 0.95,
                       beta = 0){
  #Initialize S
  old_x <- start_state[1]
  old_y <- start_state[2]
  episode_correction <- 0
  reward <- 0
  repeat{
    #print(paste0('State - (', old_x, ',', old_y,')'))
    #Choose A using policy
    action <- EpsilonGreedyPolicy(old_x,old_y,epsilon)
    #Take A and get the next state(S')
    s_prime <- transition_model(old_x, old_y, action, beta)
    new_x <- s_prime[1]
    new_y <- s_prime[2]
    #Observe the reward for moving into new state
    current_reward <- reward_map[new_x, new_y]
    #Calculate the Q-value temporal difference
    old_q <- q_table[old_x, old_y, action]
    new_q <- max(q_table[new_x, new_y, ])
    temporal_diff <- current_reward + (gamma * new_q) - old_q
    #Update Q-table
    q_table[old_x, old_y, action] <<- old_q + alpha * temporal_diff
    #Update S <- S'
    old_x <- new_x
    old_y <- new_y
    #Compute the episode correction - sum of temporal difference
    #Question: is it absolute or w/o
    episode_correction <- episode_correction + temporal_diff
    #Update the reward for the episode
    reward <- reward + current_reward
    if(current_reward == 10 | current_reward == -10)
      # End episode.
      return (c(reward,episode_correction))
  }
}
sarsa <- function(start_state, epsilon = 0.5, alpha = 0.1, gamma = 0.95,
                  beta = 0){
  #start_stae <- c(1,1)
  #Initialize S
  episode_correction <- 0
  reward <- 0
  actions <- rep(0,2)#Contains a and a'
  #Contains s, s'(after taking a) and s''(after taking a')
  states <- matrix(rep(0,4), nrow = 2)
  states[1,] <- c(start_state[1], start_state[2])
  rewards <- rep(0,1)#Contains r1(for taking a) and r2(for taking a')
  
  #Produce a, s', r
  actions[1] <- EpsilonGreedyPolicy(states[1,1],states[1,2],epsilon)
  #states[2,] <- transition_model(states[1,1],states[1,2], actions[1], beta)
  #rewards[1] <- reward_map[states[2,1], states[2,2]]
  
  repeat{
    #print(paste0('State - (', old_x, ',', old_y,')'))
    
    #Produce the next state - a', s'', r'
    states[2,] <- transition_model(states[1,1],states[1,2], actions[1], beta)
    rewards[1] <- reward_map[states[2,1], states[2,2]]
    actions[2] <- EpsilonGreedyPolicy(states[2,1],states[2,2],epsilon)
    
    #Calculate the Q-value temporal difference
    old_q <- q_table[states[1,1],states[1,2], actions[1]]
    new_q <- rewards[1] + gamma * q_table[states[2,1],states[2,2], actions[2]] - old_q
    
    #Update Q-table
    q_table[states[1,1],states[1,2], actions[1]] <<- 
      old_q + alpha * new_q
    
    #reward plus current reward
    reward <- reward + rewards[1]
    
    #Update states, rewards, actions
    states[1,] <- states[2,]#Move s' to s
    #states[2,] <- states[3,]#Move s'' to s'
    actions[1] <- actions[2]#Move a' to a
    #rewards[1] <- rewards[2]#Move r' to r
    
    if(rewards[1] == 10 | rewards[1] == -10)
      # End episode.
      return (c(reward,episode_correction))
  }
}

# Environment C (the effect of beta).
H <- 3
W <- 6
reward_map <- matrix(-1, nrow = H, ncol = W)
reward_map[1,2:5] <- -10
reward_map[1,6] <- 10
q_table <- array(0,dim = c(H,W,4))
vis_environment()


reward <- NULL
episode_correction <- NULL
for(i in 1:5000){
  foo <- q_learning(gamma = 1, beta = 0, epsilon = 0.5, alpha = 0.1,
                    start_state = c(1,1))
  reward <- c(reward,foo[1])
  episode_correction <- c(episode_correction,foo[2])
}
vis_environment(gamma = 1, beta = 0, epsilon = 0.5, alpha = 0.1,)
#q_table
plot(MovingAverage(reward,100),type = "l")


reward <- NULL
episode_correction <- NULL
for(i in 1:5000){
  foo <- sarsa(gamma = 1, beta = 0, epsilon = 0.5, alpha = 0.1,
                    start_state = c(1,1))
  reward <- c(reward,foo[1])
  episode_correction <- c(episode_correction,foo[2])
}
vis_environment(gamma = 1, beta = 0, epsilon = 0.5, alpha = 0.1,)
q_table
plot(MovingAverage(reward,100),type = "l")


