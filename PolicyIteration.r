policy_iteration <- function(gamma = 0.95) {
  # Initialize the policy arbitrarily
  policy <- matrix(sample(1:4, H * W, replace = TRUE), nrow = H, ncol = W)
  
  while (TRUE) {
    # Policy Evaluation: Evaluate the current policy
    V <- matrix(0, nrow = H, ncol = W)  # Initialize the value function arbitrarily
    
    while (TRUE) {
      delta <- 0  # Initialize the change in value function
      
      for (x in 1:H) {
        for (y in 1:W) {
          if (is_terminal(x, y)) {
            # Skip terminal states
            next
          }
          
          action <- policy[x, y]
          s_prime <- transition_model(x, y, action, beta)
          new_x <- s_prime[1]
          new_y <- s_prime[2]
          current_reward <- reward_map[new_x, new_y]
          new_value <- current_reward + gamma * V[new_x, new_y]
          
          delta <- max(delta, abs(new_value - V[x, y]))
          V[x, y] <- new_value
        }
      }
      
      # Check for convergence
      if (delta < epsilon) {
        break
      }
    }
    
    # Policy Improvement: Improve the policy using the current value function
    policy_stable <- TRUE
    
    for (x in 1:H) {
      for (y in 1:W) {
        if (is_terminal(x, y)) {
          # Skip terminal states
          next
        }
        
        old_action <- policy[x, y]
        max_action_value <- -Inf
        best_action <- 0
        
        for (action in 1:4) {  # Assuming there are four possible actions
          s_prime <- transition_model(x, y, action, beta)
          new_x <- s_prime[1]
          new_y <- s_prime[2]
          current_reward <- reward_map[new_x, new_y]
          new_value <- current_reward + gamma * V[new_x, new_y]
          
          if (new_value > max_action_value) {
            max_action_value <- new_value
            best_action <- action
          }
        }
        
        # Update the policy if a better action is found
        if (best_action != old_action) {
          policy_stable <- FALSE
          policy[x, y] <- best_action
        }
      }
    }
    
    # If the policy is stable (no change), we have found the optimal policy
    if (policy_stable) {
      break
    }
  }
  
  return(list(Policy = policy))
}
