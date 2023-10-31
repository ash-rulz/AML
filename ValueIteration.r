value_iteration <- function(epsilon = 0.0001, gamma = 0.95) {
  # Initialize the value function arbitrarily (V0)
  V <- matrix(0, nrow = H, ncol = W)  # H and W are the height and width of the grid world

  while (TRUE) {
    delta <- 0  # Initialize the change in value function

    for (x in 1:H) {
      for (y in 1:W) {
        if (is_terminal(x, y)) {
          # Skip terminal states
          next
        }

        # Calculate the new value using the Bellman equation
        max_action_value <- -Inf
        for (action in 1:4) {  # Assuming there are four possible actions
          s_prime <- transition_model(x, y, action, beta)
          new_x <- s_prime[1]
          new_y <- s_prime[2]
          current_reward <- reward_map[new_x, new_y]
          new_value <- current_reward + gamma * V[new_x, new_y]
          max_action_value <- max(max_action_value, new_value)
        }

        # Update the change in value function
        delta <- max(delta, abs(max_action_value - V[x, y]))

        # Update the value function
        V[x, y] <- max_action_value
      }
    }

    # Check for convergence
    if (delta < epsilon) {
      break
    }
  }

  # Now, V contains the optimal value function V*.

  # To extract the optimal policy π*, you can use the value function:
  policy <- matrix(0, nrow = H, ncol = W)  # Initialize the policy
  for (x in 1:H) {
    for (y in 1:W) {
      if (is_terminal(x, y)) {
        # Skip terminal states
        next
      }

      # Find the action that maximizes the expected value
      max_action_value <- -Inf
      max_action <- 0
      for (action in 1:4) {  # Assuming there are four possible actions
        s_prime <- transition_model(x, y, action, beta)
        new_x <- s_prime[1]
        new_y <- s_prime[2]
        current_reward <- reward_map[new_x, new_y]
        new_value <- current_reward + gamma * V[new_x, new_y]
        if (new_value > max_action_value) {
          max_action_value <- new_value
          max_action <- action
        }
      }

      # Set the optimal action in the policy
      policy[x, y] <- max_action
    }
  }

  return(list(ValueFunction = V, Policy = policy))
}
