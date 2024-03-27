#Creating a binomial tree representation of the Vasicek model:
#Create two function to solve for rates and probabilities on external nodes of tree
#Value of internal nodes will be the expected value of the rate based on the rate given by the "parallel" or inline node from 2 steps ago

#Function to solve for upper branch

solve_binomial_upper <- function(lower_rate, expected_value, variance) {
  # Define the equations based on expected value and variance
  equation1 <- function(upper_rate, p) {
    return(p * upper_rate + (1 - p) * lower_rate - expected_value)
  }
  
  equation2 <- function(upper_rate, p) {
    return(p * (upper_rate - expected_value)^2 + (1 - p) * (lower_rate - expected_value)^2 - variance)
  }
  
  # Define the objective function to minimize
  objective_function <- function(x) {
    upper_rate <- x[1]
    p <- x[2]
    return(equation1(upper_rate, p)^2 + equation2(upper_rate, p)^2)
  }
  
  # Initial guess for upper_rate and p
  initial_guess <- c(expected_value * 1.1, 0.5)
  
  # Use optimization to find the values that minimize the objective function
  result <- optim(initial_guess, objective_function)
  
  # Extract the optimal values
  upper_rate <- result$par[1]
  p <- result$par[2]
  
  result <- return(list(upper_rate = upper_rate, p = p))
}

#Function to calculate lower branch

solve_binomial_lower <- function(upper_rate, expected_value, variance) {
  # Define the equations based on expected value and variance
  equation1 <- function(lower_rate, q) {
    return(q * upper_rate + (1 - q) * lower_rate - expected_value)
  }
  
  equation2 <- function(lower_rate, q) {
    return(q * (upper_rate - expected_value)^2 + (1 - q) * (lower_rate - expected_value)^2 - variance)
  }
  
  # Define the objective function to minimize
  objective_function <- function(x) {
    lower_rate <- x[1]
    q <- x[2]
    return(equation1(lower_rate, q)^2 + equation2(lower_rate, q)^2)
  }
  
  # Initial guess for lower_rate and q
  initial_guess <- c(expected_value * 0.9, 0.5)
  
  # Use optimization to find the values that minimize the objective function
  result <- optim(initial_guess, objective_function)
  
  # Extract the optimal values
  lower_rate <- result$par[1]
  q <- 1 - result$par[2]
  
  return(list(lower_rate = lower_rate, q = q))
}

#Function to find expected value of rate

expected_rate <- function(delta_t, k, theta, rate) {
  exp_rate <- (k * (theta - rate) * (2/12)) + rate
  return(exp_rate)
}

#Function to calculate branches

generate_binomial_tree <- function(steps, initial_value, k, theta, delta_t, sigma) {
  tree <- matrix(0, nrow = steps + 1, ncol = steps + 1)
  tree[1, 1] <- initial_value
  tree[2, 1] <- initial_value + (k * (theta - initial_value) * delta_t) + (sigma * sqrt(delta_t))
  tree[2, 2] <- initial_value + (k * (theta - initial_value) * delta_t) - (sigma * sqrt(delta_t))
  
  #Generate a probability matrix as well
  Prob_Matrix <- matrix(0, nrow = steps + 1, ncol = 2 * (steps + 1))
  Prob_Matrix[1, 1] <- 1
  Prob_Matrix[2, 1] <- Prob_Matrix[2, 2] <- 0.5
  
  # Populate the tree with values
  for (i in 3:(steps + 1)) {
    for (j in 1:i) {
      # Calculate the value at each node
      if (j == 1) {
        result_u <- solve_binomial_upper((expected_rate((2/12), k, theta, tree[i-2, 1])), expected_rate((1/12), k, theta, tree[i-1, j]), (sigma^2 * delta_t))
        tree[i, j] <- result_u$upper_rate
        Prob_Matrix[i, j] <- result_u$p * Prob_Matrix[i-1, j]
      }
      else if (j != 1 & j != i) {
        tree[i, j] <- expected_rate(2/12, k, theta, tree[i-2, j-1])
        Prob_Matrix[i, j] <- 0.5
      }
      else if (j == i) {
        result_d <- solve_binomial_lower((expected_rate((2/12), k, theta, tree[i-2, j-2])), expected_rate((1/12), k, theta, tree[i-1, j-1]), (sigma^2 * delta_t))
        tree[i, j] <- result_d$lower_rate
        Prob_Matrix[i, j] <- result_d$q * Prob_Matrix[i-1, j-1]
      }
    }
  }
  
  # Return the binomial tree matrix and the prob matrix
  return(list(tree = tree, Prob_Matrix = Prob_Matrix))
}

#Running the function

R_0 <- 5.121
long_term_rate <- 15.339
int_sig <- 0.0126 * 100
mean_rev <- 0.025

int_rate_tree <- generate_binomial_tree(6, R_0, mean_rev, long_term_rate, 1/12, int_sig)