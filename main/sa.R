# define object function
objective_function <- function(design_matrix, rep_matrix) {
  # check if exceed
  col_sums <- colSums(design_matrix * rep_matrix)
  if (any(col_sums != total_per_env)) {
    return(-Inf)
  }
  
  # make sure that the sum of each row in the design matrix is 2
  if (any(rowSums(design_matrix) != 2)) {
    return(-Inf)
  }
  
  # make sure the number of replication matrix elements does not exceed 3
  if (any(rep_matrix > 3)) {
    return(-Inf)
  }
  
  # simulation experimental calculation results
  result <- tryCatch({
    simulate_experiment(design_matrix, rep_matrix)
  }, error = function(e) {
    message(sprintf("Error occurred in simulate_experiment: %s", e$message))
    return(0)
  })
  
  return(result)
}

# define the simulated annealing algorithm
simulated_annealing <- function(nGenos, nEnvs, max_iter = 1000, temp_start = 100, temp_end = 1, alpha = 0.99) {
  # init matrices
  current_design_matrix <- matrix(sample(c(0, 1), nGenos * nEnvs, replace = TRUE, prob = c(0.8, 0.2)), nrow = nGenos, ncol = nEnvs)
  current_rep_matrix <- matrix(sample(0:3, nGenos * nEnvs, replace = TRUE), nrow = nGenos, ncol = nEnvs)
  
  current_score <- objective_function(current_design_matrix, current_rep_matrix)
  best_score <- current_score
  best_design_matrix <- current_design_matrix
  best_rep_matrix <- current_rep_matrix
  
  temp <- temp_start
  
  for (iter in 1:max_iter) {
    # generate new solution
    new_design_matrix <- current_design_matrix
    new_rep_matrix <- current_rep_matrix
    
    # Randomly select a gene-environment combination and mutate it
    i <- sample(1:nGenos, 1)
    j <- sample(1:nEnvs, 1)
    
    # Mutate the design matrix
    new_design_matrix[i, j] <- ifelse(new_design_matrix[i, j] == 1, 0, 1)
    
    # Mutate the replication matrix
    new_rep_matrix[i, j] <- sample(0:3, 1)
    
    # calculate new fitness score
    new_score <- objective_function(new_design_matrix, new_rep_matrix)
    
    # if the new solution is better, or accept the suboptimal solution based on the temperature
    if (new_score > current_score || runif(1) < exp((new_score - current_score) / temp)) {
      current_design_matrix <- new_design_matrix
      current_rep_matrix <- new_rep_matrix
      current_score <- new_score
      
      # update best solution
      if (new_score > best_score) {
        best_score <- new_score
        best_design_matrix <- new_design_matrix
        best_rep_matrix <- new_rep_matrix
      }
    }
    
    # lower the temp
    temp <- temp * alpha
    
    # print progress
    if (iter %% 100 == 0) {
      print(paste("Iteration:", iter, "Best Score:", best_score, "Current Score:", current_score, "Temperature:", temp))
    }
  }
  
  return(list(best_design_matrix = best_design_matrix, best_rep_matrix = best_rep_matrix, best_score = best_score))
}

# run here
result <- simulated_annealing(nGenos, nEnvs)

# print results
print("Optimal Design Matrix:")
print(result$best_design_matrix)
print("Optimal Replication Matrix:")
print(result$best_rep_matrix)
print("Optimal Fitness Value:")
print(result$best_score)
