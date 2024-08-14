# 定义目标函数
objective_function <- function(design_matrix, rep_matrix) {
  # 检查是否有超出总体的情况
  col_sums <- colSums(design_matrix * rep_matrix)
  if (any(col_sums != total_per_env)) {
    return(-Inf)
  }
  
  # 确保设计矩阵每行和为2
  if (any(rowSums(design_matrix) != 2)) {
    return(-Inf)
  }
  
  # 确保重复次数矩阵元素不超过3
  if (any(rep_matrix > 3)) {
    return(-Inf)
  }
  
  # 模拟实验计算结果
  result <- tryCatch({
    simulate_experiment(design_matrix, rep_matrix)
  }, error = function(e) {
    message(sprintf("Error occurred in simulate_experiment: %s", e$message))
    return(0)
  })
  
  return(result)
}

# 定义模拟退火算法
simulated_annealing <- function(nGenos, nEnvs, max_iter = 1000, temp_start = 100, temp_end = 1, alpha = 0.99) {
  # 初始化设计矩阵和重复次数矩阵
  current_design_matrix <- matrix(sample(c(0, 1), nGenos * nEnvs, replace = TRUE, prob = c(0.8, 0.2)), nrow = nGenos, ncol = nEnvs)
  current_rep_matrix <- matrix(sample(0:3, nGenos * nEnvs, replace = TRUE), nrow = nGenos, ncol = nEnvs)
  
  current_score <- objective_function(current_design_matrix, current_rep_matrix)
  best_score <- current_score
  best_design_matrix <- current_design_matrix
  best_rep_matrix <- current_rep_matrix
  
  temp <- temp_start
  
  for (iter in 1:max_iter) {
    # 生成新解
    new_design_matrix <- current_design_matrix
    new_rep_matrix <- current_rep_matrix
    
    # 随机选择一个基因型-环境组合，进行突变
    i <- sample(1:nGenos, 1)
    j <- sample(1:nEnvs, 1)
    
    # 对设计矩阵进行变异
    new_design_matrix[i, j] <- ifelse(new_design_matrix[i, j] == 1, 0, 1)
    
    # 对重复次数矩阵进行变异
    new_rep_matrix[i, j] <- sample(0:3, 1)
    
    # 计算新解的适应度
    new_score <- objective_function(new_design_matrix, new_rep_matrix)
    
    # 如果新解更好，或者根据温度接受次优解
    if (new_score > current_score || runif(1) < exp((new_score - current_score) / temp)) {
      current_design_matrix <- new_design_matrix
      current_rep_matrix <- new_rep_matrix
      current_score <- new_score
      
      # 更新最佳解
      if (new_score > best_score) {
        best_score <- new_score
        best_design_matrix <- new_design_matrix
        best_rep_matrix <- new_rep_matrix
      }
    }
    
    # 降低温度
    temp <- temp * alpha
    
    # 打印当前迭代的进度
    if (iter %% 100 == 0) {
      print(paste("Iteration:", iter, "Best Score:", best_score, "Current Score:", current_score, "Temperature:", temp))
    }
  }
  
  return(list(best_design_matrix = best_design_matrix, best_rep_matrix = best_rep_matrix, best_score = best_score))
}

# 运行模拟退火算法
result <- simulated_annealing(nGenos, nEnvs)

# 打印最佳设计矩阵和重复次数矩阵
print("Optimal Design Matrix:")
print(result$best_design_matrix)
print("Optimal Replication Matrix:")
print(result$best_rep_matrix)
print("Optimal Fitness Value:")
print(result$best_score)
