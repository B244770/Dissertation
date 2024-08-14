library(xgboost)
library(caret)
library(ggplot2)
library(reshape2)

# function for extracting features
extract_features_from_matrices <- function(matrices) {
  # expand the matrices to vectors
  design_vector <- as.vector(matrices$design_matrix)
  rep_vector <- as.vector(matrices$rep_matrix)
  
  # 合并向量
  combined_vector <- c(design_vector, rep_vector)
  return(combined_vector)
}

# 从结果中提取数据的函数
extract_data <- function(results_all_configs) {
  results_list <- list()
  
  for (env_config in names(results_all_configs)) {
    config_results <- results_all_configs[[env_config]]
    for (run_index in seq_along(config_results)) {
      score <- as.numeric(config_results[[run_index]]$score)
      features_vector <- extract_features_from_matrices(config_results[[run_index]]$matrices)
      
      # 直接存储分数和特征向量在一个列表中
      results_list[[length(results_list) + 1]] <- list(Score = score, Features = features_vector)
    }
  }
  
  return(results_list)
}

# 假设已经加载了 results_all_configs 数据
# results_all_configs <- choose 20 envs as example
nEnvs <- 20    # 环境数量
results_list <- extract_data(results_all_configs["20"])

# # 转换数据格式
# features_matrix <- do.call(rbind, lapply(results_list, function(x) x$Features))
# scores_vector <- sapply(results_list, function(x) x$Score)
# 找到最大的向量长度
max_length <- max(sapply(results_list, function(x) length(x$Features)))

# 使用 lapply 进行填充，然后合并
features_matrix <- do.call(rbind, lapply(results_list, function(x) {
  # 填充 Features 向量到 max_length，用 NA 补充
  c(x$Features, rep(NA, max_length - length(x$Features)))
}))

scores_vector <- sapply(results_list, function(x) x$Score)

# 提取列名
colnames(features_matrix) <- paste0("V", 1:ncol(features_matrix))

# 标准化特征
scaled_features <- scale(features_matrix)

# 恢复列名
colnames(scaled_features) <- colnames(features_matrix)

# 使用XGBoost评估特征重要性
dtrain <- xgb.DMatrix(data = scaled_features, label = scores_vector)
xgb_model <- xgboost(data = dtrain, max.depth = 3, nrounds = 100, objective = "reg:squarederror", importance_type = "gain")

importance_matrix <- xgb.importance(model = xgb_model)
selected_features <- importance_matrix$Feature[importance_matrix$Gain > mean(importance_matrix$Gain)]

# 选择重要特征
reduced_features_matrix <- scaled_features[, selected_features, drop = FALSE]

#################################################

# 拟合XGBoost模型
dtrain_reduced <- xgb.DMatrix(data = scaled_features, label = scores_vector)
xgb_model_reduced <- xgboost(data = dtrain_reduced, max.depth = 3, nrounds = 500, objective = "reg:squarederror", importance_type = "gain")
print(summary(xgb_model_reduced))

#################################################

# 定义交叉验证函数
cross_validation_xgb <- function(features, scores, folds = 5) {
  set.seed(123)
  n <- nrow(features)
  fold_ids <- sample(rep(1:folds, length.out = n))
  
  cv_results <- sapply(1:folds, function(f) {
    train_features <- features[fold_ids != f, , drop = FALSE]
    test_features <- features[fold_ids == f, , drop = FALSE]
    train_scores <- scores[fold_ids != f]
    test_scores <- scores[fold_ids == f]
    
    dtrain <- xgb.DMatrix(data = train_features, label = train_scores)
    model <- xgboost(data = dtrain, max.depth = 3, nrounds = 500, objective = "reg:squarederror", importance_type = "gain")
    
    dtest <- xgb.DMatrix(data = test_features)
    predicted_scores <- predict(model, newdata = dtest)
    
    cor(predicted_scores, test_scores)
  })
  
  return(cv_results)
}

# 运行交叉验证
cv_results_xgb <- cross_validation_xgb(reduced_features_matrix, scores_vector)

# 计算平均交叉验证准确性
mean_cv_accuracy_xgb <- mean(cv_results_xgb)
print(mean_cv_accuracy_xgb)

#################################################

# 可视化重要性
importance_df <- as.data.frame(importance_matrix)
colnames(importance_df) <- c("Feature", "Gain", "Cover", "Frequency")

# 热图
# 对数据进行对数变换
heatmap_data <- data.frame(Feature = factor(importance_df$Feature), Importance = importance_df$Gain)
heatmap_data$LogImportance <- log1p(heatmap_data$Importance)


# 绘制对数变换后的热图
ggplot(heatmap_data, aes(x = Feature, y = "Importance")) +
  geom_tile(aes(fill = LogImportance), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(title = "Feature Importance Heatmap (Log Transformed)", x = "Feature", y = "Importance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#################################################

# 折线图
ggplot(importance_df, aes(x = Feature, y = Gain, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Feature Importance Line Plot", x = "Feature", y = "Gain") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#################################################

# 交叉验证结果折线图
cv_results_df <- data.frame(Fold = 1:5, Accuracy = cv_results_xgb)

ggplot(cv_results_df, aes(x = Fold, y = Accuracy)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Cross-Validation Accuracy by Fold", x = "Fold", y = "Accuracy")

#################################################

library(GA)

# 定义特征提取函数
extract_features_from_matrices <- function(matrices) {
  design_vector <- as.vector(matrices$design_matrix)
  rep_vector <- as.vector(matrices$rep_matrix)
  
  combined_vector <- c(design_vector, rep_vector)
  return(combined_vector)
}

# 定义适应度函数
fitness <- function(chromosome) {
  # 根据染色体生成设计矩阵和重复矩阵
  design_matrix <- matrix(as.integer(chromosome[1:(nEnvs * nGenos)] > 0.5), nrow = nGenos, ncol = nEnvs)
  
  rep_part <- chromosome[(nEnvs * nGenos + 1):(2 * nEnvs * nGenos)]
  scaled_rep_values <- floor(rep_part * 4)  # 乘以4以包含0到3的范围
  scaled_rep_values[scaled_rep_values > 3] <- 3  # 保证不超过3
  rep_matrix <- matrix(scaled_rep_values, nrow = nGenos, ncol = nEnvs)
  
  # # 检查列和是否符合总和约束
  # col_sums <- colSums(design_matrix * rep_matrix)
  # if (any(col_sums != 1000)) {  # 假设总和约束为1000
  #   return(-Inf)  # 如果列和不符合要求，返回负无穷
  # }
  # 
  # # 确保设计矩阵每行和为2
  # if (any(rowSums(design_matrix) != 2)) {
  #   return(-Inf)  # 如果行和不为2，返回负无穷
  # }
  # 
  # # 确保重复矩阵的值在0到3之间
  # if (any(rep_matrix > 3)) {
  #   return(-Inf)  # 如果重复矩阵的值超出范围，返回负无穷
  # }
  
  # 使用设计矩阵和重复矩阵提取特征
  features_vector <- extract_features_from_matrices(list(design_matrix = design_matrix, rep_matrix = rep_matrix))
  
  # 标准化特征
  scaled_combined_features <- scale(matrix(features_vector, nrow = 1))
  colnames(scaled_combined_features) <- colnames(features_matrix)
  
  # 只保留重要特征
  reduced_combined_features <- scaled_combined_features[, selected_features, drop = FALSE]
  
  # 使用XGBoost模型进行预测
  dtest <- xgb.DMatrix(data = scaled_combined_features)
  predicted_score <- predict(xgb_model_reduced, newdata = dtest)
  
  return(predicted_score)  # 返回预测得分作为适应度值
}

# 遗传算法设置
ga <- ga(
  type = "real-valued",
  fitness = fitness,
  nBits = 2 * nEnvs * nGenos,
  popSize = 200,
  maxiter = 10,  # 增加迭代次数以增加优化效果
  pmutation = 0.05,
  lower = c(rep(0, nGenos*nEnvs), rep(0, nGenos*nEnvs)),
  upper = c(rep(1, nGenos*nEnvs), rep(3, nGenos*nEnvs)),
  suggestions = c(as.vector(results_all_configs[["20"]][[1]]$matrices$design_matrix),as.vector(results_all_configs[["20"]][[1]]$matrices$rep_matrix))
)

# 获取最佳解的染色体
best_chromosome <- ga@solution

# 转换最佳染色体为矩阵
design_matrix <- matrix(as.integer(best_chromosome[1:(nEnvs * nGenos)] > 0.5), nrow = nGenos, ncol = nEnvs)
rep_matrix <- matrix(floor(best_chromosome[(nEnvs * nGenos + 1):(2 * nEnvs * nGenos)] * 4), nrow = nGenos, ncol = nEnvs)

# 打印结果
print("Optimal Design Matrix:")
print(design_matrix)
print("Optimal Replication Matrix:")
print(rep_matrix)
print(ga@fitnessValue)
