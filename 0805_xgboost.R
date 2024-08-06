# 加载必要的包
library(xgboost)
library(caret)
library(ggplot2)
library(reshape2)

# 提取特征向量的函数
extract_features_from_matrices <- function(matrices) {
  # 展平设计矩阵和重复矩阵
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
# results_all_configs <- your_loaded_data_here
results_list <- extract_data(results_all_configs)

# 转换数据格式
features_matrix <- do.call(rbind, lapply(results_list, function(x) x$Features))
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
dtrain_reduced <- xgb.DMatrix(data = reduced_features_matrix, label = scores_vector)
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

