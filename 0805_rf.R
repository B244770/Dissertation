# 加载必要的包
library(randomForest)
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

# 标准化特征
scaled_features <- scale(features_matrix)

# 训练初始的随机森林模型以评估特征重要性
initial_rf_model <- randomForest(x = scaled_features, y = scores_vector, ntree = 500)
importance_scores <- importance(initial_rf_model)
selected_features <- names(importance_scores)[importance_scores > mean(importance_scores)]

# 选择重要特征
reduced_features_matrix <- scaled_features[, selected_features]

#################################################

# 拟合随机森林模型
rf_model <- randomForest(x = reduced_features_matrix, y = scores_vector, ntree = 500)
print(summary(rf_model))

#################################################

# 定义交叉验证函数
cross_validation_rf <- function(features, scores, folds = 5) {
  set.seed(123)
  n <- nrow(features)
  fold_ids <- sample(rep(1:folds, length.out = n))
  
  cv_results <- sapply(1:folds, function(f) {
    train_features <- features[fold_ids != f, ]
    test_features <- features[fold_ids == f, ]
    train_scores <- scores[fold_ids != f]
    test_scores <- scores[fold_ids == f]
    
    model <- randomForest(x = train_features, y = train_scores, ntree = 500)
    
    predicted_scores <- predict(model, newdata = test_features)
    
    cor(predicted_scores, test_scores)
  })
  
  return(cv_results)
}

# 运行交叉验证
cv_results_rf <- cross_validation_rf(reduced_features_matrix, scores_vector)

# 计算平均交叉验证准确性
mean_cv_accuracy_rf <- mean(cv_results_rf)
print(mean_cv_accuracy_rf)

#################################################

# 可视化重要性
importance_df <- as.data.frame(importance(rf_model))
importance_df$Feature <- rownames(importance_df)

# 热图
heatmap_data <- data.frame(Feature = factor(rownames(importance_df)), Importance = importance_df$MeanDecreaseGini)
ggplot(heatmap_data, aes(x = Feature, y = Importance)) +
  geom_tile(aes(fill = Importance), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(title = "Feature Importance Heatmap", x = "Feature", y = "Importance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#################################################

# 折线图
ggplot(importance_df, aes(x = Feature, y = MeanDecreaseGini, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Feature Importance Line Plot", x = "Feature", y = "Mean Decrease Gini") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#################################################

# 交叉验证结果折线图
cv_results_df <- data.frame(Fold = 1:5, Accuracy = cv_results_rf)

ggplot(cv_results_df, aes(x = Fold, y = Accuracy)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Cross-Validation Accuracy by Fold", x = "Fold", y = "Accuracy")
