# 加载必要的包
library(gbm)
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

# 训练初始的GBM模型以评估特征重要性
initial_gbm_model <- gbm(y = scores_vector, x = scaled_features, distribution = "gaussian", n.trees = 100, interaction.depth = 3)
importance_scores <- summary(initial_gbm_model, plot = FALSE)
selected_features <- importance_scores$var[importance_scores$rel.inf > mean(importance_scores$rel.inf)]

# 选择重要特征
reduced_features_matrix <- scaled_features[, selected_features, drop = FALSE]


#################################################

# 拟合GBM模型
gbm_model <- gbm(y = scores_vector, x = reduced_features_matrix, distribution = "gaussian", n.trees = 500, interaction.depth = 3, shrinkage = 0.01, cv.folds = 5)
print(summary(gbm_model))

#################################################

# 定义交叉验证函数
cross_validation_gbm <- function(features, scores, folds = 5) {
  set.seed(123)
  n <- nrow(features)
  fold_ids <- sample(rep(1:folds, length.out = n))
  
  cv_results <- sapply(1:folds, function(f) {
    train_features <- features[fold_ids != f, , drop = FALSE]
    test_features <- features[fold_ids == f, , drop = FALSE]
    train_scores <- scores[fold_ids != f]
    test_scores <- scores[fold_ids == f]
    
    model <- gbm(y = train_scores, x = train_features, distribution = "gaussian", n.trees = 500, interaction.depth = 3, shrinkage = 0.01)
    
    predicted_scores <- predict(model, newdata = test_features, n.trees = 500)
    
    cor(predicted_scores, test_scores)
  })
  
  return(cv_results)
}

# 运行交叉验证
cv_results_gbm <- cross_validation_gbm(reduced_features_matrix, scores_vector)

# 计算平均交叉验证准确性
mean_cv_accuracy_gbm <- mean(cv_results_gbm)
print(mean_cv_accuracy_gbm)

#################################################

# 可视化重要性
importance_df <- data.frame(Feature = importance_scores$var, Importance = importance_scores$rel.inf)

# 热图
heatmap_data <- data.frame(Feature = factor(importance_df$Feature), Importance = importance_df$Importance)
ggplot(heatmap_data, aes(x = Feature, y = Importance)) +
  geom_tile(aes(fill = Importance), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(title = "Feature Importance Heatmap", x = "Feature", y = "Importance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#################################################

# 折线图
ggplot(importance_df, aes(x = Feature, y = Importance, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Feature Importance Line Plot", x = "Feature", y = "Importance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#################################################

# 交叉验证结果折线图
cv_results_df <- data.frame(Fold = 1:5, Accuracy = cv_results_gbm)

ggplot(cv_results_df, aes(x = Fold, y = Accuracy)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Cross-Validation Accuracy by Fold", x = "Fold", y = "Accuracy")
