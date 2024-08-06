# 加载必要的包
library(e1071)
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

# 拟合SVM模型
svm_model <- svm(x = scaled_features, y = scores_vector, kernel = "radial")
print(summary(svm_model))

# 定义交叉验证函数
cross_validation_svm <- function(features, scores, folds = 5) {
  set.seed(123)
  n <- nrow(features)
  fold_ids <- sample(rep(1:folds, length.out = n))
  
  cv_results <- sapply(1:folds, function(f) {
    train_features <- features[fold_ids != f, , drop = FALSE]
    test_features <- features[fold_ids == f, , drop = FALSE]
    train_scores <- scores[fold_ids != f]
    test_scores <- scores[fold_ids == f]
    
    model <- svm(x = train_features, y = train_scores, kernel = "radial")
    
    predicted_scores <- predict(model, newdata = test_features)
    
    cor(predicted_scores, test_scores)
  })
  
  return(cv_results)
}

# 运行交叉验证
cv_results_svm <- cross_validation_svm(scaled_features, scores_vector)

# 计算平均交叉验证准确性
mean_cv_accuracy_svm <- mean(cv_results_svm)
print(mean_cv_accuracy_svm)

# 可视化支持向量特征
support_vectors <- scaled_features[svm_model$index, ]
support_vectors_df <- as.data.frame(support_vectors)
support_vectors_df$Index <- svm_model$index

# 绘制支持向量特征
ggplot(melt(support_vectors_df, id.vars = "Index"), aes(x = variable, y = value)) +
  geom_point(aes(color = factor(Index))) +
  theme_minimal() +
  labs(title = "Support Vectors Features", x = "Feature", y = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 交叉验证结果折线图
cv_results_df <- data.frame(Fold = 1:5, Accuracy = cv_results_svm)

ggplot(cv_results_df, aes(x = Fold, y = Accuracy)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Cross-Validation Accuracy by Fold", x = "Fold", y = "Accuracy")
