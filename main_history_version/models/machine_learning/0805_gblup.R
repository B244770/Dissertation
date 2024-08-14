library(dplyr)
library(asreml)

# 提取数据并转换为适合模型使用的格式
extract_data_for_asreml <- function(results_all_configs) {
  data_list <- list()
  for (env_config in names(results_all_configs)) {
    config_results <- results_all_configs[[env_config]]
    for (run_index in seq_along(config_results)) {
      result <- config_results[[run_index]]
      design_matrix <- result$matrices$design_matrix
      rep_matrix <- result$matrices$rep_matrix
      score <- result$score
      
      design_df <- as.data.frame(design_matrix)
      colnames(design_df) <- paste0("env", 1:ncol(design_df))
      design_df$genotype <- 1:nrow(design_df)
      
      rep_df <- as.data.frame(rep_matrix)
      colnames(rep_df) <- paste0("rep_env", 1:ncol(rep_df))
      rep_df$genotype <- 1:nrow(rep_df)
      
      combined_df <- cbind(design_df, rep_df)
      combined_df$score <- score
      combined_df$env_config <- env_config
      combined_df$run <- run_index
      
      data_list[[length(data_list) + 1]] <- combined_df
    }
  }
  full_data <- bind_rows(data_list)
  return(full_data)
}

# 提取并转换数据
full_data <- extract_data_for_asreml(results_all_configs)

#################################################

# 模型1：环境数量作为固定效应，基因型ID作为随机效应
model1 <- asreml(fixed = score ~ env_config, random = ~ genotype, data = full_data)

# 模型2：环境数量和重复次数作为固定效应，基因型ID作为随机效应
model2 <- asreml(fixed = score ~ env_config + rep_env1 + rep_env2 + rep_env3 + rep_env4 + rep_env5, random = ~ genotype, data = full_data)

# 模型3：环境数量、行和列作为固定效应，基因型ID作为随机效应
model3 <- asreml(fixed = score ~ env_config + row + col, random = ~ genotype, data = full_data)

# 模型4：环境数量、行和列作为固定效应，基因型ID和重复次数作为随机效应
model4 <- asreml(fixed = score ~ env_config + row + col, random = ~ genotype + rep_env1 + rep_env2 + rep_env3 + rep_env4 + rep_env5, data = full_data)

# 提取 AIC 和 BIC 值
aic_values <- c(AIC(model1), AIC(model2), AIC(model3), AIC(model4))
bic_values <- c(BIC(model1), BIC(model2), BIC(model3), BIC(model4))

# 打印 AIC 和 BIC 值
print(aic_values)
print(bic_values)

#################################################

cross_validation_asreml <- function(model_formula, data, folds = 5) {
  set.seed(123)
  n <- nrow(data)
  fold_ids <- sample(rep(1:folds, length.out = n))
  
  cv_results <- sapply(1:folds, function(f) {
    train_data <- data[fold_ids != f, ]
    test_data <- data[fold_ids == f, ]
    
    model <- asreml(fixed = as.formula(model_formula), random = ~ genotype + rep_env1 + rep_env2 + rep_env3 + rep_env4 + rep_env5, data = train_data, na.method.X = "omit", na.method.Y = "omit")
    
    predicted_values <- predict(model, classify = "genotype")$pvals$predicted.value
    cor(predicted_values, test_data$score)
  })
  
  return(cv_results)
}

# 评估各个模型的交叉验证准确性
cv_results_model1 <- cross_validation_asreml("score ~ env_config", full_data)
cv_results_model2 <- cross_validation_asreml("score ~ env_config + rep_env1 + rep_env2 + rep_env3 + rep_env4 + rep_env5", full_data)
cv_results_model3 <- cross_validation_asreml("score ~ env_config + row + col", full_data)
cv_results_model4 <- cross_validation_asreml("score ~ env_config + row + col", full_data)

# 计算平均交叉验证准确性
mean_cv_accuracy_model1 <- mean(cv_results_model1)
mean_cv_accuracy_model2 <- mean(cv_results_model2)
mean_cv_accuracy_model3 <- mean(cv_results_model3)
mean_cv_accuracy_model4 <- mean(cv_results_model4)

print(mean_cv_accuracy_model1)
print(mean_cv_accuracy_model2)
print(mean_cv_accuracy_model3)
print(mean_cv_accuracy_model4)

#################################################

# 创建AIC和BIC值的数据框
aic_bic_df <- data.frame(
  Model = factor(c("Model 1", "Model 2", "Model 3", "Model 4")),
  AIC = aic_values,
  BIC = bic_values
)

# 绘制AIC和BIC值的柱状图
ggplot(aic_bic_df, aes(x = Model)) +
  geom_bar(aes(y = AIC, fill = "AIC"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = BIC, fill = "BIC"), stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("AIC" = "blue", "BIC" = "red")) +
  theme_minimal() +
  labs(title = "AIC and BIC Values for Different Models", x = "Model", y = "Value")

#################################################

# 创建交叉验证结果的数据框
cv_results_df <- data.frame(
  Fold = rep(1:5, times = 4),
  Model = factor(rep(c("Model 1", "Model 2", "Model 3", "Model 4"), each = 5)),
  Accuracy = c(cv_results_model1, cv_results_model2, cv_results_model3, cv_results_model4)
)

# 绘制交叉验证结果的箱线图
ggplot(cv_results_df, aes(x = Model, y = Accuracy)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, aes(color = Model)) +
  theme_minimal() +
  labs(title = "Cross-Validation Accuracy by Model", x = "Model", y = "Accuracy") +
  theme(legend.position = "none")
