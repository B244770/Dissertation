# 加载必要的包
library(keras3)
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

# 转换为矩阵格式
x_train <- as.matrix(scaled_features)
y_train <- scores_vector

#################################################

# 定义神经网络模型
model <- keras_model_sequential()
model %>%
  layer_dense(units = 64, activation = 'relu', input_shape = ncol(x_train)) %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = 32, activation = 'relu') %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = 1)

# 编译模型
model %>% compile(
  loss = 'mse',
  optimizer = optimizer_adam(),
  metrics = c('mae')
)

# 拟合模型
history <- model %>% fit(
  x_train, y_train,
  epochs = 50,
  batch_size = 32,
  validation_split = 0.2
)

#################################################

# 定义交叉验证函数
cross_validation_nn <- function(x, y, folds = 5) {
  set.seed(123)
  n <- nrow(x)
  fold_ids <- sample(rep(1:folds, length.out = n))
  
  cv_results <- sapply(1:folds, function(f) {
    x_train <- x[fold_ids != f, ]
    y_train <- y[fold_ids != f]
    x_test <- x[fold_ids == f, ]
    y_test <- y[fold_ids == f]
    
    model <- keras_model_sequential()
    model %>%
      layer_dense(units = 64, activation = 'relu', input_shape = ncol(x)) %>%
      layer_dropout(rate = 0.5) %>%
      layer_dense(units = 32, activation = 'relu') %>%
      layer_dropout(rate = 0.5) %>%
      layer_dense(units = 1)
    
    model %>% compile(
      loss = 'mse',
      optimizer = optimizer_adam(),
      metrics = c('mae')
    )
    
    model %>% fit(
      x_train, y_train,
      epochs = 50,
      batch_size = 32,
      verbose = 0
    )
    
    y_pred <- model %>% predict(x_test)
    cor(y_pred, y_test)
  })
  
  return(cv_results)
}

# 运行交叉验证
cv_results_nn <- cross_validation_nn(x_train, y_train)

# 计算平均交叉验证准确性
mean_cv_accuracy_nn <- mean(cv_results_nn)
print(mean_cv_accuracy_nn)

#################################################

# 绘制训练和验证损失
plot(history) +
  ggtitle("Training and Validation Loss") +
  xlab("Epochs") +
  ylab("Loss")

# 交叉验证结果折线图
cv_results_df <- data.frame(Fold = 1:5, Accuracy = cv_results_nn)

ggplot(cv_results_df, aes(x = Fold, y = Accuracy)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Cross-Validation Accuracy by Fold", x = "Fold", y = "Accuracy")
