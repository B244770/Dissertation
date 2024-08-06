# 加载必要的包
library(keras3)
library(caret)
library(ggplot2)
library(reshape2)

# 提取特征向量的函数
extract_features_from_matrices <- function(matrices) {
  design_vector <- as.vector(matrices$design_matrix)
  rep_vector <- as.vector(matrices$rep_matrix)
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

# 定义生成器模型
generator <- keras_model_sequential() %>%
  layer_dense(units = 128, activation = 'relu', input_shape = 100) %>%
  layer_dense(units = ncol(x_train), activation = 'linear')

# 定义判别器模型
discriminator <- keras_model_sequential() %>%
  layer_dense(units = 128, activation = 'relu', input_shape = ncol(x_train)) %>%
  layer_dense(units = 1, activation = 'sigmoid')

# 编译判别器模型
discriminator %>% compile(
  loss = 'binary_crossentropy',
  optimizer = optimizer_adam(),
  metrics = c('accuracy')
)

# 定义GAN模型
gan_input <- layer_input(shape = 100)
gan_output <- discriminator(generator(gan_input))
gan <- keras_model(gan_input, gan_output)

# 编译GAN模型
gan %>% compile(
  loss = 'binary_crossentropy',
  optimizer = optimizer_adam()
)

#################################################

# 训练GAN模型
train_gan <- function(generator, discriminator, gan, x_train, epochs = 10000, batch_size = 32) {
  for (epoch in 1:epochs) {
    noise <- matrix(runif(batch_size * 100), nrow = batch_size, ncol = 100)
    generated_data <- generator %>% predict(noise)
    real_data <- x_train[sample(nrow(x_train), batch_size), ]
    combined_data <- rbind(generated_data, real_data)
    labels <- c(rep(0, batch_size), rep(1, batch_size))
    d_loss <- discriminator %>% train_on_batch(combined_data, labels)
    noise <- matrix(runif(batch_size * 100), nrow = batch_size, ncol = 100)
    misleading_labels <- rep(1, batch_size)
    g_loss <- gan %>% train_on_batch(noise, misleading_labels)
    if (epoch %% 1000 == 0) {
      cat("Epoch:", epoch, "Discriminator Loss:", d_loss, "Generator Loss:", g_loss, "\n")
    }
  }
}

# 运行训练
train_gan(generator, discriminator, gan, x_train, epochs = 10000, batch_size = 32)

#################################################

# 生成新的数据
noise <- matrix(runif(nrow(x_train) * 100), nrow = nrow(x_train), ncol = 100)
generated_data <- generator %>% predict(noise)

# 可视化生成的数据与真实数据的比较
generated_data_df <- as.data.frame(generated_data)
real_data_df <- as.data.frame(x_train)
generated_data_df$type <- 'Generated'
real_data_df$type <- 'Real'
combined_data_df <- rbind(generated_data_df, real_data_df)

ggplot(melt(combined_data_df, id.vars = "type"), aes(x = variable, y = value, color = type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Generated vs Real Data", x = "Feature", y = "Value")

#################################################

# 定义交叉验证函数
cross_validation_gan <- function(generator, discriminator, gan, x, y, folds = 5) {
  set.seed(123)
  n <- nrow(x)
  fold_ids <- sample(rep(1:folds, length.out = n))
  
  cv_results <- sapply(1:folds, function(f) {
    x_train <- x[fold_ids != f, ]
    y_train <- y[fold_ids != f]
    x_test <- x[fold_ids == f, ]
    y_test <- y[fold_ids == f]
    
    train_gan(generator, discriminator, gan, x_train, epochs = 1000, batch_size = 32)
    
    noise <- matrix(runif(nrow(x_test) * 100), nrow = nrow(x_test), ncol = 100)
    generated_data <- generator %>% predict(noise)
    
    cor(generated_data, x_test)
  })
  
  return(cv_results)
}

# 运行交叉验证
cv_results_gan <- cross_validation_gan(generator, discriminator, gan, x_train, y_train)

# 计算平均交叉验证准确性
mean_cv_accuracy_gan <- mean(cv_results_gan)
print(mean_cv_accuracy_gan)

#################################################

# 可视化交叉验证结果
cv_results_df <- data.frame(Fold = 1:5, Accuracy = cv_results_gan)

ggplot(cv_results_df, aes(x = Fold, y = Accuracy)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Cross-Validation Accuracy by Fold", x = "Fold", y = "Accuracy")
