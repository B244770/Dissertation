library(asreml)

# 模型1：环境作为固定效应，基因型ID作为随机效应
model1 <- asreml(fixed = y.Trait1 ~ env, random = ~ id, data = pheno_df)

# 模型2：环境和重复次数作为固定效应，基因型ID作为随机效应
model2 <- asreml(fixed = y.Trait1 ~ env + rep, random = ~ id, data = pheno_df)

# 模型3：环境、行和列作为固定效应，基因型ID作为随机效应
model3 <- asreml(fixed = y.Trait1 ~ env + row + col, random = ~ id, data = pheno_df)

# 模型4：环境、行和列作为固定效应，基因型ID和重复次数作为随机效应
model4 <- asreml(fixed = y.Trait1 ~ env + row + col, random = ~ id + rep, data = pheno_df)

#################################################

# 提取AIC和BIC值
aic_values <- c(AIC(model1), AIC(model2), AIC(model3), AIC(model4))
bic_values <- c(BIC(model1), BIC(model2), BIC(model3), BIC(model4))

# 打印AIC和BIC值
aic_values
bic_values

#################################################

cross_validation_asreml <- function(model_formula, data, folds = 5) {
  set.seed(123)
  n <- nrow(data)
  fold_ids <- sample(rep(1:folds, length.out = n))
  
  cv_results <- sapply(1:folds, function(f) {
    train_data <- data[fold_ids != f, ]
    test_data <- data[fold_ids == f, ]
    
    model <- asreml(fixed = as.formula(model_formula), data = train_data, na.method.X = "omit", na.method.Y = "omit")
    
    test_data$predicted_gv <- predict(model, classify = "id")$pvals$predicted.value
    
    cor(test_data$predicted_gv, test_data$gv)
  })
  
  return(mean(cv_results))
}

# 评估各个模型的交叉验证准确性
cv_accuracy_model1 <- cross_validation_asreml("y.Trait1 ~ env", pheno_df)
cv_accuracy_model2 <- cross_validation_asreml("y.Trait1 ~ env + rep", pheno_df)
cv_accuracy_model3 <- cross_validation_asreml("y.Trait1 ~ env + row + col", pheno_df)
cv_accuracy_model4 <- cross_validation_asreml("y.Trait1 ~ env + row + col", pheno_df)

# 打印交叉验证准确性
cv_accuracy_model1
cv_accuracy_model2
cv_accuracy_model3
cv_accuracy_model4

#################################################


#################################################


#################################################