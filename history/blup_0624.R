library(lme4)
library(ggplot2)
library(dplyr)

gv_df_factor <- gv_df

gv_df_factor$id <- as.factor(gv_df_factor$id)
gv_df_factor$env <- as.factor(gv_df_factor$env)


model_data <- full_join(pheno_df, gv_df_factor, by = c("id", "env"))
model_data <- full_join(model_data, error_ls$error.df, by = c("env", "block", "col", "row"))

# check data
head(model_data)

# get specific data set, here e.g. rep=1
# model_data_unique <- model_data[!duplicated(model_data), ]
model_data_unique <- subset(model_data, rep==1)

# genetics values as random effect, others as fixed
model_blup <- lmer(y.Trait1 ~ env + block + col + row + e.Trait1 + (1|id) + (1|gv_df$gv.Trait1), data = model_data_unique)

    # get blup value
    random_effects <- ranef(model)
    print(random_effects)
    
    
    predicted_values <- predict(model, re.form = NA)
    
    # create df for plot
    data_for_plot <- data.frame(Observed = model_data_unique$y.Trait1, Predicted = predicted_values)
    
    # ggplot
    p <- ggplot(data_for_plot, aes(x = Observed, y = Predicted)) +
      geom_point(alpha = 0.6) +  # 添加点，透明度为0.6
      geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed") +  # 添加 y=x 的虚线
      labs(x = "Observed Values", y = "Predicted Values", title = "Scatter Plot of Observed vs. Predicted Values") +
      theme_minimal()  # 使用简洁主题
    print(p)
    
    # get the minimum error
    data_for_plot$Error <- abs(model_data_unique$y.Trait1 - predicted_values)
    best_predictions <- model_data_unique[which.min(data_for_plot$Error), ]
    best_predictions


# add reaction gxe
model2 <- lmer(y.Trait1 ~ env * gv_df$gv.Trait1 + block + col + row + (1|id), data = model_data_unique)

    # 计算预测值
    predictions2 <- predict(model2, re.form = NULL)
    
    # 将预测值添加到数据框
    model_data_unique$Predicted2 <- predictions2
  
    ggplot(model_data_unique, aes(x = env, y = Predicted2, color = as.factor(gv_df$gv.Trait1))) +
      geom_point(alpha = 0.6) +
      geom_line(aes(group = gv_df$gv.Trait1)) +
      labs(title = "Predicted Trait Values across Environments by Genotype",
           x = "Environment",
           y = "Predicted Trait Value") +
      theme_minimal()


# spatial self-relativity
model3 <- lmer(y.Trait1 ~ env + block + (1|row) + (1|col) + (1|id) + (1|gv_df$gv.Trait1), data = model_data_unique)

    # 提取随机效应
    random_effects3_row <- ranef(model3)$row
    random_effects3_col <- ranef(model3)$col
  
    # 转换行号为数值类型
    random_effects3_row$row_number <- as.numeric(rownames(random_effects3_row))
    random_effects3_col$col_number <- as.numeric(rownames(random_effects3_col))
    
    # 绘制行和列的随机效应
    ggplot(as.data.frame(random_effects3_row), aes(x = random_effects3_row$row_number, y = random_effects3_row$`(Intercept)`)) +
      geom_bar(stat = "identity") +
      labs(title = "Random Effects for Rows", x = "Row", y = "Effect")
    
    ggplot(as.data.frame(random_effects3_col), aes(x = random_effects3_col$col_number, y = random_effects3_col$`(Intercept)`)) +
      geom_bar(stat = "identity") +
      labs(title = "Random Effects for Columns", x = "Column", y = "Effect")


# random effect among blocks
model4 <- lmer(y.Trait1 ~ env + (1|block) + (1|col) + (1|row) + (1|id) + (1|gv_df$gv.Trait1), data = model_data_unique)

    # 提取区块的随机效应
    random_effects4_block <- ranef(model4)$block
    
    # 绘制区块的随机效应
    ggplot(as.data.frame(random_effects4_block), aes(x = rownames(random_effects4_block), y = random_effects4_block$`(Intercept)`)) +
      geom_bar(stat = "identity") +
      labs(title = "Random Effects for Blocks", x = "Block", y = "Effect")

