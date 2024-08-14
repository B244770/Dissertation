library(asreml)
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
# model_data_unique <- subset(model_data, rep==1)

# Fit the model with genetic values as random effects, others as fixed
model_blup <- asreml(fixed = y.Trait1 ~ env,
                     random = ~ id + env:block,
                     residual = ~dsum(units|env),
                     data = model_data_unique)

model_gblup <- asreml(fixed = y.Trait1 ~ env,
                     random = ~ vm(id, G) + env:block,
                     residual = ~dsum(units|env),
                     data = model_data_unique)

summary(model_blup)$varcom

# get BLUP values for random effects
blup_values <- predict(model_blup, classify = "id")
print(blup_values$pvals)

# predicted values
predicted_values <- predict(model_blup, classify = "id", only = "fixed")$pvals$predicted.value

# create df for plot
data_for_plot <- data.frame(Observed = model_data_unique$y.Trait1, Predicted = predicted_values)

# ggplot
p <- ggplot(data_for_plot, aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed") +
  labs(x = "Observed Values", y = "Predicted Values", title = "Scatter Plot of Observed vs. Predicted Values") +
  theme_minimal()
print(p)

# get the minimum error
data_for_plot$Error <- abs(model_data_unique$y.Trait1 - predicted_values)
best_predictions <- model_data_unique[which.min(data_for_plot$Error), ]
print(best_predictions)

# add reaction gxe
model2 <- asreml(fixed = y.Trait1 ~ env * gv.Trait1 + block + col + row,
                 random = ~ id,
                 data = model_data_unique)

# predicted values
predictions2 <- predict(model2, classify = "env:gv.Trait1")$pvals$predicted.value

# add predictions to data
model_data_unique$Predicted2 <- predictions2

# ggplot for interaction
ggplot(model_data_unique, aes(x = env, y = Predicted2, color = as.factor(gv.Trait1))) +
  geom_point(alpha = 0.6) +
  geom_line(aes(group = gv.Trait1)) +
  labs(title = "Predicted Trait Values across Environments by Genotype",
       x = "Environment",
       y = "Predicted Trait Value") +
  theme_minimal()

# spatial self-relativity
model3 <- asreml(fixed = y.Trait1 ~ env + block,
                 random = ~ row + col + id + gv.Trait1,
                 data = model_data_unique)

# extract random effects
random_effects3_row <- predict(model3, classify = "row")$pvals
random_effects3_col <- predict(model3, classify = "col")$pvals

# plot random effects for rows
ggplot(as.data.frame(random_effects3_row), aes(x = row, y = predicted.value)) +
  geom_bar(stat = "identity") +
  labs(title = "Random Effects for Rows", x = "Row", y = "Effect")

# plot random effects for columns
ggplot(as.data.frame(random_effects3_col), aes(x = col, y = predicted.value)) +
  geom_bar(stat = "identity") +
  labs(title = "Random Effects for Columns", x = "Column", y = "Effect")

# random effect among blocks
model4 <- asreml(fixed = y.Trait1 ~ env,
                 random = ~ block + col + row + id + gv.Trait1,
                 data = model_data_unique)

# extract random effects for blocks
random_effects4_block <- predict(model4, classify = "block")$pvals

# plot random effects for blocks
ggplot(as.data.frame(random_effects4_block), aes(x = block, y = predicted.value)) +
  geom_bar(stat = "identity") +
  labs(title = "Random Effects for Blocks", x = "Block", y = "Effect")
