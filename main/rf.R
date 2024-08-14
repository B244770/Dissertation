library(randomForest)
library(caret)
library(ggplot2)
library(reshape2)

# extract features
extract_features_from_matrices <- function(matrices) {
  design_vector <- as.vector(matrices$design_matrix)
  rep_vector <- as.vector(matrices$rep_matrix)
  combined_vector <- c(design_vector, rep_vector)
  return(combined_vector)
}

# extract data from scenarios
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
# load data from scenario 2 when nEnvs = 20
nEnvs <- 20
results_list <- extract_data(results_all_configs["20"])

# convert data formats
features_matrix <- do.call(rbind, lapply(results_list, function(x) x$Features))
scores_vector <- sapply(results_list, function(x) x$Score)

# extract column names
colnames(features_matrix) <- paste0("V", 1:ncol(features_matrix))

# standardise features
scaled_features <- scale(features_matrix)

# restore colnames
colnames(scaled_features) <- colnames(features_matrix)

# feature importance
initial_rf_model <- randomForest(x = scaled_features, y = scores_vector, ntree = 500)
importance_scores <- importance(initial_rf_model)
selected_features <- names(importance_scores)[importance_scores > mean(importance_scores)]

# select features
reduced_features_matrix <- scaled_features[, selected_features]

#################################################

# fitting
rf_model <- randomForest(x = scaled_features, y = scores_vector, ntree = 500)
print(summary(rf_model))

#################################################

# define cv
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

# run
cv_results_rf <- cross_validation_rf(scaled_features, scores_vector)

# mean accuracy
mean_cv_accuracy_rf <- mean(cv_results_rf)
print(mean_cv_accuracy_rf)

#################################################

# visualise imporatance
importance_df <- as.data.frame(importance(rf_model))
importance_df$Feature <- rownames(importance_df)

heatmap_data <- data.frame(Feature = importance_df$Feature, Importance = importance_df$IncNodePurity)
heatmap_data$LogImportance <- 0-(log(importance_df$IncNodePurity + 1e-10))
ggplot(heatmap_data, aes(x = Feature, y = LogImportance)) +
  geom_tile(aes(fill = LogImportance), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(title = "Feature Importance Heatmap", x = "Feature", y = "Importance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#################################################

# line plot
ggplot(importance_df, aes(x = Feature, y = MeanDecreaseGini, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Feature Importance Line Plot", x = "Feature", y = "Mean Decrease Gini") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#################################################

# results cv
cv_results_df <- data.frame(Fold = 1:5, Accuracy = cv_results_rf)

ggplot(cv_results_df, aes(x = Fold, y = Accuracy)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Cross-Validation Accuracy by Fold", x = "Fold", y = "Accuracy")
