library(xgboost)
library(caret)
library(ggplot2)
library(reshape2)

# function for extracting features
extract_features_from_matrices <- function(matrices) {
  # expand the matrices to vectors
  design_vector <- as.vector(matrices$design_matrix)
  rep_vector <- as.vector(matrices$rep_matrix)
  
  # combine flattened matrices to a vector
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
      
      # save score and feature matrix to a list
      results_list[[length(results_list) + 1]] <- list(Score = score, Features = features_vector)
    }
  }
  
  return(results_list)
}

# load data from scenario 2 when nEnvs = 20
nEnvs <- 20
results_list <- extract_data(results_all_configs["20"])

# # convert data format
# features_matrix <- do.call(rbind, lapply(results_list, function(x) x$Features))
# scores_vector <- sapply(results_list, function(x) x$Score)

# fill Features vector to max_length, supplemented with NA
max_length <- max(sapply(results_list, function(x) length(x$Features)))
features_matrix <- do.call(rbind, lapply(results_list, function(x) {
  c(x$Features, rep(NA, max_length - length(x$Features)))
}))

scores_vector <- sapply(results_list, function(x) x$Score)

# extract column names
colnames(features_matrix) <- paste0("V", 1:ncol(features_matrix))

# standardise features
scaled_features <- scale(features_matrix)

# restore colnames
colnames(scaled_features) <- colnames(features_matrix)

# evaluate feature importance
dtrain <- xgb.DMatrix(data = scaled_features, label = scores_vector)
xgb_model <- xgboost(data = dtrain, max.depth = 3, nrounds = 100, objective = "reg:squarederror", importance_type = "gain")

importance_matrix <- xgb.importance(model = xgb_model)
selected_features <- importance_matrix$Feature[importance_matrix$Gain > mean(importance_matrix$Gain)]

# select features by importance
reduced_features_matrix <- scaled_features[, selected_features, drop = FALSE]

#################################################

# fit xgboost
dtrain_reduced <- xgb.DMatrix(data = scaled_features, label = scores_vector)
xgb_model_reduced <- xgboost(data = dtrain_reduced, max.depth = 3, nrounds = 500, objective = "reg:squarederror", importance_type = "gain")
print(summary(xgb_model_reduced))

#################################################

# cross validation
cross_validation_xgb <- function(features, scores, folds = 5) {
  set.seed(123)
  n <- nrow(features)
  fold_ids <- sample(rep(1:folds, length.out = n))
  
  cv_results <- sapply(1:folds, function(f) {
    train_features <- features[fold_ids != f, , drop = FALSE]
    test_features <- features[fold_ids == f, , drop = FALSE]
    train_scores <- scores[fold_ids != f]
    test_scores <- scores[fold_ids == f]
    
    dtrain <- xgb.DMatrix(data = train_features, label = train_scores)
    model <- xgboost(data = dtrain, max.depth = 3, nrounds = 500, objective = "reg:squarederror", importance_type = "gain")
    
    dtest <- xgb.DMatrix(data = test_features)
    predicted_scores <- predict(model, newdata = dtest)
    
    cor(predicted_scores, test_scores)
  })
  
  return(cv_results)
}

# run cv
cv_results_xgb <- cross_validation_xgb(reduced_features_matrix, scores_vector)

# calculate mean accuracy
mean_cv_accuracy_xgb <- mean(cv_results_xgb)
print(mean_cv_accuracy_xgb)

#################################################

# visualise importance
importance_df <- as.data.frame(importance_matrix)
colnames(importance_df) <- c("Feature", "Gain", "Cover", "Frequency")

# heatmap
heatmap_data <- data.frame(Feature = factor(importance_df$Feature), Importance = importance_df$Gain)
heatmap_data$LogImportance <- log1p(heatmap_data$Importance)


# heatmap by log-transformed
ggplot(heatmap_data, aes(x = Feature, y = "Importance")) +
  geom_tile(aes(fill = LogImportance), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(title = "Feature Importance Heatmap (Log Transformed)", x = "Feature", y = "Importance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#################################################

# line plot
ggplot(importance_df, aes(x = Feature, y = Gain, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Feature Importance Line Plot", x = "Feature", y = "Gain") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#################################################

# cv result line plot
cv_results_df <- data.frame(Fold = 1:5, Accuracy = cv_results_xgb)

ggplot(cv_results_df, aes(x = Fold, y = Accuracy)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Cross-Validation Accuracy by Fold", x = "Fold", y = "Accuracy")

#################################################

library(GA)

# feature extraction
extract_features_from_matrices <- function(matrices) {
  design_vector <- as.vector(matrices$design_matrix)
  rep_vector <- as.vector(matrices$rep_matrix)
  
  combined_vector <- c(design_vector, rep_vector)
  return(combined_vector)
}

# fiteness function
fitness <- function(chromosome) {
  # extract matrices from chromosome
  design_matrix <- matrix(as.integer(chromosome[1:(nEnvs * nGenos)] > 0.5), nrow = nGenos, ncol = nEnvs)
  
  rep_part <- chromosome[(nEnvs * nGenos + 1):(2 * nEnvs * nGenos)]
  scaled_rep_values <- floor(rep_part * 4)  # multiply by 4 to include the range 0 to 3
  scaled_rep_values[scaled_rep_values > 3] <- 3  # not exceed 3
  rep_matrix <- matrix(scaled_rep_values, nrow = nGenos, ncol = nEnvs)

  # extract features
  features_vector <- extract_features_from_matrices(list(design_matrix = design_matrix, rep_matrix = rep_matrix))
  
  # standardise features
  scaled_combined_features <- scale(matrix(features_vector, nrow = 1))
  colnames(scaled_combined_features) <- colnames(features_matrix)
  
  # keep important features ONLY
  reduced_combined_features <- scaled_combined_features[, selected_features, drop = FALSE]
  
  # run xgboost
  dtest <- xgb.DMatrix(data = scaled_combined_features)
  predicted_score <- predict(xgb_model_reduced, newdata = dtest)
  
  return(predicted_score)
}

# genetic algorithm settings
ga <- ga(
  type = "real-valued",
  fitness = fitness,
  nBits = 2 * nEnvs * nGenos,
  popSize = 200,
  maxiter = 10,  # max iter
  pmutation = 0.05,
  lower = c(rep(0, nGenos*nEnvs), rep(0, nGenos*nEnvs)),
  upper = c(rep(1, nGenos*nEnvs), rep(3, nGenos*nEnvs)),
  suggestions = c(as.vector(results_all_configs[["20"]][[1]]$matrices$design_matrix),as.vector(results_all_configs[["20"]][[1]]$matrices$rep_matrix))
)

# best chromosome
best_chromosome <- ga@solution

# turn to matrices
design_matrix <- matrix(as.integer(best_chromosome[1:(nEnvs * nGenos)] > 0.5), nrow = nGenos, ncol = nEnvs)
rep_matrix <- matrix(floor(best_chromosome[(nEnvs * nGenos + 1):(2 * nEnvs * nGenos)] * 4), nrow = nGenos, ncol = nEnvs)

# print results
print("Optimal Design Matrix:")
print(design_matrix)
print("Optimal Replication Matrix:")
print(rep_matrix)
print(ga@fitnessValue)
