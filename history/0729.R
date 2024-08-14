library(FieldSimR)
library(AlphaSimR)
library(ggplot2)
library(dplyr)

# Clear the work environment
rm(list = ls())

# # Define number of envs and genos
# nEnvs <- 20
# nGenos <- 1000
# nBlocks <- 1
# total_per_env <- 1000




# overload of make_phenotypes
make_phenotypes <- function(gv.df, error.df, design_matrix, rep_matrix, randomise = FALSE, return.effects = FALSE) {
  if (!is.data.frame(gv.df)) {
    stop("'gv.df' must be a data frame")
  }
  if (!is.data.frame(error.df)) {
    stop("'error.df' must be a data frame")
  }
  
  # 确保包含必需的列
  if (any(!c("env", "rep", "id") %in% colnames(gv.df))) {
    stop("'gv.df' must contain the columns 'env', 'rep', 'id', and the genetic values for each trait")
  }
  if (any(!c("env", "block", "col", "row") %in% colnames(error.df))) {
    stop("'error.df' must contain the columns 'env', 'block', 'col', 'row', and the plot errors for each trait")
  }
  
  # 根据设计矩阵和重复次数矩阵重新构建 gv.df 和 error.df
  # 这里假设 error.df 是已经按照设计矩阵和重复次数矩阵预先准备好的
  # 合并 gv.df 和 error.df 数据
  gv.df$env <- factor(as.numeric(as.character(gv.df$env)))
  gv.df$rep <- factor(as.numeric(as.character(gv.df$rep)))
  gv.df$id <- factor(as.numeric(as.character(gv.df$id)))
  
  error.df$env <- factor(as.numeric(as.character(error.df$env)))
  error.df$block <- factor(as.numeric(as.character(error.df$block)))
  error.df$col <- factor(as.numeric(as.character(error.df$col)))
  error.df$row <- factor(as.numeric(as.character(error.df$row)))
  
  # 排序并唯一化
  gv.df <- gv.df[order(gv.df$env, gv.df$rep), ]
  gv.df <- unique(gv.df)
  error.df <- error.df[order(error.df$env, error.df$block), ]
  error.df <- unique(error.df)
  
  # check rows
  if (nrow(gv.df) != nrow(error.df)) 
    stop("number of rows in 'gv.df' and 'error.df' must match")
  
  # randomise
  if (randomise) {
    gv.df$ord <- sample(nrow(gv.df))
    gv.df <- gv.df[order(gv.df$ord), ]
    gv.df$ord <- NULL # 去除序号列
  }
  
  # generate pheno data
  # ntraits <- ncol(gv.df) - 3  # 减去 'env', 'rep', 'id'
  ntraits <- 1
  y <- error.df[, !(colnames(error.df) %in% c("env", "block", "col", "row"))] + 
    gv.df[, !(colnames(gv.df) %in% c("env", "rep", "id"))]
  pheno_df <- data.frame(env = gv.df$env, block = error.df$block, col = error.df$col, row = error.df$row, id = gv.df$id, y)
  
  # rename pheno column
  colnames(pheno_df)[6:(5+ntraits)] <- paste0("y.Trait", 1:ntraits)
  
  # 如果需要返回基因效应和环境误差
  if (return.effects) {
    effects_df <- list()
    for (i in 1:ntraits) {
      trait_name <- paste0("Trait", i)
      effects_df[[trait_name]] <- data.frame(
        env = gv.df$env,
        id = gv.df$id,
        gen_effect = gv.df[[trait_name]],
        env_error = error.df[[paste0("error.Trait", i)]]
      )
    }
    return(list(pheno_df = pheno_df, effects_df = effects_df))
  }
  
  return(pheno_df)
}

simulate_experiment <- function(design_matrix, rep_matrix) {
  # Population parameters
  nFounders <- 20
  nChr <- 21
  nQTN <- 300
  nSNP <- 300
  nSegSites <- (nQTN + nSNP)
  
  mu <- rep(4, nEnvs) # environment/trial means, in tonnes/ha
  sigma2 <- rep(0.2, nEnvs) # environment/trial genetic variances
  H2 <- 0.3 # plot level heritability 
  sigma2e <- sigma2 / H2 - sigma2
  
  # population simulation
  # generate a random positive-definite matrix
  # set.seed(nBlocks)  # seeds
  random_matrix <- matrix(rnorm(nEnvs^2), nrow = nEnvs, ncol = nEnvs)
  cov_matrix <- abs(random_matrix %*% t(random_matrix))  # create covariance matrix
  # eigen(cov_matrix)$values # add while loop to run again if negative eigenvalues
  
  # 将协方差矩阵转换为相关矩阵
  diag_inv_sqrt <- diag(1 / sqrt(diag(cov_matrix)), nrow = nEnvs)
  cor_matrix <- diag_inv_sqrt %*% cov_matrix %*% diag_inv_sqrt
  # mean(cor_matrix[upper.tri(cor_matrix)])
  
  # create founder group
  founders <- quickHaplo(nInd = nFounders, nChr = nChr, segSites = nSegSites)
  
  # simulation params
  SP <- SimParam$new(founders)
  SP$addTraitA(nQtlPerChr = nQTN,
               mean = mu, 
               var = sigma2,
               corA = cor_matrix)
  founders <- newPop(founders)
  
  SP$addSnpChip(nSNP)
  
  # generate f1 and DH groups
  f1s <- randCross(pop = founders, nCrosses = 40, nProgeny = 25)
  # M <- pullSnpGeno(f1s)
  # unique(c(M))
  # dim(M)
  DHs <- makeDH(pop = f1s, nDH = 1)
  # cov2cor(var(DHs@gv))
  300*21
  M <- pullSnpGeno(DHs)
  unique(c(M))
  dim(M)
  M <- scale(M, center = T)
  G <- M %*% t(M)
  G <- G/mean(diag(G)) + diag(0.0000001, nrow = nrow(G))
  mean(diag(G))
  dim(G)
  # cols and rows 
  # # nrows * ncols = nGenos * nEnvs = total plotS
  # nplots <- nGenos*nBlocks
  nrows <- 50
  ncols <- total_per_env/nrows
  
  # generate field trial error
  error_ls <- field_trial_error(
    ntraits = 1,
    nenvs = nEnvs,
    nblocks = 2, # using envs to simulate blocks
    block.dir = "col", # block layouts ## !!bug here col is not working
    ncols = ncols,
    nrows = nrows,
    varR = sigma2e,
    spatial.model = "AR1",
    col.cor = 0.5,
    row.cor = 0.7,
    prop.spatial = 0.5,
    return.effects = TRUE
  )
  
  # get gv
  dhs_gvs_df <- data.frame(id = DHs@id)
  
  # DHs@gv 是一个矩阵或数据框，其中列名为 Trait1, Trait2, ..., Traitn
  # 将这些性状值合并到 dhs_gvs_df 数据框中
  dhs_gvs_df <- cbind(dhs_gvs_df, DHs@gv)
  
  
  # sampling and parsing as interger
  selected_gvs <- sample_n(dhs_gvs_df, nGenos)
  selected_gvs$id <- as.integer(selected_gvs$id)
  
  
  # generate gv data frame using sampling
  # init gv_df
  gv_df <- data.frame(env = integer(), rep = integer(), id = integer(), gv = numeric())
  
  # go through each env and genotype
  for (env in 1:ncol(design_matrix)) {
    for (geno in 1:nrow(design_matrix)) {
      if (design_matrix[geno, env] == 1) {
        # this geno IS planted in this env
        num_reps <- rep_matrix[geno, env]  # get rep
        geno_id <- selected_gvs$id[geno]   # get id
        trait_name <- paste("Trait", env, sep = "")  # create Trait n
        geno_value <- selected_gvs[[trait_name]][geno]  # get gv
        
        # 为每个重复生成一行数据
        for (rep in 1:num_reps) {
          gv_df <- rbind(gv_df, data.frame(env = env, rep = rep, id = geno_id, gv = geno_value))
        }
      }
    }
  }
  
  # create empty df to fill
  rows_to_add <- total_per_env * nEnvs - nrow(gv_df)
  if (rows_to_add > 0) {
    empty_rows <- data.frame(env = integer(rows_to_add), 
                             rep = integer(rows_to_add), 
                             id = integer(rows_to_add), 
                             gv = numeric(rows_to_add))
    # add 0
    empty_rows[] <- 0
    
    # add empty rows to gv_df
    gv_df <- rbind(gv_df, empty_rows)
  }
  # # check
  # print(head(gv_df))
  
  
  # generate phenotype data frame
  # View(gv_df)
  # View(error_ls$error.df)
  pheno_df <- make_phenotypes(
    gv.df = gv_df,
    error.df = error_ls$error.df,
    randomise = TRUE
  )
  
  # calculate pearson's coefficient
  result <- cor(
    with(pheno_df, tapply(y.Trait1, env, mean)),
    with(gv_df, tapply(gv, env, mean))
  )
  return(result)
}

# monte carlo

generate_matrices <- function(nGenos, nEnvs, total_per_env) {
  # init design matrix 0
  design_matrix <- matrix(0, nrow = nGenos, ncol = nEnvs, 
                          dimnames = list(paste("Geno", 1:nGenos, sep = ""), 
                                          paste("Env", 1:nEnvs, sep = "")))
  
  # init rep matrix 0
  rep_matrix <- matrix(0, nrow = nGenos, ncol = nEnvs, 
                       dimnames = list(paste("Geno", 1:nGenos, sep = ""), 
                                       paste("Env", 1:nEnvs, sep = "")))
  
  # geno and rep allocation, total guaranteed
  for (env in 1:nEnvs) {
    remaining_plants <- total_per_env
    while (remaining_plants > 0) {
      # randomly choose a genotype
      selected_genotype <- sample(nGenos, 1)
      # check current reps, make sure less than 3
      current_reps <- rep_matrix[selected_genotype, env]
      if (current_reps < 3) {
        max_reps <- min(remaining_plants, sample(1:3, 1), 3 - current_reps)  # 添加限制条件
        if (design_matrix[selected_genotype, env] == 0) {
          design_matrix[selected_genotype, env] <- 1
        }
        rep_matrix[selected_genotype, env] <- current_reps + max_reps
        remaining_plants <- remaining_plants - max_reps
      }
    }
  }
  return(list(design_matrix = design_matrix, rep_matrix = rep_matrix))
}

# evaluation
# function for matrix generation and evaluation
generate_and_evaluate <- function(nGenos, nEnvs, total_per_env) {
  matrices <- generate_matrices(nGenos, nEnvs, total_per_env)
  if (is.null(matrices)) {
    return(NULL)  # fail generation return NULL
  }
  score <- simulate_experiment(matrices$design_matrix, matrices$rep_matrix)
  # return(list(score = score, matrices = matrices))
  return(score)
}


# main

# Define configurations
env_configs <- c(4, 5, 10, 20, 40, 50)
nRep <- 10
total_per_env <- 1000
results_all_configs <- list()  # This will store all results

# Loop over each configuration
for (nEnvs in env_configs) {
  nGenos <- 1000
  nBlocks <- 1
  
  # Population parameters
  nFounders <- 20
  nChr <- 21
  nQTN <- 300
  nSNP <- 300
  nSegSites <- (nQTN + nSNP)
  
  mu <- rep(4, nEnvs) # environment/trial means, in tonnes/ha
  sigma2 <- rep(0.2, nEnvs) # environment/trial genetic variances
  H2 <- 0.3 # plot level heritability 
  sigma2e <- sigma2 / H2 - sigma2
  
  # Initialize the storage for results of this config
  results_for_current_config <- vector("list", nRep)
  
  # Simulation and evaluation per configuration
  for (i in 1:nRep) {
    # Random matrix for GxE interaction
    random_matrix <- matrix(rnorm(nEnvs^2), nrow = nEnvs, ncol = nEnvs)
    cov_matrix <- abs(random_matrix %*% t(random_matrix))
    diag_inv_sqrt <- diag(1 / sqrt(diag(cov_matrix)), nrow = nEnvs)
    cor_matrix <- diag_inv_sqrt %*% cov_matrix %*% diag_inv_sqrt
    
    # Founder and trait simulation
    founders <- quickHaplo(nInd = nFounders, nChr = nChr, segSites = nSegSites)
    SP <- SimParam$new(founders)
    SP$addTraitA(nQtlPerChr = nQTN, mean = mu, var = sigma2, corA = cor_matrix)
    founders <- newPop(founders)
    SP$addSnpChip(nSNP)
    f1s <- randCross(pop = founders, nCrosses = 40, nProgeny = 25)
    DHs <- makeDH(pop = f1s, nDH = 1)
    
    # Generate matrices and simulate experiments
    matrices <- generate_matrices(nGenos, nEnvs, total_per_env)
    simulation_result <- simulate_experiment(matrices$design_matrix, matrices$rep_matrix)
    
    # Store result
    results_for_current_config[[i]] <- list(score = simulation_result, matrices = matrices)
  }
  
  # Store results for current nEnvs
  results_all_configs[[as.character(nEnvs)]] <- results_for_current_config
}

# Initialize an empty data frame to store the results
results_df <- data.frame(Scenario = integer(), Accuracy = numeric())

# Loop through the results for all configurations
for (env_config in names(results_all_configs)) {
  # Extract scores for the current configuration
  current_results <- results_all_configs[[env_config]]
  current_scores <- sapply(current_results, function(x) x$score)
  
  # Create a temporary data frame for the current scenario
  temp_df <- data.frame(
    Scenario = factor(rep(env_config, length(current_scores))),
    Accuracy = current_scores
  )
  
  # Combine with the main results data frame
  results_df <- rbind(results_df, temp_df)
}

# Convert Scenario to a factor with more meaningful labels if necessary
results_df$Scenario <- factor(results_df$Scenario, labels = c("4 Envs", "5 Envs", "10 Envs", "20 Envs", "40 Envs", "50 Envs"))
# results_df$Scenario <- factor(results_df$Scenario, labels = c("4 Envs", "5 Envs", "10 Envs", "20 Envs"))


# Create a boxplot of Accuracy by Scenario
ggplot(results_df, aes(x = Scenario, y = Accuracy)) +
  geom_boxplot() +
  labs(title = "Accuracy Across Different Environmental Scenarios",
       x = "Scenario",
       y = "Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust text angle for better readability
#############

# 假设results_df是你的数据框，包含列Scenario和Accuracy
results_df$Scenario <- as.factor(results_df$Scenario)  # 确保Scenario是因子类型

# 如果Scenario是类别型，转化为数值型
results_df$ScenarioNumeric <- as.numeric(results_df$Scenario)

# 构建线性回归模型
model <- lm(Accuracy ~ ScenarioNumeric, data = results_df)

# 查看模型摘要
summary(model)

# 绘制模型的诊断图
par(mfrow = c(2, 2))
plot(model)

plot(model_quad)
# 对新数据进行预测，假设新数据也有Scenario
new_data <- data.frame(ScenarioNumeric = seq(min(results_df$ScenarioNumeric), max(results_df$ScenarioNumeric), length.out = 10))
predictions <- predict(model, newdata = new_data)

# 查看预测结果
predictions
plot(predictions)
