library(FieldSimR)
library(AlphaSimR)
library(ggplot2)
library(dplyr)
library(GA)


# Clear the work environment
rm(list = ls())

# 定义环境和基因型的数量
nEnvs <- 5
nGenos <- 1000
nBlocks <- 1
total_per_env <- 1000
# Population parameters
nFounders <- 20
nChr <- 21
nSegSites <- 300
mu <- rep(4, nEnvs) # environment/trial means, in tonnes/ha
sigma2 <- rep(0.2, nEnvs) # environment/trial genetic variances
H2 <- 0.3 # plot level heritability 
sigma2e <- sigma2 / H2 - sigma2

# 种群模拟
# 生成一个随机正定矩阵
# set.seed(nBlocks)  # 设置随机种子以便结果可重复
random_matrix <- matrix(rnorm(nEnvs^2), nrow = nEnvs, ncol = nEnvs)
cov_matrix <- abs(random_matrix %*% t(random_matrix))  # 生成协方差矩阵
# eigen(cov_matrix)$values # add while loop to run again if negative eigenvalues

# 将协方差矩阵转换为相关矩阵
diag_inv_sqrt <- diag(1 / sqrt(diag(cov_matrix)), nrow = nEnvs)
cor_matrix <- diag_inv_sqrt %*% cov_matrix %*% diag_inv_sqrt
# mean(cor_matrix[upper.tri(cor_matrix)])

# create founder group
founders <- quickHaplo(nInd = nFounders, nChr = nChr, segSites = nSegSites)

# simulation params
SP <- SimParam$new(founders)
SP$addTraitA(nQtlPerChr = nSegSites,
             mean = mu, 
             var = sigma2,
             corA = cor_matrix)
founders <- newPop(founders)

# generate f1 and DH groups
f1s <- randCross(pop = founders, nCrosses = 40, nProgeny = 25)
DHs <- makeDH(pop = f1s, nDH = 1)
# cov2cor(var(DHs@gv))

# cols and rows 
# # nrows * ncols = nGenos * nEnvs = total plotS
# nplots <- nGenos*nBlocks
nrows <- 50
ncols <- total_per_env/nrows


# generate field trial error
error_ls <- field_trial_error(
  ntraits = 1,
  nenvs = nEnvs,
  nblocks = nEnvs, # using envs to simulate blocks
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

# 获取id和性状值
dhs_gvs_df <- data.frame(id = DHs@id)

# DHs@gv 是一个矩阵或数据框，其中列名为 Trait1, Trait2, ..., Traitn
# 将这些性状值合并到 dhs_gvs_df 数据框中
dhs_gvs_df <- cbind(dhs_gvs_df, DHs@gv)


# sampling and parsing as interger
selected_gvs <- sample_n(dhs_gvs_df, nGenos)
selected_gvs$id <- as.integer(selected_gvs$id)

# reload of make_phenotypes
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
  
  # 检查行数是否匹配
  if (nrow(gv.df) != nrow(error.df)) 
    stop("number of rows in 'gv.df' and 'error.df' must match")
  
  # 随机化处理
  if (randomise) {
    gv.df$ord <- sample(nrow(gv.df))
    gv.df <- gv.df[order(gv.df$ord), ]
    gv.df$ord <- NULL # 去除序号列
  }
  
  # 生成表型数据
  # ntraits <- ncol(gv.df) - 3  # 减去 'env', 'rep', 'id'
  ntraits <- 1
  y <- error.df[, !(colnames(error.df) %in% c("env", "block", "col", "row"))] + 
    gv.df[, !(colnames(gv.df) %in% c("env", "rep", "id"))]
  pheno_df <- data.frame(env = gv.df$env, block = error.df$block, col = error.df$col, row = error.df$row, id = gv.df$id, y)
  
  # 重命名表型列
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
  # generate gv data frame using sampling
  # 初始化gv_df数据框
  gv_df <- data.frame(env = integer(), rep = integer(), id = integer(), gv = numeric())
  
  # 遍历每个环境和基因型
  for (env in 1:ncol(design_matrix)) {
    for (geno in 1:nrow(design_matrix)) {
      if (design_matrix[geno, env] == 1) {
        # 这个基因型在这个环境中有种植
        num_reps <- rep_matrix[geno, env]  # 获取重复次数
        geno_id <- selected_gvs$id[geno]   # 获取基因型ID
        trait_name <- paste("Trait", env, sep = "")  # 构建Trait列名，例如'Trait1'
        geno_value <- selected_gvs[[trait_name]][geno]  # 获取基因值
        
        # 为每个重复生成一行数据
        for (rep in 1:num_reps) {
          gv_df <- rbind(gv_df, data.frame(env = env, rep = rep, id = geno_id, gv = geno_value))
        }
      }
    }
  }
  
  # # 查看结果
  # print(head(gv_df))
  
  
  # generate phenotype data frame
  pheno_df <- make_phenotypes(
    gv.df = gv_df,
    error.df = error_ls$error.df,
    randomise = TRUE
  )
  
  # 计算皮尔森相关系数
  result <- cor(
    with(pheno_df, tapply(y.Trait1, env, mean)),
    with(gv_df, tapply(gv, env, mean))
  )
  return(result)
}

# monte carlo
# 调整函数
adjust_env <- function(design_matrix, rep_matrix, nGenos, env, total_per_env) {
  # 随机选择一个基因型并尝试减少重复次数
  selected_genotype <- sample(nGenos, 1)
  if (rep_matrix[selected_genotype, env] > 0) {
    change <- sample(-1:-min(1, rep_matrix[selected_genotype, env]), 1)  # 尝试减少1到最小剩余数量
    rep_matrix[selected_genotype, env] <- rep_matrix[selected_genotype, env] + change
    # 验证列总和是否符合要求
    if (sum(rep_matrix[, env] * design_matrix[, env]) != total_per_env) {
      # 如果不符合，撤销更改
      rep_matrix[selected_genotype, env] <- rep_matrix[selected_genotype, env] - change
    }
  }
}
# 填充函数
fill_env <- function(design_matrix, rep_matrix, nGenos, env, total_per_env) {
  remaining_plants <- total_per_env
  while (remaining_plants > 0) {
    selected_genotype <- sample(nGenos, 1)
    # 给未被完全种植的基因型添加植物
    if (rep_matrix[selected_genotype, env] < 3) {
      max_reps <- min(remaining_plants, sample(1:3, 1), 3 - rep_matrix[selected_genotype, env])
      rep_matrix[selected_genotype, env] <- rep_matrix[selected_genotype, env] + max_reps
      design_matrix[selected_genotype, env] <- 1
      remaining_plants <- remaining_plants - max_reps
    }
  }
}

generate_matrices <- function(nGenos, nEnvs, total_per_env, nIterations = 100) {
  # 设计矩阵和重复次数矩阵的初始化
  design_matrix <- matrix(0, nrow = nGenos, ncol = nEnvs, 
                          dimnames = list(paste("Geno", 1:nGenos, sep = ""), 
                                          paste("Env", 1:nEnvs, sep = "")))
  
  rep_matrix <- matrix(0, nrow = nGenos, ncol = nEnvs, 
                       dimnames = list(paste("Geno", 1:nGenos, sep = ""), 
                                       paste("Env", 1:nEnvs, sep = "")))
  
  # 初始填充
  for (env in 1:nEnvs) {
    fill_env(design_matrix, rep_matrix, nGenos, env, total_per_env)
  }
  
  # SGD样式的迭代改进
  for (iter in 1:nIterations) {
    for (env in 1:nEnvs) {
      # 尝试改变当前环境中的分配
      adjust_env(design_matrix, rep_matrix, nGenos, env, total_per_env)
    }
  }
  print(colSums(design_matrix*rep_matrix))
  return(list(design_matrix = design_matrix, rep_matrix = rep_matrix))
}

# evaluation
# 定义生成矩阵和评估结果的函数
generate_and_evaluate <- function(nGenos, nEnvs, total_per_env) {
  matrices <- generate_matrices(nGenos, nEnvs, total_per_env)
  if (is.null(matrices)) {
    return(NULL)  # 如果生成失败，返回NULL
  }
  score <- simulate_experiment(matrices$design_matrix, matrices$rep_matrix)
  return(list(score = score, matrices = matrices))
}

# 进行多次实验
nIterations <- 10
results <- vector("list", nIterations)  # 创建一个列表来存储结果

for (i in 1:nIterations) {
  result <- generate_and_evaluate(1000, 5, 1000)
  results[[i]] <- result
}

# 过滤掉任何失败的实验
valid_results <- Filter(Negate(is.null), results)

# 查看有效结果的统计信息
scores <- sapply(valid_results, function(x) x$score)
max_score <- max(scores)  # 找到最高分
best_result <- valid_results[[which.max(scores)]]  # 获取最高分的结果



