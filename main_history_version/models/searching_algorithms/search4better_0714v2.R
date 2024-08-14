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
set.seed(nBlocks)  # 设置随机种子以便结果可重复
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
  # 假设design_matrix, rep_matrix, 和 selected_gvs 已经定义并初始化
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
  View(gv_df)
  View(error_ls$error.df)
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

# 设计矩阵初始化，初始都不种植（即设为0）
design_matrix <- matrix(0, nrow = nGenos, ncol = nEnvs, 
                        dimnames = list(paste("Geno", 1:nGenos, sep = ""), 
                                        paste("Env", 1:nEnvs, sep = "")))

# 重复次数矩阵初始化，初始都为0次
rep_matrix <- matrix(0, nrow = nGenos, ncol = nEnvs, 
                     dimnames = list(paste("Geno", 1:nGenos, sep = ""), 
                                     paste("Env", 1:nEnvs, sep = "")))

# 定义适应度函数
fitness <- function(chromosome) {
  design_matrix <- matrix(as.integer(chromosome[1:(nGenos*nEnvs)] > 0.5), nrow = nGenos, ncol = nEnvs)
  
  # 从染色体中获取部分，假设值范围是[0, 1]
  rep_part = chromosome[(nGenos * nEnvs + 1):(2 * nGenos * nEnvs)]
  # 缩放到[0, 3]范围并取整
  scaled_rep_values = floor(rep_part * 4)  # 乘以4，因为需要包括3这个端点
  # 保证不超过3
  scaled_rep_values[scaled_rep_values > 3] <- 3
  # 转换成矩阵
  rep_matrix <- matrix(scaled_rep_values, nrow = nGenos, ncol = nEnvs)
  
  print(colSums(design_matrix*rep_matrix))
  print(total_per_env)
  # 转换为数据框，以便于操作和展示
  design_df <- as.data.frame(design_matrix)
  rep_df <- as.data.frame(rep_matrix)
  
  # 计算设计矩阵的列和的惩罚
  if(any(colSums(design_matrix*rep_matrix) - total_per_env)!=0){
    return(-Inf)
  }
  
  # 确保设计矩阵每行和为2
  if(!(rowSums(design_matrix) - 2)){
    return(-Inf)
  }

  
  # 确保重复次数矩阵元素不超过3
  if(sum(rep_matrix[rep_matrix > 3] - 3) > 0){
    return(-Inf)
  }
  
  # 模拟实验计算结果
  result <- tryCatch({
    # 尝试运行 simulate_experiment 函数
    simulate_experiment(design_matrix, rep_matrix)
  }, error = function(e) {
    # 如果发生错误，打印错误消息并返回0
    message("Error occurred: ", e$message)
    return(0)
  }
  )
  
  return(result)
}


# 重构生成染色体的函数
generate_valid_chromosome <- function(nGenos, nEnvs, total_per_env) {
  design_matrix <- matrix(0, nrow = nGenos, ncol = nEnvs)
  rep_matrix <- matrix(0, nrow = nGenos, ncol = nEnvs)
  
  for (env in 1:nEnvs) {
    remaining_plants = total_per_env
    attempts = 0
    
    while (remaining_plants > 0 && attempts < total_per_env) {
      # 随机选择一个基因型
      geno <- sample(nGenos, 1)
      
      if (sum(design_matrix[geno, ]) < 2) {  # 确保每个基因型在所有环境中最多被选择两次
        # 确定可以添加的最大数量
        max_add = min(3 - rep_matrix[geno, env], remaining_plants)  # 最多可增加到3，同时不能超过剩余需要的量
        if (design_matrix[geno, env] == 0 && max_add > 0) {
          design_matrix[geno, env] <- 1
          rep_matrix[geno, env] <- rep_matrix[geno, env] + max_add
          remaining_plants <- remaining_plants - max_add
        }
      }
      attempts <- attempts + 1
    }
    
    # 如果尝试次数过多仍未满足条件，则简单均匀分配剩余量
    if (remaining_plants > 0) {
      eligible_indices <- which(design_matrix[, env] == 0 & rowSums(design_matrix) < 2)
      num_eligible <- length(eligible_indices)
      if (num_eligible > 0) {
        evenly_dist <- floor(remaining_plants / num_eligible)
        extra <- remaining_plants %% num_eligible
        
        rep_matrix[eligible_indices, env] <- evenly_dist + (1:num_eligible <= extra)
        design_matrix[eligible_indices, env] <- 1
      }
    }
  }
  
  chromosome <- c(as.vector(design_matrix), as.vector(rep_matrix))
  return(chromosome)
}

# 遗传算法设置
ga <- ga(type = "real-valued", 
         fitness = fitness, 
         nBits = 2*nGenos*nEnvs, 
         popSize = 200, 
         maxiter = 500, 
         pmutation = 0.05,
         lower = c(rep(0, nGenos*nEnvs), rep(0, nGenos*nEnvs)),
         upper = c(rep(1, nGenos*nEnvs), rep(3, nGenos*nEnvs)),
         popInit = generate_valid_chromosome(nGenos, nEnvs, total_per_env)
)


# 获取最佳解的染色体
best_chromosome <- ga@solution

# 转换染色体为矩阵
design_matrix <- matrix(as.integer(best_chromosome[1:(nGenos*nEnvs)] > 0.5), nrow = nGenos, ncol = nEnvs)
rep_matrix <- matrix(as.integer(best_chromosome[(nGenos*nEnvs+1):(2*nGenos*nEnvs)]), nrow = nGenos, ncol = nEnvs)

# 打印矩阵
print("Optimal Design Matrix:")
print(design_matrix)
print("Optimal Replication Matrix:")
print(rep_matrix)
print("Optimal Fitness Value:", ga@fitnessValue)
