library(FieldSimR)
library(AlphaSimR)
library(ggplot2)
library(dplyr)

# Clear the work environment
rm(list = ls())

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
  }
  
  # 生成表型数据
  ntraits <- ncol(gv.df) - 3  # 减去 'env', 'rep', 'id'
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


# 定义环境和基因型的数量
nEnvs <- 5
nGenos <- 1000
nBlocks <- 1
total_per_env <- 1000

# 设计矩阵初始化，初始都不种植（即设为0）
design_matrix <- matrix(0, nrow = nGenos, ncol = nEnvs, 
                        dimnames = list(paste("Geno", 1:nGenos, sep = ""), 
                                        paste("Env", 1:nEnvs, sep = "")))

# 重复次数矩阵初始化，初始都为0次
rep_matrix <- matrix(0, nrow = nGenos, ncol = nEnvs, 
                     dimnames = list(paste("Geno", 1:nGenos, sep = ""), 
                                     paste("Env", 1:nEnvs, sep = "")))

# 设置随机种子以保证可重复性
set.seed(123)

# 为每个环境随机分配基因型和重复次数，以保证总数为2000
for (env in 1:nEnvs) {
  remaining_plants <- total_per_env
  while (remaining_plants > 0) {
    # 随机选择一个基因型
    selected_genotype <- sample(nGenos, 1)
    # 为了避免单个基因型重复次数过多，限制单次分配的最大重复次数
    max_reps <- min(remaining_plants, sample(1:10, 1))
    
    # 如果这个基因型已经被选择过，就更新重复次数
    if (design_matrix[selected_genotype, env] == 1) {
      added_reps <- min(max_reps, remaining_plants)
      rep_matrix[selected_genotype, env] <- rep_matrix[selected_genotype, env] + added_reps
      remaining_plants <- remaining_plants - added_reps
    } else {
      # 标记这个基因型为已种植
      design_matrix[selected_genotype, env] <- 1
      added_reps <- min(max_reps, remaining_plants)
      rep_matrix[selected_genotype, env] <- added_reps
      remaining_plants <- remaining_plants - added_reps
    }
  }
}

# 转换为数据框，以便于操作和展示
design_df <- as.data.frame(design_matrix)
rep_df <- as.data.frame(rep_matrix)

# 显示部分结果
print("Design Matrix:")
print(design_df[1:10, 1:5])
print("Replication Matrix:")
print(rep_df[1:10, 1:5])

# 计算每个环境的种植总数
colSums(design_matrix)
rowSums(design_matrix)
total_plants_per_env <- colSums(design_matrix * rep_matrix)
# 打印结果
print(total_plants_per_env)

# Population parameters
nFounders <- 20
nChr <- 21
nSegSites <- 300
mu <- rep(4, nEnvs) # environment/trial means, in tonnes/ha
sigma2 <- rep(0.2, nEnvs) # environment/trial genetic variances
H2 <- 0.3 # plot level heritability 
sigma2e <- sigma2 / H2 - sigma2

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

# generate gv data frame using sampling
# 总行数

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

# 查看结果
print(head(gv_df))


# generate phenotype data frame
pheno_df <- make_phenotypes(
  gv.df = gv_df,
  error.df = error_ls$error.df,
  randomise = TRUE
)

# Calculate Pearson's coefficient
# 假设已有环境数 nEnvs
results <- list()  # 使用列表存储每个环境的结果

# 计算每个环境的相关系数
for (env in unique(gv_df$env)) {
  pheno_env <- pheno_df[pheno_df$env == env, ]  # 筛选当前环境的表型数据
  gv_env <- gv_df[gv_df$env == env, ]  # 筛选当前环境的基因型数据
  
  # 计算皮尔森相关系数
  results[[as.character(env)]] <- cor(
    with(pheno_env, tapply(y.Trait1, env, mean)),  # 表型数据
    with(gv_env, tapply(gv, env, mean))  # 基因型数据
  )
}

# 将结果输出，每个环境一个相关系数
for (env in names(results)) {
  cat(sprintf("环境 %s 的皮尔森相关系数为: %f\n", env, results[[env]]))
}



# show results
print(results)
results_df <- data.frame(
  Scenario = factor(rep(1:nScenarios, each = nReps)),
  Accuracy = c(results)
)
head(results_df)

# plot results
ggplot(results_df, aes(x = Scenario, y = Accuracy)) + 
  geom_boxplot() + 
  theme_minimal() + 
  labs(
    title = "Prediction Accuracy Across Scenarios",
    x = "Scenario",
    y = "Accuracy"
  )

#https://cran.r-project.org/web/packages/AlphaSimR/AlphaSimR.pdf