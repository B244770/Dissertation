library(FieldSimR)
library(AlphaSimR)
library(ggplot2)
library(dplyr)

# Clear the work environment
rm(list = ls())


# 定义环境和基因型的数量
nEnvs <- 5
nGenos <- 1000
nBlocks <- 1
total_per_env <- 2000

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
total_plants_per_env <- colSums(design_matrix * rep_matrix)
# 打印结果
print(total_plants_per_env)


# initial results matrix
results <- matrix(NA, nrow = nReps, ncol = nScenarios)

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
# nrows * ncols = nGenos * nEnvs = total plotS
nplots <- nGenos*nBlocks
nrows <- 50
ncols <- nplots/nrows


# generate field trial error
error_ls <- field_trial_error(
  ntraits = 1,
  nenvs = nEnvs,
  nblocks = nBlocks,
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
# 初始化数据框，只包含环境、重复块和ID信息
gv_df <- data.frame(
  env = rep(1:nEnvs, each = nplots),
  rep = rep(rep(1:nBlocks, each = nplots / nBlocks), times = nEnvs),
  id = rep(rep(selected_gvs$id, times = nBlocks), times = nEnvs)
)

with(gv_df, table(id, env))
# with(gv_df, table(id, env, rep))
# 初始化gv.Trait列，先设置为NA
gv_df$gv.Trait1 <- NA

# 根据环境编号动态选择Trait
for (i in 1:nrow(gv_df)) {
  env_num <- gv_df$env[i]  # 当前行的环境编号
  trait_name <- paste("Trait", env_num, sep = "")  # 构造Trait名称，如Trait2, Trait3等
  if (trait_name %in% names(selected_gvs)) {
    # 查找相应id和trait_name的值
    trait_value <- selected_gvs[selected_gvs$id == gv_df$id[i], trait_name]
    if (length(trait_value) > 0) {
      gv_df$gv.Trait1[i] <- trait_value
    }
  }
}

# 查看生成的数据框
# head(gv_df)

# generate phenotype data frame
pheno_df <- make_phenotypes(
  gv.df = gv_df,
  error.df = error_ls$error.df,
  randomise = TRUE
)

# Calculate Pearson's coefficient
results[rep, scenario] <- cor(
  with(pheno_df, tapply(y.Trait1, id, mean)),
  with(gv_df, tapply(gv.Trait1, id, mean))
)


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