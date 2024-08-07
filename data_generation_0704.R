library(FieldSimR)
library(AlphaSimR)
library(ggplot2)
library(dplyr)

# Clear the work environment
rm(list = ls())


# 定义场景数和重复数
nScenarios <- 12
nReps <- 10

# scenarios params
scenarios <- list(
  list(nGenos = 1000, nEnvs = 1, nBlocks = 1),
  list(nGenos = 500, nEnvs = 2, nBlocks = 1),
  list(nGenos = 250, nEnvs = 4, nBlocks = 1),
  list(nGenos = 200, nEnvs = 5, nBlocks = 1),
  list(nGenos = 1000, nEnvs = 1, nBlocks = 2),
  list(nGenos = 500, nEnvs = 2, nBlocks = 2),
  list(nGenos = 250, nEnvs = 4, nBlocks = 2), ## !!bug here
  list(nGenos = 200, nEnvs = 5, nBlocks = 2),
  list(nGenos = 1000, nEnvs = 1, nBlocks = 3),
  list(nGenos = 500, nEnvs = 2, nBlocks = 3),
  list(nGenos = 250, nEnvs = 4, nBlocks = 3),
  list(nGenos = 200, nEnvs = 5, nBlocks = 3)
)

# initial results matrix
results <- matrix(NA, nrow = nReps, ncol = nScenarios)

# Replication is implemented in both way of duplicated blocks and repeated simulation
for (scenario in 1:nScenarios) { # scenario <- 2; rep <- 1
  for (rep in 1:nReps) {
    
    # params assignment
    params <- scenarios[[scenario]]
    nGenos <- params$nGenos
    nEnvs <- params$nEnvs
    nBlocks <- params$nBlocks
    
    # Population parameters
    nFounders <- 20
    nChr <- 21
    nSegSites <- 300
    mu <- rep(4, nEnvs) # environment/trial means, in tonnes/ha
    sigma2 <- rep(0.2, nEnvs) # environment/trial genetic variances
    H2 <- 0.3 # plot level heritability 
    sigma2e <- sigma2 / H2 - sigma2
    
    # 生成一个随机正定矩阵
    set.seed(123)  # 设置随机种子以便结果可重复
    random_matrix <- matrix(rnorm(nEnvs^2), nrow = nEnvs, ncol = nEnvs)
    cov_matrix <- random_matrix %*% t(random_matrix)  # 生成协方差矩阵
    
    # 将协方差矩阵转换为相关矩阵
    diag_inv_sqrt <- diag(1 / sqrt(diag(cov_matrix)))
    cor_matrix <- diag_inv_sqrt %*% cov_matrix %*% diag_inv_sqrt
    
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
    
    # 检查并添加最多100个性状
    for (i in 1:100) {
      trait_name <- paste("Trait", i, sep = "")  # 构造性状名，如 Trait1, Trait2, ..., Trait100
      if (trait_name %in% names(selected_gvs)) {  # 如果该性状存在于selected_gvs中
        gv_df[[paste("gv", trait_name, sep = ".")]] <- rep(rep(selected_gvs[[trait_name]], times = nBlocks), times = nEnvs)
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
  }
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