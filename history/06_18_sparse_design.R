library(FieldSimR)
library(AlphaSimR)
library(ggplot2)

# Clear the work environment
rm(list = ls())

# Population parameters
nFounders <- 20
nChr <- 21
nSegSites <- 300
mu <- 4
sigma2 <- 0.2
H2 <- 0.3
sigma2e <- sigma2 / H2 - sigma2

# create founder group
founders <- quickHaplo(nInd = nFounders, nChr = nChr, segSites = nSegSites)

# simulation params
SP <- SimParam$new(founders)
SP$addTraitA(nQtlPerChr = nSegSites, mean = mu, var = sigma2)
founders <- newPop(founders)

# generate f1 and DH groups
f1s <- randCross(pop = founders, nCrosses = 40, nProgeny = 25)
DHs <- makeDH(pop = f1s, nDH = 1)

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
for (scenario in 1:nScenarios) {
  for (rep in 1:nReps) {
    # params assignment
    params <- scenarios[[scenario]]
    nGenos <- params$nGenos
    nEnvs <- params$nEnvs
    nBlocks <- params$nBlocks
    
    # cols and rows 
    # nrows * ncols = nGenos * nEnvs = total plotS
    ncols <- 20
    nrows <- 50
    
    # generate field trial error
    error_ls <- field_trial_error(
      ntraits = 1,
      nenvs = nEnvs,
      nblocks = nBlocks,
      block.dir = "col", # block layouts ## !!bug here col is not working
      ncols = rep(ncols*nBlocks/nEnvs, nEnvs),
      nrows = nrows,
      varR = sigma2e,
      spatial.model = "AR1",
      col.cor = 0.5,
      row.cor = 0.7,
      prop.spatial = 0.5,
      return.effects = TRUE
    )
    
    # sample gene values from DH group
    selected_gvs <- sample(DHs@gv, nGenos)
    
    # generate gv data frame using sampling
    # 重新计算每个环境应有的行数
    rows_per_env <- nBlocks * ncols * nrows / nEnvs
    
    # generate gv data frame using sampling
    # 总行数
    total_rows_per_env <- nBlocks * ncols * nrows / nEnvs
    
    # generate gv data frame using sampling
    gv_df <- data.frame(
      env = rep(1:nEnvs, each = total_rows_per_env), # 每个环境重复它应有的行数
      rep = rep(rep(1:nBlocks, each = ncols * nrows), times = nEnvs), # 每个块在每个环境中重复 ncols * nrows 次
      id = rep(rep(1:nGenos, each = total_rows_per_env / nGenos), times = nEnvs), # 每个基因型按照需要的总行数均匀分布
      gv.Trait1 = runif(n = nEnvs * total_rows_per_env) # 这里假设你需要随机生成Trait1的值，你可以按需要调整
    )
    
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

