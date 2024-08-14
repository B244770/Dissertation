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
    corA <- matrix(0.5, ncol = nEnvs, nrow = nEnvs)
    diag(corA) <- 1
    
    # create founder group
    founders <- quickHaplo(nInd = nFounders, nChr = nChr, segSites = nSegSites)
    
    # simulation params
    SP <- SimParam$new(founders)
    SP$addTraitA(nQtlPerChr = nSegSites,
                 mean = mu, 
                 var = sigma2,
                 corA = corA)
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
    
    # sample gene values from DH group
    dhs_gvs_df <- data.frame(
      id = DHs@id,         # 获取id信息
      Trait1 = DHs@gv      # 获取对应的遗传值
    )
    # sampling and parsing as interger
    selected_gvs <- sample_n(dhs_gvs_df, nGenos)
    selected_gvs$id <- as.integer(selected_gvs$id)
    names(selected_gvs)[names(selected_gvs) == "Trait1"] <- "gv.Trait1"
    
    # generate gv data frame using sampling
    # 总行数
    
    # generate gv data frame using sampling
    gv_df <- data.frame(
      env = rep(1:nEnvs, each = nplots), # 每个环境重复它应有的行数
      rep = rep(rep(1:nBlocks, each = nplots/nBlocks), times = nEnvs), # 每个块在每个环境中重复 ncols * nrows 次
      id = rep(rep(selected_gvs$id, times = nBlocks), times = nEnvs), # 每个基因型按照需要的总行数均匀分布
      gv.Trait1 = rep(rep(selected_gvs$Trait1.Trait1, times = nBlocks), times = nEnvs) # 这里假设你需要随机生成Trait1的值，你可以按需要调整
    )
    # 使用 left_join 合并数据，并直接用 selected_gvs 中的 gv.Trait1 替换 gv_df 中的 gv.Trait1
    gv_df <- gv_df %>%
      left_join(selected_gvs, by = "id") %>%
      mutate(gv.Trait1 = coalesce(gv.Trait1.y, gv.Trait1.x)) %>%
      select(-gv.Trait1.x, -gv.Trait1.y)  # 删除辅助列
    
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