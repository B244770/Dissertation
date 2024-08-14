library(FieldSimR)
library(AlphaSimR)
library(ggplot2)

# 清空工作环境
rm(list = ls())

# 种群参数
nFounders <- 20
nChr <- 21
nSegSites <- 300
mu <- 4
sigma2 <- 0.2
H2 <- 0.3
sigma2e <- sigma2 / H2 - sigma2

# 创建创始群体
founders <- quickHaplo(nInd = nFounders, nChr = nChr, segSites = nSegSites)

# 设置模拟参数
SP <- SimParam$new(founders)
SP$addTraitA(nQtlPerChr = nSegSites, mean = mu, var = sigma2)
founders <- newPop(founders)

# 生成F1和DH群体
f1s <- randCross(pop = founders, nCrosses = 40, nProgeny = 25)
DHs <- makeDH(pop = f1s, nDH = 1)

# 定义场景数和重复数
nScenarios <- 8
nReps <- 10

# 示例场景参数设置
scenarios <- list(
  list(nGenos = 1000, nEnvs = 1, nBlocks = 2),
  list(nGenos = 500, nEnvs = 2, nBlocks = 2),
  list(nGenos = 250, nEnvs = 4, nBlocks = 2),
  list(nGenos = 1000, nEnvs = 1, nBlocks = 2),
  list(nGenos = 1000, nEnvs = 1, nBlocks = 1),
  list(nGenos = 1000, nEnvs = 1, nBlocks = 2),
  list(nGenos = 750, nEnvs = 2, nBlocks = 1.5),
  list(nGenos = 500, nEnvs = 2, nBlocks = 2)
)

# 初始化gv数据框
original_gv <- DHs@gv  # 假设这是已定义的原始基因值数组

# 最大参数设
nGenos_max <- 1000
nEnvs_max <- 4
nBlocks_max <- 2

# 创建基础数据框
gv_df <- data.frame()

# 初始化结果矩阵
results <- matrix(NA, ncol = nScenarios, nrow = nReps)

for (scenario in 1:nScenarios) {
  for (rep in 1:nReps) {
    # 场景参数
    params <- scenarios[[scenario]]
    nGenos <- params$nGenos
    nEnvs <- params$nEnvs
    nBlocks <- params$nBlocks
    # 计算每个环境的行数和列数
    total_plots <- 2000  # 已根据公倍数确定
    ncols <- 40
    nrows <- total_plots / (ncols * nEnvs)
    
    # 生成田间试验误差
    error_ls <- field_trial_error(
      ntraits = 1,
      nenvs = nEnvs,
      nblocks = nBlocks,
      block.dir = "col",
      ncols = ncols,
      nrows = nrows,
      varR = sigma2e,
      spatial.model = "AR1",
      col.cor = 0.5,
      row.cor = 0.7,
      prop.spatial = 0.5,
      return.effects = TRUE
    )
    # 特殊处理
    if (scenario == 1) {
      # 场景1: nGenos = 1000, nEnvs = 1, nBlocks = 2
      nGenos <- 1000
      nEnvs <- 1
      nBlocks <- 2
      gv_df <- data.frame(
        env = rep(1:nEnvs, each = nGenos),
        rep = rep(1:nBlocks, each = nGenos),
        id = rep(1:nGenos, times = nEnvs * nBlocks),
        gv.Trait1 = rep(original_gv, times = nEnvs * nBlocks)
      )
    } else if (scenario == 2) {
      # 场景2: nGenos = 500, nEnvs = 2, nBlocks = 2
      nGenos <- 500
      nEnvs <- 2
      nBlocks <- 2
      genos_env1 <- sample(1:1000, nGenos)
      full_df <- data.frame(
        env = rep(1:2, each = 1000),
        rep = rep(1:nBlocks, each = 1000),
        id = rep(1:1000, times = 2 * nBlocks),
        gv.Trait1 = rep(original_gv, times = 2 * nBlocks)
      )
      gv_df <- full_df[(full_df$env == 1 & full_df$id %in% genos_env1) | 
                         (full_df$env == 2 & !full_df$id %in% genos_env1), ]
    } else if (scenario == 3) {
      # 场景3: nGenos = 250, nEnvs = 4, nBlocks = 2
      nGenos <- 250
      nEnvs <- 4
      nBlocks <- 2
      genos_envs <- split(1:1000, rep(1:nEnvs, each = 250))
      gv_df <- data.frame(
        env = rep(1:nEnvs, each = nGenos),
        rep = rep(1:nBlocks, each = nGenos),
        id = rep(unlist(genos_envs), times = nBlocks),
        gv.Trait1 = rep(original_gv[unlist(genos_envs)], times = nBlocks)
      )
    } else if (scenario == 4 || scenario == 6 || scenario == 8) {
      # 场景4, 6, 8: nGenos = 1000, nEnvs = 1 (或2), nBlocks = 2
      nGenos <- 1000
      nEnvs <- ifelse(scenario == 8, 2, 1)
      nBlocks <- 2
      gv_df <- data.frame(
        env = rep(1:nEnvs, each = nGenos),
        rep = rep(1:nBlocks, each = nGenos),
        id = rep(1:nGenos, times = nEnvs * nBlocks),
        gv.Trait1 = rep(original_gv, times = nEnvs * nBlocks)
      )
    } else if (scenario == 5) {
      # 场景5: nGenos = 1000, nEnvs = 1, nBlocks = 1
      nGenos <- 1000
      nEnvs <- 1
      nBlocks <- 1
      gv_df <- data.frame(
        env = rep(1:nEnvs, each = nGenos),
        rep = rep(1:nBlocks, each = nGenos),
        id = rep(1:nGenos, times = nEnvs * nBlocks),
        gv.Trait1 = rep(original_gv, times = nEnvs * nBlocks)
      )
    } 
    # else if (scenario == 7) {
    #   # 场景7: nGenos = 750, nEnvs = 2, nBlocks = 2
    #   nGenos <- 750
    #   nEnvs <- 2
    #   nBlocks <- 2
    #   genos_env1 <- sample(1:1000, 375)
    #   full_df <- data.frame(
    #     env = rep(1:2, each = 1000),
    #     rep = rep(1:nBlocks, each = 1000),
    #     id = rep(1:1000, times = 2 * nBlocks),
    #     gv.Trait1 = rep(original_gv, times = 2 * nBlocks)
    #   )
    #   gv_df <- full_df[(full_df$env == 1 & full_df$id %in% genos_env1) | 
    #                      (full_df$env == 2 & !full_df$id %in% genos_env1), ]
    # }
    
    
    # # 构建基因值数据框
    # gv_df <- data.frame(
    #   env = rep(1:nEnvs, each = nGenos),
    #   rep = rep(1:nBlocks, each = nGenos),
    #   id = rep(1:nGenos, each = nEnvs * nBlocks),
    #   gv.Trait1 = rep(DHs@gv, each = nBlocks * nEnvs)
    # )
    
    # # 处理特定场景的基因型样本
    # if (scenario == 1) {
    #   # 场景1: nGenos = 1000, nEnvs = 1, nBlocks = 2
    #   # 没有特殊处理
    # } else if (scenario == 2) {
    #   # 场景2: nGenos = 500, nEnvs = 2, nBlocks = 2
    #   genos_env1 <- sample(1:1000, nGenos)
    #   gv_df <- gv_df[(gv_df$env == 1 & gv_df$id %in% genos_env1) | 
    #                    (gv_df$env == 2 & !gv_df$id %in% genos_env1), ]
    # } else if (scenario == 3) {
    #   # 场景3: nGenos = 250, nEnvs = 4, nBlocks = 2
    #   # 使用交错基因型进行环境的划分
    #   genos_envs <- split(1:1000, 1:4)
    #   gv_df <- gv_df[gv_df$env == rep(1:4, each = 250), ]
    # } else if (scenario == 4) {
    #   # 场景4: nGenos = 1000, nEnvs = 1, nBlocks = 2
    #   # 没有特殊处理
    # } else if (scenario == 5) {
    #   # 场景5: nGenos = 1000, nEnvs = 1, nBlocks = 1
    #   # 只有一个block
    #   gv_df$rep <- 1
    # } else if (scenario == 6) {
    #   # 场景6: nGenos = 1000, nEnvs = 1, nBlocks = 2
    #   # 没有特殊处理
    # } else if (scenario == 7) {
    #   # 场景7: nGenos = 750, nEnvs = 2, nBlocks = 1.5
    #   genos_env1 <- sample(1:750, 375)
    #   gv_df <- gv_df[(gv_df$env == 1 & gv_df$id %in% genos_env1) | 
    #                    (gv_df$env == 2 & !gv_df$id %in% genos_env1), ]
    # } else if (scenario == 8) {
    #   # 场景8: nGenos = 500, nEnvs = 2, nBlocks = 2
    #   # 没有特殊处理
    # }
    
    # # 假设原始 gv_df 已经正确构建并包含了基因值 'gv.Trait1'
    # original_gv_df <- gv_df  # 保存原始的 gv_df，其中包含基因值
    # 
    # # 使用 expand.grid 生成所有可能的组合
    # expanded_gv_df <- expand.grid(env = unique(error_df$env),
    #                               block = unique(error_df$block),
    #                               col = unique(error_df$col),
    #                               row = unique(error_df$row),
    #                               id = unique(original_gv_df$id))
    # 
    # # 将原始的基因值合并回来
    # gv_df <- merge(expanded_gv_df, original_gv_df, by = c("env", "block", "col", "row", "id"), all.x = TRUE)
    
    # 生成表型数据框
    print(scenario)
    print(nrow(gv_df))
    print(nrow(error_ls$error.df))
    pheno_df <- make_phenotypes(
      gv.df = gv_df,
      error.df = error_ls$error.df,
      randomise = TRUE
    )
    
    # 计算准确度并存储结果
    results[rep, scenario] <- cor(
      with(pheno_df, tapply(y.Trait1, id, mean)),
      with(gv_df, tapply(gv.Trait1, id, mean))
    )
  }
}

# 显示结果
print(results)
results_df <- data.frame(
  Scenario = factor(rep(1:nScenarios, each = nReps)),
  Accuracy = c(results)
)
head(results_df)

# 绘制结果图
ggplot(results_df, aes(x = Scenario, y = Accuracy)) + 
  geom_boxplot() + 
  theme_minimal() + 
  labs(
    title = "Prediction Accuracy Across Scenarios",
    x = "Scenario",
    y = "Accuracy"
  )

