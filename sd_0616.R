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
nScenarios <- 8
nReps <- 10

# scenarios params
scenarios <- list(
  list(nGenos = 1000, nEnvs = 1, nBlocks = 1),
  list(nGenos = 500, nEnvs = 2, nBlocks = 1),
  list(nGenos = 250, nEnvs = 4, nBlocks = 1),
  list(nGenos = 200, nEnvs = 5, nBlocks = 1),
  list(nGenos = 1000, nEnvs = 1, nBlocks = 2),
  list(nGenos = 500, nEnvs = 2, nBlocks = 2),
  list(nGenos = 250, nEnvs = 4, nBlocks = 2),
  list(nGenos = 200, nEnvs = 5, nBlocks = 2)
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
      # total_plots <- 2000 
      ncols <- 20
      nrows <- 50
      
      # generate field trial error
      error_ls <- field_trial_error(
        ntraits = 1,
        nenvs = nEnvs,
        nblocks = nBlocks,
        block.dir = "side-by-side", # layouts
        ncols = ncols/nEnvs,
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
      gv_df <- data.frame(
        env = rep(1:nEnvs, each = nGenos*nBlocks),
        rep = rep(nBlocks, each = nGenos*nBlocks),
        id = rep(rep(1:nGenos, each = nBlocks), times = nEnvs),
        gv.Trait1 = rep(rep(selected_gvs, each = nBlocks), nEnvs)
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

