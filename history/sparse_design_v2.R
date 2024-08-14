library(FieldSimR)
library(AlphaSimR)
library(ggplot2)

# 清空工作环境
rm(list = ls())


nScenarios <- 3
nReps <- 10
results <- matrix(NA, ncol = nScenarios, nrow = nReps)

for(scenario in 1:nScenarios){ # scenario <- 2

  for (rep in 1:nReps) { # rep <- 1
    
# Params
nGenos <- c(1000, 500, 250)[scenario]
nEnvs <- c(1, 2, 4)[scenario]
nFounders <- 20
nChr <- 21
nSegSites <- 300
mu <- 4
sigma2 <- 0.2
H2 <- 0.3
sigma2e <- sigma2 / H2 - sigma2

# Founders
founders <- quickHaplo(nInd = nFounders, nChr = nChr, segSites = nSegSites) #runMacs(nInd = nFounders, nChr = nChr, segSites = nSegSites, inbred = TRUE, species = "WHEAT", nThreads = 2)

# Sim Params
SP <- SimParam$new(founders)
SP$addTraitA(nQtlPerChr = nSegSites, mean = mu, var = sigma2)
founders <- newPop(founders)

# F1 genos & DHs
f1s <- randCross(pop = founders, nCrosses = 40, nProgeny = 25)
DHs <- makeDH(pop = f1s, nDH = 1)

# Count rows and columns for each env
print(paste0("Scenario: ",scenario, ", Rep: ", rep))
total_plots <- 2000 / nEnvs
ncols <- 40
(nrows <- ceiling(total_plots/ncols))
  
# field_trial_error parameters
error_df <- field_trial_error(
  ntraits = 1,
  nenvs = nEnvs,
  nblocks = 2,
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

# plot_effects(error_df$Trait1[error_df$Trait1$env == 1,], effect = "e.spat")

# now construct phenotypes...
gv.df <- data.frame(env = rep(1:nEnvs, each = 1000 * 2),
                    rep = rep(1:2, each = 1000),
                    id = 1:1000,
                    gv.Trait1 = c(DHs@gv))
nrow(gv.df)
if(scenario == 2){genos_env1 <- sample(1:1000, nGenos)
  gv.df <- gv.df[(gv.df$env == 1 & gv.df$id %in% genos_env1) | (gv.df$env == 2 & !gv.df$id %in% genos_env1),]
  }
rowSums(with(gv.df, table(id, env)))
pheno_df <- make_phenotypes(gv.df = gv.df,
                            error.df = error_df,
                            randomise = TRUE)
  # plot_effects(pheno_df, effect = "y.Trait1") # visualise phenotypes
  
#plot(gv.df$gv.Trait1[gv.df$rep == 1],
      # with(pheno_df, tapply(y.Trait1, id, mean)))
  # TO DO: confusion 
# select top 100
# results[rep, scenario] <- cor(with(pheno_df, tapply(y.Trait1, id, mean)),
#                               gv.df$gv.Trait1[gv.df$rep == 1])
results[rep, scenario] <- cor(with(pheno_df, tapply(y.Trait1, id, mean)),
                              with(gv.df, tapply(gv.Trait1, id, mean)))
}
}
 
results
results_df <- data.frame(Scenario = factor(rep(1:nScenarios, each = nReps)),
                         Accuracy = c(results))
head(results_df)
ggplot(results_df, aes(x = Scenario, y = Accuracy)) + geom_boxplot()



