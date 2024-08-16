library(FieldSimR)
library(AlphaSimR)
library(ggplot2)
library(dplyr)

# Clear the work environment
rm(list = ls())

# Population parameters
nFounders <- 50
nChr <- 21
nSegSites <- 1000
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
f1s <- randCross(pop = founders, nCrosses = 250, nProgeny = 60)
DHs <- makeDH(pop = f1s, nDH = 1)

# define number of scenarios and replications
nScenarios <- 52
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
  list(nGenos = 200, nEnvs = 5, nBlocks = 2),
  list(nGenos = 1000, nEnvs = 1, nBlocks = 3),
  list(nGenos = 500, nEnvs = 2, nBlocks = 3),
  list(nGenos = 250, nEnvs = 4, nBlocks = 3),
  list(nGenos = 200, nEnvs = 5, nBlocks = 3),
  list(nGenos = 1000, nEnvs = 1, nBlocks = 4),
  list(nGenos = 500, nEnvs = 2, nBlocks = 4),
  list(nGenos = 250, nEnvs = 4, nBlocks = 4),
  list(nGenos = 200, nEnvs = 5, nBlocks = 4),
  
  list(nGenos = 2000, nEnvs = 1, nBlocks = 1),
  list(nGenos = 1000, nEnvs = 2, nBlocks = 1),
  list(nGenos = 500, nEnvs = 4, nBlocks = 1),
  list(nGenos = 250, nEnvs = 8, nBlocks = 1),
  list(nGenos = 2000, nEnvs = 1, nBlocks = 2),
  list(nGenos = 1000, nEnvs = 2, nBlocks = 2),
  list(nGenos = 500, nEnvs = 4, nBlocks = 2),
  list(nGenos = 250, nEnvs = 8, nBlocks = 2),
  list(nGenos = 2000, nEnvs = 1, nBlocks = 3),
  list(nGenos = 1000, nEnvs = 2, nBlocks = 3),
  list(nGenos = 500, nEnvs = 4, nBlocks = 3),
  list(nGenos = 250, nEnvs = 8, nBlocks = 3),
  list(nGenos = 2000, nEnvs = 1, nBlocks = 4),
  list(nGenos = 1000, nEnvs = 2, nBlocks = 4),
  list(nGenos = 500, nEnvs = 4, nBlocks = 4),
  list(nGenos = 250, nEnvs = 8, nBlocks = 4),
  
  list(nGenos = 3000, nEnvs = 1, nBlocks = 1),
  list(nGenos = 1500, nEnvs = 2, nBlocks = 1),
  list(nGenos = 750, nEnvs = 4, nBlocks = 1),
  list(nGenos = 600, nEnvs = 5, nBlocks = 1),
  list(nGenos = 3000, nEnvs = 1, nBlocks = 2),
  list(nGenos = 1500, nEnvs = 2, nBlocks = 2),
  list(nGenos = 750, nEnvs = 4, nBlocks = 2),
  list(nGenos = 600, nEnvs = 5, nBlocks = 2),
  list(nGenos = 3000, nEnvs = 1, nBlocks = 3),
  list(nGenos = 1500, nEnvs = 2, nBlocks = 3),
  list(nGenos = 750, nEnvs = 4, nBlocks = 3),
  list(nGenos = 600, nEnvs = 5, nBlocks = 3),
  list(nGenos = 3000, nEnvs = 1, nBlocks = 4),
  list(nGenos = 1500, nEnvs = 2, nBlocks = 4),
  list(nGenos = 750, nEnvs = 4, nBlocks = 4),
  list(nGenos = 600, nEnvs = 5, nBlocks = 4),
  
  list(nGenos = 600, nEnvs = 10, nBlocks = 4),
  list(nGenos = 600, nEnvs = 15, nBlocks = 4),
  list(nGenos = 600, nEnvs = 20, nBlocks = 4),
  list(nGenos = 600, nEnvs = 25, nBlocks = 4)
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
      id = DHs@id,
      Trait1 = DHs@gv
    )
    # sampling and parsing as interger
    selected_gvs <- sample_n(dhs_gvs_df, nGenos)
    selected_gvs$id <- as.integer(selected_gvs$id)
    names(selected_gvs)[names(selected_gvs) == "Trait1"] <- "gv.Trait1"
    
    # generate gv data frame using sampling
    gv_df <- data.frame(
      env = rep(1:nEnvs, each = nplots),
      rep = rep(rep(1:nBlocks, each = nplots/nBlocks), times = nEnvs),
      id = rep(rep(selected_gvs$id, times = nBlocks), times = nEnvs),
      gv.Trait1 = runif(n = nplots * nEnvs) 
    )
    # merge
    gv_df <- gv_df %>%
      left_join(selected_gvs, by = "id") %>%
      mutate(gv.Trait1 = coalesce(gv.Trait1.y, gv.Trait1.x)) %>%
      select(-gv.Trait1.x, -gv.Trait1.y)
    
    # sort gv df, to make them grown by order
    sorted_gv <- gv_df %>%
      arrange(env, rep, id)
    
    
    # generate phenotype data frame
    pheno_df <- make_phenotypes(
      gv.df = sorted_gv,
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

