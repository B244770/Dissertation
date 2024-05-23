library(AlphaSimR)
library(rsm)
library(dplyr)

# Set population parameters
nFounders <- 20
nChr <- 21
nSegSites <- 300
mu <- 4
sigma2 <- 0.2
H2 <- 0.3
sigma2e <- sigma2 / H2 - sigma2
nGenos <- 1000  # Total number of genotypes remains constant

# Create founder population
founders <- quickHaplo(nInd = nFounders, nChr = nChr, segSites = nSegSites)

# Setup simulation parameters
SP <- SimParam$new(founders)
SP$addTraitA(nQtlPerChr = nSegSites, mean = mu, var = sigma2)
founders <- newPop(founders)

# Generate F1 and DH populations
f1s <- randCross(pop = founders, nCrosses = 40, nProgeny = 25)
DHs <- makeDH(pop = f1s, nDH = 1)

# Define a function to simulate and analyze data
simulate_experiment <- function(nEnvs, nReps, plotsPerEnv) {
  results <- numeric(nEnvs)
  
  for (i in 1:nEnvs) {
    SP$setReplications(nReps)
    pop <- sampleInd(DHs, nInd = plotsPerEnv * nEnvs)
    simData <- runSim(pop, SP, nRep = nReps)
    results[i] <- mean(simData$pheno[,'Trait1'])
  }
  
  # Use the variance as a response to minimize
  return(var(results))
}

# Define ranges for parameters
envs_range <- 2:5        # Environment range from 2 to 5
reps_range <- 1:3        # Repetitions per environment from 1 to 3
plots_range <- c(200, 250, 300)  # Plots per environment options

# Generate a design data frame with all combinations
design <- expand.grid(nEnvs = envs_range, nReps = reps_range, plotsPerEnv = plots_range)

# Evaluate all combinations
design$Y <- apply(design, 1, function(x) simulate_experiment(x['nEnvs'], x['nReps'], x['plotsPerEnv']))

# Fit the response surface model
model <- rsm(Y ~ FO(nEnvs, nReps, plotsPerEnv), data = design)

# Perform optimization using steepest ascent
optimal_params <- steepest(model, path = TRUE)

# View the results of the optimization
print(summary(optimal_params))
