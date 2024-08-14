library(AlphaSimR)
library(GA)
library(FieldSimR)

# Set population parameters
nFounders <- 20
nChr <- 21
nSegSites <- 300
mu <- 4
sigma2 <- 0.2
H2 <- 0.3
sigma2e <- sigma2 / H2 - sigma2

# Create founder population
founders <- quickHaplo(nInd = nFounders, nChr = nChr, segSites = nSegSites)

# Setup simulation parameters
SP <- SimParam$new(founders)
SP$addTraitA(nQtlPerChr = nSegSites, mean = mu, var = sigma2)
founders <- newPop(founders)

# Generate F1 and DH populations
f1s <- randCross(pop = founders, nCrosses = 40, nProgeny = 25)
DHs <- makeDH(pop = f1s, nDH = 1)

# Define the fitness function for the genetic algorithm
fitness_function <- function(x) {
  nrows <- round(x[1])
  ncols <- round(x[2])
  nblocks <- round(x[3])
  
  # Ensure nblocks is at least 1
  nblocks <- max(nblocks, 1)
  
  # Generate field trial errors
  error_ls <- field_trial_error(
    ntraits = 1,
    nenvs = 1,
    nblocks = nblocks,
    block.dir = "col",
    ncols = ncols,
    nrows = nrows,
    varR = sigma2e,
    spatial.model = "AR1",
    col.cor = 0.5,
    row.cor = 0.7,
    prop.spatial = 0.5,
    return.effects = FALSE
  )
  
  # Simulate phenotypic data
  SP$setPlotDims(nrow = nrows, ncol = ncols)
  phenotypes <- runTrial(pop = DHs, SP = SP, errorDF = error_ls)
  
  # The goal is to minimize the variance of the trait
  return(-var(phenotypes$pheno[, "Trait1"]))  # GA minimizes, so negate variance
}

# Set up and run the genetic algorithm
ga <- ga(
  type = "real-valued",
  fitness = fitness_function,
  lower = c(10, 10, 1),  # Lower bounds for nrows, ncols, and nblocks
  upper = c(50, 50, 5),  # Upper bounds for nrows, ncols, and nblocks
  popSize = 50,       # Population size
  maxiter = 100,      # Maximum number of iterations
  run = 100           # Maximum number of generations without improvement
)

# Print the results
cat("Optimal rows, columns, and number of blocks:\n")
print(round(ga@solution))
cat("Objective function value (negative variance):\n")
print(ga@fitnessValue)
