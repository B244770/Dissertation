library(FieldSimR)
library(AlphaSimR)
library(ggplot2)
library(dplyr)
library(asreml)

# library(devtools)
# install_github("crWerner/fieldsimr") # install newest version

nTPE <- 50
nEnvs_all <- c(1, 2, 4, 5)
nInds <- 1000
nPlots <- 2000
(nPlots_per_Ind <- nPlots/nInds)
(nPlots_per_Env <- floor(nPlots/nEnvs_all))

# range of correlations between environments
min_cor <- 0
max_cor <- 1
  
nReps <- 20
met_tpe <- gblup_acc <- blup_acc <- pheno_acc <- matrix(0, ncol = length(nEnvs_all), nrow = nReps)

for(j in 1:nReps){
  
for(i in 1:length(nEnvs_all)){ # i <- 1

nEnvs <- nEnvs_all[i]

# Population parameters
nFounders <- 20
nChr <- 21
nQTN <- 300
nSNP <- 300
nSegSites <- (nQTN + nSNP)

mu <- rep(4, nTPE) # environment/trial means, in tonnes/ha
sigma2 <- rep(0.2, nTPE) # environment/trial genetic variances
H2 <- 0.3 # plot level heritability 
sigma2e <- (sigma2 / H2 - sigma2)[1:nEnvs]

# correlation matrix between environments
# set.seed(123)
cor_matrix <- FieldSimR::rand_cor_mat(n = nTPE, min.cor = min_cor, max.cor = max_cor, pos.def = T)

# create founders
founders <- quickHaplo(nInd = nFounders, nChr = nChr, segSites = nSegSites)

# simulation params
SP <- SimParam$new(founders)
SP$addTraitA(nQtlPerChr = nQTN,
             mean = mu, 
             var = sigma2,
             corA = cor_matrix)
founders <- newPop(founders)

SP$addSnpChip(nSNP)

# generate f1 and DH groups
f1s <- randCross(pop = founders, nCrosses = 40, nProgeny = 25)
DHs <- makeDH(pop = f1s, nDH = 1)

M <- pullSnpGeno(DHs)
M <- scale(M, center = T)
G <- M %*% t(M)
G <- G/mean(diag(G)) + diag(0.0000001, nrow = nrow(G))

# genetic values

gv_df <- data.frame(env = factor(rep(1:nTPE, each = nInds)),
                    id = factor(as.numeric(DHs@id)),
                    rep = factor(1),
                    gv.Trait1 = c(DHs@gv))

# sample nEnvs from gv_df
gv_df_sample <- droplevels(gv_df[gv_df$env %in% 1:nEnvs_all[i],])
met_tpe[j,i] <- cor(with(gv_df_sample, tapply(gv.Trait1, id, mean)), with(gv_df, tapply(gv.Trait1, id, mean)))

# create an experimental design
source("make_design.R")
if(nEnvs == 1){design_df <- data.frame(env = factor(1), 
                                       id = levels(gv_df_sample$id),
                                       nreps = 2)}
if(nEnvs > 1){design_df <- make_design(nenvs = nEnvs,
                                       ninds = nInds,
                                       nplots = nPlots_per_Ind,
                                       prep = 0.25)
design_df$id <- factor(design_df$id, labels = levels(gv_df_sample$id))
}


# if(nEnvs == 1){design_df <- as.table(matrix(2, ncol = nEnvs, nrow = nInds))}
# if(nEnvs == 2){
#   # lets say we want 250 individuals to have two reps in a single environment
#   # and 750 to have one rep in each environment
#   design_df <- as.table(matrix(1, ncol = 2, nrow = nInds))
#   design_df[sample(1:nInds, 250),1] <- 2
#   design_df[design_df[,1] == 2, 2] <- 0
#   design_df[sample(which(design_df[,1] != 2), 250), 2] <- 2
#   design_df[design_df[,2] == 2, 1] <- 0
#   # colSums(design_df) # 1000 1000
#   # rowSums(design_df) # 2
#   }
# rownames(design_df) <- levels(gv_df$id)
# colnames(design_df) <- levels(gv_df$env)[1:i]


# now create errors
nRows <- 50

nCols <- with(design_df, tapply(nreps, env, sum))/nRows

error_df <- field_trial_error(nenvs = nEnvs, 
                              ncols = nCols,
                              nrows = nRows,
                              nblocks = 2)
nrow(error_df)
plot_effects(error_df[error_df$env == 1,], effect = "e.Trait1")

# construct phenotypes, based on the design above
pheno_df <- make_phenotypes(gv.df = gv_df_sample,
                            error.df = error_df,
                            design.df = design_df,
                            randomise = TRUE)
plot_effects(pheno_df[pheno_df$env == 1,], effect = "y.Trait1")

# Obtain prediction 
# i) Genotype effects across environments
# obtain accuracy from phenotypes
(pheno_acc[j,i] <- cor(with(pheno_df, tapply(y.Trait1, id, mean)), with(gv_df, tapply(gv.Trait1, id, mean))))

# obtain accuracy from blups
asr1 <- asreml(y.Trait1 ~ env,
               random = ~id,
               # residual = ~dsum(~ar1(col):ar1(row)| env),
               data = pheno_df)
(blup_acc[j,i] <- cor(asr1$coefficients$random, with(gv_df, tapply(gv.Trait1, id, mean))))

# obtain accuracy from genomic blups
asr2 <- asreml(y.Trait1 ~ env,
               random = ~ vm(id, G),
               # residual = ~dsum(~ar1(col):ar1(row)| env),
               data = pheno_df)
(gblup_acc[j,i] <- cor(asr2$coefficients$random, with(gv_df, tapply(gv.Trait1, id, mean))))

}

}

apply(gblup_acc, 2, mean)
apply(gblup_acc, 2, min)
apply(gblup_acc, 2, max)

summary(gblup_acc)
