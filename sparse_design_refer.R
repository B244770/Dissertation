library(AlphaSimR)
library(FieldSimR)
library(ggplot2)


# Clear the work environment
rm(list = ls())


ntraits <- 2 # Number of traits
nenvs <- 3 # Number of environments
nblocks <- c(2, 2, 3) # Number of blocks per environment
block_dir <- "col" # Arrangement of blocks ("side-by-side")
ncols <- c(10, 10, 15) # Number of columns per environment
nrows <- 20 # Number of rows per environment
plot_length <- 8 # Plot length; here in meters (column direction)
plot_width <- 2 # Plot width; here in meters (row direction)


H2 <- c(0.3, 0.3, 0.3, 0.5, 0.5, 0.5) # c(Yld:E1, Yld:E2, Yld:E3, Pht:E1, Pht:E2, Pht:E3)
var <- c(0.086, 0.12, 0.06, 15.1, 8.5, 11.7) # c(Yld:E1, Yld:E2, Yld:E3, Pht:E1, Pht:E2, Pht:E3)

# Calculation of error variances based on the genetic variance and target heritability vectors.
calc_varR <- function(var, H2) {
  varR <- (var / H2) - var
  return(varR)
}

varR <- calc_varR(var, H2)
round(varR, 2) # Vector of error variances: c(Yld:E1, Yld:E2, Yld:E3, Pht:E1, Pht:E2, Pht:E3)

# spatical error
spatial_model <- "Bivariate" # Spatial error model.
prop_spatial <- 0.4 # Proportion of spatial trend.

ScorR <- rand_cor_mat(ntraits, min.cor = 0, max.cor = 0.5, pos.def = TRUE)
round(ScorR, 2)

# extraneous error
ext_ord <- "zig-zag"
ext_dir <- "row"
prop_ext <- 0.2

EcorR <- rand_cor_mat(ntraits, min.cor = 0, max.cor = 0.5, pos.def = TRUE)
round(EcorR, 2)

error_ls <- field_trial_error(
  ntraits = ntraits,
  nenvs = nenvs,
  nblocks = nblocks,
  block.dir = block_dir,
  ncols = ncols,
  nrows = nrows,
  plot.length = plot_length,
  plot.width = plot_width,
  varR = varR,
  ScorR = ScorR,
  EcorR = EcorR,
  RcorR = NULL,
  spatial.model = spatial_model,
  prop.spatial = prop_spatial,
  ext.ord = ext_ord,
  ext.dir = ext_dir,
  prop.ext = prop_ext,
  return.effects = TRUE
)

e_total_env1 <- error_ls$error.df[error_ls$error.df$env == 1, ]
e_terms_env1 <- error_ls$Trait1[error_ls$Trait1$env == 1, ]

# Total plot error
plot_effects(e_total_env1, effect = "e.Trait1", labels = TRUE)


# Spatial error simulated using bivariate interpolation
plot_effects(e_terms_env1, effect = "e.spat", labels = TRUE)


# Random error
plot_effects(e_terms_env1, effect = "e.rand", labels = TRUE)


# Extraneous variation in the row direction
plot_effects(e_terms_env1, effect = "e.ext.row")


gv_df <- gv_df_unstr

pheno_df <- make_phenotypes(
  gv_df,
  error_ls$error.df,
  randomise = TRUE
)

pheno_env1 <- pheno_df[pheno_df$env == 1, ] # Extract phenotypes in environment 1.
