
make_design <- function(nenvs = 4, 
                        ninds = 1000,
                        nplots = 2,
                        prep = 0.25
                        ){
  if (!(is.atomic(nenvs) && length(nenvs) == 1L)) stop("'nenvs' must be a scalar")
  if (nenvs < 1 || nenvs %% 1 != 0) stop("'nenvs' must be a positive integer")
  if (!(is.atomic(ninds) && length(ninds) == 1L)) stop("'ninds' must be a scalar")
  if (ninds < 1 || ninds %% 1 != 0) stop("'ninds' must be a positive integer")
  if (!(is.atomic(nplots) && length(nplots) == 1L)) stop("'nplots' must be a scalar")
  if (nplots < 1 || nplots %% 1 != 0) stop("'nplots' must be a positive integer")
  if (!(is.atomic(prep) && length(prep) == 1L)) stop("'prep' must be a scalar")
  if (prep < 0 || prep > 1) stop("'prep' must be a positive value <= 1")
  
  # general checks of dimension compatibility
  
  nplots_total <- nplots * ninds
  nplots_env <- floor(nplots_total/nenvs)
  # check leftovers 
  # nplots_total/nenvs - nplots_env
  
  double_rep <- nplots_env - floor(nplots_env/(prep + 1))
  # check leftovers 
  # prep * nplots_env - double_rep
  
  design_df1 <- kronecker(diag(nenvs), matrix(2, ncol = 1, nrow = double_rep))
  
  env_pairs <- t(combn(seq_len(nenvs), 2))
  single_rep <- floor((ninds - double_rep*nenvs) / nrow(env_pairs))
  # floor(2 * (nplots_env - double_rep*2) / nrow(env_pairs))
  # check leftovers
  
  design_df2 <- matrix(0, ncol = nenvs, nrow = single_rep * nrow(env_pairs))
  for(i in 1:nrow(env_pairs)){
    design_df2[((i-1)*single_rep + 1):(single_rep*i),env_pairs[i,]] <- 1
  }
  
  design_df <- rbind(design_df1, design_df2)
  design_df <- design_df[sample(1:nrow(design_df)), ]
  
  # print actual p-rep level, only if target not achieved
  # check this is <= 1
  prep_observed <- round(double_rep/(nplots_env - double_rep), 2)
  if(prep_observed != prep) {
    warning(paste0("p-rep level of ", prep_observed," obtained, istead of ", prep))
  }
  # print dropped number of indivduals / plots
  # print observed no. individuals per environment
  ninds_env <- nplots_env - double_rep # check this is <= ninds
  if(ninds_env != ninds) {
    warning(paste0(ninds_env, " individuals per environment out of ", ninds, " individuals in total (", round(ninds_env/ninds, 2), ")"))
  }
  # total plots
  plots_total_obs <- sum(design_df)
  if(plots_total_obs != ninds * nplots) {
    warning(paste0(plots_total_obs, " plots out of ", ninds * nplots, " used (", plots_total_obs/nenvs ," per environment)"))
  }
  # colSums(design_df)
  # rowSums(design_df)
  design_df <- as.table(design_df)
  colnames(design_df) <- 1:nenvs
  rownames(design_df) <- 1:ninds
  design_df <- as.data.frame(t(design_df))
  colnames(design_df) <- c("env", "id", "nreps")
  design_df$env <- factor(as.numeric(as.character(design_df$env)))
  design_df$id <- factor(as.numeric(as.character(design_df$id)))
  design_df <- design_df[order(design_df$env, design_df$id),]
  return(design_df)
}
