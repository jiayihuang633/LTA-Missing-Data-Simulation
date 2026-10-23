library(dplyr)

# ------------------------------------------------------------
# 1) Output folder
# ------------------------------------------------------------

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

# ------------------------------------------------------------
# 2) Helper function for label alignment
# ------------------------------------------------------------

align_fit <- function(fit) {
  
  # Get P(Y = 1 | state)
  psi_y1 <- fit$Psi[2, , , drop = FALSE]
  
  # Average across indicators
  state_means <- apply(psi_y1, 2, mean)
  
  if (state_means[1] < state_means[2]) {
    fit$Psi <- fit$Psi[, c(2, 1), , drop = FALSE]
    fit$Pi  <- fit$Pi[c(2, 1), c(2, 1), , drop = FALSE]
    fit$Piv <- fit$Piv[c(2, 1)]
  }
  
  return(fit)
}

# ------------------------------------------------------------
# 3) Function to combine one method's results
# ------------------------------------------------------------

combine_results <- function(root_dir, output_file) {
  
  files <- list.files(
    root_dir,
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  cat("Processing", root_dir, "-", length(files), "files found\n")
  
  results_list <- vector("list", length(files))
  
  for (i in seq_along(files)) {
    
    file <- files[i]
    fit  <- readRDS(file)
    
    # Extract replication number from filename
    rep_num <- as.integer(
      gsub("rep|\\.rds", "", basename(file))
    )
    
    # Failed model fit
    if (is.na(fit$Pi)[1]) {
      
      results_list[[i]] <- data.frame(
        rep = rep_num,
        n = fit$n,
        T = fit$T,
        p = fit$p,
        miss_mech = fit$miss_mech,
        miss_rate = fit$miss_rate,
        converged = FALSE
      )
      
    } else {
      
      # Align latent-class labels
      fit <- align_fit(fit)
      
      # Extract parameters
      piv <- fit$Piv
      Pi  <- fit$Pi[, , 2]
      Psi <- fit$Psi
      
      results_list[[i]] <- data.frame(
        rep = rep_num,
        n = fit$n,
        T = fit$T,
        p = fit$p,
        miss_mech = fit$miss_mech,
        miss_rate = fit$miss_rate,
        
        # Initial probabilities
        piv1 = piv[1],
        piv2 = piv[2],
        
        # Transition probabilities
        pi11 = Pi[1, 1],
        pi12 = Pi[1, 2],
        pi21 = Pi[2, 1],
        pi22 = Pi[2, 2],
        
        # Item-response probabilities
        psi1 = mean(Psi[2, 1, ]),
        psi2 = mean(Psi[2, 2, ]),
        
        logLik = fit$logLik,
        BIC = fit$BIC,
        AIC = fit$AIC,
        converged = TRUE
      )
    }
  }
  
  results_df <- bind_rows(results_list)
  
  saveRDS(
    results_df,
    file.path(results_dir, output_file)
  )
  
  cat("Saved:", file.path(results_dir, output_file), "\n")
}

# ------------------------------------------------------------
# 4) Combine results for all three methods
# ------------------------------------------------------------

combine_results(
  root_dir = "analysis",
  output_file = "all_results_FIML.rds"
)

combine_results(
  root_dir = "analysis_MI",
  output_file = "all_results_MI.rds"
)

combine_results(
  root_dir = "analysis_SI",
  output_file = "all_results_SI.rds"
)

cat("Finished combining all results.\n")