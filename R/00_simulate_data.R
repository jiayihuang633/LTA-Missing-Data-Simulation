############################################################
## Basic LTA Simulation (2 classes, binary indicators)
############################################################

set.seed(123)

# ------------------------------------------------------------
# 1) Helper functions
# ------------------------------------------------------------

# Emission probabilities (the probability that Indicator_i = 1, given the specific class)
# rows = classes, cols = indicators
make_emission_probs <- function(k = 2, p,
                                high = 0.80,
                                low = 0.30,
                                jitter = 0.10,
                                seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  stopifnot(k == 2)
  
  offset <- runif(p, -jitter, jitter)
  
  class1 <- high + offset
  class2 <- low  + offset
  
  class1 <- pmin(pmax(class1, 0.05), 0.95)
  class2 <- pmin(pmax(class2, 0.05), 0.95)
  
  rbind(class1, class2)
}

# ------------------------------------------------------------
# 2.1) MCAR Missingness
# ------------------------------------------------------------
  
apply_MCAR <- function(dat, miss_rate) {
  
  dat_miss <- dat
  
  # Only apply to observed indicator variables
  y_cols <- grep("^y", names(dat), value = TRUE)
  
  # To match MAR, make sure it is not possible to have t1 indicators missing
  y_cols <- y_cols[!grepl("_t1$", y_cols)]
  
  for (col in y_cols) {
    miss_index <- runif(nrow(dat)) < miss_rate
    dat_miss[miss_index, col] <- NA
  }
  
  return(dat_miss)
}

# ------------------------------------------------------------
# 2.2) MAR Missingness (depends on same indicator at t-1)
# ------------------------------------------------------------
apply_MAR <- function(dat, miss_rate, T, p, gamma = 1) {
  
  dat_miss <- dat
  
  for (t in 2:T) {
    
    for (j in 1:p) {
      
      current_col  <- paste0("y", j, "_t", t)
      previous_col <- paste0("y", j, "_t", t-1)
      
      prev_values <- dat[[previous_col]]
      
      # Linear predictor without intercept
      lin_no_int <- gamma * prev_values
      
      # Solve for intercept alpha so that:
      # mean(plogis(alpha + lin_no_int)) = miss_rate
      
      f <- function(alpha) {
        mean(plogis(alpha + lin_no_int)) - miss_rate
      }
      
      alpha_hat <- uniroot(f, interval = c(-10, 10))$root
      
      miss_prob <- plogis(alpha_hat + lin_no_int)
      
      miss_index <- runif(nrow(dat)) < miss_prob
      
      dat_miss[miss_index, current_col] <- NA
    }
  }
  
  return(dat_miss)
}


# ------------------------------------------------------------
# 3) single_run(): one simulation run
# ------------------------------------------------------------
single_run <- function(n, T, p, k = 2,
                       pi_k,
                       trans_mat,
                       miss_mech = c("none", "MCAR", "MAR"),
                       miss_rate = 0,
                       gamma = 1,
                       seed = NULL) {
  
  miss_mech <- match.arg(miss_mech)
  if (!is.null(seed)) set.seed(seed)
  
  # Emission probs must be k x p
  emission_probs <- make_emission_probs(k = k, p = p)
  
  # Generate latent class trajectories
  Z <- matrix(NA_integer_, nrow = n, ncol = T)
  Z[, 1] <- sample(1:k, n, replace = TRUE, prob = pi_k)
  
  for (t in 2:T) {
    for (i in 1:n) {
      Z[i, t] <- sample(1:k, 1, prob = trans_mat[Z[i, t-1], ])
    }
  }
  
  # Generate binary indicators (y{j}_t{t})
  dat <- data.frame(id = 1:n)
  
  for (t in 1:T) {
    for (j in 1:p) {
      prob_vec <- emission_probs[Z[, t], j]
      dat[[paste0("y", j, "_t", t)]] <- rbinom(n, size = 1, prob = prob_vec)
    }
  }
  
  # Save true latent states (class_t{t})
  for (t in 1:T) {
    dat[[paste0("class_t", t)]] <- Z[, t]
  }
  
  # Apply missingness
  dat_obs <- dat
  if (miss_mech == "MCAR" && miss_rate > 0) {
    dat_obs <- apply_MCAR(dat_obs, miss_rate)
  }
  if (miss_mech == "MAR" && miss_rate > 0) {
    dat_obs <- apply_MAR(dat_obs, miss_rate, T = T, p = p, gamma = gamma)
  }
  
  # Ensure miss_rate = 0 when miss_mech == "none"
  dat_obs$N <- n
  dat_obs$TT <- T
  dat_obs$P <- p
  dat_obs$miss_mech <- miss_mech
  dat_obs$miss_rate <- ifelse(miss_mech == "none", 0, miss_rate)
  
  dat_obs
}

# ------------------------------------------------------------
# 4) Loop over conditions and save outputs
# ------------------------------------------------------------

# Choose vectors
n_vec <- c(100, 200, 500, 1000)
T_vec <- c(3, 5, 10)
p_vec <- c(1, 5, 10)
miss_rates <- c(0.10, 0.20, 0.40)  # non-zero rates
R <- 100                     # repetitions per condition

# Fixed model parameters
pi_k <- c(0.5, 0.5)
trans_mat <- matrix(c(
  0.80, 0.20,
  0.25, 0.75
), nrow = 2, byrow = TRUE)

# Parameter grid
grid <- expand.grid(
  n = n_vec,
  T = T_vec,
  p = p_vec,
  miss_mech = c("none", "MCAR", "MAR"),
  miss_rate = c(0, miss_rates),
  rep = 1:R,
  KEEP.OUT.ATTRS = FALSE
)

# Keep only valid pairs:
# - none must have rate 0
# - MCAR/MAR must have rate > 0
grid <- subset(
  grid,
  (miss_mech == "none" & miss_rate == 0) |
    (miss_mech != "none" & miss_rate > 0)
)

results <- vector("list", nrow(grid))

# Run and save
for (r in seq_len(nrow(grid))) {
  g <- grid[r, ]
  
  dat_run <- single_run(
    n = g$n, T = g$T, p = g$p, k = 2,
    pi_k = pi_k,
    trans_mat = trans_mat,
    miss_mech = as.character(g$miss_mech),
    miss_rate = g$miss_rate,
    gamma = 1,
    seed = 10000 + r
  )
  
  results[[r]] <- dat_run
}

# ------------------------------------------------------------
# 5) Run simulations and save individual datasets
# ------------------------------------------------------------

data_dir <- "data"
dir.create(data_dir, showWarnings = FALSE)

for (r in seq_len(nrow(grid))) {
  
  g <- grid[r, ]
  
  dat_run <- single_run(
    n = g$n,
    T = g$T,
    p = g$p,
    k = 2,
    pi_k = pi_k,
    trans_mat = trans_mat,
    miss_mech = as.character(g$miss_mech),
    miss_rate = g$miss_rate,
    gamma = 1,
    seed = 10000 + r
  )
  
  condition_dir <- file.path(
    data_dir,
    paste0(
      "n", g$n,
      "_T", g$T,
      "_p", g$p,
      "_", g$miss_mech,
      "_", g$miss_rate
    )
  )
  
  dir.create(condition_dir, recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(
    dat_run,
    file.path(condition_dir, paste0("rep", g$rep, ".rds"))
  )
}

cat("Finished generating simulation datasets.\n")

