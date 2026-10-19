library(LMest)
library(tidyr)
library(dplyr)

# ------------------------------------------------------------
# 1) Get folder path from command line
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No folder path supplied")
}

folder <- args[1]

cat("Processing folder:", folder, "\n")

# ------------------------------------------------------------
# 2) List all .rds files in this folder
# ------------------------------------------------------------

files <- list.files(folder, pattern = "\\.rds$", full.names = TRUE)

cat("Found", length(files), "files\n")

# ------------------------------------------------------------
# 3) Loop over datasets in this folder
# ------------------------------------------------------------

for (i in seq_along(files)) {
  
  file <- files[i]
  cat(sprintf("[%d/%d] %s\n", i, length(files), file))
  
  # ----------------------------
  # Read data
  # ----------------------------
  dat <- readRDS(file)
  
  # ----------------------------
  # Remove true classes
  # ----------------------------
  dat_no_class <- dat %>%
    select(-starts_with("class_"))
  
  # ----------------------------
  # Reshape
  # ----------------------------
  dat_long <- dat_no_class %>%
    pivot_longer(
      cols = starts_with("y"),
      names_to = c("indicator", "time"),
      names_pattern = "y([0-9]+)_t([0-9]+)"
    ) %>%
    mutate(
      indicator = paste0("y", indicator),
      time = as.integer(time)
    ) %>%
    arrange(id, time)
  
  dat_long <- dat_long %>%
    pivot_wider(
      names_from = indicator,
      values_from = value
    ) %>%
    arrange(id, time)
  
  dat_long <- as.data.frame(dat_long)
  
  cat("---- CHECK DATA ----\n")
cat("dim:", dim(dat_long), "\n")
cat("unique ids:", length(unique(dat_long$id)), "\n")
cat("time table:\n")
print(table(dat_long$time))
  
  # ----------------------------
  # Build formula
  # ----------------------------
  y_cols <- grep("^y[0-9]+$", names(dat_long), value = TRUE)
 
     dat_long_fit <- dat_long %>%
      select(id, time, all_of(y_cols)) %>%
      as.data.frame()
  
  cat("Total NA in data:", sum(is.na(dat_long[, y_cols])), "\n")
  # ----------------------------
  # Fit model
  # ----------------------------
  fit <- tryCatch({
  
  LMest::lmest(
    responsesFormula = NULL,
    data = dat_long_fit,
    index = c("id", "time"),
    k = 2,
    modBasic = 1,
    start = 1, 
    ntry = 5 
  )
  
}, error = function(e) {
  message("Error in file: ", file)
  message("Actual error: ", e$message)
  
  return(list(
    Pi = NA,
    Piv = NA,
    Psi = NA,
    logLik = NA,
    BIC = NA,
    AIC = NA,
    n = unique(dat$N),
    T = unique(dat$TT),
    p = unique(dat$P),
    miss_mech = unique(dat$miss_mech),
    miss_rate = unique(dat$miss_rate)
  ))
})

  # ------------------------------------------------------------
  # 4) Save result
  # ------------------------------------------------------------
  
  # Extract quantities of interest
res <- list(
  Pi      = fit$Pi,      # transition matrix
  Piv     = fit$piv,     # initial probs
  Psi     = fit$Psi,     # emission probs
  logLik  = fit$lk,
  BIC    = fit$bic,
  AIC    = fit$aic,
  
  # Optional: condition metadata (if present in dataset)
  n       = unique(dat$N),
  T       = unique(dat$TT),
  p       = unique(dat$P),
  miss_mech  = unique(dat$miss_mech),
  miss_rate  = unique(dat$miss_rate)
)

  out_file <- sub("^data", "analysis", file)
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(res, file = out_file)
}

cat("Finished folder:", folder, "\n")