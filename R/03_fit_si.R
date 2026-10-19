library(LMest)
library(tidyr)
library(dplyr)
library(mice)

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

# ------------------------------------------------------------
# Helper Function for Label Switching
# ------------------------------------------------------------

align_fit <- function(fit) {
  
  if (is.null(fit)) return(NULL)
  
  # Get P(Y=1 | state) safely for any p
  psi_y1 <- fit$Psi[2, , , drop = FALSE]
  
  # Collapse over items if multiple
  state_means <- apply(psi_y1, 2, mean)
  
  if (state_means[1] < state_means[2]) {
    fit$Psi <- fit$Psi[, c(2,1), , drop=FALSE]
    fit$Pi  <- fit$Pi[c(2,1), c(2,1), , drop=FALSE]
    fit$piv <- fit$piv[c(2,1)]
  }
  
  return(fit)
}

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
  # SI Step (Wide format)
  # ----------------------------
  
  y_wide <- dat_no_class %>%
    select(starts_with("y"))%>%
    mutate(across(everything(), as.factor))
  
  meth <- make.method(y_wide)
  meth[] <- "logreg"
  
  pred <- make.predictorMatrix(y_wide)
  
M <- 1

imp <- mice(
  y_wide,
  m = M,
  method = meth,
  predictorMatrix = pred,
  maxit = 5,
  printFlag = FALSE,
  seed = 10000 + i
)

# ----------------------------
# Single Stochastic Imputation
# ----------------------------

dat_complete <- complete(imp, 1)

dat_complete$id <- dat_no_class$id

dat_complete <- dat_complete %>%
  mutate(across(starts_with("y"), ~ as.numeric(as.character(.))))

dat_long_SI <- dat_complete %>%
  pivot_longer(
    cols = starts_with("y"),
    names_to = c("indicator", "time"),
    names_pattern = "y([0-9]+)_t([0-9]+)"
  ) %>%
  mutate(
    indicator = paste0("y", indicator),
    time = as.integer(time)
  ) %>%
  arrange(id, time) %>%
  pivot_wider(
    names_from = indicator,
    values_from = value
  ) %>%
  arrange(id, time) %>%
  select(id, time, starts_with("y")) %>%
  as.data.frame()

# ----------------------------
# Fit LTA to single imputed dataset
# ----------------------------

fit_SI <- tryCatch({
  LMest::lmest(
    responsesFormula = NULL,
    data = dat_long_SI,
    index = c("id", "time"),
    k = 2,
    modBasic = 1,
    start = 1,
    ntry = 5
  )
}, error = function(e) NULL)

# Align latent class labels
fit_SI <- align_fit(fit_SI)
  
# ----------------------------
# Save SI estimates
# ----------------------------

if (is.null(fit_SI)) {
  
  res <- list(
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
  )
  
} else {
  
  res <- list(
    Pi = fit_SI$Pi,
    Piv = fit_SI$piv,
    Psi = fit_SI$Psi,
    logLik = fit_SI$lk,
    BIC = fit_SI$bic,
    AIC = fit_SI$aic,
    n = unique(dat$N),
    T = unique(dat$TT),
    p = unique(dat$P),
    miss_mech = unique(dat$miss_mech),
    miss_rate = unique(dat$miss_rate)
  )
}

# ----------------------------
# Save result
# ----------------------------

out_file <- sub("^data", "analysis_SI", file)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

saveRDS(res, file = out_file)

}  # closes for (i in seq_along(files))

cat("Finished folder:", folder, "\n")