# ============================================================
# LTA result processing + revised paper plots
# FIML vs MI vs Single Stochastic Imputation (SI)
# ALL parameters plotted using absolute relative bias
# Uses ggplot2 facet_grid for stable saving (no ggh4x dependency)
# ============================================================

# Set project-relative paths for results and figures.

library(dplyr)
library(tidyr)
library(ggplot2)

results_dir <- "results"
figures_dir <- "figures"

dir.create(results_dir, showWarnings = FALSE)
dir.create(figures_dir, showWarnings = FALSE)

# ------------------------------------------------------------
# 0) Read result files
# ------------------------------------------------------------

all_results_FIML <- readRDS(file.path(results_dir, "all_results_FIML.rds"))
all_results_MI   <- readRDS(file.path(results_dir, "all_results_MI.rds"))
all_results_SI   <- readRDS(file.path(results_dir, "all_results_SI.rds"))

# ------------------------------------------------------------
# 1) True parameter values
# ------------------------------------------------------------

true_values <- c(
  piv1 = 0.50,
  piv2 = 0.50,
  pi11 = 0.80,
  pi12 = 0.20,
  pi21 = 0.25,
  pi22 = 0.75,
  psi1 = 0.80,
  psi2 = 0.30
)

param_cols <- names(true_values)
required_cols <- c("n", "T", "p", "miss_mech", "miss_rate", param_cols)

check_result_file <- function(dat, method_name) {
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop(
      method_name,
      " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

check_result_file(all_results_FIML, "FIML")
check_result_file(all_results_MI,   "MI")
check_result_file(all_results_SI,   "SI")

# ------------------------------------------------------------
# 2) Completeness check
# Expected: 100 rows/replications per condition
# ------------------------------------------------------------

condition_counts <- bind_rows(
  all_results_FIML %>% mutate(method = "FIML"),
  all_results_MI   %>% mutate(method = "MI"),
  all_results_SI   %>% mutate(method = "SI")
) %>%
  count(method, n, T, p, miss_mech, miss_rate, name = "n_rows")

cat("\nConditions not containing exactly 100 result rows:\n")
print(condition_counts %>% filter(n_rows != 100), n = Inf)

write.csv(
  condition_counts,
  file.path(results_dir, "condition_counts_FIML_MI_SI.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 3) Relative bias for each replication
# ------------------------------------------------------------

make_iteration_bias <- function(dat, method_name) {
  dat %>%
    mutate(
      iteration_id = row_number(),
      method = method_name
    ) %>%
    pivot_longer(
      cols = all_of(param_cols),
      names_to = "parameter",
      values_to = "estimate"
    ) %>%
    mutate(
      true_value = unname(true_values[parameter]),
      rel_bias = (estimate - true_value) / true_value
    )
}

rel_bias_each_iteration_FIML <- make_iteration_bias(all_results_FIML, "FIML")
rel_bias_each_iteration_MI   <- make_iteration_bias(all_results_MI,   "MI")
rel_bias_each_iteration_SI   <- make_iteration_bias(all_results_SI,   "SI")

rel_bias_each_iteration <- bind_rows(
  rel_bias_each_iteration_FIML,
  rel_bias_each_iteration_MI,
  rel_bias_each_iteration_SI
)

# ------------------------------------------------------------
# 4) Mean relative bias within each simulation condition
# ------------------------------------------------------------

rel_bias_by_condition <- rel_bias_each_iteration %>%
  group_by(method, n, T, p, miss_mech, miss_rate, parameter) %>%
  summarise(
    mean_rel_bias = mean(rel_bias, na.rm = TRUE),
    sd_rel_bias   = sd(rel_bias, na.rm = TRUE),
    n_iterations  = sum(!is.na(rel_bias)),
    .groups = "drop"
  )

saveRDS(
  rel_bias_each_iteration,
  file.path(results_dir, "rel_bias_each_iteration_FIML_MI_SI.rds")
)

saveRDS(
  rel_bias_by_condition,
  file.path(results_dir, "rel_bias_by_condition_FIML_MI_SI.rds")
)

write.csv(
  rel_bias_by_condition,
  file.path(results_dir, "rel_bias_by_condition_FIML_MI_SI.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5) Prepare data for all revised plots
# ------------------------------------------------------------

all_plot_params <- c(
  "piv1", "piv2",
  "pi11", "pi12", "pi21", "pi22",
  "psi1", "psi2"
)

# All parameters are plotted as the absolute value of the
# condition-specific mean relative bias.
abs_bias_params <- all_plot_params

plot_data_all <- rel_bias_by_condition %>%
  filter(
    parameter %in% all_plot_params,
    miss_mech %in% c("MAR", "MCAR"),
    miss_rate %in% c(0.1, 0.2, 0.4)
  ) %>%
  mutate(
    plot_bias = abs(mean_rel_bias),
    n = factor(n, levels = c(100, 200, 500, 1000)),
    T = factor(T, levels = c(3, 5, 10)),
    p = factor(p, levels = c(1, 5, 10)),
    miss_rate = factor(
      miss_rate,
      levels = c(0.1, 0.2, 0.4),
      labels = c("10%", "20%", "40%")
    ),
    miss_mech = factor(
      miss_mech,
      levels = c("MAR", "MCAR")
    ),
    method = factor(
      method,
      levels = c("FIML", "MI", "SI")
    ),
    parameter = factor(
      parameter,
      levels = all_plot_params
    )
  )

# ------------------------------------------------------------
# 6) MAR-only paper plots
# Rows: missingness rate (10%, 20%, 40%)
# Major columns: sample size (100, 200, 500, 1000)
# Subcolumns: number of indicators (1, 5, 10)
# Within each panel: T = 3, 5, 10; lines = FIML, MI, SI
# ------------------------------------------------------------

make_mar_plot <- function(param_name) {
  metric_label <- "Absolute Relative Bias"

  this_data <- plot_data_all %>%
    filter(
      parameter == param_name,
      miss_mech == "MAR"
    )

  ggplot(
    this_data,
    aes(
      x = T,
      y = plot_bias,
      color = method,
      group = method
    )
  ) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.7) +
    facet_grid(
      rows = vars(miss_rate),
      cols = vars(n, p),
      scales = "fixed"
    ) +
    labs(
      title = paste(metric_label, "of", param_name, "under MAR"),
      x = "Number of Measurement Occasions",
      y = metric_label,
      color = "Method"
    ) +
    theme_minimal() +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(size = 8, face = "bold"),
      strip.text.y = element_text(size = 9, face = "bold"),
      panel.spacing.x = unit(0.35, "lines"),
      panel.spacing.y = unit(0.6, "lines"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray85"),
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
}

# ------------------------------------------------------------
# 7) MAR + MCAR comparison plots
# Major rows: missingness rate (10%, 20%, 40%)
# Subrows within each rate: MAR and MCAR
# Major columns: sample size; subcolumns: indicators
# ------------------------------------------------------------

make_mar_mcar_plot <- function(param_name) {
  metric_label <- "Absolute Relative Bias"

  this_data <- plot_data_all %>%
    filter(
      parameter == param_name,
      miss_mech %in% c("MAR", "MCAR")
    )

  ggplot(
    this_data,
    aes(
      x = T,
      y = plot_bias,
      color = method,
      group = method
    )
  ) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.7) +
    facet_grid(
      rows = vars(miss_rate, miss_mech),
      cols = vars(n, p),
      scales = "fixed"
    ) +
    labs(
      title = paste(metric_label, "of", param_name, "under MAR and MCAR"),
      x = "Number of Measurement Occasions",
      y = metric_label,
      color = "Method"
    ) +
    theme_minimal() +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(size = 8, face = "bold"),
      strip.text.y = element_text(size = 8, face = "bold"),
      panel.spacing.x = unit(0.35, "lines"),
      panel.spacing.y = unit(0.45, "lines"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray85"),
      axis.text.x = element_text(size = 6.5),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
}

# ------------------------------------------------------------
# 8) Create and save every MAR-only and MAR+MCAR plot directly
# This avoids storing ggplot objects in lists, which can cause
# compatibility problems in some ggplot2/ggh4x setups.
# ------------------------------------------------------------

for (param_name in all_plot_params) {

  # MAR-only plot
  p_mar <- make_mar_plot(param_name)

  ggsave(
    filename = file.path(
      figures_dir,
      paste0("MAR_", param_name, ".png")
    ),
    plot = p_mar,
    width = 18,
    height = 8,
    dpi = 300
  )

  # MAR + MCAR plot
  p_both <- make_mar_mcar_plot(param_name)

  ggsave(
    filename = file.path(
      figures_dir,
      paste0("BOTH_", param_name, ".png")
    ),
    plot = p_both,
    width = 18,
    height = 12,
    dpi = 300
  )

  cat("Saved:", param_name, "\n")
}



# ------------------------------------------------------------
# 9) Main-paper MAR vs MCAR comparison for psi2
# Collapse across methods; lines compare MAR vs MCAR
# Rows: missingness rate (10%, 20%, 40%)
# Major columns: sample size; subcolumns: indicators
# Within each panel: T = 3, 5, 10
# ------------------------------------------------------------

psi2_mar_mcar_data <- plot_data_all %>%
  filter(parameter == "psi2") %>%
  group_by(n, T, p, miss_rate, miss_mech) %>%
  summarise(
    plot_bias = mean(plot_bias, na.rm = TRUE),
    .groups = "drop"
  )

p_psi2_mar_mcar <- ggplot(
  psi2_mar_mcar_data,
  aes(
    x = T,
    y = plot_bias,
    color = miss_mech,
    group = miss_mech
  )
) +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.9) +
  facet_grid(
    rows = vars(miss_rate),
    cols = vars(n, p),
    scales = "fixed"
  ) +
  labs(
    title = expression(paste("Absolute Relative Bias of ", psi[2], " under MAR and MCAR")),
    x = "Number of Measurement Occasions",
    y = "Absolute Relative Bias",
    color = "Missingness Mechanism"
  ) +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 8, face = "bold"),
    strip.text.y = element_text(size = 9, face = "bold"),
    panel.spacing.x = unit(0.35, "lines"),
    panel.spacing.y = unit(0.6, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray85"),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

print(p_psi2_mar_mcar)

ggsave(
  file.path(figures_dir, "MAR_MCAR_psi2_collapsed_methods.png"),
  p_psi2_mar_mcar,
  width = 18,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# 10) No-missingness baseline plot (kept from prior workflow)
# ------------------------------------------------------------

baseline_params <- c("piv1", "pi22", "psi2")

none_data <- rel_bias_by_condition %>%
  filter(
    method == "FIML",
    miss_mech == "none",
    parameter %in% baseline_params
  ) %>%
  mutate(
    plot_bias = abs(mean_rel_bias),
    n = factor(n, levels = c(100, 200, 500, 1000)),
    T = factor(T, levels = c(3, 5, 10)),
    p = factor(p, levels = c(1, 5, 10)),
    parameter = factor(parameter, levels = baseline_params)
  )

p_none <- ggplot(
  none_data,
  aes(
    x = T,
    y = plot_bias,
    group = 1
  )
) +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  facet_grid(
    rows = vars(parameter),
    cols = vars(n, p),
    scales = "free_y"
  ) +
  labs(
    title = "Absolute Relative Bias under No Missingness",
    x = "Number of Measurement Occasions",
    y = "Absolute Relative Bias"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 8, face = "bold"),
    strip.text.y = element_text(size = 8, face = "bold"),
    panel.border = element_blank(),
    panel.spacing.x = unit(0.5, "lines"),
    panel.spacing.y = unit(2, "lines"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_none)

ggsave(
  file.path(figures_dir, "RB_none_FIML_baseline.png"),
  p_none,
  width = 14,
  height = 10,
  dpi = 300
)

cat("\nDone.\n")
cat("All MAR-only and MAR+MCAR plots were saved in: ",
    normalizePath(figures_dir), "\n")
