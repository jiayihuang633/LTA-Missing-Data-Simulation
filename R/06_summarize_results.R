library(dplyr)

rel_bias_by_condition <- read.csv(
  file.path("results", "rel_bias_by_condition_FIML_MI_SI.csv")
)

# Summary of absolute relative bias by method and parameter
# Using absolute value

overall_method_summary <- rel_bias_by_condition %>%
  filter(miss_mech %in% c("MCAR", "MAR")) %>%
  mutate(abs_rel_bias = abs(mean_rel_bias)) %>%
  group_by(method, parameter) %>%
  summarise(
    mean_abs_rel_bias = mean(abs_rel_bias, na.rm = TRUE),
    median_abs_rel_bias = median(abs_rel_bias, na.rm = TRUE),
    max_abs_rel_bias = max(abs_rel_bias, na.rm = TRUE),
    .groups = "drop"
  )

overall_method_summary

write.csv(
  overall_method_summary,
  file.path("results", "summary_by_method_parameter.csv"),
  row.names = FALSE
)

# Summary by missing-data rate

missing_rate_summary <- rel_bias_by_condition %>%
  filter(miss_mech %in% c("MCAR", "MAR")) %>%
  mutate(abs_rel_bias = abs(mean_rel_bias)) %>%
  group_by(method, parameter, miss_mech, miss_rate) %>%
  summarise(
    mean_abs_rel_bias = mean(abs_rel_bias, na.rm = TRUE),
    median_abs_rel_bias = median(abs_rel_bias, na.rm = TRUE),
    .groups = "drop"
  )

missing_rate_summary

write.csv(
  missing_rate_summary,
  file.path("results", "summary_by_missing_rate.csv"),
  row.names = FALSE
)

# Summary by design factors
# by sample size, number of indicators, and measurement occasions

design_summary <- rel_bias_by_condition %>%
  filter(miss_mech %in% c("MCAR", "MAR")) %>%
  mutate(abs_rel_bias = abs(mean_rel_bias)) %>%
  group_by(method, parameter, n, T, p) %>%
  summarise(
    mean_abs_rel_bias = mean(abs_rel_bias, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  design_summary,
  file.path("results", "summary_by_design.csv"),
  row.names = FALSE
)

design_summary

