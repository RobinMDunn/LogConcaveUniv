# Create table of average time to run simulations

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(tidyr))
suppressMessages(library(xtable))

# Read in data
perm_test <- fread("sim_data/tab01_perm_test.csv")

full_oracle_ddim <- fread("sim_data/tab01_full_oracle_ddim.csv")

partial_oracle_ddim <- fread("sim_data/tab01_partial_oracle_ddim.csv")

partial_oracle_axis <- fread("sim_data/tab01_partial_oracle_axis.csv")

fully_NP_axis <- fread("sim_data/tab01_fully_NP_axis.csv")

partial_oracle_randproj <- fread("sim_data/tab01_partial_oracle_randproj.csv")

fully_NP_randproj <- fread("sim_data/tab01_fully_NP_randproj.csv")

# Combine data
all_data <- rbind(perm_test, full_oracle_ddim, partial_oracle_ddim,
                  partial_oracle_axis, fully_NP_axis, partial_oracle_randproj,
                  fully_NP_randproj, fill = TRUE)

all_data <- all_data %>% 
  mutate(Method = factor(Method, 
                         levels = c("Partial oracle, random projections",
                                    "Fully nonparametric, random projections",
                                    "Permutation test",
                                    "Partial oracle, axis-aligned projections",
                                    "Fully nonparametric, axis-aligned projections",
                                    "Partial oracle, d-dim",
                                    "Full oracle, d-dim"),
                         labels = c("Partial oracle, random projections",
                                    "Fully nonparametric, random projections",
                                    "Permutation test",
                                    "Partial oracle, axis-aligned projections",
                                    "Fully nonparametric, axis-aligned projections",
                                    "Partial oracle, d-dim",
                                    "Full oracle, d-dim")))

# Get average time
time_df <- all_data %>% 
  group_by(Method, d, mu_norm) %>% 
  dplyr::summarise(avg_time = mean(time_sec))

# Create LaTeX table
time_df %>% 
  pivot_wider(id_cols = c("Method", "d"), names_from = "mu_norm",
              values_from = "avg_time", names_prefix = "mu_norm_") %>% 
  mutate(avg_vec = paste0("(", signif(mu_norm_0, 2), ", ", 
                          signif(mu_norm_5, 2), ", ", 
                          signif(mu_norm_10, 2), ")")) %>% 
  dplyr::select(Method, d, avg_vec) %>% 
  pivot_wider(id_cols = "Method", names_from = "d", values_from = "avg_vec",
              names_prefix = "d = ") %>% 
  xtable()
