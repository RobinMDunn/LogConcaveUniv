# Create Figure XX in main paper.
# Plot results for tests of
# H_0: log-concave versus H_1: not log-concave.
# Density is 0.5*N(0, I_2) + 0.5*N(0, sigma^2 I_2), where sigma = sqrt(3).
# Density is not LC, but 1-d projections are LC.

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

# Create theme
paper_theme <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        panel.spacing = unit(1.2, "lines"))

# Read in data
full_oracle_ddim <- fread("sim_data/new_example/full_oracle_ddim.csv") %>% 
  dplyr::mutate(Method = "Full oracle, d-dim")

partial_oracle_ddim <- fread("sim_data/new_example/partial_oracle_ddim.csv") %>% 
  dplyr::mutate(Method = "Partial oracle, d-dim")

partial_oracle_randproj <- fread("sim_data/new_example/partial_oracle_randproj.csv") %>% 
  dplyr::mutate(Method = "Partial oracle, random projections")

fully_NP_randproj <- fread("sim_data/new_example/fully_NP_randproj.csv") %>% 
  dplyr::mutate(Method = "Fully nonparametric, random projections")

# Combine results
results <- rbind(full_oracle_ddim, partial_oracle_ddim, 
                 partial_oracle_randproj, fully_NP_randproj,
                 fill = TRUE)

# Get rejection proportion at each combination
reject_df <- results %>% 
  dplyr::mutate(reject_prop = n_reject / n_sim)

# Check parameters
stopifnot(unique(full_oracle_ddim$B) == 100,
          unique(partial_oracle_ddim$B) == 100,
          unique(partial_oracle_randproj$B) == 100,
          unique(fully_NP_randproj$B) == 100,
          unique(partial_oracle_randproj$n_proj == 100),
          unique(fully_NP_randproj$n_proj == 100),
          unique(reject_df$n_sim) == 200)

######################
##### Plot power #####
######################

# Plot rejection proportions
reject_props <- reject_df %>% 
  dplyr::mutate(Method = factor(Method,
                         levels = c("Partial oracle, random projections",
                                    "Fully nonparametric, random projections",
                                    "Full oracle, d-dim",
                                    "Partial oracle, d-dim"),
                         labels = c("Partial oracle, random projections",
                                    "Fully nonparametric, random projections",
                                    "Full oracle, d-dim",
                                    "Partial oracle, d-dim"))) %>%
  ggplot(aes(x = n_obs, y = reject_prop)) +
  geom_line(aes(col = Method), alpha = 0.7) +
  geom_point(aes(col = Method), alpha = 0.3) +
  geom_hline(yintercept = 0.10, lty = "dashed", col = "darkgrey") +
  labs(x = "Number of observations", 
       y = "Rejection proportion",
       col = "Method",
       title = expression("Tests for H"[0]*
                            ": Log-concave vs H"[1]*": Not log-concave"),
       subtitle = expression("Normal mixture 0.5*(N(0,"~I[2]*") +"*
                             "N(0,"~sigma**2~I[2]*")), where"~sigma~"="~sqrt(3))) +
  scale_color_manual(values = c("#66a3ff", "#ff8080", 
                                "#5200cc", "#b380ff")) +
  paper_theme +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  scale_x_continuous(breaks = seq(200, 1600, by = 200)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1))
  
######################
##### Save plots #####
######################

ggsave(plot = reject_props,
       filename = "sim_plots/new_example.pdf",
       width = 7, height = 5)
