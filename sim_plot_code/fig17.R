# Plot results from various tests of
# H_0: log-concave versus H_1: not log-concave.
# Use beta distribution. Density is log-concave iff alpha_param, beta_param >= 1.

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
perm_test <- fread("sim_data/fig17_perm_test.csv") %>% 
  mutate(Method = "Permutation test (B = 99 shuffles)")

full_oracle <- fread("sim_data/fig17_full_oracle.csv") %>% 
  mutate(Method = "Full oracle (B = 100)")

partial_oracle <- fread("sim_data/fig17_partial_oracle.csv") %>% 
  mutate(Method = "Partial oracle (B = 100)")

fully_NP <- fread("sim_data/fig17_fully_NP.csv") %>% 
  mutate(Method = "Fully nonparametric (B = 100)")

# Combine results
results <- rbind(perm_test, full_oracle, partial_oracle, fully_NP, fill = TRUE)

# Get rejection proportion at each (d, mu_norm) combination
reject_df <- results %>% 
  group_by(n_obs, d, beta_param, alpha_param, B, Method) %>% 
  dplyr::summarise(reject_prop = mean(reject),
                   sim_count = n())

# Check parameters
stopifnot(unique(perm_test$B) == 99,
          unique(full_oracle$B) == 100,
          unique(partial_oracle$B) == 100,
          unique(fully_NP$B) == 100,
          unique(reject_df$n_obs) == 100,
          unique(reject_df$sim_count) == 200)

######################
##### Plot power #####
######################

# Plot rejection proportions of permutation and log-concave mixture methods
reject_plot_beta <- reject_df %>% 
  mutate(log_concave = case_when(alpha_param >= 1 & beta_param >= 1 ~ "Log-concave",
                                 alpha_param <  1 | beta_param < 1~ "Not log-concave"),
         Method = 
           factor(Method,
                  levels = c("Permutation test (B = 99 shuffles)",
                             "Full oracle (B = 100)",
                             "Partial oracle (B = 100)",
                             "Fully nonparametric (B = 100)"),
                  labels = c("Permutation test (B = 99 shuffles)",
                             "Full oracle (B = 100)",
                             "Partial oracle (B = 100)",
                             "Fully nonparametric (B = 100)")),
         beta_param = factor(beta_param,
                             levels = c(0.5, 1, 2),
                             labels = c("Shape~parameter~beta==0.5", 
                                        "Shape~parameter~beta==1", 
                                        "Shape~parameter~beta==2"))) %>%
  ggplot(aes(x = alpha_param, y = reject_prop)) +
  geom_line(aes(col = Method), alpha = 0.6) +
  geom_point(aes(col = Method, shape = log_concave), size = 2, alpha = 0.6) +
  geom_hline(yintercept = 0.10, lty = "dashed", col = "darkgrey") +
  facet_grid(. ~ beta_param, labeller = label_parsed) +
  labs(x = expression("Shape parameter"~alpha), 
       y = "Rejection proportion",
       col = "Method",
       shape = "Log-concave?",
       title = expression("Permutation test and universal tests for H"[0]*
                            ": Log-concave vs H"[1]*": Not log-concave"),
       subtitle = expression("Beta distribution. 200 sims. n = 100.")) +
  scale_shape_manual(values = c(16, 4)) +
  scale_color_manual(values = c("black", "red", "orange", "blue")) +
  paper_theme +
  theme(legend.position = "bottom") +
  guides(shape = guide_legend(nrow = 2),
         colour = guide_legend(nrow = 4))

#####################
##### Save plot #####
#####################

ggsave(plot = reject_plot_beta,
       filename = "sim_plots/figure_17.pdf",
       width = 10, height = 4.8)
