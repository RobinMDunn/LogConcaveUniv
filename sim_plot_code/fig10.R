# Create Figure 10 in appendix.
# Plot log-concave MLE densities for normal mixtures at n = 500 and d = 2.

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
density_values <- fread("sim_data/fig10_densities.csv")

# Check parameters
stopifnot(unique(density_values$d) == 2,
          unique(density_values$n_obs) == 500)

# Extract (1/n)*loglik values
loglik_df <- density_values %>% 
  group_by(mu_norm) %>% 
  slice(1) %>% 
  pivot_longer(cols = c("mean_loglik_true_dens", "mean_loglik_LogConcDEAD",
                        "mean_loglik_logcondens")) %>% 
  mutate(name = factor(name, 
                       levels = c("mean_loglik_true_dens", 
                                  "mean_loglik_LogConcDEAD",
                                  "mean_loglik_logcondens"),
                       labels = c("True density", 
                                  "Log-concave MLE (LogConcDEAD)",
                                  "Log-concave MLE (logcondens)")),
         mu_norm = factor(mu_norm, levels = c(0, 2, 4),
                          labels = c("||u|| = 0", "||u|| = 2", "||u|| = 4")))

# Plot contours at n = 500 and d = 2
logconc_contours_n500_d2 <- density_values %>% 
  dplyr::select(-logcondens_density, -mean_loglik_logcondens) %>% 
  pivot_longer(cols = c("true_density", "LogConcDEAD_density")) %>% 
  mutate(name = factor(name, levels = c("true_density", "LogConcDEAD_density"),
                       labels = c("True density", 
                                  "Log-concave MLE (LogConcDEAD)")),
         mu_norm = factor(mu_norm, levels = c(0, 2, 4),
                          labels = c("||u|| = 0", "||u|| = 2", "||u|| = 4"))) %>% 
  ggplot(aes(x = x, y = y, z = value)) + 
  stat_contour(aes(colour = ..level..), binwidth = 0.01) +
  facet_grid(mu_norm ~ name, scales = "free") +
  scale_color_gradient(low = "#035efc", high = "red") +
  geom_text(aes(x = -4.5, y = 3.2, 
                label = paste("\n", "(1/n)*loglik =", round(value, 2))), 
            data = loglik_df %>% 
              filter(name != "Log-concave MLE (logcondens)")) +
  paper_theme +
  labs(color = "Density",
       title = "True density and log-concave MLE estimates, n = 500, d = 2")

#####################
##### Save plot #####
#####################

ggsave(plot = logconc_contours_n500_d2,
       filename = "sim_plots/figure_10.pdf",
       width = 8.5, height = 5)
