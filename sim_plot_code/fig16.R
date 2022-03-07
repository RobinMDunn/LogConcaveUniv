# Create Figure 16 in appendix.
# Plot beta densities under variety of alpha and beta parameters.

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

# Evaluate beta densities 

x_seq <- seq(0.01, 0.99, by = 0.01)

beta_df <- data.table(expand.grid(x = x_seq,
                                  alpha_param = c(0.1, 0.5, 1, 1.5, 2),
                                  beta_param = c(0.5, 1, 2)),
                      density_eval = NA_real_)

beta_df[, density_eval := dbeta(x = x, shape1 = alpha_param, shape2 = beta_param),
         by = seq_len(nrow(beta_df))]

# Plot beta densities

beta_densities <- beta_df %>% 
  mutate(log_concave = case_when(alpha_param <  1 | beta_param < 1 ~ 
                                   "Not log-concave",
                                 alpha_param >= 1 & beta_param >= 1 ~
                                   "Log-concave")) %>% 
  mutate(log_concave = factor(log_concave, 
                              levels = c("Not log-concave", "Log-concave"),
                              labels = c("Not log-concave", "Log-concave"))) %>% 
  ggplot(aes(x = x, y = density_eval, col = log_concave)) +
  facet_wrap(. ~ beta_param + alpha_param,
             labeller = label_bquote(cols = alpha==.(alpha_param)*","~
                                       beta==.(beta_param)), nrow = 3,
             scales = "free") +
  geom_line() +
  labs(y = "Density",
       col = "Log-concave?",
       title = expression("Beta("*alpha*","~beta*") densities with varying"~alpha~
                            "and"~beta~"parameters")) +
  scale_color_manual(values = c("red", "black")) +
  paper_theme +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10))

#####################
##### Save plot #####
#####################

ggsave(plot = beta_densities,
       filename = "sim_plots/figure_16.pdf",
       width = 10, height = 5.5)
