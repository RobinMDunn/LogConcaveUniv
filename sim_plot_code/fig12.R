# Create Figure 12 in appendix.
# Plot results from Cule et al. (2010) permutation test of
# H_0: log-concave versus H_1: not log-concave.
# Use two-component normal location model.

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
results <- fread("sim_data/fig12_perm_test.csv")

# Get rejection proportion at each (d, mu_norm) combination
reject_df <- results %>% 
  group_by(n_obs, d, mu_norm, B) %>% 
  dplyr::summarise(reject_prop = mean(reject),
                   sim_count = n())

# Check parameters
stopifnot(unique(results$equal_space_mu) == 0,
          unique(results$n_obs) == 100,
          sort(unique(results$mu_norm)) == 0:5,
          unique(reject_df$sim_count) == 200)

# Plot rejection proportions at varying numbers of shuffles (B), n = 100.
perm_test_reject_B_vary <- reject_df %>% 
  mutate(d = factor(d, levels = 1:5, 
                    labels = c("d = 1", "d = 2", "d = 3", "d = 4", "d = 5")),
         log_concave = 
           factor(mu_norm, levels = 0:5, 
                  labels = c("Log-concave", "Log-concave", "Log-concave", 
                             "Not log-concave", "Not log-concave", 
                             "Not log-concave")),
         B = factor(B, levels = c(100, 200, 300, 400, 500),
                    labels = c("B = 100", "B = 200", "B = 300", 
                               "B = 400", "B = 500"))) %>% 
  ggplot(aes(x = mu_norm, y = reject_prop, col = log_concave)) +
  facet_grid(d ~ B) +
  geom_line(color = "darkgrey") +
  geom_point() +
  geom_hline(yintercept = 0.10, lty = "dashed") +
  labs(x = expression("||"*mu*"|| in normal location family"), 
       y = "Rejection proportion",
       col = "",
       title = expression("Permutation test for H"[0]*
                            ": Log-concave vs H"[1]*": Not log-concave"),
       subtitle = expression("Normal location family f(x) = 0.5"*phi[d]*"(x)"~
                               "+ 0.5"*phi[d]*"(x -"~mu*"). Varying B. n = 100.")) +
  scale_color_manual(values = c("blue", "red")) +
  paper_theme +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10))

#####################
##### Save plot #####
#####################

ggsave(plot = perm_test_reject_B_vary,
       filename = "sim_plots/figure_12.pdf",
       width = 8, height = 7)
