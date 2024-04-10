# Create Figure 3 in main paper.
# Plot results from full oracle Gaussian mixture/universal approach and 
# Cule et al. (2010) permutation test of
# H_0: log-concave versus H_1: not log-concave.

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(gtable))
suppressMessages(library(cowplot))

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
perm_test <- fread("sim_data/fig03_perm_test.csv") %>% 
  mutate(Method = "Permutation test (B = 99 shuffles)")

full_oracle <- fread("sim_data/fig03_full_oracle_ddim.csv") %>% 
  mutate(Method = "Full oracle, d-dim (B = 100)")

# Combine results
results <- rbind(perm_test, full_oracle, fill = TRUE)

# Get rejection proportion at each (d, mu_norm) combination
reject_df <- results %>% 
  group_by(n_obs, d, mu_norm, B, Method) %>% 
  dplyr::summarise(reject_prop = mean(reject),
                   sim_count = n())

# Check parameters
stopifnot(unique(perm_test$B) == 99,
          unique(full_oracle$B) == 100,
          unique(results$equal_space_mu) == 0,
          unique(results$n_obs) == 100,
          sort(unique(results$d)) == 1:4,
          unique(reject_df$sim_count) == 200)

######################
##### Plot power #####
######################

# Plot rejection proportions of permutation and log-concave mixture methods
reject_perm_true_n100 <- reject_df %>% 
  mutate(Method = factor(Method,
                         levels = c("Permutation test (B = 99 shuffles)",
                                    "Full oracle, d-dim (B = 100)"),
                         labels = c("Permutation test (B = 99 shuffles)",
                                    "Full oracle, d-dim (B = 100)")),
         d = factor(d, levels = 1:4, 
                    labels = c("d = 1", "d = 2", "d = 3", "d = 4"))) %>%
  ggplot(aes(x = mu_norm, y = reject_prop)) +
  geom_line(aes(col = Method, lty = Method), alpha = 0.7) +
  geom_point(aes(col = Method), alpha = 0.3) +
  geom_hline(yintercept = 0.10, lty = "dashed", col = "darkgrey") +
  geom_vline(xintercept = 2, lty = "dashed", col = "darkgrey") +
  annotate(geom = "text", x = 1, y = 0.7, label = "LC", angle = 90) +
  annotate(geom = "text", x = 3, y = 0.7, label = "Not LC", angle = 90) +
  facet_wrap(. ~ d) +
  labs(x = expression("||"*mu*"|| in normal location family"), 
       y = "Rejection proportion",
       col = "Method",
       lty = "Method",
       shape = "Log-concave?",
       title = expression("Permutation test and true density universal test for H"[0]*
                            ": Log-concave vs H"[1]*": Not log-concave"),
       subtitle = expression("Normal location family f(x) = 0.5"*phi[d]*"(x)"~
                               "+ 0.5"*phi[d]*"(x -"~mu*"). n = 100 obs. 200 sims.")) +
  scale_shape_manual(values = c(16, 4)) +
  scale_color_manual(values = c("black", "#5200cc")) +
  paper_theme +
  theme(legend.position = "bottom", 
        legend.key.width = unit(2, "line")) +
  guides(colour = guide_legend(nrow = 1))

#####################
##### Save plot #####
#####################

ggsave(plot = reject_perm_true_n100,
       filename = "sim_plots/figure_03.pdf",
       width = 10, height = 6)
