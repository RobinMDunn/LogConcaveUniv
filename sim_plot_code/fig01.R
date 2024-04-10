# Create Figure 1 in main paper.
# Plot results from Cule et al. (2010) permutation test
# and fully nonparametric random projection universal test of
# H_0: log-concave versus H_1: not log-concave.
# Use two-component normal location model.

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
perm_test <- fread("sim_data/fig01_perm_test.csv") %>% 
  mutate(Method = "Permutation test")

fully_NP_randproj <- fread("sim_data/fig01_fully_NP_randproj.csv") %>% 
  mutate(Method = "Universal test (fully nonparametric, random projections)")

results <- rbind(perm_test, fully_NP_randproj, fill = TRUE)

# Get rejection proportion at each (d, mu_norm) combination
reject_df <- results %>% 
  group_by(Method, n_obs, d, mu_norm, B) %>% 
  dplyr::summarise(reject_prop = mean(reject),
                   sim_count = n())

# Check parameters
stopifnot(unique(results$equal_space_mu) == 0,
          unique(perm_test$B) == 99,
          unique(fully_NP_randproj$B) == 100,
          unique(results$n_obs) == 100,
          sort(unique(results$mu_norm)) == 0:10,
          unique(reject_df$sim_count) == 200)

# Plot rejection proportions
perm_randproj_tests <- reject_df %>% 
  mutate(Method = factor(Method,
                         levels = c("Permutation test",
                                    "Universal test (fully nonparametric, random projections)"),
                         labels = c("Permutation test",
                                    "Universal test (fully nonparametric, random projections)")),
         d = factor(d, levels = 1:5, 
                    labels = c("d = 1", "d = 2", "d = 3", 
                               "d = 4", "d = 5"))) %>%
  ggplot(aes(x = mu_norm, y = reject_prop)) +
  facet_wrap(. ~ d) +
  geom_hline(yintercept = 0.10, lty = "dashed", col = "darkgrey") +
  geom_vline(xintercept = 2, lty = "dashed", col = "darkgrey") +
  geom_line(aes(col = Method, lty = Method), alpha = 0.7) +
  geom_point(aes(col = Method), alpha = 0.7) +
  annotate(geom = "text", x = 1.5, y = 0.7, label = "LC", angle = 90) +
  annotate(geom = "text", x = 2.5, y = 0.7, label = "Not LC", angle = 90) +
  labs(x = expression("||"*mu*"|| in normal location family"), 
       y = "Rejection proportion",
       col = "",
       lty = "",
       title = expression("Tests for H"[0]*
                            ": Log-concave vs H"[1]*": Not log-concave"),
       subtitle = expression("Normal location family f(x) = 0.5"*phi[d]*"(x)"~
                               "+ 0.5"*phi[d]*"(x -"~mu*"). n = 100.")) +
  scale_color_manual(values = c("black", "red")) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  paper_theme +
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "line"))

#####################
##### Save plot #####
#####################

ggsave(plot = perm_randproj_tests,
       filename = "sim_plots/figure_01.pdf",
       width = 9, height = 5.5)
