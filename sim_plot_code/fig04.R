# Create Figure 4 in main paper.
# Plot results for tests of
# H_0: log-concave versus H_1: not log-concave.
# Use normal mixture where mean of second component is 
# mu = -(||mu||, 0, ..., 0).

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
perm_test <- fread("sim_data/fig04_perm_test.csv") %>% 
  mutate(Method = "Permutation test")

full_oracle_ddim <- fread("sim_data/fig04_full_oracle_ddim.csv") %>% 
  mutate(Method = "Full oracle, d-dim")

partial_oracle_ddim <- fread("sim_data/fig04_partial_oracle_ddim.csv") %>% 
  mutate(Method = "Partial oracle, d-dim")

partial_oracle_axis <- fread("sim_data/fig04_partial_oracle_axis.csv") %>% 
  mutate(Method = "Partial oracle, axis-aligned projections")

fully_NP_axis <- fread("sim_data/fig04_fully_NP_axis.csv") %>% 
  mutate(Method = "Fully nonparametric, axis-aligned projections")

partial_oracle_randproj <- fread("sim_data/fig04_partial_oracle_randproj.csv") %>% 
  mutate(Method = "Partial oracle, random projections")

fully_NP_randproj <- fread("sim_data/fig04_fully_NP_randproj.csv") %>% 
  mutate(Method = "Fully nonparametric, random projections")

# Combine results
results <- rbind(perm_test, full_oracle_ddim, partial_oracle_ddim, 
                 partial_oracle_axis, fully_NP_axis, partial_oracle_randproj, 
                 fully_NP_randproj,
                 fill = TRUE)

# Get rejection proportion at each combination
reject_df <- results %>% 
  group_by(d, mu_norm, Method) %>% 
  dplyr::summarise(reject_prop = mean(reject),
                   sim_count = n())

# Check parameters
stopifnot(unique(perm_test$B) == 99,
          unique(full_oracle_ddim$B) == 100,
          unique(partial_oracle_ddim$B) == 100,
          unique(partial_oracle_axis$B) == 100,
          unique(fully_NP_axis$B) == 100,
          unique(partial_oracle_randproj$B) == 100,
          unique(fully_NP_randproj$B) == 100,
          unique(partial_oracle_randproj$n_proj == 100),
          unique(fully_NP_randproj$n_proj == 100),
          unique(results$n_obs == 100),
          unique(results$equal_space_mu) == 0,
          unique(reject_df$sim_count) == 200)

######################
##### Plot power #####
######################

# Plot rejection proportions at (||mu||, 0, ..., 0)
reject_props_unequal_spaced <- reject_df %>% 
  mutate(Method = factor(Method,
                         levels = c("Partial oracle, axis-aligned projections",
                                    "Partial oracle, random projections",
                                    "Fully nonparametric, axis-aligned projections",
                                    "Fully nonparametric, random projections",
                                    "Full oracle, d-dim",
                                    "Partial oracle, d-dim",
                                    "Permutation test"),
                         labels = c("Partial oracle, axis-aligned projections",
                                    "Partial oracle, random projections",
                                    "Fully nonparametric, axis-aligned projections",
                                    "Fully nonparametric, random projections",
                                    "Full oracle, d-dim",
                                    "Partial oracle, d-dim",
                                    "Permutation test")),
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
  labs(x = expression("||"*mu*"|| in normal location family, where"~
                        mu=="-(||"*mu*"||, 0, ..., 0)"), 
       y = "Rejection proportion",
       col = "Method",
       lty = "Method",
       title = expression("Tests for H"[0]*
                            ": Log-concave vs H"[1]*": Not log-concave"),
       subtitle = expression("Normal location family f(x) = 0.5"*phi[d]*"(x)"~
                               "+ 0.5"*phi[d]*"(x -"~mu*"). n = 100 obs. 200 sims.")) +
  scale_color_manual(values = c("#0047b3", "#66a3ff", "#b30000", "#ff8080", 
                                "#5200cc", "#b380ff", "black")) +
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid",
                                   "longdash", "longdash", "dashed")) +
  paper_theme +
  theme(legend.position = "bottom", 
        legend.key.width = unit(2.5, "line")) +
  guides(colour = guide_legend(nrow = 2))

# Plot rejection proportions at (||mu||, 0, ..., 0) for single dim approaches,
# zooming in further on ||mu|| \in [4, 8]
reject_props_unequal_zoom_onerow <- reject_df %>% 
  filter(mu_norm >= 4, mu_norm <= 8,
         Method %in% c("Partial oracle, axis-aligned projections",
                       "Partial oracle, random projections",
                       "Fully nonparametric, axis-aligned projections",
                       "Fully nonparametric, random projections")) %>% 
  mutate(Method = factor(Method,
                         levels = c("Partial oracle, axis-aligned projections",
                                    "Partial oracle, random projections",
                                    "Fully nonparametric, axis-aligned projections",
                                    "Fully nonparametric, random projections"),
                         labels = c("Partial oracle, axis-aligned projections",
                                    "Partial oracle, random projections",
                                    "Fully nonparametric, axis-aligned projections",
                                    "Fully nonparametric, random projections")),
         d = factor(d, levels = 1:4, 
                    labels = c("d = 1", "d = 2", "d = 3", "d = 4"))) %>%
  ggplot(aes(x = mu_norm, y = reject_prop)) +
  geom_line(aes(col = Method), alpha = 0.7) +
  geom_point(aes(col = Method, shape = Method), alpha = 0.8) +
  geom_hline(yintercept = 0.10, lty = "dashed", col = "darkgrey") +
  facet_grid(. ~ d) +
  labs(x = expression("||"*mu*"|| in normal location family, where"~
                        mu=="-(||"*mu*"||, 0, ..., 0)"), 
       y = "Rejection proportion",
       col = "Method",
       shape = "Method",
       title = NULL,
       subtitle = NULL) +
  scale_color_manual(values = c("#0047b3", "#66a3ff", "#b30000", "#ff8080")) +
  paper_theme +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2))

######################
##### Save plots #####
######################

ggsave(plot = reject_props_unequal_spaced,
       filename = "sim_plots/figure_04a.pdf",
       width = 13, height = 7)

ggsave(plot = reject_props_unequal_zoom_onerow,
       filename = "sim_plots/figure_04b.pdf",
       width = 12, height = 3.3)
