# Create Figure 13 in appendix.
# Plot quantiles of permutation test from Cule et al. (2010).
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
results <- fread("sim_data/fig13_perm_test_stats.csv")

# Data frame of quantiles and original test stats
quantile_df <- results %>% 
  group_by(d, mu_norm, B, orig_test_stat) %>% 
  dplyr::summarise(q90 = quantile(shuffle_test_stat, 0.90),
                   q95 = quantile(shuffle_test_stat, 0.95),
                   q99 = quantile(shuffle_test_stat, 0.99),
                   sim_count = n()) %>% 
  mutate(details = 
           factor(paste(d, mu_norm, B),
                  levels = c("1 0 100", "5 0 100", "1 2 100", "5 2 100",
                             "1 0 500", "5 0 500", "1 2 500", "5 2 500"),
                  labels = c("d==1 ~~~~ ll*mu*ll==0 ~~~~ B==100", 
                             "d==5 ~~~~ ll*mu*ll==0 ~~~~ B==100",
                             "d==1 ~~~~ ll*mu*ll==2 ~~~~ B==100", 
                             "d==5 ~~~~ ll*mu*ll==2 ~~~~ B==100",
                             "d==1 ~~~~ ll*mu*ll==0 ~~~~ B==500", 
                             "d==5 ~~~~ ll*mu*ll==0 ~~~~ B==500",
                             "d==1 ~~~~ ll*mu*ll==2 ~~~~ B==500",
                             "d==5 ~~~~ ll*mu*ll==2 ~~~~ B==500")))

# Check parameters
stopifnot(quantile_df$B == quantile_df$sim_count)

# Plot shuffled permutation test statistics, quantiles, and orig test statistics 
perm_test_quantiles <- results %>% 
  mutate(details = 
           factor(paste(d, mu_norm, B),
                  levels = c("1 0 100", "5 0 100", "1 2 100", "5 2 100",
                             "1 0 500", "5 0 500", "1 2 500", "5 2 500"),
                  labels = c("d==1 ~~~~ ll*mu*ll==0 ~~~~ B==100", 
                             "d==5 ~~~~ ll*mu*ll==0 ~~~~ B==100",
                             "d==1 ~~~~ ll*mu*ll==2 ~~~~ B==100", 
                             "d==5 ~~~~ ll*mu*ll==2 ~~~~ B==100",
                             "d==1 ~~~~ ll*mu*ll==0 ~~~~ B==500", 
                             "d==5 ~~~~ ll*mu*ll==0 ~~~~ B==500",
                             "d==1 ~~~~ ll*mu*ll==2 ~~~~ B==500",
                             "d==5 ~~~~ ll*mu*ll==2 ~~~~ B==500"))) %>% 
  ggplot(aes(x = shuffle_test_stat, y = ..density..)) +
  geom_histogram(fill = "lightskyblue", bins = 30) +
  facet_wrap(~ details, nrow = 4, labeller = label_parsed) +
  geom_vline(aes(xintercept = orig_test_stat), data = quantile_df,
             col = "black") + 
  geom_text(aes(x = orig_test_stat - 0.01, y = 15, label = "Test stat"), 
            data = quantile_df, angle = 90) +
  geom_vline(aes(xintercept = q90), data = quantile_df, 
             col = "blue", lty = "dashed", alpha = 0.5) +
  geom_text(aes(x = q90, y = 20, label = "q['0.90']"), data = quantile_df,
            parse = T, angle = 0) +
  geom_vline(aes(xintercept = q95), data = quantile_df, 
             col = "blue", lty = "dashed", alpha = 0.5) +
  geom_text(aes(x = q95, y = 13, label = "q[0.95]"), data = quantile_df,
            parse = T, angle = 0) +
  geom_vline(aes(xintercept = q99), data = quantile_df, 
             col = "blue", lty = "dashed", alpha = 0.5) +
  geom_text(aes(x = q99, y = 6, label = "q[0.99]"), data = quantile_df,
            parse = T, angle = 0) +
  labs(x = "Test statistic on shuffled data",
       y = "Count",
       title = "Distribution of shuffled data test statistics in eight simulations.",
       subtitle =  expression(atop("Permutation test for H"[0]*
                                ": Log-concave vs H"[1]*": Not log-concave.", 
                                "Normal location family f(x) = 0.5"*phi[d]*"(x)"~
                                  "+ 0.5"*phi[d]*"(x -"~mu*")."))) +
  paper_theme

#####################
##### Save plot #####
#####################

ggsave(plot = perm_test_quantiles,
       filename = "sim_plots/figure_13.pdf",
       width = 8, height = 7)
