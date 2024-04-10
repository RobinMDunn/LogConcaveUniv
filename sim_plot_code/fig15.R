# Test statistics of partial oracle tests in d = 2 dimensions or
# projected to d = 1 dimension.
# True density is (1/2)N((0,0), I_2) + (1/2)N((-6,0), I_2).

library(tidyverse)
library(data.table)

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
test_stats_d1_proj <- fread(file = "sim_data/fig15_test_stats_d1_proj.csv") %>% 
  dplyr::mutate(test_dim = "Projection to d = 1")

test_stats_d2 <- fread(file = "sim_data/fig15_test_stats_d2.csv") %>% 
  dplyr::mutate(test_dim = "Full dimension (d = 2)")

# Combine test stat data
test_stats <- rbind(test_stats_d1_proj, test_stats_d2, fill = TRUE)

# Distribution of log-concave log test statistics
d1_d2_test_stats <- ggplot(test_stats, aes(x = log(test_stat))) +
  geom_histogram(bins = 60, fill = "lightblue", col = "black") +
  facet_grid(test_dim ~ .) +
  geom_vline(aes(xintercept = log(10)), lty = "dashed", col = "blue") +
  annotate(geom = "text", x = log(10) + 8, y = 162, 
           label = "log(1/0.1)", col = "blue", size = 6) +
  paper_theme +
  labs(x = "Log of test statistic",
       y = "Count",
       title = expression("Partial oracle test statistics for H"[0]*
                            ": Log-concave vs H"[1]*": Not log-concave"),
       subtitle = expression("Normal location family f(x) = 0.5"*phi[2]*"(x)"~
                               "+ 0.5"*phi[2]*"(x + (6,0)"**"T"*")."~
                               "n = 100 obs. 1000 sims."))

# Rejection % of two-dimensional test: 4.5%
mean(test_stats_d2$test_stat > log(1/0.1)) * 100

# Rejection % of one-dimensional projection test: 44.3%
mean(test_stats_d1_proj$test_stat > log(1/0.1)) * 100

# 19% of projection test stats are above 1000
mean(test_stats_d1_proj$test_stat > 1000) * 100

# 65% of the projection test statistics with abs(random_vector_1) > 0.9
# are above 1000
mean(test_stats_d1_proj[abs(random_vector_1) > 0.9]$test_stat > 1000) * 100

# Probability that first component of 2-d random projection exceeds 0.9,
# when drawing projection vectors uniformly from boundary of the unit circle.
# Probability is 0.29.
(2/pi)*acos(.9)

#####################
##### Save plot #####
#####################

ggsave(plot = d1_d2_test_stats, 
       filename = "sim_plots/figure_15.pdf",
       width = 8, height = 6)
       