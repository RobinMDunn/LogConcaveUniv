# Create Figure 14 in appendix.
# Plot results from Gaussian mixture/universal approach and 
# Cule et al. (2010) permutation test of
# H_0: log-concave versus H_1: not log-concave.
# Use true density as numerator.

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
full_oracle <- fread("sim_data/fig14_full_oracle.csv")

# Get rejection proportion at each (d, mu_norm) combination
reject_props <- full_oracle %>% 
  group_by(n_obs, d, mu_norm, B) %>% 
  dplyr::summarise(reject_prop = mean(reject),
                   sim_count = n())

# Check correct parameters
stopifnot(unique(full_oracle$B) == 100, 
          unique(full_oracle$equal_space_mu) == 0,
          unique(full_oracle$n_obs) == 100,
          unique(reject_props$sim_count) == 200)

######################################################################
##### Plot mu_norm necessary for power of approx 0.90, varying d #####
######################################################################

power_pt90_df_vary_d <- reject_props %>% 
  group_by(d) %>% 
  filter(abs(reject_prop - 0.90) == min(abs(reject_prop - 0.90))) %>% 
  slice(1) %>% 
  ungroup()

power_vary_d <- power_pt90_df_vary_d %>% 
  mutate(approx_power = cut(reject_prop, 
                            breaks = c(0.88, 0.89, 0.90, 0.91, 0.92),
                            include.lowest = TRUE)) %>%
  mutate(slope = lm(mu_norm ~ I(exp(1)^d))$coefs[2]) %>% 
  ggplot(aes(x = d, y = mu_norm)) +
  scale_y_continuous(limits = c(0, 90)) +
  geom_line(col = "grey") +
  geom_point(aes(col = approx_power), size = 3) +  
  scale_color_brewer(palette = "RdBu", drop = FALSE) +
  stat_smooth(geom = "line", method = "lm", formula = y ~ I(exp(1)^x), 
              col = "blue", se = F,
              lty = "dashed", alpha = 0.5, size = 1) +
  labs(y = expression("||"*mu*"||"),
       col = "Approx power",
       title = expression("How does ||"*mu*"|| need to grow with dimension"~
                            "for power of 0.90?"),
       subtitle = expression("With best fit"~
                               hat("||"*mu*"||")~"= a + b*"*exp(d)~"curve."~
                               "n = 100 obs. 200 sims.")) +
  paper_theme

#####################
##### Save plot #####
#####################

ggsave(plot = power_vary_d,
       filename = "sim_plots/figure_14.pdf",
       width = 10, height = 4)
