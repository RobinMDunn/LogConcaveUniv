# Create Figures 16 and 17 in appendix.
# Get densities to plot appearance of log-concave MLEs from
# non-log-concave beta densities.

suppressMessages(library(logcondens))
suppressMessages(library(MASS))
suppressMessages(library(data.table))
suppressMessages(library(progress))
suppressMessages(library(tidyverse))
suppressMessages(library(ReIns))

# Create data frame to store results
results <- data.table(d = 1,
                      x = rep(seq(0.001, 0.999, length.out = 999), times = 2),
                      alpha = c(rep(0.5, 999), rep(0.5, 999)),
                      beta = c(rep(0.5, 999), rep(1, 999)),
                      n_obs = 100000,
                      true_density = NA_real_,
                      logcondens_density = NA_real_)

results[alpha == 0.5 & beta == 0.5, sim := 1]
results[alpha == 0.5 & beta == 1, sim := 2]

# Run simulations to check whether to reject H_0
for(sim_val in 1:2) {
  
  # Print sim_val
  print(sim_val)
  
  # Set seed
  set.seed(sim_val)
  
  # Extract n_obs, alpha, beta
  n_obs <- results[sim == sim_val, n_obs][1]
  alpha_param <- results[sim == sim_val, alpha][1]
  beta_param <- results[sim == sim_val, beta][1]
  
  # Generate sample from beta density
  true_sample <- rbeta(n = n_obs, shape1 = alpha_param, shape2 = beta_param)
  
  # Get MLE log-concave estimate using logcondens package
  mle_logcondens <- logConDens(x = true_sample)
  
  # Extract x and y coords
  x_coord <- results[sim == sim_val, x]
  
  # Evaluate true beta density
  results[sim == sim_val,
          true_density := dbeta(x = x_coord, 
                                shape1 = alpha_param, shape2 = beta_param)]
  
  # Evaluate MLE log-concave estimate using logcondens package
  results[sim == sim_val,
          logcondens_density := evaluateLogConDens(
            xs = x_coord, res = mle_logcondens, which = 2)[, 3]]

}

# Estimated truncated exponential density
results[sim == 2, guess := dtexp(x = x, rate = 2.18, endpoint = 1)]

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

# Plot true and logcondens densities
true_and_logcondens <- results %>% 
  pivot_longer(cols = c("true_density", "logcondens_density")) %>% 
  mutate(name = factor(name, levels = c("true_density", "logcondens_density"),
                       labels = c("True~density", 
                                  "Log-concave~MLE~(logcondens)")),
         sim = factor(sim, levels = c(1, 2),
                      labels = c("alpha==0.5~~~beta==0.5",
                                 "alpha==0.5~~~beta==1"))) %>% 
  ggplot(aes(x = x, y = value)) + 
  geom_line() +
  facet_wrap(sim ~ name, labeller = label_parsed, 
             ncol = 2, scales = "free") +
  expand_limits(y = c(0.9, 1.1)) +
  paper_theme +
  labs(y = "Density",
       title = expression("True densities and log-concave MLE estimates for"~
                            "Beta("*alpha*","~beta*") densities"))

# Plot logcondens and estimate for Beta(0.5, 1) sim
truncated_exp <- results %>% 
  filter(alpha == 0.5, beta == 1) %>% 
  ggplot(aes(x = x, y = logcondens_density)) + 
  geom_line() +
  geom_line(aes(x = x, y = guess), lty = "dashed", col = "red") +
  paper_theme +
  labs(y = "Density",
       title = "Log-concave MLE estimate for Beta(0.5, 1)",
       subtitle = "Black curve: logcondens. Dashed red curve: truncated exponential.")

######################
##### Save plots #####
######################

ggsave(plot = true_and_logcondens,
       filename = "sim_plots/figure_16.pdf",
       width = 8, height = 6)

ggsave(plot = truncated_exp,
       filename = "sim_plots/figure_17.pdf",
       width = 7, height = 4)
