# Create Figure 6 in appendix.
# Plot log-concave MLE densities for normal mixtures at n = 50 and d = 1.

# Read in library
library(LogConcaveUniv)

# Create theme
paper_theme <- ggplot2::theme_bw() +
  ggplot2::theme(plot.title    = ggplot2::element_text(hjust = 0.5, size = 16),
                 plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 14),
                 legend.title  = ggplot2::element_text(size = 14),
                 axis.title    = ggplot2::element_text(size = 14),
                 legend.text   = ggplot2::element_text(size = 12),
                 axis.text     = ggplot2::element_text(size = 12),
                 strip.text    = ggplot2::element_text(size = 12),
                 panel.spacing = ggplot2::unit(1.2, "lines"))

# Read in data
density_points <- data.table::fread("sim_data/fig06_points.csv")

density_values <- data.table::fread("sim_data/fig06_densities.csv")

# Check parameters
stopifnot(unique(density_values$d) == 1,
          unique(density_values$n_obs) == 50)

# Extract (1/n)*loglik values
loglik_df <- density_values %>%
  dplyr::group_by(mu_norm) %>%
  dplyr::slice(1) %>%
  tidyr::pivot_longer(cols = c("mean_loglik_true_dens",
                               "mean_loglik_LogConcDEAD",
                               "mean_loglik_logcondens")) %>%
  dplyr::mutate(name = factor(name,
                              levels = c("mean_loglik_true_dens",
                                         "mean_loglik_LogConcDEAD",
                                         "mean_loglik_logcondens"),
                              labels = c("True density",
                                         "Log-concave MLE (LogConcDEAD)",
                                         "Log-concave MLE (logcondens)")),
                mu_norm = factor(mu_norm, levels = c(0, 2, 4),
                                 labels = c("||u|| = 0", "||u|| = 2", "||u|| = 4")))

# Plot densities
logconc_densities_n50_d1 <- density_values %>%
  tidyr::pivot_longer(cols = c("true_density", "LogConcDEAD_density",
                               "logcondens_density")) %>%
  dplyr::mutate(name = factor(name, levels = c("true_density",
                                               "LogConcDEAD_density",
                                               "logcondens_density"),
                              labels = c("True density",
                                         "Log-concave MLE (LogConcDEAD)",
                                         "Log-concave MLE (logcondens)")),
                mu_norm = factor(mu_norm, levels = c(0, 2, 4),
                                 labels = c("||u|| = 0", "||u|| = 2", "||u|| = 4"))) %>%
  ggplot2::ggplot(ggplot2::aes(x = x, y = value)) +
  ggplot2::geom_line() +
  ggplot2::facet_grid(mu_norm ~ name, scales = "free") +
  ggplot2::geom_rug(ggplot2::aes(x = x), sides = "b", inherit.aes = F,
                    data = density_points) +
  ggplot2::geom_text(ggplot2::aes(x = -4.5, y = 0.35,
                                  label = paste("\n", "(1/n)*loglik =",
                                                round(value, 2))),
                     data = loglik_df %>% dplyr::filter(d == 1, n_obs == 50)) +
  paper_theme +
  ggplot2::labs(y = "Density",
                title = "True density and log-concave MLE estimates, n = 50, d = 1")

#####################
##### Save plot #####
#####################

ggplot2::ggsave(plot = logconc_densities_n50_d1,
                filename = "sim_plots/figure_06.pdf",
                width = 9.7, height = 5)
