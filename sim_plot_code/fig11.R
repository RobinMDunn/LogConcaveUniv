# Create Figure 11 in appendix.
# Plot results from Cule et al. (2010) permutation test of
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

# Function to position legend into wasted space.
# Source: https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2
shift_legend <- function(p){
  
  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }
  
  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")
  
  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")
  
  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")
  
  return(gp)
}

# Read in data
results <- fread("sim_data/fig11_perm_test.csv")

# Get rejection proportion at each (d, mu_norm) combination
reject_df <- results %>% 
  filter(equal_space_mu == 0) %>% 
  group_by(n_obs, d, mu_norm, B) %>% 
  dplyr::summarise(reject_prop = mean(reject),
                   sim_count = n())

# Check parameters
stopifnot(unique(results$equal_space_mu) == 0,
          unique(results$B) == 99,
          unique(results$n_obs) == 250,
          sort(unique(results$mu_norm)) == 0:5,
          unique(reject_df$sim_count) == 200)

# Plot rejection proportions at B = 99 shuffles, n = 250.
perm_test_reject_n250 <- reject_df %>% 
  mutate(d = factor(d, levels = 1:5, 
                    labels = c("d = 1", "d = 2", "d = 3", "d = 4", "d = 5")),
         log_concave = 
           factor(mu_norm, levels = 0:5, 
                  labels = c("Log-concave", "Log-concave", "Log-concave", 
                             "Not log-concave", "Not log-concave", 
                             "Not log-concave"))) %>% 
  ggplot(aes(x = mu_norm, y = reject_prop, col = log_concave)) +
  facet_wrap(. ~ d) +
  geom_line(color = "darkgrey") +
  geom_point() +
  geom_hline(yintercept = 0.10, lty = "dashed") +
  labs(x = expression("||"*mu*"|| in normal location family"), 
       y = "Rejection proportion",
       col = "",
       title = expression("Permutation test for H"[0]*
                            ": Log-concave vs H"[1]*": Not log-concave"),
       subtitle = expression("Normal location family f(x) = 0.5"*phi[d]*"(x)"~
                               "+ 0.5"*phi[d]*"(x -"~mu*"). B = 99, n = 250.")) +
  scale_color_manual(values = c("blue", "red")) +
  paper_theme

#####################
##### Save plot #####
#####################

ggsave(plot = shift_legend(perm_test_reject_n250),
       filename = "sim_plots/figure_11.pdf",
       width = 8, height = 5)
