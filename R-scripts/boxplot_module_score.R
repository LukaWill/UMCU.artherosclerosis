setwd("~/Graduation/Rstudio")

library(ggplot2)
module_score_cluster <- read.csv("module_score_networks.csv", sep = "\t", header = T)
module_score_cluster <- module_score_cluster[, -1]
module_score_cluster <- data.frame(module_score_cluster)

module_score_cluster$cluster <- factor(module_score_cluster$cluster, levels = c("0", "1", "2", "3", "4"))

box.plot = "~/Graduation/Rstudio/violin_plot_module_score.pdf"
pdf(file = box.plot)

generate_plot <- function(column_index, title) {
  ggplot(module_score_cluster, aes(x = cluster, y = module_score_cluster[, column_index], fill = cluster)) +
    geom_violin() +
    geom_jitter(color = "black", size = 0.1, alpha = 0.9) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(), 
      plot.title = element_text(size = 10),
      legend.text = element_text(size = 10), legend.title = element_text(size=10)
    ) +
    labs(title = title, y = 'modulescore', x = 'cluster')
}

plots <- list(
  list(column_index = 2, title = "N39"),
  list(column_index = 3, title = "N195"),
  list(column_index = 4, title = "N177"),
  list(column_index = 5, title = "N24"),
  list(column_index = 6, title = "N35"),
  list(column_index = 7, title = "N74"),
  list(column_index = 8, title = "N109"),
  list(column_index = 9, title = "N172"),
  list(column_index = 10, title = "N147")
)

plots <- lapply(plots, function(plot) {
  generate_plot(plot$column_index, plot$title)
})

# Plot the generated plots

gridExtra::grid.arrange(grobs = plots, ncol = 3)

dev.off()