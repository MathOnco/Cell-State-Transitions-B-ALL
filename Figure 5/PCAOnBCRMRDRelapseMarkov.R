library(ggplot2)
library(GGally)
library(ggfortify) 
library(factoextra)
library("corrplot")
library(MASS)
library(caret)
library(ggrepel)


create_pca_plot <- function(data, labels, title_prefix, color_values) {
  # PCA
  pca_result <- prcomp(data[,2:17], scale. = TRUE)
  
  var_data <- as.data.frame(pca_result$rotation[, 1:2])
  var_data$Variable <- rownames(var_data)
  ind_data <- as.data.frame(pca_result$x[, 1:2])
  ind_data$Group <- factor(labels, levels = c(0, 1, 2), 
                          labels = c(paste0(title_prefix, "Negative"), 
                                   paste0(title_prefix, "Positive"), 
                                   "Remission"))
  

  scatter_range <- max(abs(ind_data$PC1), abs(ind_data$PC2)) 
  loading_range <- max(abs(var_data$PC1), abs(var_data$PC2)) 
  scaling_factor <- scatter_range / loading_range
  
 
  var_data$PC1 <- var_data$PC1 * 0.6 * scaling_factor
  var_data$PC2 <- var_data$PC2 * 0.6 * scaling_factor
  

  ggplot() +
    geom_segment(data = var_data, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.3, "cm")), 
                 color = "black", 
                 size = 0.7,
                 alpha = 0.8) +
    geom_text_repel(data = var_data, aes(x = PC1, y = PC2, label = Variable),
                    size = 3.5, color = "black", 
                    max.overlaps = Inf,
                    box.padding = 0.5,
                    point.padding = 0.3) +
    geom_point(data = ind_data, aes(x = PC1, y = PC2, color = Group), 
               size = 3, alpha = 0.8) +
    scale_color_manual(values = color_values) +
    labs(title = "",
         x = paste("PC1 (", round(summary(pca_result)$importance[2,1] * 100, 2), "%)", sep = ""),
         y = paste("PC2 (", round(summary(pca_result)$importance[2,2] * 100, 2), "%)", sep = "")) +
    scale_x_continuous(breaks = seq(-6, 6, by = 1), limits = c(-6, 6)) +
    scale_y_continuous(breaks = seq(-6, 6, by = 1), limits = c(-6, 6)) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_line(size = 0.5, color = "gray80"),
      panel.grid.minor = element_line(size = 0.25, color = "gray90"),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
}

# BCR-ABL1 data
bcr_data <- read.csv("BCRABL1Markov.csv")
bcr_plot <- create_pca_plot(bcr_data, bcr_data$BCR_ABL1, "BCRABL",
                           c("BCRABLNegative" = "red",
                             "BCRABLPositive" = "blue",
                             "Remission" = "gray"))
print(bcr_plot)

#  MRD data
mrd_data <- read.csv("MRDMarkov.csv")
mrd_plot <- create_pca_plot(mrd_data, mrd_data$MRD, "MRD",
                           c("MRDNegative" = "#A0C1C2",
                             "MRDPositive" = "#CC9933",
                             "Remission" = "gray"))
print(mrd_plot)

#  Relapse data
relapse_data <- read.csv("RelapseMarkov.csv")
relapse_plot <- create_pca_plot(relapse_data, relapse_data$Relapse, "Relapse",
                               c("RelapseNegative" = "#bc4749",
                                 "RelapsePositive" = "#386641",
                                 "Remission" = "gray"))
print(relapse_plot)

