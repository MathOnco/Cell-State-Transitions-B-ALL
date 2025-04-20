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
  pca_result <- prcomp(data, scale. = TRUE)  # Removed [,2:17] since we'll pass the correct columns
  
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
setwd("/Users/4477547/Documents/GitHub/Cell-State-Transitions-B-ALL/Data")
Data <- read.csv("markov.csv")
# Transpose the data and set the first row as column names
transposed_data <- as.data.frame(t(Data))
colnames(transposed_data) <- transposed_data[1,]
transposed_data <- transposed_data[-1,]

# For MRD(flow)
transposed_data$MRD_new <- sapply(1:nrow(transposed_data), function(i) {
  if (transposed_data$"DiseaseProgression"[i] == "Remission") {
    return(2)
  } else if (grepl("Positive \\(> 10\\^-3\\)", transposed_data$"MRD(flow)"[i])) {
    return(1)
  } else {
    return(0)  # For Negative and Positive (< 10^-3)
  }
})

# For BCR::ABL1
transposed_data$BCR_new <- sapply(1:nrow(transposed_data), function(i) {
  if (transposed_data$"DiseaseProgression"[i] == "Remission") {
    return(2)
  } else if (grepl("BCR::ABL1|BCR::ABL1-like", transposed_data$"BCR::ABL1/like"[i])) {
    return(1)
  } else {
    return(0)  # For Negative
  }
})

# For DiseaseProgression
transposed_data$Disease_new <- sapply(transposed_data$"DiseaseProgression", function(x) {
  if (x == "Relapse") {
    return(1)
  } else if (x == "Remission") {
    return(2)
  } else {
    return(0)
  }
})



m_columns <- grep("^M[1-4][1-4]$", colnames(transposed_data), value = TRUE)
bcr_data <- transposed_data[, c(m_columns,"BCR_new")]
bcr_data[] <- lapply(bcr_data, function(x) as.numeric(as.character(x)))

mrd_data <- transposed_data[, c(m_columns, "MRD_new")]
mrd_data[] <- lapply(mrd_data, function(x) as.numeric(as.character(x)))

relapse_data <- transposed_data[, c(m_columns, "Disease_new")]
relapse_data[] <- lapply(relapse_data, function(x) as.numeric(as.character(x)))

# Create plots with new labels
bcr_plot <- create_pca_plot(bcr_data[,1:16], transposed_data$BCR_new, "BCRABL",
                            c("BCRABLNegative" = "red",
                              "BCRABLPositive" = "blue",
                              "Remission" = "gray"))
print(bcr_plot)

mrd_plot <- create_pca_plot(mrd_data[,1:16], transposed_data$MRD_new, "MRD",
                            c("MRDNegative" = "#A0C1C2",
                              "MRDPositive" = "#CC9933",
                              "Remission" = "gray"))
print(mrd_plot)

relapse_plot <- create_pca_plot(relapse_data[,1:16], transposed_data$Disease_new, "Relapse",
                                c("RelapseNegative" = "#bc4749",
                                  "RelapsePositive" = "#386641",
                                  "Remission" = "gray"))
print(relapse_plot)

pdf("BCR_ABL1_PCA.pdf", width = 6, height = 6)
print(bcr_plot)
dev.off()

pdf("MRD_PCA.pdf", width = 6, height = 6)
print(mrd_plot)
dev.off()

pdf("Relapse_PCA.pdf", width = 6, height = 6)
print(relapse_plot)
dev.off()
