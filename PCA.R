library(ggplot2)

# Set input file paths
input <- "C:/Users/59609/Desktop/huatu/最终/pca/PCA_out.eigenvec"
x <- 1  # PC1
y <- 2  # PC2
pop_file <- "C:/Users/59609/Desktop/huatu/最终/pca/sample.pop"
outpre <- "output_prefix1"

# Read data
vec <- read.table(input, header = FALSE, row.names = 1, sep = " ")
pop <- read.table(pop_file, header = FALSE, row.names = 1, sep = "\t")
population <- factor(pop[rownames(vec), 1])  # Convert to factor

# Custom colors
custom_colors <- c(
  "SA" = "#56B4E9",
  "MM+WYN" = "#8cd0c3",
  "SEA+EYN" = "#CC79A7",
  "HN" = "#E69F00"
)

# Calculate PCA variance explained (assuming .eigenval file exists)
eigenval_file <- sub(".eigenvec", ".eigenval", input)
if (file.exists(eigenval_file)) {
  eigenval <- scan(eigenval_file)
  total_var <- sum(eigenval)
  pc1_var <- round(eigenval[x] / total_var * 100, 1)
  pc2_var <- round(eigenval[y] / total_var * 100, 1)
} else {
  warning(".eigenval file not found. Variance explained will not be displayed.")
  pc1_var <- pc2_var <- NA
}


p1 <- ggplot(mapping = aes(
  x = vec[, x + 1],  # Assuming column 1 is sample ID, PC1 starts at column 2
  y = vec[, y + 1],
  colour = population
)) +
  geom_point(size = 3, shape = 16) +  # shape = 16 for solid circles
  stat_ellipse(linewidth = 0.8, linetype = "solid") +
  theme_bw() +
  xlab(ifelse(is.na(pc1_var), 
              paste0("PC", x), 
              paste0("PC", x, " (", pc1_var, "%)"))) +
  ylab(ifelse(is.na(pc2_var), 
              paste0("PC", y), 
              paste0("PC", y, " (", pc2_var, "%)"))) +
  scale_color_manual(values = custom_colors, name = "Population") +
  labs(title = "PCA Plot") +
  theme(
    text = element_text(family = "Arial"),  # Global font
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.border = element_rect(linewidth = 1.3, color = "black"),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    # Optional: Adjust legend position
    legend.position = "right",
    # Optional: Increase legend key size
    legend.key.size = unit(1.2, "lines")
  )


# Save plots
ggsave(
  paste0(outpre, "_PCA.png"), 
  p1, 
  width = 10, 
  height = 8, 
  dpi = 300
)

ggsave("ggplot.pdf", plot = p1, width = 6, height = 4, device = cairo_pdf)