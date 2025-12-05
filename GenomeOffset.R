setwd('C:/Users/59609/Desktop/plotting/Final/gf/regional_analysis/') # Update this path
require(raster)
require(geosphere)
require(gdm)
require(foreach)
require(parallel)
require(doParallel)
require(gradientForest)
require(fields)
library(sf)
library(raster)
library(corrplot)
library(ggrepel)
library(ggplot2)

# 1. Data Preparation
# Note: Update the path below to your actual data file
gf.data <- read.csv("D:/xwechat_files/wxid_7276162773812_3671/msg/file/2025-07/output.csv", header = T, row.names = 1)

# Filter SNPs
SNP_vector <- apply(gf.data, 2, function(x) length(unique(x)))
df <- data.frame(ColumnIndex = seq_along(SNP_vector), UniqueCount = SNP_vector)
df <- df[df$UniqueCount > 5, ]
gf.alle <- gf.data[, df$ColumnIndex]

# Load environmental layers
BIO_files <- list.files('D:/study/environmental layers', pattern = 'tif', full.names = T)
BIO <- stack(BIO_files)

# Load location data
loc <- read.csv("D:/xwechat_files/wxid_7276162773812_3671/msg/file/2025-07/loc.csv", header = T)
loc.env <- data.frame(extract(BIO, loc[,1:2]))
colnames(loc.env) <- c("BIO1", paste0("BIO", 10:19), paste0("BIO", 2:9))
loc.env <- loc.env[, paste0("BIO", 1:19)]

# Gradient Forest Setup
maxLevel <- log2(0.368 * nrow(gf.alle) / 2) # Corrected '\' to '/'
gf.env <- loc.env[,c("BIO2","BIO3","BIO5","BIO13","BIO19")]
gf1.env <- loc.env
set.seed(123)

gfMod <- gradientForest(cbind(gf1.env, gf.alle),
                        predictor.vars=colnames(gf.env),
                        response.vars=colnames(gf.alle), ntree=1000,
                        maxLevel=maxLevel, trace=T, corr.threshold=0.50)

# Plot Split Density
pdf(file="Splitsdensityplots.pdf",width = 6,height = 4)
plot(gfMod, plot.type="S", imp.vars= colnames(gf.env), leg.posn="topright",cex.legend=0.4, cex.axis=0.6, cex.lab=0.7, line.ylab=0.9,par.args=list(mgp=c(1.5, 0.5, 0), mar=c(3.1,1.5,0.1,1)))
dev.off()

# Plot Response Curves
pdf(file="Response_curve.pdf",width = 4,height = 4)
plot(gfMod, plot.type="C", imp.vars= colnames(gf.env), show.species=F,common.scale=T, cex.axis=0.6, cex.lab=0.7, line.ylab=0.9,par.args=list(mgp=c(1.5, 0.5, 0), mar=c(2.5,1.0,0.1,0.5), omi=c(0,0.3,0,0)))
dev.off()

# Quantify Environmental Gradients
greengrid <- gf.env # Load environmental data for all points in the study area
SDM_point <- read.csv("D:/xwechat_files/wxid_7276162773812_3671/msg/file/2025-07/loc.csv", header = TRUE) # Load lat/long for all points in study area

# Quantify environmental gradients based on Gradient Forest
all_tgrid = cbind(SDM_point[,c("Long","Lat")], predict(gfMod, greengrid)) 

n <- sum(is.na(all_tgrid)) # Check for missing values in dataframe
Trns_grid <- na.omit(all_tgrid) # Remove missing values
n <- sum(is.na(Trns_grid)) # Re-check for missing values

# PCA Analysis based on environmental data in the study area
all_PCs <- prcomp(Trns_grid[,c('BIO2','BIO3','BIO5','BIO13','BIO19')],center=TRUE, scale.=FALSE)
summary(all_PCs)

a1 <- all_PCs$x[,1] # PC1
a2 <- all_PCs$x[,2] # PC2
a3 <- all_PCs$x[,3] # PC3

# RGB Color Calculation
r <- a1+a2
g <- -a2
b <- a3+a2-a1

# Normalize RGB values (Corrected '\' to '/')
r <- (r-min(r)) / (max(r)-min(r)) * 255
g <- (g-min(g)) / (max(g)-min(g)) * 255
b <- (b-min(b)) / (max(b)-min(b)) * 255

grid <- SDM_point[,c("Long","Lat")] # Lat/Long of study area
grid$R=r
grid$G=g
grid$B=b

nvs <- dim(all_PCs$rotation)[1] # Check data dimensions
vec <- c('BIO2','BIO3','BIO5','BIO13','BIO19')
lv <- length(vec)
vind <- rownames(all_PCs$rotation) %in% vec
scal <- 60

# Corrected '\' to '/' for ranges
xrng <- range(all_PCs$x[,1], all_PCs$rotation[,1]/scal)*1.1
yrng <- range(all_PCs$x[,2], all_PCs$rotation[,2]/scal)*1.1

# PCA Plot of Gradient Forest Transformed Climate Variables
pdf(file="all_PCplot.pdf") 
plot((all_PCs$x[,1:2]), xlim=xrng, ylim=yrng, pch=".", cex=7, col=rgb(r,g,b, max = 255), asp=1)
# Corrected '\' to '/' for arrows and text
arrows(rep(0,lv), rep(0,lv), all_PCs$rotation[,1]/scal*15, all_PCs$rotation[,2]/scal*15, length = 0.1)
jit <- 0.0015
text(all_PCs$rotation[,1]/scal*5+jit*sign(all_PCs$rotation[,1]), all_PCs$rotation[,2]/scal*5+jit*sign(all_PCs$rotation[,2]), labels = vec)
dev.off()

# Future Climate Data Loading
ssp126_2050 <- stack("D:/study/CMIP6(BCC-CSM2-MR)_2.5/ssp126_2.5/wc2.1_2.5m_bioc_BCC-CSM2-MR_ssp126_2041-2060.tif")
ssp585_2050 <- stack("D:/study/CMIP6(BCC-CSM2-MR)_2.5/ssp585_2.5/wc2.1_2.5m_bioc_BCC-CSM2-MR_ssp585_2041-2060.tif")

SDM <- raster("Sampling_points_current_fitness.asc") # Renamed from Chinese
plot(SDM)
SDM.points <- rasterToPoints(SDM)[rasterToPoints(SDM)[,3] > 1, ]
#SDM.points <- rasterToPoints(SDM)

SDM.ENV <- data.frame(extract(BIO,SDM.points[,c("x","y")]))
SDM.ENV_clean <- na.omit(SDM.ENV)
SDM.points_clean <- SDM.points[!rowSums(is.na(SDM.ENV)), ]

SDM.ENV <- data.frame(extract(BIO,SDM.points_clean[,c("x","y")]))
colnames(SDM.ENV) <- c("BIO1", paste0("BIO", 10:19), paste0("BIO", 2:9))
SDM.ENV <- SDM.ENV[,c("BIO2","BIO3","BIO5","BIO13","BIO19")]

ssp585_2050.SDM <- data.frame(extract(ssp585_2050,SDM.points_clean[,c("x","y")]))
ssp126_2050.SDM <- data.frame(extract(ssp126_2050,SDM.points_clean[,c("x","y")]))

# Set column names
colnames(ssp585_2050.SDM) <- paste0("BIO", 1:19)
colnames(ssp126_2050.SDM) <- paste0("BIO", 1:19)

# Select relevant bioclim variables
ssp585_2050.SDM <- ssp585_2050.SDM[, c("BIO2","BIO3","BIO5","BIO13","BIO19")]
ssp126_2050.SDM <- ssp126_2050.SDM[, c("BIO2","BIO3","BIO5","BIO13","BIO19")]

SDM.ENV
SDM.ENV.2 <- SDM.ENV[complete.cases(ssp585_2050.SDM), ]
ssp585_2050.SDM.2 <- ssp585_2050.SDM[complete.cases(ssp585_2050.SDM), ]
ssp126_2050.SDM.2 <- ssp126_2050.SDM[complete.cases(ssp585_2050.SDM), ]
SDM.points_clean.2 <- SDM.points_clean[complete.cases(ssp585_2050.SDM), ]

all_tgrid <- cbind(SDM.points_clean.2[,c("x","y")], predict(gfMod, SDM.ENV.2)) 

# Predict Future Genomic Vulnerability
future_585 <- cbind(SDM.points_clean.2[,c("x","y")],
                    predict(gfMod,ssp585_2050.SDM.2))
colSums(is.na(SDM.ENV)) 

future_126 <- cbind(SDM.points_clean.2[,c("x","y")],
                    predict(gfMod,ssp126_2050.SDM.2))

# Calculate Euclidean distance for columns 3 to 7 (Genetic Offset)
genOffsetAll <- sqrt(rowSums((future_585[, 3:7] - all_tgrid[, 3:7])^2))
Offset_ssp585 <- cbind(future_585[,c("x","y")],genOffsetAll)
Offset_ssp126 <- cbind(future_126[,c("x","y")],genOffsetAll)
colnames(Offset_ssp585)[3]<-"offset"
colnames(Offset_ssp126)[3]<-"offset"
max(Offset_ssp126[,3])
max(Offset_ssp585[,3])

# Extract offset data
offset_values <- Offset_ssp585[, 3]

# Generate 14 quantile breaks (corresponding to 13 groups, probabilities from 0 to 1)
# probs = seq(0, 1, length.out = 14) means 13 equal probability quantiles
breaks <- quantile(offset_values, probs = seq(0, 1, length.out = 14), na.rm = TRUE, names = FALSE)

# Handle potential duplicate breaks (if there are many duplicate values, quantiles might be identical)
breaks <- unique(breaks)
# If breaks are fewer than 14 after deduplication, automatically fill to 14 (avoid insufficient groups)
if (length(breaks) < 14) {
  breaks <- seq(min(breaks), max(breaks), length.out = 14)
}

# Generate group names (corresponding to quantile breaks)
groups <- paste0(
  round(breaks[-length(breaks)], 4),  # Start point (4 decimal places)
  "-", 
  round(breaks[-1], 4)                # End point
)
groups[length(groups)] <- paste0(round(breaks[length(groups)-1], 4), "+")  # Last group uses "+"

# Color vector (Keep 13 colors, order consistent with groups)
offset_color <- viridis::viridis(13)
names(offset_color) <- groups  # Strictly match color names with group names

# Grouping Function (Force level to be an ordered factor, specifying order)
creategroup <- function(gf) {
  colnames(gf) <- c("x", "y", "offset")
  # Use cut to group as character, then convert to factor (specifying levels order)
  level_char <- cut(
    gf$offset,
    breaks = breaks,
    labels = groups,
    right = FALSE,
    include.lowest = TRUE
  )
  # Key step: Convert level to factor, levels order consistent with groups (ensuring colors sort correctly)
  gf$level <- factor(level_char, levels = groups)
  return(gf)
}

# Generate grouped data (level is an ordered factor)
resdf_585 <- creategroup(Offset_ssp585)
resdf_126 <- creategroup(Offset_ssp126)

# Check factor level order (should match groups exactly)
levels(resdf_585$level)

# Check color vector order (should match groups)
names(offset_color)

# =======================
# Load Packages
# =======================
library(ggplot2)
library(sf)
library(grid)
library(showtext)

# =======================
# Enable showtext and add Arial font
# =======================
showtext_auto()
font_add("Arial", regular = "C:/Windows/Fonts/arial.ttf")  # Windows system path

# =======================
# Read Country Boundary Shapefile
# =======================
# Update path: Sampling_Points_Country/Sampling_Points_Country_Complete.shp
country_shp <- sf::st_read("Sampling_Points_Country/Sampling_Points_Country_Complete.shp")

# =======================
# Plot SSP585 Map
# =======================
pdf(file = "SSP585.pdf", height = 5, width = 8)
ggplot() +
  geom_tile(data = resdf_585, aes(x = x, y = y, fill = level), alpha = 1) +
  geom_sf(data = country_shp, fill = "transparent", color = "black", size = 0.5) +
  scale_fill_manual(values = offset_color) +
  labs(x = "Longitude", y = "Latitude", title = "SSP585 Genomic Vulnerability") +
  theme_bw() +
  theme(
    text = element_text(family = "Arial"),         # Use Arial font
    axis.ticks.length = unit(0.25, "lines"),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    legend.key.size = unit(0.07, 'inch'),         # Shrink legend
    legend.title = element_blank(),
    legend.position = c(0.95, 0.85),              # Legend at top right
    legend.background = element_rect(colour = NA, fill = NA),  # Remove legend border
    legend.text = element_text(size = 6, color = "black"),
    panel.grid.major = element_line(colour = NA),
    plot.background = element_rect(fill = "white"),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "mm")
  ) +
  # Add grey background for title
  annotation_custom(grob = rectGrob(gp = gpar(fill = "grey", col = NA)),
                    xmin = -Inf, xmax = Inf, ymin = Inf, ymax = Inf) + 
  theme(plot.title = element_text(size = 18, hjust = 0.5)) 
dev.off()

# =======================
# Plot SSP126 Map
# =======================
pdf(file = "SSP126.pdf", height = 5, width = 8)
ggplot() +
  geom_tile(data = resdf_126, aes(x = x, y = y, fill = level), alpha = 1) +
  geom_sf(data = country_shp, fill = "transparent", color = "black", size = 0.5) +
  scale_fill_manual(values = offset_color) +
  labs(x = "Longitude", y = "Latitude", title = "SSP126 Genomic Vulnerability") +
  theme_bw() +
  theme(
    text = element_text(family = "Arial"),         # Use Arial font
    axis.ticks.length = unit(0.25, "lines"),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    legend.key.size = unit(0.07, 'inch'),         # Shrink legend
    legend.title = element_blank(),
    legend.position = c(0.95, 0.85),              # Legend at top right
    legend.background = element_rect(colour = NA, fill = NA),  # Remove legend border
    legend.text = element_text(size = 6, color = "black"),
    panel.grid.major = element_line(colour = NA),
    plot.background = element_rect(fill = "white"),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "mm")
  ) +
  # Add grey background for title
  annotation_custom(grob = rectGrob(gp = gpar(fill = "grey", col = NA)),
                    xmin = -Inf, xmax = Inf, ymin = Inf, ymax = Inf) + 
  theme(plot.title = element_text(size = 18, hjust = 0.5)) 
dev.off()