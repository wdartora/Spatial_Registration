# 01_preprocessing_visium.R

rm(list = ls())

####################################################################################################
# STEP 1: Load and preprocess Visium spatial transcriptomics data
####################################################################################################

# Required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(magrittr)
library(EBImage)
library(gridExtra)
library(dplyr)

# Define project root directory and set working directory
project_root <- "C:/Users/darto/Documents/William/Cornell/Project/Spatial_Registration/Proj02/MPI_Registration/"
setwd(project_root)

# Load custom function to handle multi-slice Visium objects
source(file.path("scripts", "utils", "backup_functions.R"))

# User-defined input: slide with two samples separated by underscore
slide_name <- "ABCA7-3-a1-8_BIN1-6-a3-8"
slice_names <- strsplit(slide_name, "_")[[1]]
sliceA <- unname(slice_names[1])
sliceB <- unname(slice_names[2])

# Output directories
output_fig_dir <- file.path(project_root, "results", "plots")
output_r_dir <- file.path(project_root, "results", "r_output")

# Increase memory for Seurat
options(future.globals.maxSize = 891289600)

# Define path to Space Ranger 'outs' folder
visium_base_dir <- file.path(project_root, "data")
data_dir <- file.path(visium_base_dir, slide_name, "outs")

# Create output directories if they do not exist
if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE)
if (!dir.exists(output_r_dir)) dir.create(output_r_dir, recursive = TRUE)

####################################################################################################
# Helper functions
####################################################################################################

# Fix missing spatial orientation metadata required by SpatialDimPlot in some Seurat versions
fix_spatial_orientation <- function(obj) {
  img_name <- names(obj@images)[1]
  img <- obj@images[[img_name]]
  
  if (length(img@coords_x_orientation) == 0) {
    img@coords_x_orientation <- "left-to-right"
  }
  
  obj@images[[img_name]] <- img
  return(obj)
}

# Save one Seurat object per slice
save_slice_object <- function(se_obj, sample_name, output_dir) {
  se_obj$slice <- sample_name
  saveRDS(se_obj, file = file.path(output_dir, paste0(sample_name, ".rds")))
}

####################################################################################################
# Load full Visium object
####################################################################################################

se <- Load10X_Spatial_2(data.dir = data_dir, assay = "Spatial")
se <- fix_spatial_orientation(se)

# Add spatial coordinates (x and y) to metadata
se$x <- se@images$slice1@boundaries$centroids@coords[, 1]
se$y <- se@images$slice1@boundaries$centroids@coords[, 2]

####################################################################################################
# Visualize full coordinate space and define cutoff for manual splitting
####################################################################################################

cutoff_x <- 8400

p_coords <- ggplot(se@meta.data, aes(x = x, y = y)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_y_reverse() +
  scale_x_continuous(breaks = seq(0, max(se$x), by = 1000)) +
  theme_minimal() +
  labs(
    title = "Spot coordinates (x vs y)",
    x = "x",
    y = "y"
  ) +
  geom_vline(xintercept = cutoff_x, linetype = "dashed", color = "red")

p_coords

ggsave(
  filename = file.path(output_fig_dir, paste0(slide_name, "_spot_coordinates_split_check.png")),
  plot = p_coords,
  height = 6,
  width = 10
)

####################################################################################################
# Split full object into two slices using x cutoff
####################################################################################################

se1 <- subset(se, subset = x < cutoff_x)
se2 <- subset(se, subset = x >= cutoff_x)

# Fix orientation again after subset
se1 <- fix_spatial_orientation(se1)
se2 <- fix_spatial_orientation(se2)

# Update Seurat objects to synchronize spatial image and cell data
se1 <- UpdateSeuratObject(se1)
se2 <- UpdateSeuratObject(se2)

# Rename image slots to match slice names
names(se1@images) <- sliceA
names(se2@images) <- sliceB

# Add slice label to metadata
se1$slice <- sliceA
se2$slice <- sliceB

####################################################################################################
# Quick visualization to verify split
####################################################################################################

p_se1 <- SpatialDimPlot(se1, images = sliceA) + NoLegend() + ggtitle(sliceA)
p_se2 <- SpatialDimPlot(se2, images = sliceB) + NoLegend() + ggtitle(sliceB)


p_se1
ggsave(
  filename = file.path(output_fig_dir, paste0(sliceA, "_SpatialDimPlot.png")),
  plot = p_se1,
  height = 6,
  width = 6
)

p_se2
ggsave(
  filename = file.path(output_fig_dir, paste0(sliceB, "_SpatialDimPlot.png")),
  plot = p_se2,
  height = 6,
  width = 6
)

####################################################################################################
# Save split Seurat objects
####################################################################################################

save_slice_object(se1, sliceA, output_r_dir)
save_slice_object(se2, sliceB, output_r_dir)

####################################################################################################
# Reload and merge individual slices
####################################################################################################

paths <- unname(slice_names)

se_merged <- readRDS(file.path(output_r_dir, paste0(paths[1], ".rds")))

if (length(paths) > 1) {
  for (path in paths[-1]) {
    print(path)
    se_tmp <- readRDS(file.path(output_r_dir, paste0(path, ".rds")))
    se_merged <- merge(se_merged, se_tmp)
  }
}

####################################################################################################
# Standard Seurat workflow
####################################################################################################

se_merged <- SCTransform(
  se_merged,
  assay = "Spatial",
  vars.to.regress = "nCount_Spatial",
  verbose = FALSE
)

DefaultAssay(se_merged) <- "SCT"

se_merged <- RunPCA(se_merged, assay = "SCT", verbose = FALSE)
se_merged <- FindNeighbors(se_merged, reduction = "pca", dims = 1:30)
se_merged <- FindClusters(se_merged, verbose = FALSE)
se_merged <- RunUMAP(se_merged, reduction = "pca", dims = 1:30)

# Make sure image names remain linked to the correct slice labels
names(se_merged@images) <- paths

####################################################################################################
# Plot spatial clustering per slice
####################################################################################################

spatial_plot_list <- lapply(paths, function(img_name) {
  SpatialDimPlot(se_merged, images = img_name) +
    ggtitle(img_name) +
    NoLegend() +
    theme(aspect.ratio = 3000 / 5000)
})

combined_spatial_plot <- gridExtra::grid.arrange(grobs = spatial_plot_list, nrow = 2)

ggsave(
  filename = file.path(output_fig_dir, "All_Spatial_DimPlot.png"),
  plot = combined_spatial_plot,
  height = 8,
  width = 6
)

####################################################################################################
# Save merged Seurat object
####################################################################################################

saveRDS(se_merged, file = file.path(output_r_dir, "Visium_Seurat_object.rds"))



####################################################################################################
# STEP 2: Extract hi-res image and scaled coordinates for use in registration
####################################################################################################

library(purrr)
library(dplyr)
library(EBImage)
library(ggplot2)

# This step assumes the following objects already exist from STEP 1:
# se_merged, paths, sliceA, sliceB, slide_name, output_fig_dir, output_r_dir, visium_base_dir

####################################################################################################
# Extract spatial coordinates and scale factor for each slice
####################################################################################################

coords <- map_dfr(paths, function(img_name) {
  current_coords <- data.frame(se_merged@images[[img_name]]@boundaries$centroids@coords)
  rownames(current_coords) <- se_merged@images[[img_name]]@boundaries$centroids@cells
  current_coords$slice <- img_name
  current_coords$hires_sf <- se_merged@images[[img_name]]@scale.factors$hires
  current_coords %>% relocate(slice)
})

write.csv(
  coords,
  file = file.path(output_fig_dir, "Coords_slices.csv"),
  row.names = TRUE
)

####################################################################################################
# Read the high-resolution image used for spot mapping
####################################################################################################

hires_img_path <- file.path(
  visium_base_dir,
  slide_name,
  "outs",
  "spatial",
  "tissue_hires_image.png"
)

img <- readImage(hires_img_path)

####################################################################################################
# Convert spot coordinates to hires image coordinate system
####################################################################################################

coords_sub <- coords

# Keep this convention consistent with your original workflow
coords_sub$imagerow <- coords_sub$x * coords_sub$hires_sf
coords_sub$imagecol <- coords_sub$y * coords_sub$hires_sf

sf <- unique(coords_sub$hires_sf)
print(sf)

# Quick checks by slice
print(max(coords_sub$imagerow[coords_sub$slice == sliceA]))
print(min(coords_sub$imagerow[coords_sub$slice == sliceB]))

####################################################################################################
# Visualize spots scaled to hires image
####################################################################################################

scaled_plot <- ggplot(coords_sub, aes(x = imagecol, y = imagerow, color = slice)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_y_reverse() +
  coord_fixed() +
  theme_light() +
  labs(
    title = "Spots scaled to hires image",
    x = "imagecol (X)",
    y = "imagerow (Y)"
  ) +
  theme(legend.position = "bottom")

scaled_plot

ggsave(
  filename = file.path(output_fig_dir, paste0("Spots_scaled_to_hires_image_", slide_name, ".png")),
  plot = scaled_plot,
  width = 7,
  height = 6,
  dpi = 300
)

####################################################################################################
# Define manual cutoffs for splitting the hires image
####################################################################################################

# Check these values visually for each slide
bottom_cutoff <- 1800
top_cutoff <- 1850

scaled_plot +
  geom_vline(xintercept = bottom_cutoff, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = top_cutoff, color = "red", linetype = "dashed")

# Split hires image into left and right portions
img1 <- img[, 1:bottom_cutoff, ]
img2 <- img[, top_cutoff:ncol(img), ]

####################################################################################################
# Save split hires images
####################################################################################################

writeImage(img1, file.path(output_fig_dir, paste0(sliceA, "_hires_image.png")), quality = 100)
writeImage(img2, file.path(output_fig_dir, paste0(sliceB, "_hires_image.png")), quality = 100)
####################################################################################################
# Adjust spot coordinates to match the split hires images
####################################################################################################

# Because the second image starts at 'top_cutoff' in the original image columns,
# subtract this offset from the horizontal coordinate of sliceB
coords_spots <- coords_sub %>%
  mutate(
    bottom_cutoff = bottom_cutoff,
    top_cutoff = top_cutoff,
    imagecol = ifelse(slice == sliceB, imagecol - top_cutoff, imagecol)
  )

write.csv(
  coords_spots,
  file = file.path(output_r_dir, "spots_coords_hi_res.csv"),
  quote = FALSE,
  row.names = FALSE
)

####################################################################################################
# Final message
####################################################################################################

message("Preprocessing completed for slide: ", slide_name)
message("Now open: 02_generate_landmarks.R")
