# 01_preprocessing_visium.R

rm(list=ls())

####################################################################################################
# STEP 1: Load and preprocess Visium spatial transcriptomics data
####################################################################################################

# --- Required packages ---
library(Seurat)
library(ggplot2)
library(patchwork)
library(magrittr)
library(EBImage)
library(gridExtra)
library(dplyr)

# --- Define project root directory and set working directory ---
project_root <- "W:/Project/MPI_Registration/Proj_Spatial_Registration"
setwd(project_root)

# --- Load custom function to handle multi-slice Visium objects ---
source(file.path("scripts", "utils", "backup_functions.R"))

# --- User-defined input: slide with two samples separated by underscore ---
# Set the folder's name from VISIUM data (slices)
slide_name <- "ABCA7-5-a1-7_BIN1-5-a3-8"
slice_names <- strsplit(slide_name, "_")[[1]]
sliceA <- unname(slice_names[1])
sliceB <- unname(slice_names[2])



# --- Output directories ---
output_fig_dir <- "results/plots/"
output_r_dir <- "results/r_output/"

# --- Increase memory for Seurat ---
options(future.globals.maxSize = 891289600)

# --- Define path to Space Ranger 'outs' folder (short, local path recommended) ---
visium_base_dir <- "W:/Data_base_VISIUM/MPI/data/"
data_dir <- file.path(visium_base_dir, slide_name, "outs")

# --- Create output directories if they don't exist ---
if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE)
if (!dir.exists(output_r_dir)) dir.create(output_r_dir, recursive = TRUE)

# --- Load Visium object and plot to visualize boundary between slices ---
se <- Load10X_Spatial_2(data.dir = data_dir, assay = "Spatial")
SpatialDimPlot(se) + NoLegend()

# --- Add spatial coordinates (X and Y) to metadata ---
se$x <- se@images$slice1@boundaries$centroids@coords[,1]
se$y <- se@images$slice1@boundaries$centroids@coords[,2]

# --- Visualize full tissue image with spot coordinates ---
ggplot(se@meta.data, aes(x = x, y = y)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_y_reverse() +
  scale_x_continuous(breaks = seq(0, max(se$x), by = 1000)) +
  theme_minimal() +
  labs(title = "Coordenadas dos spots (x vs y)", x = "x", y = "y") +
  geom_vline(xintercept = 9800, linetype = "dashed", color = "red")

# --- Split full object into two slices manually using x cutoff ---
se1 <- subset(se,  x < 9400, invert =F)
se2 <- subset(se,  x < 10000, invert =T)

# --- Quick visualization to verify split ---
SpatialDimPlot(se1) + NoLegend()
SpatialDimPlot(se2) + NoLegend()



# Standard Seurat workflow to process data, normalize, dimension reductionality, cluster,
# not strictly needed this early
sct_transform <- function(se, sample_name){
  
  # se <- SCTransform(se, assay = "Spatial", verbose = FALSE)
  # se <- RunPCA(se, assay = "SCT", verbose = FALSE)
  # se <- FindNeighbors(se, reduction = "pca", dims = 1:30)
  # se <- FindClusters(se, verbose = FALSE)
  # se <- RunUMAP(se, reduction = "pca", dims = 1:30)
  # plt <- SpatialDimPlot(se,label=T)+ theme(aspect.ratio = 3000/5000)
  # ggsave(paste0(output_fig_dir,sample_name,"_Spatial_DimPlot.png"),plt, height=10,width=11)
  se$slice <- sample_name
  # --- Save individual slices as RDS files ---
  saveRDS(se, file = paste0(output_r_dir, sample_name,".rds"))
  return()
}


sct_transform(se1, sliceA)
sct_transform(se2, sliceB)

# --- Reload and merge individual slices ---
se <- readRDS(paste0(output_r_dir, sliceA, ".rds"))
paths <- unname(slice_names)

read_file <- function(path) {
  print(path)
  se1 <- readRDS(paste0(output_r_dir, path, ".rds"))
  se <<- merge(se, se1)
}
sapply(paths[-1], read_file)

# --- Apply standard Seurat workflow ---
se <- SCTransform(se, assay = "Spatial", vars.to.regress = "nCount_Spatial", verbose = FALSE)
DefaultAssay(se) <- "SCT"
se <- RunPCA(se, assay = "SCT", verbose = FALSE)
se <- FindNeighbors(se, reduction = "pca", dims = 1:30)
se <- FindClusters(se, verbose = FALSE)
se <- RunUMAP(se, reduction = "pca", dims = 1:30)

names(se@images) <- paths

# --- Plot spatial clustering per slice ---
sdp <- lapply(paths, function(x)
  SpatialDimPlot(se, images = x) +
    theme(aspect.ratio = 3000/5000) +
    ggtitle(x) + NoLegend()
)
sdp <- gridExtra::grid.arrange(grobs = sdp, nrow = 2)
ggsave(paste0(output_fig_dir, "All_Spatial_DimPlot.png"), sdp, height = 8, width = 6)

# --- Save merged Seurat object ---
saveRDS(se, paste0(output_r_dir, "Visium_Seurat_object.rds"))


####################################################################################################
# STEP 2: Extract hi-res image and scaled coordinates for use in registration
####################################################################################################

library(purrr)
library(dplyr)

# --- Extract spatial coordinates and scale factor for each slice ---
coords <- map_dfr(paths, function(x) {
  coords <- data.frame(se@images[[x]]@boundaries$centroids@coords)
  rownames(coords) <- se@images[[x]]@boundaries$centroids@cells
  coords$slice <- x
  coords$hires_sf <- se@images[[x]]@scale.factors$hires
  return(coords %>% relocate(slice))
})
coords %>% write.csv(paste0(output_fig_dir, "Coords_slices.csv"))

# --- Read the high-resolution image used for spot mapping ---
hires_img_path <- file.path(visium_base_dir, slide_name, "outs", "spatial", "tissue_hires_image.png")
img <- readImage(hires_img_path)

# --- Convert to image-based coordinates using scale factor ---
coords_sub <- coords
coords_sub$imagerow = coords_sub$x * coords_sub$hires_sf
coords_sub$imagecol = coords_sub$y * coords_sub$hires_sf
sf <- coords_sub %>% pull(hires_sf) %>% unique()

max(coords_sub$imagerow[grepl(sliceA, coords_sub$slice)])
min(coords_sub$imagerow[grepl(sliceB, coords_sub$slice)])

# --- Visualize spots scaled to image resolution ---
scaled_plot <- ggplot(coords_sub, aes(x = imagecol, y = imagerow, color = slice)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_y_reverse() +
  coord_fixed() +
  theme_light() +
  labs(title = "Spots scaled to hires image",
       x = "imagecol (X)", y = "imagerow (Y)") +
  theme(legend.position = "bottom")
ggsave(filename = file.path(output_fig_dir, paste0("Spots_scaled_to_hires_image_", slide_name, ".png")),
       plot = scaled_plot, width = 7, height = 6, dpi = 300)

# --- Define manual cutoffs for splitting the image ---
bottom_cutoff <- 2100
top_cutoff <- 2000

img1 <- img[, 1:bottom_cutoff, ]
img2 <- img[, top_cutoff:ncol(img), ]

# --- Save top half image ---
png(paste0(output_fig_dir, sliceA, "_hires_image.png"), height = dim(img1)[2], width = dim(img1)[1], res = 800)
display(img1, method = "raster")
dev.off()

# --- Save bottom half image ---
png(paste0(output_fig_dir, sliceB, "_hires_image.png"), height = dim(img2)[2], width = dim(img2)[1], res = 800)
display(img2, method = "raster")
dev.off()

# --- Adjust coordinates for lower slice to match split image ---
coords_spots <- coords_sub %>%
  mutate(bottom_cutoff = bottom_cutoff) %>%
  mutate(top_cutoff = top_cutoff) %>%
  mutate(imagerow = ifelse(slice == sliceB, imagerow - top_cutoff, imagerow))

# Save coordinates
coords_spots %>% write.csv(paste0(output_r_dir, "spots_coords_hi_res.csv"), quote = FALSE, row.names = FALSE)

# --- Launch landmark script ---
message("Preprocessing completed for slide: ", slide_name)
message("Now open: 02_generate_landmarks.R")
