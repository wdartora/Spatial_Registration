# 04_annotate_seurat.R

# This script imports the Python output and integrates it into the Seurat objects
library(Seurat)
library(tidyverse)

# Import Visium data object and region file from Python script
output_r_dir <- "results/r_output/"

se <- readRDS( paste0(output_r_dir, "Visium_Seurat_object.rds"))
se$cellid <- paste0(se$slice,"#",colnames(se))
py_output_dir <- "results/py_output/"
region_spots <- read.csv(paste0(py_output_dir, "spots_coords_regions.csv")) %>%
  select(slice, x,y, region) %>%
  mutate(region = ifelse(grepl("Serie",region),"unidentified",region))

# Add region information
se$region <- region_spots$region

# Visualize to ensure done properly
sdp <- lapply(paths, function(x)
  SpatialDimPlot(se, images = x, group.by="region",label=T,repel=T, label.size=2)+
    theme(aspect.ratio=3000/5000)+ggtitle(x)+NoLegend()
  
)
sdp <- gridExtra::grid.arrange(grobs = sdp,nrow=2)
ggsave(paste0(output_fig_dir,"All_Spatial_DimPlot_Regions.png"),sdp, height=6,width=8)

saveRDS(se, paste0(output_r_dir, "Visium_Seurat_object.rds"))

message("Todos os objetos Seurat foram anotados com as regiões espaciais!")
