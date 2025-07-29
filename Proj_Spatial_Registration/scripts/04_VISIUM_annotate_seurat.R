# 04_annotate_seurat.R

# This script imports the Python output and integrates it into the Seurat objects
library(Seurat)
library(tidyverse)

# Import Visium data object and region file from Python script
output_r_dir <- "results/r_output/"
py_output_dir <- "results/py_output/"

# --- Load Seurat object containing Visium slices ---
se <- readRDS( paste0(output_r_dir, "Visium_Seurat_object.rds"))

# --- Build a unique cell ID for each spot: <slice>#<barcode> ---
se$cellid <- paste0(se$slice,"#",colnames(se))

# --- Load Python output (spot regions) ---
region_spots <- read.csv(paste0(py_output_dir, "spots_coords_regions.csv")) %>%
  select(slice, x,y, region) %>%
  mutate(region = ifelse(grepl("Serie",region),"unidentified",region))

# --- Add region information to Seurat metadata ---
# This assumes same order of barcodes — safe if using same pipeline throughout
se$region <- region_spots$region


# --- Plot regions for each Visium slice ---
sdp <- lapply(paths, function(x)
  SpatialDimPlot(se, images = x, group.by="region",label=T,repel=T, label.size=2)+
    theme(aspect.ratio=3000/5000)+ggtitle(x)+NoLegend()
  
)

# --- Combine and save all plots ---
sdp <- gridExtra::grid.arrange(grobs = sdp,nrow=2)
ggsave(paste0(output_fig_dir,"All_Spatial_DimPlot_Regions.png"),sdp, height=6,width=8)


# --- Save the updated Seurat object ---
saveRDS(se, paste0(output_r_dir, "Visium_Seurat_object.rds"))


# --- Save IDs ---

# Define diretório de saída
project_root <- "W:/Project/MPI_Registration/Proj02/"
setwd(project_root)

py_output_dir <- "results/py_output"

# Cria dataframe com IDs
sample_IDs <- data.frame(
  sample_ID = c("PICALM-10-a1-7", "ABCA7-3-a1-8", "BIN1-6-a3-8", "ABCA7-8-a3",
                "BIN1-5-a1", "PICALM-5-a1-", "WTP-5-a7-8", "WTP-5-a2-8"),
  desi = c("10a1_PICALM", "3a1_ABCA7", "6a3_BIN1", "8a3_ABCA7",
           "5a1_BIN1", "5a1_PICALM", "5a7_WTP", "5a2_WTP")
)

# Salva arquivo CSV corretamente
write.csv(sample_IDs, file = file.path(py_output_dir, "convert_IDs.csv"), row.names = FALSE, quote = FALSE)



# --- Final message ---
message("Todos os objetos Seurat foram anotados com as regiões espaciais!")
