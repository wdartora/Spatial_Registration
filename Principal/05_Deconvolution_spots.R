
# 1) ------------------->>>>>>>>>>>>>>>>>>>>>>
#Clear all
rm(list = ls())
###load libraries###
require(spacexr)
require(Matrix)
require(ggplot2)
require(hdf5r)
require(Seurat)


# Load Visium data directly from Space Ranger output directory
#VisiumData<-read.VisiumSpatialRNA("../outs/")
VisiumData<-read.VisiumSpatialRNA("W:/Data_base_VISIUM/MPI/data/ABCA7-5-a1-7_BIN1-5-a3-8/outs/")

# Create a list of barcodes from the column names of the count matrix
barcodes <- colnames(VisiumData@counts)

# Load single cell data from Cell Ranger output directory
Counts<-Read10X_h5("W:/Data_base_VISIUM/MPI/data/ABCA7-5-a1-7_BIN1-5-a3-8/outs/filtered_feature_bc_matrix.h5")

# Create new Seurat Object from the count input
count_SeuratObject<-CreateSeuratObject(Counts)

# Create a count matrix object and take a look at the matrix
#sc_counts <- count_SeuratObject@assays$RNA@counts
sc_counts <- GetAssayData(count_SeuratObject, assay = "RNA", slot = "counts")
sc_counts[1:5,1:5]


# Load reference
refdata=readRDS("C:/Users/wjd4002/Documents/William/Project/VISION/Proj03/data/allen_mop_2020.Rds")
cell_types=as.character(refdata@meta.data$class)
cell_types[which(cell_types=="Non-Neuronal")]=refdata@meta.data$subclass[which(cell_types=="Non-Neuronal")]


# Convert to factor data type
cell_types <- as.factor(cell_types) 
names(cell_types)=colnames(refdata)
head(cell_types)


###single-cell reference requires all cell types to have at least 25 members, so remove any cell types that have fewer than that
celltypetable=table(cell_types)
keepclasses=names(celltypetable)[which(celltypetable>=25)]
keepcells=which(cell_types %in% keepclasses)
newcell_types=as.factor(as.character(cell_types[keepcells]))
names(newcell_types)=names(cell_types)[keepcells]

# Create single cell reference object
SCreference <- Reference(refdata@assays$RNA@counts[,keepcells], newcell_types)#, sc_umis)

# runs faster with more cores
Sys.setenv("OPENBLAS_NUM_THREADS"=2)

# Create and run RCTD algorithm
myRCTD <- create.RCTD(VisiumData, SCreference, max_cores = 2)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")


#Results 
# Create the output directory in your working directory
resultsdir <- "W:/Data_base_VISIUM/Projects/Proj02/deconvolution/ABCA7-5-a1-7_BIN1-5-a3-8"
dir.create(resultsdir) 

# Create variables from the myRCTD object to plot results
barcodes <- colnames(myRCTD@spatialRNA@counts) # list of spatial barcodes
weights <- myRCTD@results$weights # Weights for each cell type per barcode

# Normalize per spot weights so cell type probabilities sum to 1 for each spot
norm_weights <- normalize_weights(weights) 
cell_type_names<-colnames(norm_weights) # List of cell types
dim(norm_weights)

# Look at cell type normalized weights 
subset_df <- as.data.frame(t(as.data.frame(norm_weights[1:2,])))
subset_df$celltypes <- rownames(subset_df); rownames(subset_df) <- NULL
subset_df[order(subset_df$AAACAGAGCGACTCCT-1, decreasing=T),]

# 1. Acessar os metadados de posição espacial:
coords <- myRCTD@spatialRNA@coords
head(coords)

# Plot simples pra ver onde cortar
coords$barcode <- rownames(coords)
ggplot(coords, aes(x = x, y = y)) +
  geom_point(size = 0.5) +
  scale_y_reverse() +
  coord_fixed() +
  theme_minimal() +
  ggtitle("Visualização para split do slide")


# 2. Identificar os dois cortes manualmente:
# Adiciona coluna com os barcodes e divide as regiões
coords$barcode <- rownames(coords)
coords$sample <- ifelse(coords$y < 115, "ABCA7-5-a1-7", "BIN1-5-a3-8")

# 3. Separar os pesos normalizados por corte:
norm_weights_df <- as.data.frame(norm_weights)
norm_weights_df$barcode <- rownames(norm_weights_df)
norm_weights_df <- merge(norm_weights_df, coords[, c("barcode", "x", "y", "sample")], by = "barcode")

# 4. Plotar separadamente por amostra:
library(dplyr)

for (cell_type in cell_type_names) {
  for (sample_name in unique(norm_weights_df$sample)) {
    df_sub <- norm_weights_df %>%
      filter(sample == sample_name)
    
    # Ivert plot
    df_sub$y <- -df_sub$y
    
    p <- ggplot(df_sub, aes(x = x, y = y, color = !!sym(cell_type))) +
      geom_point(size = 1.2) +
      scale_color_viridis_c(option = "C") +
      coord_fixed() +
      scale_y_reverse() +
      theme_void() +
      ggtitle(paste("Deconvolution:", cell_type, "-", sample_name))
    
    ggsave(filename = file.path(resultsdir, paste0("split_", sample_name, "_", cell_type, ".jpg")),
           plot = p, height = 5, width = 8, units = "in", dpi = 300)
  }
}