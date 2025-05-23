rm(list=ls())
# --- Pacotes necessários ---
library(Seurat)
library(ggplot2)
library(patchwork)
library(magrittr)
library(EBImage)
library(gridExtra)
library(dplyr)

# --- Diretório de trabalho onde estão scripts e saídas --- /////////////
setwd("Y:/William/Data_base_VISIUM/Projects/Proj02")

# --- Diretórios de saída para resultados ---
output_fig_dir <- "r_output6/"
output_r_dir <- "r_output6/"

# --- Carregar função Load10X_Spatial_2 ---
source("script/backup_functions.R")

# --- Aumentar limite de memória futura ---
options(future.globals.maxSize = 891289600)

########################
# STEP 1: Load Visium data
########################

# Caminho dos dados Visium do Space Ranger (pasta "outs") /////////////
visium_path <- "Y:/William/Data_base_VISIUM/MPI/data/ABCA7-5-a1-7_BIN1-5-a3-8/outs"

# Carregar objeto Seurat usando função customizada
se <- Load10X_Spatial_2(visium_path, assay = "Spatial")
SpatialDimPlot(se) + NoLegend()
# note, two sections, needs to be split apart

# add spot coordinates to meta informations
se$x <-  se@images$slice1@boundaries$centroids@coords[,1]
se$y <-  se@images$slice1@boundaries$centroids@coords[,2]


#-------------------------------------------->>
# --- Pré-visualização da imagem completa com eixos e pontos dos spots ---
# Pacote necessário

# Plotar coordenadas x e y dos spots com eixo real
ggplot(se@meta.data, aes(x = x, y = y)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_y_reverse() +
  scale_x_continuous(breaks = seq(0, max(se$x), by = 1000)) +  # Mostrar de 1000 em 1000
  theme_minimal() +
  labs(title = "Coordenadas dos spots (x vs y)", x = "x", y = "y") +
  geom_vline(xintercept = 8500, linetype = "dashed", color = "red")

#-------------------------------------------->>



#plot multiple times to find boundary of imagerow between both slices,
# split object into two slices /////////////
se1 <- subset(se,  x < 9400, invert =F)
se2 <- subset(se,  x < 10000, invert =T)

# Ensure no overlap
SpatialDimPlot(se1)+NoLegend()
SpatialDimPlot(se2)+NoLegend()

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
  saveRDS(se, file = paste0(output_r_dir, sample_name,".rds"))
  
  return()
  
}
# Criar diretórios de saída se não existirem
if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE)
if (!dir.exists(output_r_dir)) dir.create(output_r_dir, recursive = TRUE)

# Salvar os objetos separados /////////////
sct_transform(se1, "ABCA7-5-a1-7")
sct_transform(se2, "BIN1-5-a3-8")

#####Merge all objects together and re-run pipeline all together

se <- readRDS(paste0(output_r_dir, "ABCA7-5-a1-7",".rds"))
paths <- c("ABCA7-5-a1-7","BIN1-5-a3-8")

read_file <- function(path){
  print(path)
  se1 <- readRDS(paste0(output_r_dir,path,".rds"))
  se <<- merge(se, se1)
}
sapply(paths[-1], read_file)

# Now perform standard seurat pipeline, normalize, filter, dimension reduction
se <- SCTransform(se, assay = "Spatial",vars.to.regress="nCount_Spatial", verbose = FALSE)
DefaultAssay(se) <- "SCT"
se <- RunPCA(se, assay = "SCT", verbose = FALSE)
se <- FindNeighbors(se, reduction = "pca", dims = 1:30)
se <- FindClusters(se, verbose = FALSE)
se <- RunUMAP(se, reduction = "pca", dims = 1:30)

names(se@images) <- paths

# Visualize default clusters across tissue 
sdp <- lapply(paths, function(x)
  SpatialDimPlot(se, images = x)+theme(aspect.ratio=3000/5000)+ggtitle(x)+NoLegend()
  
)
sdp <- gridExtra::grid.arrange(grobs = sdp,nrow=2)
ggsave(paste0(output_fig_dir,"All_Spatial_DimPlot.png"),sdp, height=8,width=6)

# Save Visium data object
saveRDS(se, paste0(output_r_dir, "Visium_Seurat_object.rds"))

#######################################


##### STEP 2: Extract hi-res images and spot coordinates for Seurat 5 objects


#########################################
library(purrr) # //////////////////////
library(dplyr)  # também necessário para %>% e relocate() ////////////////////

# Compile spatial spot center coordinates as well as additional metainformation, like scale factor
coords <- map_dfr(paths, function(x){
  coords <- data.frame(se@images[[x]]@boundaries$centroids@coords)
  rownames(coords) <- se@images[[x]]@boundaries$centroids@cells
  coords$slice <- x
  coords$hires_sf <- se@images[[x]]@scale.factors$hires
  return(coords %>% relocate(slice))
}
)
coords %>% write.csv(paste0(output_fig_dir,"Coords_slices.csv"))

#Single hires image contains both tissue slices, needs to be split 
# spot coordinates also different scale, also needs to be adjusted
#img <- readImage("raw_Visium/spatial/tissue_hires_image.png")

# Read the hires (high Resolution) image /////
img <- readImage("Y:/William/Data_base_VISIUM/MPI/data/ABCA7-5-a1-7_BIN1-5-a3-8/outs/spatial/tissue_hires_image.png")

coords_sub <- coords
coords_sub$imagerow = coords_sub$x*coords_sub$hires_sf
coords_sub$imagecol = coords_sub$y*coords_sub$hires_sf
sf <- coords_sub %>% pull(hires_sf) %>% unique

max(coords_sub$imagerow[grepl("ABCA7-5-a1-7",coords_sub$slice)])
min(coords_sub$imagerow[grepl("BIN1-5-a3-8",coords_sub$slice)])

# Split the image literally  ////////////////
ggplot(coords_sub, aes(x = imagecol, y = imagerow, color = slice)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_y_reverse() +  # Necessário para alinhar com a imagem (origem no topo)
  coord_fixed() +
  theme_light() +
  labs(title = "Spots scaled to hires image",
       x = "imagecol (X)", y = "imagerow (Y)") +
  theme(legend.position = "bottom")
ggsave("r_output6/Spots_scaled_to_hires_image_ABCA7-5-a1-7_BIN1-5-a3-8.png", width = 7, height = 6, dpi = 300)


bottom_cutoff <- 2100 # 975
top_cutoff <- 2000 # NOTE DIFFERENT NUMBER 990

img1 <- img[,1:bottom_cutoff,]
img2 <- img[,top_cutoff:ncol(img),]

png(paste0(output_fig_dir,"ABCA7-5-a1-7_hires_image.png"), height=dim(img1)[2],width=dim(img1)[1],res=800)
display(img1, method="raster")
coords_sub2 <- coords_sub %>% filter(slice  == "ABCA7-5-a1-7")
#points(coords_sub2$imagecol,coords_sub2$imagerow,pch=3,cex=.1)
dev.off()

png(paste0(output_fig_dir,"BIN1-5-a3-8_hires_image.png"), height=dim(img2)[2],width=dim(img2)[1],res=800)
display(img2, method="raster")
coords_sub2 <- coords_sub %>% filter(slice == "BIN1-5-a3-8")
#points(coords_sub2$imagecol,coords_sub2$imagerow-top_cutoff,pch=3,cex=.1)
dev.off()

# Re-adjust coordinate system for the bottom half image
coords_spots <- coords_sub %>%
  mutate(bottom_cutoff = bottom_cutoff)%>%
  mutate(top_cutoff = top_cutoff) %>%
  mutate(imagerow = ifelse(slice == "BIN1-5-a3-8", imagerow-top_cutoff, imagerow)) 

coords_spots  %>% write.csv(paste0(output_r_dir,"spots_coords_hi_res.csv"),quote=F, row.names=F)