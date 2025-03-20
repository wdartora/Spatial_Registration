#.libPaths(c("/library", "/rpackages"))

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(magrittr)
library(EBImage)
setwd("VizDezI")

output_fig_dir <- "r_output/"
output_r_dir <- "r_output/"
source("scripts/backup_functions.R")
options(future.globals.maxSize= 891289600)

###################


######## STEP 1: Loading Visium data using Seurat v5


###################

se <- Load10X_Spatial_2("raw_Visium", assay = "Spatial")
SpatialDimPlot(se)+NoLegend()
# note, two sections, needs to be split apart

# add spot coordinates to meta informations
se$x <-  se@images$slice1@boundaries$centroids@coords[,1]
se$y <-  se@images$slice1@boundaries$centroids@coords[,2]

#plot multiple times to find boundary of imagerow between both slices,
# split object into two slices
se1 <- subset(se,  x < 5000, invert =F)
se2 <- subset(se,  x < 5000, invert =T)

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

sct_transform(se1, "12mo_270a8")
sct_transform(se2, "12mo_299a3")

#####Merge all objects together and re-run pipeline all together

se <- readRDS(paste0(output_r_dir, "12mo_270a8",".rds"))
paths <- c("12mo_270a8","12mo_299a3")

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
img <- readImage("raw_Visium/spatial/tissue_hires_image.png")
coords_sub <- coords
coords_sub$imagerow = coords_sub$x*coords_sub$hires_sf
coords_sub$imagecol = coords_sub$y*coords_sub$hires_sf
sf <- coords_sub %>% pull(hires_sf) %>% unique

max(coords_sub$imagerow[grepl("12mo_270a8",coords_sub$slice)])
min(coords_sub$imagerow[grepl("12mo_299a3",coords_sub$slice)])

# Split the image literally 
bottom_cutoff <- 975
top_cutoff <- 990 # NOTE DIFFERENT NUMBER

img1 <- img[,1:bottom_cutoff,]
img2 <- img[,top_cutoff:ncol(img),]

png(paste0(output_fig_dir,"12mo_270a8_hires_image.png"), height=dim(img1)[2],width=dim(img1)[1],res=800)
display(img1, method="raster")
coords_sub2 <- coords_sub %>% filter(slice  == "12mo_270a8")
#points(coords_sub2$imagecol,coords_sub2$imagerow,pch=3,cex=.1)
dev.off()

png(paste0(output_fig_dir,"12mo_299a3_hires_image.png"), height=dim(img2)[2],width=dim(img2)[1],res=800)
display(img2, method="raster")
coords_sub2 <- coords_sub %>% filter(slice == "12mo_299a3")
#points(coords_sub2$imagecol,coords_sub2$imagerow-top_cutoff,pch=3,cex=.1)
dev.off()

# Re-adjust coordinate system for the bottom half image
coords_spots <- coords_sub %>%
  mutate(bottom_cutoff = bottom_cutoff)%>%
  mutate(top_cutoff = top_cutoff) %>%
  mutate(imagerow = ifelse(slice == "12mo_299a3", imagerow-top_cutoff, imagerow)) 

coords_spots  %>% write.csv(paste0(output_r_dir,"spots_coords_hi_res.csv"),quote=F, row.names=F)

########################


#### STEP 3 : LANDMARKS - DO THIS IN R NOT R STUDIO


###################

library(tidyverse)
library(EBImage)
setwd("VizDezI")

output_fig_dir <- "r_output/"
output_r_dir <- "r_output/"

files <- list.files(output_r_dir, full.names=T,pattern="_hires_image.png$")


# This next step is manual landmarking, ensure the # of landmarks is consistent
# across all images, see align_images/Landmarks_12.png
# This next step is manual! You will need to click on the image once it pops up, and
# click stop in order to go to the next image
y <- list()
for(i in 1:length(files)){
  img <- readImage(files[i])
  display(img, method="raster")
  y[[i]] <- locator() # Click the pre-determined landmarks in order
}

# Creates data frame with two rows per brain slice, for the X & Y coordinates of the landmatks
y_df <- map_dfr(1:length(y), function(i)
  data.frame(X=y[[i]]["x"], Y = y[[i]]["y"]) %>% t() %>% as.data.frame())
rownames(y_df) <- paste(rep((1:length(y)),each=2),rep(c("Q_X","Q_Y"),length(y)), sep="")

# Add some important meta data information
to_csv <- function(landmarks, file_name){
  landmarks <- as.data.frame(landmarks)
  landmarks$Sample <- as.numeric(str_extract(rownames(landmarks),"[0-9]+"))
  landmarks$Coordinates <- str_extract(rownames(landmarks),"X|Y")
  landmarks$sample_ID <- rep(str_remove_all(files,paste0(output_fig_dir,"|_hires_image.png|[A-Z]+[0-9]+_")),each=2)
  landmarks %>% 
    relocate(Sample,Coordinates,sample_ID) %>%
    arrange(Sample, Coordinates,sample_ID)%>%
    write.csv(file_name)
}

to_csv(y_df, paste0(output_r_dir,"Query_Reference_Landmarks_for_Visium_Data.csv"))

# Do the same now for the reference image
z <- list()
img <- readImage("reference/40x_Wt-317-a4_6mo_(1)_Visium_2A_nl.png")
display(img, method="raster")
z[[1]] <- locator()

z_df <- data.frame(X=z[[1]]$x, Y = z[[1]]$y) %>% t() %>% as.data.frame()
rownames(z_df) <- c("0R_X","0R_Y") 
z_df$Sample <- 0
z_df$Coordinates <- str_extract(rownames(z_df),"X|Y")
z_df$sample_ID <- "reference"
z_df %>% 
  relocate(Sample,Coordinates,sample_ID) %>%
  arrange(Sample, Coordinates,sample_ID)%>%
  write_csv(paste0(output_r_dir,"Reference_Landmarks_for_Spatial_Registration.csv"))

################



#########  STEP 4 - PROCEED AFTER PYTHON REGISTRATION



##############

##############


## STEP 5: MERGED REGISTRATION OUTPUT TO SEURAT OBJECT


################
library(Seurat)
library(tidyverse)

# Import Visium data object and region file from Python script

se <- readRDS( paste0(output_r_dir, "Visium_Seurat_object.rds"))
se$cellid <- paste0(se$slice,"#",colnames(se))
region_spots <- read.csv(paste0("py_output/spots_coords_regions.csv")) %>%
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
