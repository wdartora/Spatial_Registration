# =============================
# Description: Normalize DESI files, convert to .mtx, generate grayscale PNGs
# =============================


library(tidyverse)
library(stringr)
library(MBA)
library(EBImage)

setwd("W:/Project/MPI_Registration/Proj_Spatial_Registration/results/desi_output")
output_r_dir <- "W:/Project/MPI_Registration/Proj_Spatial_Registration/results/r_output/"
output_fig_dir <- "W:/Project/MPI_Registration/Proj_Spatial_Registration/results/r_output/"

############


######## STEP 1 - Code below formats data into producing into a black and white image.


#############

mat_files <- list.files(".","_mat.csv$", recursive=T)
list_files <- list.files(".","_ions_list.csv$", recursive=T)

for(i in 1:length(mat_files)){
  
  file <- mat_files[i]
  ion_file <- list_files[i]
  
  #### Reads files, one is matrix data, other is ion list
  mat <- read.csv(file)[,-1]
  ion_list <- read.csv(ion_file) %>% drop_na(mass)
  
  # +2 in the line below comes from the fact the matrix data should have two more columns for X,Y coordinates
  if(ncol(mat)==(nrow(ion_list)+2)){
    mat <- mat[,c(1,2,2+ion_list$X)]
    colnames(mat)[3:ncol(mat)] <- paste0("m_", paste0(ion_list$mass,":",ion_list$rt))
    
    #Normalizing coordinates into integers with minimum value of 0
    mat[[1]] <- ((mat[[1]] -.02)/.06)+17
    mat[[1]] <- mat[[1]] - min(mat[[1]])
    mat[[2]] <- ((mat[[2]] +.02)/.06) 
    mat[[2]] <- mat[[2]] - min(mat[[2]])
    
    mat %>% as.matrix() %>% write.table(str_replace(mat_files[i],"_mat.csv",
                                                    "_matrix_formatted.mtx"), quote=F)
  }
}




#######################


#### STEP 2 - create images of different ions that best reflect cortical structure


############

files <- list.files("W:/Project/MPI_Registration/Proj_Spatial_Registration/results/desi_output","_matrix_formatted.mtx",recursive=T)

for(file in files){
  mat <- read.table(paste0("W:/Project/MPI_Registration/Proj_Spatial_Registration/results/desi_output/",file))
  # Create images from matrix data for first couple of ions, using min-max normalization
  for(i in 3:4){
    df <- as.data.frame(mat[,c(1,2, i+2)])
    df[,3] <- (df[,3]*255)/max(df[,3])
    
    # Converts three column data (X,Y, ion) to 2-D image matrix
    y <- df %>% arrange(desc(X)) %>% 
      pivot_wider(names_from = "X", values_from = contains("m")) %>%
      dplyr::select(-Y) %>% as.matrix()
    
    # Creates the image
    png(paste0(output_fig_dir,str_replace(file, "matrix_formatted.mtx",paste0("to_image_",i,".png"))),
        width = 2*max(mat$Y), height=2*max(mat$X) )  #### note we double size of image
    par(mar = rep(0, 4))
    image(y, axes = FALSE, col = grey(seq(0, 1, length = 256)))
    dev.off()
  }
}

# look through the positive and negative ions separately to find the ion that best outlines brain in desi images
# select one representative positive and negative ion respectively



