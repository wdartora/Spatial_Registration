library(tidyverse)
library(stringr)
library(MBA)
library(EBImage)

setwd("VizDezI")
output_r_dir <- "r_output/"
output_fig_dir <- "r_output/"

############


######## STEP 1 - Code below formats data into producing into a black and white image.


#############

mat_files <- list.files(".","HDMS_TL_mat.csv$", recursive=T)
list_files <- list.files(".","HDMS_TL_ions_list.csv$", recursive=T)

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

files <- list.files("raw_DESI","TL_matrix_formatted.mtx",recursive=T)

for(file in files){
  mat <- read.table(paste0("raw_DESI/",file))
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


#############


### STEP 3 - LANDMARK, pretty much the same as that done for Visium
# DO IN BASE R, not R STUDIO


##############

library(tidyverse)
library(EBImage)

setwd("VizDezI")
output_r_dir <- "r_output/"
output_fig_dir <- "r_output/"


files <- c(list.files(output_r_dir, full.names=T,pattern="pos.+_4"),
           list.files(output_r_dir, full.names=T,pattern="neg.+_3"))

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
  landmarks$image_file <- rep(str_remove(files,paste0(output_fig_dir)),each=2)
  landmarks$sample_ID <- rep(c("12mo_270a8","12mo_299a3","12mo_299a3","12mo_270a8"),each=2)
  landmarks$charge <- rep(str_extract(files,"pos|neg"),each=2)
  landmarks %>% 
    relocate(Sample,Coordinates,sample_ID,charge, image_file) %>%
    arrange(Sample, Coordinates,sample_ID,charge)%>%
    write.csv(file_name)
}

to_csv(y_df, paste0(output_r_dir,"Query_Reference_Landmarks_for_DESI_Data.csv"))




##############


#############


### STEP 4- PROCEED AFTER PYTHON ALIGNMENT


##############

###############


##### STEP 5 - LINK ION DATA TO VISIUM SPOTS


################

# For positive ion data, import matrix data and add some meta data
files = list.files(".","pos.+matrix_formatted.mtx",recursive=T)
all_images = map_dfr(files, function(file) 
  mat <- read.table(file) %>% as.data.frame() %>% mutate(sample_ID = parse_number(str_extract(file, "Laura_[0-9]{3,4}")),
                                                         direction = str_extract(file,"pos|neg"))
)

all_images <- all_images %>% relocate(sample_ID,direction) %>% mutate(sample_ID=as.character(sample_ID)) %>%
  mutate(X = round(X), Y= round(Y))
all_images$sample_ID[all_images$sample_ID==270] <- "12mo_270a8"
all_images$sample_ID[all_images$sample_ID=="299"] <- "12mo_299a3"
all_images[is.na(all_images)] <- 0


# Import spots-to-pixels links file, round up transformed pixel coordinates and then
# join these linked pixels to Visium spots 
spots_pixel <- read.csv("py_output/spots_coords_ion_pixels_pos.csv") %>%
  mutate(transformed_row = round(transformed_row/2),
         transformed_col = round(transformed_col/2)) %>%
  mutate(sample_ID = as.character(slice))
spots_pixel <- spots_pixel %>% left_join(all_images %>% mutate(sample_ID = sample_ID),
                                         by = c("sample_ID"="sample_ID","transformed_col"="Y","transformed_row"="X")) 
spots_pixel[is.na(spots_pixel)] <- 0
spots_pixel %>% write.csv(paste0(output_r_dir,"spots_coords_ion_intensities_pos.csv"))



# For positive ion data, import matrix data and add some meta data

files = list.files(".","neg.+matrix_formatted.mtx",recursive=T)
all_images = map_dfr(files, function(file) 
  mat <- read.table(file) %>% as.data.frame() %>% mutate(sample_ID = parse_number(str_extract(file, "Laura_[0-9]{3,4}")),
                                                         direction = str_extract(file,"pos|neg"))
)
all_images <- all_images %>% relocate(sample_ID,direction) %>% mutate(sample_ID=as.character(sample_ID)) %>%
  mutate(X = round(X), Y= round(Y))
all_images$sample_ID[all_images$sample_ID==270] <- "12mo_270a8"
all_images$sample_ID[all_images$sample_ID=="299"] <- "12mo_299a3"
all_images[is.na(all_images)] <- 0

# Import spots-to-pixels links file, round up transformed pixel coordinates and then
# join these linked pixels to Visium spots 

spots_pixel <- read.csv("py_output/spots_coords_ion_pixels_neg.csv") %>%
  mutate(transformed_row = round(transformed_row/2),
         transformed_col = round(transformed_col/2)) %>%
  mutate(sample_ID = as.character(slice))
spots_pixel <- spots_pixel %>% left_join(all_images %>% mutate(sample_ID = sample_ID),
                                         by = c("sample_ID"="sample_ID","transformed_col"="Y","transformed_row"="X")) 
spots_pixel[is.na(spots_pixel)] <- 0
spots_pixel %>% write.csv(paste0(output_r_dir,"spots_coords_ion_intensities_neg.csv"))

rm(spots_pixel)

############


############ STEP 6 - NORMALIZE/TRANSFORM ION DATA


###########

library(readxl)

#  Import linked ion data and filtered ion list, we will also rename ions to indicate
#  positive or negative charge

spots_p <- read.csv(paste0(output_r_dir,"spots_coords_ion_intensities_pos.csv"))
colnames(spots_p) <- str_replace_all(colnames(spots_p),"m_","p_")
spots_n <- read.csv(paste0(output_r_dir,"spots_coords_ion_intensities_neg.csv"))
colnames(spots_n) <- str_replace_all(colnames(spots_n),"m_","n_")

spots_data <- cbind(spots_p,spots_n)
# Combine metabolomics data and filter using pre-selected list (based on noise) OR

one_sd <- read_xlsx("reference/one_sd_filtered_ions_no_noise_POS.xlsx")
colnames(one_sd) <- c("X","X0")
one_sd <- one_sd %>% mutate(X0 = parse_number(X0)) %>% as.data.frame()
two_sd <- read.csv("reference/two_sd_filtered_ions_no_noise_pos.csv")

pattern_input <- paste0(paste0("p_",c(one_sd$X0, two_sd$X0)),
                        collapse="|")
pos_pattern_input <- pattern_input

one_sd <- read_xlsx("reference/one_sd_filtered_ions_no_noise_NEG.xlsx")
colnames(one_sd) <- c("X","X0")
one_sd <- one_sd %>% mutate(X0 = parse_number(X0)) %>% as.data.frame()
two_sd <- read.csv("reference/two_sd_filtered_ions_no_noise_neg.csv")

pattern_input <- paste0(paste0("n_",c(one_sd$X0, two_sd$X0)),
                        collapse="|")
neg_pattern_input <- pattern_input
pattern_input <- paste0(pos_pattern_input,"|",neg_pattern_input)
pattern_input <- str_replace_all(pattern_input,"\\.","\\\\\\.")

spots_data <- cbind(spots_p %>% select(contains("p_")),
                    spots_n %>% select(contains("n_")))
spots_data <- spots_data[,grepl(pattern_input, colnames(spots_data))]

### To filter by considering highest variable ions
cv_n <- spots_n[,grepl("slice|^n_", colnames(spots_n))] %>% group_by(slice) %>%
  summarise_if(is.numeric, function(col) mean(col,na.rm=T)/sd(col,na.rm=T)) %>%
  ungroup() %>%
  summarise_if(is.numeric, function(var) mean(var,na.rm=T))

cv_p <- spots_p[,grepl("slice|p_", colnames(spots_p))] %>% group_by(slice) %>%
  summarise_if(is.numeric, function(col) mean(col,na.rm=T)/sd(col,na.rm=T)) %>%
  ungroup() %>%
  summarise_if(is.numeric, function(var) mean(var,na.rm=T))

cv_n <- t(cv_n)
cv_n <- names(na.omit(cv_n[cv_n[,1]>= quantile(cv_n[,1],.6,na.rm=T),]))
cv_p <- t(cv_p)
cv_p <- names(na.omit(cv_p[cv_p[,1]>= quantile(cv_p[,1],.6,na.rm=T),]))

#### Perform normalization and then combine
spots_n1 <- map_dfr(unique(spots_n$slice), function(x) {
  spots_n %>% select(any_of(c("slice",cv_n))) %>% filter(slice==x) %>%
    mutate_if(is.numeric, function(col) scale(col/max(col))) 
  
})
spots_p1 <- map_dfr(unique(spots_p$slice), function(x) {
  spots_p %>% select(any_of(c("slice",cv_p))) %>% filter(slice==x) %>%
    mutate_if(is.numeric, function(col) scale(col/max(col))) 
  
})

spots_data <- cbind(spots_p1,spots_n1)
############


############ STEP 7 - ADD ION DATA TO RNA OBJECT


###########
library(viridis)
se <- readRDS(paste0(output_r_dir, "Visium_Seurat_object.rds"))

# Add normalized ion data to metadata
se@meta.data <- cbind(se@meta.data, spots_data)

# Plot ion data from Visium tissue section
sdp <- lapply(paths, function(x)
  SpatialFeaturePlot(se, feature ="p_132.0766.32.03",images = x,pt.size.factor = 2.85)+theme(aspect.ratio=3000/5000)+ggtitle(x)+
    theme(legend.position="none")+viridis::scale_fill_viridis()

)
sdp <- gridExtra::grid.arrange(grobs = sdp,nrow=2)

