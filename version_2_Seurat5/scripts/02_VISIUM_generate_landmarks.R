# 02_generate_landmarks.R

# Run this script OUTSIDE RStudio (via terminal)

library(tidyverse)
library(EBImage)

output_fig_dir <- "results/plots/"
output_r_dir <- "results/r_output/"

# List hires images saved in previous step
files <- list.files(output_fig_dir, full.names = TRUE, pattern = "_hires_image.png$")



# --- Manual landmarking with visual feedback ---
message("Opening images for manual landmark marking...")
y <- list()
for(i in 1:length(files)){
  img <- readImage(files[i])
  display(img, method="raster")
  y[[i]] <- locator() # Click the pre-determined landmarks in order
  n_points <- length(y[[i]]$x)
  cat("✔️  You clicked", n_points, "point(s).\n\n")
}


# Save dataframe with landmarks
message("Saving clicked landmarks...")
y_df <- map_dfr(1:length(y), function(i)
  data.frame(X=y[[i]]["x"], Y = y[[i]]["y"]) %>% t() %>% as.data.frame())

rownames(y_df) <- paste(rep((1:length(y)),each=2),rep(c("Q_X","Q_Y"),length(y)), sep="")
# Add some important meta data information
landmarks <- y_df
to_csv <- function(landmarks, file_name){
  landmarks <- as.data.frame(landmarks)
  landmarks$Sample <- as.numeric(str_extract(rownames(landmarks),"[0-9]+"))
  landmarks$Coordinates <- str_extract(rownames(landmarks),"X|Y")
  #landmarks$image_file <- rep(files,each=2)
  landmarks$image_file <- rep(basename(files), each = 2)
  #landmarks$sample_ID <- rep(str_remove(files,"_hires_image.png"),each=2)
  landmarks$sample_ID <- rep(basename(str_remove(files, "_hires_image.png")), each = 2)
  landmarks %>% 
    relocate(Sample,Coordinates,sample_ID,#,charge,
             image_file) %>%
    arrange(Sample, Coordinates,sample_ID#,charge
    )%>%
    write.csv(file_name)
}


to_csv(y_df, paste0(output_r_dir,"Query_Reference_Landmarks_for_Visium_Data.csv"))

#  --- Reference image landmarking ---
message("Opening reference image for landmarks...")
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



message("\n Landmarks successfully saved!")
message(" Now manually open the following notebook in Jupyter or VSCode:")
message("   scripts/03_image_registration.ipynb")
message(" Visually inspect the landmarks and run the notebook to the end to generate:")
message("   - PDF comparing the transformed images")
message("   - CSV file with registered spot coordinates (py_output/spots_coords_regions.csv)")
message("\nAfter that, return to R and run 04_annotate_seurat.R to complete the pipeline.")
