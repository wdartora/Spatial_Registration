
# =============================
### STEP 3 - LANDMARK, pretty much the same as that done for Visium
# DO IN BASE R, not R STUDIO


##############

library(tidyverse)
library(EBImage)

## Creating the list conver_IDs.csv
sample_IDs <- data.frame(
  sample_ID = c(
    "PICALM-10-a1-7", "ABCA7-3-a1-8", "BIN1-6-a3-8", "ABCA7-8-a3",
    "BIN1-5-a1", "BIN1-5-a3-8", "ABCA7-5-a1-7", "PICALM-5-a1", 
    "WTP-5-a7-8", "WTP-5-a2-8"
  ),
  desi = c(
    "10a1_PICALM", "3a1_ABCA7", "6a3_BIN1", "8a3_ABCA7",
    "5a1_BIN1", "5a3_BIN1", "5a1_ABCA7", "5a1_PICALM",
    "5a7_WTP", "5a2_WTP"
  )
)

write.csv(sample_IDs, file = paste0(output_r_dir, "convert_IDs.csv"), row.names = FALSE, quote = FALSE)



setwd("W:/Project/MPI_Registration/Proj_Spatial_Registration/results/desi_output")
output_r_dir <- "W:/Project/MPI_Registration/Proj_Spatial_Registration/results/r_output/"
output_fig_dir <- "W:/Project/MPI_Registration/Proj_Spatial_Registration/results/r_output/"


files <- c(list.files(output_r_dir, full.names=T,pattern="pos.+_4"),
           list.files(output_r_dir, full.names=T,pattern="neg.+_3"))

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

table(desi_ids %in% convert_IDs$desi)


# Add some important meta data information
# Add some important meta data information
landmarks <- y_df
to_csv <- function(landmarks, file_name){
  landmarks <- as.data.frame(landmarks)
  landmarks$Sample <- as.numeric(str_extract(rownames(landmarks),"[0-9]+"))
  landmarks$Coordinates <- str_extract(rownames(landmarks),"X|Y")
  landmarks$image_file <- rep(basename(files), each = 2)
  landmarks$sample_ID <- rep(str_extract(basename(files),"^[0-9]+a[0-9]_[a-zA-Z0-9]+"), each = 2)
  landmarks$charge <- rep(str_extract(files,"pos|neg1|neg2"),each=2)
  landmarks %>% 
    relocate(Sample,Coordinates,sample_ID,charge,
             image_file) %>%
    arrange(Sample, Coordinates,
            sample_ID,charge
    )%>%
    write.csv(file_name)
}

to_csv(landmarks, paste0(output_r_dir,"Query_Reference_Landmarks_for_DESI_Data.csv"))




