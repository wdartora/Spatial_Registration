# Install packages
#install.packages("qs")

# Load Libraries
library(qs)

# Load Seurat obj (rds)
objeto <- readRDS("C:/Users/wjd4002/Documents/William/GitHub/Spatial_Registration/Principal/allen_mop_2020.rds")

# Save in QS format
qsave(objeto, file = "C:/Users/wjd4002/Documents/William/GitHub/Spatial_Registration/Principal/allen_mop_2020_archive.qs", preset = "archive")



# Load QS format
library(qs)
objeto2 <- qread("C:/Users/wjd4002/Documents/William/GitHub/Spatial_Registration/Principal/allen_mop_2020_archive.qs")


