# =============================
# STEP 3 - LANDMARK
# Pretty much the same as that done for Visium
# DO IN BASE R, not RStudio
# =============================

library(tidyverse)
library(EBImage)
library(stringr)

# -----------------------------
# Project paths
# -----------------------------
project_root <- "C:/Users/darto/Documents/William/Cornell/Project/Spatial_Registration/Proj02/MPI_Registration"
output_r_dir <- file.path(project_root, "results", "r_output")
output_fig_dir <- file.path(project_root, "results", "plots", "desi")

dir.create(output_r_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Create convert_IDs.csv
# Full table for all current samples
# -----------------------------
sample_IDs <- data.frame(
  sample_ID = c(
    "ABCA7-3-a1-8",
    "BIN1-6-a3-8",
    "ABCA7-5-a1-7",
    "BIN1-5-a3-8",
    "ABCA7-8-a3",
    "BIN1-5-a1",
    "WTP-2-a1-8",
    "PICALM-10-a1-7",
    "WTP-5-a2-8",
    "PICALM-10-a2-8",
    "WTP-5-a7-8",
    "PICALM-5-a1-8"
  ),
  desi = c(
    "3a1_ABCA7",
    "6a3_BIN1",
    "5a1_ABCA7",
    "5a3_BIN1",
    "8a3_ABCA7",
    "5a1_BIN1",
    "2a1_WTP",
    "10a1_PICALM",
    "5a2_WTP",
    "10a2_PICALM",
    "5a7_WTP",
    "5a1_PICALM"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  sample_IDs,
  file = file.path(output_r_dir, "convert_IDs.csv"),
  row.names = FALSE,
  quote = FALSE
)

# -----------------------------
# Optional: define which files to landmark
# Keep this style close to your old script
# -----------------------------
# Example choices:
# files <- c(
#   list.files(output_fig_dir, full.names = TRUE, recursive = TRUE, pattern = "pos.*_to_image_4\\.png$"),
#   list.files(output_fig_dir, full.names = TRUE, recursive = TRUE, pattern = "neg1.*_to_image_3\\.png$")
# )

# For your current use, if you want only POS and image 3:
files <- list.files(
  output_fig_dir,
  full.names = TRUE,
  recursive = TRUE,
  pattern = "pos.*_to_image_3\\.png$"
)

# If you want only one pair, uncomment and edit:
# files <- files[grepl("3a1_ABCA7|6a3_BIN1", basename(files))]

files <- unique(files)

if (length(files) == 0) {
  stop("No DESI PNG files found for landmarking.")
}

print(files)

# -----------------------------
# Manual landmarking with visual feedback
# -----------------------------
message("Opening images for manual landmark marking...")
y <- list()

for (i in seq_along(files)) {
  img <- readImage(files[i])
  display(img, method = "raster")
  y[[i]] <- locator()  # Click the pre-determined landmarks in order
  n_points <- length(y[[i]]$x)
  cat("✔ You clicked", n_points, "point(s).\n\n")
}

# -----------------------------
# Save dataframe with landmarks
# -----------------------------
message("Saving clicked landmarks...")

y_df <- map_dfr(
  seq_along(y),
  function(i) {
    data.frame(X = y[[i]][["x"]], Y = y[[i]][["y"]]) %>%
      t() %>%
      as.data.frame()
  }
)

rownames(y_df) <- paste(
  rep(seq_along(y), each = 2),
  rep(c("Q_X", "Q_Y"), length(y)),
  sep = ""
)

# -----------------------------
# Build landmarks CSV
# -----------------------------
landmarks <- y_df

to_csv <- function(landmarks, file_name) {
  landmarks <- as.data.frame(landmarks)
  landmarks$Sample <- as.numeric(str_extract(rownames(landmarks), "[0-9]+"))
  landmarks$Coordinates <- str_extract(rownames(landmarks), "X|Y")
  landmarks$image_file <- rep(basename(files), each = 2)
  landmarks$sample_ID <- rep(str_extract(basename(files), "[0-9]+a[0-9]+_[A-Za-z0-9]+"), each = 2)
  landmarks$charge <- rep(str_extract(basename(files), "pos|neg1|neg2"), each = 2)
  
  landmarks %>%
    relocate(Sample, Coordinates, sample_ID, charge, image_file) %>%
    arrange(Sample, Coordinates, sample_ID, charge) %>%
    write.csv(file_name, row.names = FALSE)
}

to_csv(
  landmarks,
  file.path(output_r_dir, "Query_Reference_Landmarks_for_DESI_Data.csv")
)

message("DESI landmarks successfully saved!")
message("Output file: results/r_output/Query_Reference_Landmarks_for_DESI_Data.csv")

