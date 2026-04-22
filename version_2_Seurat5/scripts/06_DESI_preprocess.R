# =============================
# Description: Normalize DESI files, convert to .mtx, generate grayscale PNGs
# =============================
# 06_desi_format_and_images.R

library(tidyverse)
library(stringr)
library(MBA)
library(EBImage)

project_root <- "C:/Users/darto/Documents/William/Cornell/Project/Spatial_Registration/Proj02/MPI_Registration"
desi_dir <- file.path(project_root, "results", "desi_output")
output_r_dir <- file.path(project_root, "results", "r_output")
output_fig_dir <- file.path(project_root, "results", "plots", "desi")

dir.create(output_r_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

setwd(desi_dir)

# =============================
# STEP 1 - Format DESI matrix files into .mtx
# =============================

mat_files <- list.files(
  desi_dir,
  pattern = "_mat\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

list_files <- list.files(
  desi_dir,
  pattern = "_ions_list\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

# Match each matrix file to its corresponding ion list by filename stem
mat_info <- tibble(mat_file = mat_files) %>%
  mutate(
    rel_path = str_remove(mat_file, paste0("^", fixed(desi_dir), "[/\\\\]?")),
    sample_stub = str_replace(rel_path, "_mat\\.csv$", "")
  )

list_info <- tibble(ion_file = list_files) %>%
  mutate(
    rel_path = str_remove(ion_file, paste0("^", fixed(desi_dir), "[/\\\\]?")),
    sample_stub = str_replace(rel_path, "_ions_list\\.csv$", "")
  )

file_pairs <- inner_join(mat_info, list_info, by = "sample_stub")

message("Number of matched DESI matrix/ion-list pairs: ", nrow(file_pairs))

for (i in seq_len(nrow(file_pairs))) {
  file <- file_pairs$mat_file[i]
  ion_file <- file_pairs$ion_file[i]
  
  message("Processing matrix pair: ", basename(file))
  
  # Read files
  mat <- read.csv(file, check.names = FALSE)
  ion_list <- read.csv(ion_file, check.names = FALSE) %>%
    drop_na(mass)
  
  # Remove possible index column added by pandas/R export
  if (ncol(mat) > 0 && names(mat)[1] %in% c("", "X", "Unnamed: 0")) {
    # only remove if first column is clearly an index and not the real X coordinate
    if (!("Y" %in% names(mat)[1:2])) {
      mat <- mat[, -1, drop = FALSE]
    }
  }
  
  # Ensure first two columns are X and Y
  if (!all(c("X", "Y") %in% names(mat)[1:2])) {
    warning("First two columns are not X and Y in: ", basename(file))
  }
  
  # Number of ion columns in the matrix
  n_ions_mat <- ncol(mat) - 2
  
  if (n_ions_mat <= 0) {
    warning("No ion columns found in: ", basename(file))
    next
  }
  
  if (nrow(ion_list) < n_ions_mat) {
    warning(
      "Ion list shorter than matrix ion columns: ", basename(file),
      " | n_ions_mat = ", n_ions_mat,
      " | nrow(ion_list) = ", nrow(ion_list)
    )
    next
  }
  
  # Keep only as many ions as exist in the matrix
  ion_list_use <- ion_list[1:n_ions_mat, , drop = FALSE]
  
  # Rename ion columns using mass and retention time
  colnames(mat)[3:ncol(mat)] <- paste0("m_", ion_list_use$mass, ":", ion_list_use$rt)
  
  # Normalize coordinates into integers with minimum value of 0
  mat[[1]] <- ((mat[[1]] - 0.02) / 0.06) + 17
  mat[[1]] <- mat[[1]] - min(mat[[1]], na.rm = TRUE)
  
  mat[[2]] <- ((mat[[2]] + 0.02) / 0.06)
  mat[[2]] <- mat[[2]] - min(mat[[2]], na.rm = TRUE)
  
  # Save .mtx next to original file
  out_mtx <- str_replace(file, "_mat\\.csv$", "_matrix_formatted.mtx")
  
  write.table(
    mat,
    file = out_mtx,
    quote = FALSE,
    row.names = FALSE
  )
}

# Check created files
files_mtx <- list.files(
  desi_dir,
  pattern = "_matrix_formatted\\.mtx$",
  recursive = TRUE,
  full.names = TRUE
)

message("Number of formatted matrix files created: ", length(files_mtx))

if (length(files_mtx) > 0) {
  print(files_mtx[1:min(10, length(files_mtx))])
}

# =============================
# STEP 2 - Create grayscale PNGs for selected ions
# =============================

files <- list.files(
  desi_dir,
  pattern = "_matrix_formatted\\.mtx$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No _matrix_formatted.mtx files were found. STEP 1 did not generate outputs.")
}

for (file in files) {
  message("Creating ion images for: ", basename(file))
  
  mat <- read.table(file, header = TRUE, check.names = FALSE)
  
  # Create output subdirectory preserving polarity structure
  rel_path <- str_remove(file, paste0("^", fixed(desi_dir), "[/\\\\]?"))
  rel_dir <- dirname(rel_path)
  fig_subdir <- file.path(output_fig_dir, rel_dir)
  dir.create(fig_subdir, recursive = TRUE, showWarnings = FALSE)
  
  # Use the first two ion columns after X and Y for quick inspection
  ion_cols <- 3:min(4, ncol(mat))
  
  for (j in ion_cols) {
    df <- as.data.frame(mat[, c(1, 2, j), drop = FALSE])
    colnames(df) <- c("X", "Y", "intensity")
    
    # Skip empty ion images
    if (all(is.na(df$intensity)) || max(df$intensity, na.rm = TRUE) == 0) {
      next
    }
    
    # Min-max style normalization to grayscale
    df$intensity <- (df$intensity * 255) / max(df$intensity, na.rm = TRUE)
    
    # Convert three-column data (X, Y, ion) to 2D image matrix
    y <- df %>%
      arrange(desc(X)) %>%
      pivot_wider(names_from = "X", values_from = "intensity") %>%
      dplyr::select(-Y) %>%
      as.matrix()
    
    out_png <- file.path(
      fig_subdir,
      str_replace(
        basename(file),
        "_matrix_formatted\\.mtx$",
        paste0("_to_image_", j, ".png")
      )
    )
    
    png(
      filename = out_png,
      width = 2 * max(mat$Y, na.rm = TRUE),
      height = 2 * max(mat$X, na.rm = TRUE)
    )
    
    par(mar = rep(0, 4))
    image(y, axes = FALSE, col = grey(seq(0, 1, length = 256)))
    dev.off()
  }
}

message("DESI matrix formatting and ion image generation completed.")