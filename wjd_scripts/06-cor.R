# --- Configuração inicial ---
rm(list = ls())

library(Seurat)
library(dplyr)
library(readr)
library(tidyr)
library(matrixStats)
library(openxlsx)

# --- Função de limpeza de barcode ---
clean_barcode <- function(x) {
  x <- gsub("-1$", "", x)
  x <- gsub("_\\d+(?:_\\d+)*$", "", x)
  return(x)
}

# --- Caminhos base ---
seurat_dir <- "W:/Data_base_VISIUM/Projects/Proj01/data_seurat"
ion_file <- "W:/Data_base_VISIUM/Projects/Proj01/integration_outputs/visium_desi_files/spots_ion_data_0509.csv"
output_dir <- "W:/Data_base_VISIUM/Projects/Proj01/correlation_per_slice"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- Lista de arquivos .rds ---
seurat_files <- list.files(seurat_dir, pattern = "\\.rds$", full.names = TRUE)

# --- Ler dados de íons uma única vez ---
ion_data_raw <- read_csv(ion_file, show_col_types = FALSE)

# Renomear colunas do ion_data_raw
if ("...1" %in% colnames(ion_data_raw)) {
  colnames(ion_data_raw)[which(colnames(ion_data_raw) == "...1")] <- "barcode_ion"
}
slice_cols <- grep("^slice", colnames(ion_data_raw), value = TRUE)
if (length(slice_cols) > 0) {
  colnames(ion_data_raw)[which(colnames(ion_data_raw) == slice_cols[1])] <- "slice_ion"
}
ion_cols <- grep("^(p_|n1_|n2_)", colnames(ion_data_raw), value = TRUE)

# --- Loop para cada slice ---
for (seurat_path in seurat_files) {
  slice_id <- gsub("\\.rds$", "", basename(seurat_path))
  cat("\n🔁 Processando:", slice_id, "\n")
  
  # Load Seurat
  seurat_obj <- readRDS(seurat_path)
  barcodes <- rownames(seurat_obj@meta.data)
  expr <- seurat_obj@assays$Spatial@layers$counts
  colnames(expr) <- barcodes[1:ncol(expr)]
  rownames(expr) <- rownames(seurat_obj@assays$Spatial)
  
  # Filtrar para slice atual
  ion_data <- ion_data_raw %>% filter(slice_ion == slice_id)
  ion_data$barcode_ion <- clean_barcode(ion_data$barcode_ion)
  
  # Interseção de barcodes
  common_barcodes <- intersect(colnames(expr), ion_data$barcode_ion)
  cat("  • Barcodes em comum:", length(common_barcodes), "\n")
  if (length(common_barcodes) < 3) {
    cat("  ⚠️ Pulando por poucos barcodes.\n")
    next
  }
  
  # Subset
  expr <- expr[, common_barcodes, drop = FALSE]
  ion_mtx <- ion_data[match(common_barcodes, ion_data$barcode_ion), ion_cols, drop = FALSE]
  
  # Filtragem por variância
  expr <- expr[rowVars(as.matrix(expr)) > 1e-9, , drop = FALSE]
  ion_mtx <- ion_mtx[, apply(ion_mtx, 2, var, na.rm = TRUE) > 1e-9, drop = FALSE]
  if (nrow(expr) == 0 || ncol(ion_mtx) == 0) {
    cat("  ⚠️ Sem variância suficiente.\n")
    next
  }
  
  # Correlação
  expr_mat <- as.matrix(expr)
  ion_mat <- as.matrix(ion_mtx)
  cor_mat <- cor(t(expr_mat), ion_mat, method = "spearman")
  
  # p-valor
  pvals <- matrix(NA, nrow = nrow(expr_mat), ncol = ncol(ion_mat),
                  dimnames = list(rownames(expr_mat), colnames(ion_mat)))
  for (i in 1:nrow(expr_mat)) {
    for (j in 1:ncol(ion_mat)) {
      x <- as.numeric(expr_mat[i, ])
      y <- as.numeric(ion_mat[, j])
      if (sum(complete.cases(x, y)) >= 3) {
        pvals[i, j] <- suppressWarnings(cor.test(x, y, method = "spearman")$p.value)
      }
    }
  }
  
  # FDR
  fdr <- matrix(p.adjust(as.vector(pvals), method = "BH"),
                nrow = nrow(pvals), dimnames = dimnames(pvals))
  
  # DataFrame final
  df <- as.data.frame(as.table(cor_mat))
  colnames(df) <- c("Gene", "Ion", "Spearman")
  df$p_value <- as.vector(pvals)
  df$FDR <- as.vector(fdr)
  df <- df %>% filter(!is.na(FDR))
  
  # Salvar
  base_name <- gsub("[^a-zA-Z0-9]", "_", slice_id)
  out_csv <- file.path(output_dir, paste0("correlation_", base_name, ".csv"))
  out_xlsx <- file.path(output_dir, paste0("correlation_", base_name, ".xlsx"))
  write.csv(df, out_csv, row.names = FALSE, quote = TRUE)
  write.xlsx(df, out_xlsx)
  
  cat("  ✅ Salvos:", out_csv, "e", out_xlsx, "\n")
}
