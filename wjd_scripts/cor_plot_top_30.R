library(readxl)
library(dplyr)
library(ggplot2)

# --- Diretório com os arquivos .xlsx ---
dir_xlsx <- "D:/William/Data_base_VISIUM/Projects/Proj01/correlation_per_slice"
xlsx_files <- list.files(dir_xlsx, pattern = "\\.xlsx$", full.names = TRUE)

# --- Novo diretório de saída para os gráficos ---
plot_dir <- "D:/William/Data_base_VISIUM/Projects/Proj01/correlation_per_slice/plot_top30"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# --- Loop por arquivo ---
for (file in xlsx_files) {
  
  # Extrair nome base
  base_name <- gsub("^correlation_|\\.xlsx$", "", basename(file))
  
  # Ler dados
  df <- read_xlsx(file)
  
  # Selecionar top 30 por correlação absoluta
  top30 <- df %>%
    mutate(abs_spearman = abs(Spearman)) %>%
    arrange(desc(abs_spearman)) %>%
    slice_head(n = 30)
  
  # Criar gráfico
  p <- ggplot(top30, aes(x = reorder(paste(Gene, Ion, sep = "_"), Spearman), y = Spearman, fill = Spearman)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    labs(title = paste0("Top 30 Gene-Ion Spearman Correlations - ", base_name),
         x = "Gene_Ion",
         y = "Spearman Correlation") +
    theme_light(base_size = 12)
  
  # Salvar imagem
  ggsave(filename = file.path(plot_dir, paste0("top30_", base_name, ".png")),
         plot = p, width = 10, height = 6, dpi = 300)
  
  cat("??? Saved plot for:", base_name, "\n")
}
