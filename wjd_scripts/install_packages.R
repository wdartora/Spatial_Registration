# R script to install necessary packages for correlation analysis

# Update package list and install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

# Install common dependencies and Seurat using BiocManager and CRAN
# Using update=FALSE, ask=FALSE for non-interactive installation
BiocManager::install("multtest", update=FALSE, ask=FALSE, force = TRUE)

# Install other packages from CRAN
packages_to_install <- c("Seurat", "dplyr", "tidyr", "readr", "pheatmap", "remotes")

for (pkg in packages_to_install) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
}

# Special handling for SeuratObject if Seurat v5 or later (often installed with Seurat)
# but good to ensure if specific functions are used directly.
# This check might be redundant if Seurat installation handles it well.
tryCatch({
    if (packageVersion("Seurat") >= "5.0.0") {
        if (!requireNamespace("SeuratObject", quietly = TRUE)) {
            install.packages("SeuratObject", repos = "https://cloud.r-project.org/")
        }
    }
}, error = function(e) {
    print(paste("Could not check Seurat version or install SeuratObject:", e$message))
    print("Proceeding, assuming Seurat installation was complete.")
})

print("Package installation script finished.")

