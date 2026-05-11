devtools::install_github("campbio/decontX")
library(decontX)
load("merge_seurat_cycle_for.rdata")
# Extract raw count matrix from the Seurat object
counts <- GetAssayData(object = merge_seurat,assay = "RNA",slot = "counts")
# Run DecontX to estimate ambient RNA contamination for each cell
decontX_results <- decontX(counts)
# Add DecontX-estimated contamination scores to the Seurat metadata
merge_seurat$Contamination <- decontX_results$contamination
# Remove cells with high estimated ambient RNA contamination
# Here, cells with contamination scores >= 0.25 are excluded
contam_cutoff <- 0.25
merge_seurat_rm_ambient <- merge_seurat[,merge_seurat$Contamination < contam_cutoff]
# Save the filtered Seurat object
save(merge_seurat_rm_ambient,file = "../merge_seurat_rm_ambient.rdata")
