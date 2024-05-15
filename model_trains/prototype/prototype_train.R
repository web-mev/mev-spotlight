library(ggplot2)
library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)

# Download mouse spatial data
library(TENxVisiumData)
spe <- MouseKidneyCoronal()
# Use symbols instead of Ensembl IDs as feature names
rownames(spe) <- rowData(spe)$symbol

# Download mouse single cell atlas data
library(TabulaMurisSenisData)
sce <- TabulaMurisSenisDroplet(tissues = "Kidney")$Kidney

# Keep cells from 18m mice
sce <- sce[, sce$age == "18m"]
# Keep cells with clear cell type annotations
sce <- sce[, !sce$free_annotation %in% c("nan", "CD45")]

# Normalize the data
sce <- logNormCounts(sce)

# Get vector indicating which genes are neither ribosomal or mitochondrial
genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(sce))

# Model the variance of the log-expression profiles for each gene, decomposing 
# it into technical and biological components based on a fitted mean-variance trend.
dec <- modelGeneVar(sce, subset.row = genes)

# Get the top 3000 genes.
hvg <- getTopHVGs(dec, n = 3000)

# Overwrite the single cell labels with cell types
colLabels(sce) <- colData(sce)$free_annotation

# Compute marker genes
mgs <- scoreMarkers(sce, subset.row = genes)

mgs_fil <- lapply(names(mgs), function(i) {
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.8
    x <- x[x$mean.AUC > 0.8, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)

# split cell indices by identity
idx <- split(seq(ncol(sce)), sce$free_annotation)
# downsample to at most 20 per identity & subset
# We are using 5 here to speed up the process but set to 75-100 for your real
# life analysis
n_cells <- 5
cs_keep <- lapply(idx, function(i) {
    n <- length(i)
    if (n < n_cells)
        n_cells <- n
    sample(i, n_cells)
})
sce <- sce[, unlist(cs_keep)]

# Two step train and run process
# Train
mod_ls <- trainNMF(
    x = sce, 
    y = spe, 
    groups = as.character(sce$free_annotation), 
    mgs = mgs_df, 
    weight_id = "mean.AUC", 
    group_id = "cluster", 
    gene_id = "gene"
)

# saving the model to text and re-loading
save(mod_ls, ascii=T, file="model_store.txt")
rm(mod_ls)
load("model_store.txt")
