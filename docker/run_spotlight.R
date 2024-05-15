suppressMessages(library(SPOTlight))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(SpatialExperiment))
suppressMessages(library(NMF))
suppressMessages(library(Matrix))
suppressMessages(library(sparseMatrixStats))
source('/usr/local/bin/prep_stlist.R')

suppressMessages(library(optparse))

# our own scaling function. Note that in the spotlight code
# https://github.com/MarcElosua/SPOTlight/blob/devel/R/utils.R
# their scaling introduces NAs which cause downstream issues. 
# The NAs originate from a 0/0 which results from the introduction
# of zero-count genes which are needed to get the model and count
# matrix to load appropriately.
.scale_uv <- function(x) {
    sds <- sparseMatrixStats::rowSds(x, na.rm = TRUE)
    # Scale by gene (each row by its sd) for unit variance
    x / sds
}

# args from command line:
args <- commandArgs(TRUE)

option_list <- list(
    make_option(
        c('-f', '--input_file'),
        help='Path to the count matrix input.'
    ),
    make_option(
        c('-c', '--coordinates_file'),
        help='Path to the barcode spatial coordinates input.'
    ),
    make_option(
        c('-s', '--sample_name'),
        help='Sample name'
    ),
    make_option(
        c('-n', '--normalization'),
        help='Normalization method of `log` or `sct`'
    ),
    make_option(
        c('-t', '--model'),
        help='The set of genes to test under SThet.'
    ),
    make_option(
        c('-m', '--map_file'),
        help='TSV file for mapping gene identifiers'
    ),
    make_option(
        c('-i', '--gene_ids'),
        help='The gene identifier system used.'
    ),
    make_option(
        c('-o', '--output'),
        help='The name of the output file'
    )
)

opt <- parse_args(OptionParser(option_list=option_list))

# Sanity checks on the inputs:
# Check that the file was provided:
if (is.null(opt$input_file)){
    message('Need to provide a count matrix with the -f/--input_file arg.')
    quit(status=1)
}

if (is.null(opt$coordinates_file)){
    message('Need to provide a count matrix with the -c/--coordinates_file arg.')
    quit(status=1)
}

if (is.null(opt$map_file)){
    message('Need to provide a gene ID mapping file with the -m/--map_file arg.')
    quit(status=1)
}

if (is.null(opt$model)){
    message('Need to provide a model file with the -t/--model arg.')
    quit(status=1)
}

if (is.null(opt$gene_ids)){
    message('Need to provide the gene identifier with the -i/--gene_ids arg.')
    quit(status=1)
} else {
    gene_ids <- toupper(opt$gene_ids)
}

# transform the name of the normalization scheme:
if (is.null(opt$normalization)){
    message('Need to provide a normalization scheme with the -n/--normalization arg.')
    quit(status=1)
} else if(tolower(opt$normalization) == 'sctransform'){
    norm_scheme <- 'sct'
} else if(tolower(opt$normalization) == 'log'){
    norm_scheme <- 'log'
} else {
    message('We only accept `log` or `SCTransform` for the normalization scheme.')
    quit(status=1)
}

# Load the trained model
load(opt$model)

tryCatch({
    load(opt$model)
}, error = function(x){
        message('Failed when loading the model file.')
        quit(status=1)
    }
)

# change the working directory to co-locate with the counts file:
working_dir <- dirname(opt$input_file)
setwd(working_dir)

# dataframe which maps from one system to another:
gene_mapping_df <- read.table(opt$map_file, sep='\t', header=T)

# prepare an STList instance. Note that we are potentially re-mapping
# the original gene identifiers to symbols such that they work with
# the model file:
spat_list <- prep_stlist(opt$input_file, 
                         opt$coordinates_file,
                         opt$sample_name,
                         gene_mapping_df,
                         gene_ids,
                         'SYMBOL')
spat <- spat_list$spat

# normalize
spat <- transform_data(spat, method=norm_scheme)

save(spat, file='spat.RData')

W <- NMF::basis(mod_ls$mod)
model_genes <- rownames(W)

counts <- spat@counts[[opt$sample_name]]
mtx_genes <- rownames(counts)
diff_set <- setdiff(model_genes, mtx_genes)

# create a matrix of zeros, make it sparse, and add it to 
# the existing counts:
m <- matrix(0, length(diff_set), dim(counts)[2])
rownames(m) <- diff_set
colnames(m) <- colnames(counts)
ms <- Matrix(m, sparse=TRUE)
final <- rbind(counts, ms)

# scale to avoid the NAs if we call runDeconvolution
# with scale=TRUE (the default) below:
final <- .scale_uv(final)
final[is.na(final)] <- 0

# Run the NMF model against input data
res <- runDeconvolution(
    x = final,
    mod = mod_ls$mod,
    ref = mod_ls$topic,
    scale=FALSE
)

# Write output
write.table(
    t(res$mat),
    opt$output,
    sep="\t",
    quote=F,
    row.names=T
)