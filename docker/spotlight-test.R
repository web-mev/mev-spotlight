suppressMessages(library(SPOTlight))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(SpatialExperiment))
suppressMessages(library(scater))
suppressMessages(library(scran))

source('/usr/local/bin/prep_stlist.R')

suppressMessages(library(optparse))

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
        c('-m', '--model'),
        help='The set of genes to test under SThet.'
    ),
    make_option(
        c('-o', '--output'),
        help='The name of the output file'
    )
)

# prepare an STList instance. Note that we are potentially re-mapping
# the original gene identifiers to symbols such that they work with
# the MSigDB files:
spat_list <- prep_stlist(opt$input_file, 
                         opt$coordinates_file,
                         opt$sample_name,
                         gene_mapping_df,
                         gene_ids,
                         'SYMBOL')
spat <- spat_list$spat

# normalize
spat <- transform_data(spat, method=opt$normalization)

# Load the trained model
load(opt$model)

# Run the NMF model against input data
res <- runDeconvolution(
    x = spat,
    mode = mod
)

# Write output
write.table(
    t(res$mat),
    opt$output,
    sep="\t",
    quote=F,
    row.names=T
)