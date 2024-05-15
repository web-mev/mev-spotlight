#!/bin/bash

RAW_COUNTS=$1
COORDS=$2
SAMPLE_ID=$3
NORMALIZATION_METHOD=$4
GENE_IDS=$5
MODEL=$6
OUTPUT_FILENAME=$7

# TODO: handle multiple organisms and gene mapping more generally
declare -A MODELS_MAP
MODELS_MAP["Mouse:Kidney"]="https://webmev-public.s3.us-east-2.amazonaws.com/spotlight/kidney.txt"

MODEL_FILE=${MODELS_MAP[${MODEL}]}
GENE_MAP_FILE=/opt/resources/mouse_genes.tsv

if [ ! -z "$MODEL_FILE" ]
then
    curl $MODEL_FILE -o /tmp/model.txt
    Rscript /usr/local/bin/run_spotlight.R \
        -f $RAW_COUNTS \
        -c $COORDS \
        -s $SAMPLE_ID \
        -n $NORMALIZATION_METHOD \
        -t /tmp/model.txt \
        -m $GENE_MAP_FILE \
        -i $GENE_IDS \
        -o $OUTPUT_FILENAME
else
    echo "The choice of \"$MODEL\" was not a valid choice."
    exit 1;
fi