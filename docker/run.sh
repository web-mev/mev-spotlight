#!/bin/bash

RAW_COUNTS=$1
COORDS=$2
SAMPLE_ID=$3
NORMALIZATION_METHOD=$4
GENE_IDS=$5
MODEL=$6
OUTPUT_FILENAME=$7

declare -A HUMAN_MODELS_MAP
HUMAN_MODELS_MAP["Kidney"]=""

MODEL_FILE=${HUMAN_MODELS_MAP[${MODEL}]}
GENE_MAP_FILE=/opt/resources/human_genes.tsv

if [ ! -z "$MODEL_FILE" ]
then
    Rscript /usr/local/bin/stenrich.R \
        -f $RAW_COUNTS \
        -c $COORDS \
        -s $SAMPLE_ID \
        -n $NORMALIZATION_METHOD \
        -t $MODEL_FILE \
        -m $GENE_MAP_FILE \
        -i $GENE_IDS \
        -o $OUTPUT_FILENAME
else
    echo "The choice of \"$MODEL\" was not a valid choice."
    exit 1;
fi