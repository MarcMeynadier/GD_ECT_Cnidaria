#!/bin/bash

SPECIES=$1 # Genus
PAPER=$2 # autor_date
GEO_NUMBER=$3 # GEO accession number
GSM_NUMBER=$4 # Sample number

PATH_SAMPLE=/scratch2/genomes/mmeynadier/species/$SPECIES/raw/$PAPER/$GEO_NUMBER/$GSM_NUMBER
cd $PATH_SAMPLE
shopt -s nullglob
array=(*/)
shopt -u nullglob
SAMPLE=$(echo "${array[@]}")
SAMPLE_LIST=$(sed "s/\/ /,/g" <<< "$SAMPLE")
SAMPLE_LIST=$(sed "s/\///g" <<< "$SAMPLE_LIST")

cd /scratch2/genomes/mmeynadier/species/$SPECIES/analysis/CellRanger/$PAPER/$GEO_NUMBER

cellranger count --id=$GSM_NUMBER \
   --fastqs=$PATH_SAMPLE \
   --sample=$SAMPLE_LIST \
   --transcriptome=/scratch2/genomes/mmeynadier/species/Nematostella/raw/NvMkref
