#!/bin/bash

SPECIES=$1
OUTPUT=/scratch2/genomes/mmeynadier/species/$SPECIES/analysis/interproscan
cd /backup2/genomes/mmeynadier/species/$SPECIES/raw

sed 's/\*//g' longest_isoform*.fasta > longest_isoform_${SPECIES}ProteinsCured.fasta

mkdir -p $OUTPUT

/scratch1/db/interproscan-5.53-87.0/interproscan.sh -i /backup2/genomes/mmeynadier/species/$SPECIES/raw/*ProteinsCured.fasta -f tsv -o ${OUTPUT}/${SPECIES}_interproscan.tsv --iprlookup --goterms --pathways --cpu 16 -dp
