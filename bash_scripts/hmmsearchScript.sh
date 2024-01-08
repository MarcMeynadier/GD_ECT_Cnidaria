#!/bin/bash

SPECIES=$1

PFAM_PATH=/scratch1/db/Pfam-A.hmm
FASTA_INPUT=/backup2/genomes/mmeynadier/species/$SPECIES/${SPECIES}Proteins.fasta

outputDir=/scratch2/genomes/mmeynadier/species/$SPECIES/analysis/hmmsearch/
mkdir -p $outputDir
cd $outputDir
hmmsearch --tblout pfam_results_${SPECIES}.table --cpu 16 --cut_ga $PFAM_PATH $FASTA_INPUT
