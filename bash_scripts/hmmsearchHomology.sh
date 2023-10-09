#!/bin/bash

SPECIES=$1
AUTOR=$2
GEO=$3
GSM=$4

PFAM_PATH=/scratch1/db/Pfam-A.hmm
FASTA_INPUT=/scratch2/genomes/mmeynadier/species/$SPECIES/analysis/trimFastaBiomarkers/${SPECIES}_${AUTOR}_${GEO}_${GSM}.fasta

cd /scratch2/genomes/mmeynadier/species/$SPECIES/analysis/HMMER/
hmmsearch --tblout pfam_results_${SPECIES}_${AUTOR}_${GEO}_${GSM}.table --cut_ga $PFAM_PATH $FASTA_INPUT
