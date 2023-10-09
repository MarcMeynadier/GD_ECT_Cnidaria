#!/bin/bash

SPECIES=$1
PAPER=$2
GEO_NUMBER=$3
pathSRA=/scratch2/genomes/mmeynadier/scripts/SRA_files/
pathSpecies=/scratch2/genomes/mmeynadier/species
outputDir=$pathSpecies/$SPECIES/raw/$PAPER/$GEO_NUMBER

cd $outputDir
shopt -s globstar
gsmArray=(GSM*); echo ${gsmArray[@]}
for i in ${gsmArray[@]}; do
  cd $i
  shopt -s globstar
  sraArray=(SRR*); echo ${sraArray[@]}
    for j in ${sraArray[@]}; do
      cd $j
      parallel-fastq-dump --sra-id $j --split-files --threads 16 --outdir $outputDir/$i/$j --tmpdir $pathSpecies --gzip
      echo "$j finished"
      cd $outputDir/$i
    done
  echo "$i finished"
  cd $outputDir
done
