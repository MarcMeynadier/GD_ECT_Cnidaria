#!/bin/bash

SPECIES1=$1
SPECIES2=$2

pathDb=/scratch1/db/uniprot_sprot

pathEOO=/scratch2/genomes/mmeynadier/species/crossSpecies/${SPECIES1}_${SPECIES2}/EOO/

cd $pathEOO
for i in $(ls -d */); do
  pwd
  cd ../psiblast/
  mkdir $i
  outputPath=/scratch2/genomes/mmeynadier/species/crossSpecies/${SPECIES1}_${SPECIES2}/psiblast/${i}
  cd ../EOO/$i
  for j in $(ls); do
    psiblast -db $pathDb -query $j -out ${j}_psiblast.tsv -outfmt "6 qseqid sseqid evalue bitscore sgi sacc staxids ssciname scomnames stitle" -max_target_seqs 1 -num_threads 12  
    mv ${j}_psiblast.tsv $outputPath
    done
  cd ../
done 
