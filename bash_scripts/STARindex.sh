#!/bin/bash

SPECIES=$1 # Genus
genomeFastaFiles=/scratch2/genomes/mmeynadier/species/$SPECIES/raw/*Genome.fasta
sjdbGTFfile=/scratch2/genomes/mmeynadier/species/$SPECIES/raw/*.gtf
genomeDir=/scratch2/genomes/mmeynadier/species/$SPECIES/analysis/STARalign

STAR --runThreadN 16 --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $genomeFastaFiles --sjdbGTFfile $sjdbGTFfile --genomeSAindexNbases 13 
