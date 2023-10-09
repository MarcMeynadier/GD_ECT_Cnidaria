#!/bin/bash

SPECIES=$1 # Genus
genomeFastaFiles=/backup2/genomes/mmeynadier/species/$SPECIES/*Genome.fasta
sjdbGTFfile=/backup2/genomes/mmeynadier/species/$SPECIES/*.gtf
genomeDir=/scratch2/genomes/mmeynadier/species/$SPECIES/analysis/STARalign

STAR --runThreadN 16 --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $genomeFastaFiles --sjdbGTFfile $sjdbGTFfile --genomeSAindexNbases 13 
