#!/bin/bash

SPECIES=$1
pathProteom=/backup2/genomes/mmeynadier/species/${SPECIES}/${SPECIES}Proteins.fasta
pathDb=/scratch2/genomes/mmeynadier/scripts/tools/refseq89o/mitoProtDb/mitoProtDb.fasta
pathOutput=/backup2/genomes/mmeynadier/species/${SPECIES}/

cd $pathOutput
blastp -query $pathProteom -db $pathDb -outfmt "6 qseqid sacc qlen slen length nident pident evalue stitle" -evalue 1e-10 1>${SPECIES}MitoGenes.txt -num_threads 16  
