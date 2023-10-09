#!/bin/bash

<<Block_comment 
kallistoQuant : Run the Kallisto tool on the data specified as arguments when the script is launched. 
These arguments are the type of organism (juvenile or adult), the type of reference transcriptome (A, LA, LJA), 
the type of sequences (paired or single), the name of dataset directory. Also, if the sequence type is single end, 
the average length of the RNA fragments (length) as well as the standard deviation of this metric (SD) must be specified. 
They can be easily obtained using the Python script averageLength_SD_bioanalyser.

Marc Meynadier
Block_comment

SPECIES=$1 
AUTOR=$2
INDEX=/scratch2/genomes/mmeynadier/species/$SPECIES/analysis/kallisto/*.idx
INPUT=/backup2/genomes/mmeynadier/species/$SPECIES/raw/$AUTOR
OUTPUT=/scratch2/genomes/mmeynadier/species/$SPECIES/analysis/kallisto/

for i in ${INPUT}/*R1.fq.gz
  do
    SPL_NAME=${i//_R1/}
    SAMPLE_NAME=$(echo $SPL_NAME | cut -d . -f 1 | sed 's|.*/||')
    j=${i/_R1/_R2}
    echo $i $j $SAMPLE_NAME
    kallisto quant -i $INDEX -o $OUTPUT -b 100 -t 16 $i $j
    mv ${OUTPUT}/abundance.h5 ${OUTPUT}/abundance_${SAMPLE_NAME}.h5
    mv ${OUTPUT}/abundance.tsv ${OUTPUT}/abundance_${SAMPLE_NAME}.tsv
    mv ${OUTPUT}/run_info.json ${OUTPUT}/run_info_${SAMPLE_NAME}.json
done

