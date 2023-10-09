#!/bin/bash

SPECIES=$1
PAPER=$2
GEO_NUMBER=$3
pathSRA=/scratch2/genomes/mmeynadier/scripts/SRA_files/
pathSpecies=/scratch2/genomes/mmeynadier/species/


outputDir=$pathSpecies/$SPECIES/raw/$PAPER/$GEO_NUMBER
cd $pathSRA/$SPECIES/$PAPER/$GEO_NUMBER
shopt -s globstar
gsmArray=(GSM*); echo ${gsmArray[@]}
for ((i=0; i<"${#gsmArray[@]}"; i++)); do
  gsmArray[i]=$(echo ${gsmArray[i]} | sed 's/.txt//g')
done
cd $pathSpecies 
limitFolder=$outputDir/${gsmArray[-1]}
while [ ! -d ${limitFolder} ];
do
  cd $pathSpecies
  if [ -d "$SPECIES" ]; then
    cd $SPECIES
    if [ -d "raw" ]; then
      cd raw
      if [ -d "$PAPER" ]; then
        cd $PAPER
        if [ -d "$GEO_NUMBER" ]; then
          cd $GEO_NUMBER
          for i in ${gsmArray[@]}; do
            mkdir $i
          done
        else 
          mkdir $GEO_NUMBER
        fi
      else 
        mkdir $PAPER
      fi
    else 
       mkdir raw
    fi
  else 
     mkdir $SPECIES
  fi 
done


cd $outputDir
for i in ${gsmArray[@]}; do
  cd $i
  prefetch --max-size 50G --option-file $pathSRA/$SPECIES/$PAPER/$GEO_NUMBER/${i}.txt
  cd $outputDir
done  




