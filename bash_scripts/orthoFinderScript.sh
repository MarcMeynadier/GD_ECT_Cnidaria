#!/bin/bash

SPECIES=$speciesFile

cd /scratch2/genomes/mmeynadier/species/
if [ -d "crossSpecies" ]; then
  cd crossSpecies
else
  mkdir crossSpecies
  cd crossSpecies
fi

if [ -d "${SPECIES1}_${SPECIES2}" ]; then
  cd ${SPECIES1}_${SPECIES2}
else
  mkdir ${SPECIES1}_${SPECIES2}
fi


cp /scratch2/genomes/mmeynadier/species/$SPECIES1/raw/*Proteins.fasta /scratch2/genomes/mmeynadier/species/crossSpecies/${SPECIES1}_${SPECIES2}/
cp /scratch2/genomes/mmeynadier/species/$SPECIES2/raw/*Proteins.fasta /scratch2/genomes/mmeynadier/species/crossSpecies/${SPECIES1}_${SPECIES2}/

orthofinder -f ${SPECIES1}_${SPECIES2}/ -t 12

mv /scratch2/genomes/mmeynadier/species/crossSpecies/${SPECIES1}_${SPECIES2}/OrthoFinder/Results_*/* /scratch2/genomes/mmeynadier/species/crossSpecies/${SPECIES1}_${SPECIES2}/OrthoFinder/
mv /scratch2/genomes/mmeynadier/species/crossSpecies/${SPECIES1}_${SPECIES2}/*.fasta /scratch2/genomes/mmeynadier/species/crossSpecies/${SPECIES1}_${SPECIES2}/OrthoFinder/
rmdir /scratch2/genomes/mmeynadier/species/crossSpecies/${SPECIES1}_${SPECIES2}/OrthoFinder/Results_*
