#!/bin/bash

SPECIES=$1
PAPER=$2
GEO_NUMBER=$3
inputDir=cachalot:/scratch2/genomes/mmeynadier/species/$SPECIES/analysis/STARmapping/$PAPER/$GEO_NUMBER/
outputDir=/Users/mmeynadier/Documents/PhD/species/$SPECIES/analysis/STARmapping/$PAPER/$GEO_NUMBER/

rsync -rP $inputDir $outputDir
