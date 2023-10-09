#!/bin/bash

inputDir=/Users/mmeynadier/Documents/PhD/scripts/SRA_files/
outputDir=cachalot:/scratch2/genomes/mmeynadier/scripts/SRA_files/

rsync -rP $inputDir $outputDir
