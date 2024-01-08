#!/bin/bash

SPECIES=$1 # Genus
PAPER=$2 # autor_date
GEO_NUMBER=$3 # GEO accession number
GSM_NUMBER=$4 # Sample number
SOLOTYPE=$5 # Droplet / CB_UMI_Simple / CB_UMI_Complex
TECH=$6 #10X_V2 / 10X_V3 / BD_rhapsody_V1 / DropSeq
soloCBwhitelistPath=/scratch2/genomes/mmeynadier/scripts/whitelist/
genomeDir=/scratch2/genomes/mmeynadier/species/$SPECIES/analysis/STARalign
outputDir=/scratch2/genomes/mmeynadier/species/$SPECIES/analysis/STARmapping/$PAPER/$GEO_NUMBER/$GSM_NUMBER
sampleDir=/backup2/genomes/mmeynadier/species/$SPECIES/$PAPER/$GEO_NUMBER/$GSM_NUMBER

#soloUMIlen=0
#if [ $TECH=="10X_V2" ]; then
  #soloCBwhitelist=${soloCBwhitelistPath}${TECH}.txt
  #soloCBstart=1 
  #soloCBlen=16 
  #soloUMIstart=17 
  #soloUMIlen=10
#fi
#if [ $TECH=="10X_V3" ]; then
#  soloCBwhitelist=${soloCBwhitelistPath}${TECH}.txt
#  soloCBstart=1 
#  soloCBlen=16 
#  soloUMIstart=17 
#  soloUMIlen=12
#fi
#if [ $TECH=="BD_rhapsody_V1" ]; then
  #soloCBwhitelist="${soloCBwhitelistPath}${TECH}_1.txt ${soloCBwhitelistPath}${TECH}_2.txt ${soloCBwhitelistPath}${TECH}_3.txt"
  #soloType='CB_UMI_Complex'
  #soloUMIlen=8
  #soloCBposition=0_0_0_8 0_21_0_29 0_43_0_51
  #soloUMIposition=0_52_0_59
#fi
if [ $TECH=="DropSeq" ]; then
  soloCBwhitelist="None"
  soloCBstart=1
  soloCBlen=12
  soloUMIstart=13
  soloUMIlen=8
fi



cd $sampleDir
shopt -s globstar
sraArray=(SRR*); echo ${sraArray[@]}

mkdir -p $outputDir

cd $outputDir

readFilesIn1=""
readFilesIn2=""
for i in ${sraArray[@]}; do
  cd $sampleDir/$i
  mv *_1.* *_2.* ../
  cd ../
  R1=$( find ${i}*_1.* )
  R2=$( find ${i}*_2.* )
  readFilesIn1+=${sampleDir}/${R1},
  readFilesIn2+=${sampleDir}/${R2},
done
readFilesIn1=${readFilesIn1%?}
readFilesIn2=${readFilesIn2%?}
readFilesIn="${readFilesIn2} ${readFilesIn1}"
echo $readFilesIn

cd $outputDir
#STAR --runThreadN 16 --soloType $SOLOTYPE --soloCBstart $soloCBstart --soloCBlen $soloCBlen --soloUMIstart $soloUMIstart --soloUMIlen $soloUMIlen --soloCellFilter EmptyDrops_CR --readFilesCommand zcat --genomeDir $genomeDir --readFilesIn $readFilesIn --soloCBwhitelist None
STAR --runThreadN 16 --soloCBstart $soloCBstart --soloCBlen $soloCBlen --soloUMIstart $soloUMIstart --soloUMIlen $soloUMIlen --soloCellFilter EmptyDrops_CR --readFilesCommand zcat --genomeDir $genomeDir --readFilesIn $readFilesIn --soloType $SOLOTYPE --soloCBwhitelist $soloCBwhitelist
#STAR --runThreadN 16 --soloCBmatchWLtype Exact --soloUMIlen $soloUMIlen --soloCBposition 0_0_0_8 0_21_0_29 0_43_0_51 --soloUMIposition 0_52_0_59 --soloCellFilter EmptyDrops_CR --genomeDir $genomeDir --readFilesIn $readFilesIn --soloType $SOLOTYPE --soloCBwhitelist $soloCBwhitelist

cd $sampleDir
for i in ${sraArray[@]}; do
  mv ${i}*_1.* ${i}*_2.* $i
done

cd $outputDir
find . -type d -print -exec chmod 777 {} \;

echo "STAR mapping done for $GSM_NUMBER sample"
