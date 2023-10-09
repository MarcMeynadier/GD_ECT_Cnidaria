#!/bin/bash

SPECIES=$1 # Genus
PAPER=$2 # autor_date
GEO_NUMBER=$3 # GEO accession number
GSM_NUMBER=$4 # Sample number
SOLOTYPE=$5 # Droplet / CB_UMI_Simple / CB_UMI_Complex
TECH=$6 #10X_V2 / 10X_V3
soloCBwhitelist=/scratch2/genomes/mmeynadier/scripts/whitelist/${TECH}.txt
genomeDir=/scratch2/genomes/mmeynadier/species/$SPECIES/analysis/STARalign
outputDir=/scratch2/genomes/mmeynadier/species/$SPECIES/analysis/STARmapping/$PAPER/$GEO_NUMBER/$GSM_NUMBER
sampleDir=/scratch2/genomes/mmeynadier/species/$SPECIES/raw/$PAPER/$GEO_NUMBER/$GSM_NUMBER

soloBarcodeReadLength=0
soloUMIlen=0
if [ $TECH == "10X_V2" ]; then
  soloUMIlen=10
fi
if [ $TECH == "10X_V3" ]; then
  soloUMIlen=12
fi

cd $sampleDir
shopt -s globstar
sraArray=(SRR*); echo ${sraArray[@]}
cd $genomeDir/../../../ 
limitFolder=$outputDir/${sraArray[-1]}
while [ ! -d ${limitFolder} ];
do
  cd $genomeDir/../../../
  if [ -d "$SPECIES" ]; then
    cd $SPECIES
    if [ -d "analysis" ]; then
      cd analysis
      if [ -d "STARmapping" ]; then
        cd STARmapping
        if [ -d "$PAPER" ]; then
          cd $PAPER
          if [ -d "$GEO_NUMBER" ]; then
            cd $GEO_NUMBER
            if [ -d "$GSM_NUMBER" ]; then
              cd $GSM_NUMBER
              for i in ${sraArray[@]}; do
                mkdir $i
              done
            else
              mkdir $GSM_NUMBER
            fi
          else 
            mkdir $GEO_NUMBER
          fi
        else 
          mkdir $PAPER
        fi
      else
        mkdir STARmapping
      fi
    else
      mkdir analysis
    fi
  else 
     mkdir $SPECIES
  fi 
done

cd $outputDir ; rmdir *

readFilesIn1=""
readFilesIn2=""
for i in ${sraArray[@]}; do
  cd $sampleDir/$i
  mv *_R1_* *_R2_* ../
  cd ../
  R1=$( find ${i}*_R1_* )
  R2=$( find ${i}*_R2_* )
  readFilesIn1+=${sampleDir}/${R1},
  readFilesIn2+=${sampleDir}/${R2},
done
readFilesIn1=${readFilesIn1%?}
readFilesIn2=${readFilesIn2%?}
readFilesIn="${readFilesIn2} ${readFilesIn1}"

cd $outputDir
STAR --runThreadN 16 --soloUMIlen $soloUMIlen --soloBarcodeReadLength $soloBarcodeReadLength --soloCellFilter EmptyDrops_CR --readFilesCommand zcat --genomeDir $genomeDir --readFilesIn $readFilesIn --soloType $SOLOTYPE --soloCBwhitelist $soloCBwhitelist

cd $sampleDir
for i in ${sraArray[@]}; do
  mv ${i}*_R1_* ${i}*_R2_* $i
done

cd $outputDir
find . -type d -print -exec chmod 777 {} \;

echo "STAR mapping done for $GSM_NUMBER sample"
