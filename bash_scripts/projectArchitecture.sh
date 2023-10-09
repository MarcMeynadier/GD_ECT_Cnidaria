#!/bin/bash

mkdir analysis
cd analysis

mkdir species scripts
cd scripts

mkdir SRA_files whitelist
git clone https://github.com/MarcMeynadier/GD_ECT_Cnidaria.git

cd ../species

species=("Nematostella" "Clytia" "Pelagia" "Hydra" "Xenia" "Stylophora")

for i in "${species[@]}"
do
  mkdir $i
  mkdir {$i}/{genome}
done 