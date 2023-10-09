
# Preliminary samples - 2016 dataset
# Packages and dependence
packageCheckClassic <- function(x){
  # 
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE )
      require( i , character.only = TRUE )
    }
  }
}
packageCheckClassic(c('DESeq2','devtools','BiocManager','ggplot2','ggrepel','markdown','pheatmap','RColorBrewer','genefilter','gplots','vegan','dplyr'))
#BiocManager::install('tximport', force = TRUE)
#BiocManager::install('apeglm')
#BiocManager::install('ashr')
#BiocManager::install("EnhancedVolcano")
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library('ggvenn')
library('tximport')
library('apeglm')
library('ashr')
library('EnhancedVolcano')
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")


# Working environment 
scriptPath<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(scriptPath)

species = 'Astroides'
autor = 'Teixido'
dataset = list('adult','juvenile')
subdataset = list('nov2016','june2017','sept2017','may2018','sept2018','juv')
dataPath = paste('/Users/mmeynadier/Documents/PhD/',species,sep='')
dataPath = paste(dataPath,'/analysis/STARmapping/',sep='')
dataPath = paste(dataPath,autor,sep='')

filesList = list()
for (i in dataset) {
  dataPath = paste(dataPath,i,sep='/')
  for (j in subdataset) {
    dataPath = paste(dataPath,j,sep='/')
    
  }
}

#Make an object listing the path and name of each of the gene counts files (output of STAR). This is going to save us time because instead of loading all tables individually, they will be extracted all directly from the path.
#ff <- list.files( path = "/Users/mmeynadier/Documents/PhD/species/Astroides/analysis/STARmapping/teixido/", pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )


#Import all the count information from all filles in table format, skipping the first 4 lines (information not needed). The object is a list of tables.
counts.files <- lapply(ff, read.table, skip = 4 )


#Make a table only including the forth column (counts for 2nd read)
raw_counts <- as.data.frame(sapply(counts.files, function(x) x[ , 4]))


#To be able to add the column names --> Simplify the name of each file by removing the first part (path) and the end _ReadsPerGene.out.tab
ff <- gsub( "Users/mmeynadier/Documents/PhD/species/Astroides/analysis/STARmapping/teixido/adult/nov2016", "", ff )
ff <- gsub( "_ReadsPerGene.out.tab", "", ff )
ff <- gsub( "/", "", ff )
ff <- gsub( "1_", "", ff )
ff <- gsub( "2_", "", ff )
ff <- gsub( "3_", "", ff )


#Add coulmn names
colnames(raw_counts) <- ff


#Add row names = first column of each count file (they are ordered the same way)
row.names(raw_counts) <- counts.files[[1]]$V1

