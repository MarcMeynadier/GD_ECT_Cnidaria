# sc-RNA - Seurat analysis

# Packages and dependence
packageCheckClassic <- function(x){
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      #  Load package after installing
      require( i , character.only = TRUE )
    }
  }
}

packageCheckClassic(c('Seurat','dplyr','patchwork','R.utils'))

# Setup environment variables
species <- 'Nematostella'
paper <- 'steger2022'
GEO <- 'GSE200198'
GSM <- 'GSM6021902_planula2d'

nFeature_RNA_Var = 0
nCount_RNA_Var = 0
nfeatures = 0
pcsCompute = 0
nDimsElbow = 0
dimensions = 0
resolution = 0
randomSeed = 0
nNeighbours = 0
spread = 0
seedUse = 0

if (GSM == 'GSM6021899_gastrula18h') {
  nFeature_RNA_Var = 300 ; nCount_RNA_Var = 100000 ; nfeatures = 2000 ; pcsCompute = 50 ; nDimsElbow = 30
  dimensions = c(1:10) ; resolution = 0.2 ; randomSeed = 0 ; nNeighbours = 30 ; spread = 1 ; seedUse = 1
}

if (GSM == 'GSM6021900_gastrula3') {
  nFeature_RNA_Var = 250 ; nCount_RNA_Var = 10000 ; nfeatures = 2000 ; pcsCompute = 50 ; nDimsElbow = 50
  dimensions = c(1:10) ; resolution = 0.35 ; randomSeed = 0 ; nNeighbours = 30 ; spread = 1 ; seedUse = 5
}

if (GSM == 'GSM6021901_Gastrula2') {
  nFeature_RNA_Var = 250 ; nCount_RNA_Var = 30000 ; nfeatures = 2000 ; pcsCompute = 50 ; nDimsElbow = 20
  dimensions = c(1:8) ; resolution = 0.04 ; randomSeed = 0 ; nNeighbours = 50 ; spread = 0.2 ; seedUse = 0
}

if (GSM == 'GSM6021902_planula2d') {
  nFeature_RNA_Var = 300 ; nCount_RNA_Var = 100000 ; nfeatures = 2000 ; pcsCompute = 50 ; nDimsElbow = 30
  dimensions = c(1:20) ; resolution = 0.5 ; randomSeed = 0 ; nNeighbours = 30 ; spread = 0.75 ; seedUse = 1
}

if (GSM == 'GSM6021903_planula3d') {
  nFeature_RNA_Var = 250 ; nCount_RNA_Var = 25000 ; nfeatures = 2000 ; pcsCompute = 50 ; nDimsElbow = 50
  dimensions = c(1:10) ; resolution = 0.7 ; randomSeed = 0 ; nNeighbours = 25 ; spread = 0.5 ; seedUse = 1
}

if (GSM == 'GSM6021904_planula4d') {
  nFeature_RNA_Var = 300 ; nCount_RNA_Var = 10000 ; nfeatures = 2000 ; pcsCompute = 50 ; nDimsElbow = 50
  dimensions = c(1:15) ; resolution = 0.5 ; randomSeed = 0 ; nNeighbours = 30 ; spread = 0.5 ; seedUse = 1
}

if (GSM == 'GSM6021904_planula4d') {
  nFeature_RNA_Var = 300 ; nCount_RNA_Var = 10000 ; nfeatures = 2000 ; pcsCompute = 50 ; nDimsElbow = 50
  dimensions = c(1:15) ; resolution = 0.5 ; randomSeed = 0 ; nNeighbours = 30 ; spread = 0.5 ; seedUse = 1
}

if (GSM == 'GSM6021905_planula4d2c') {
  nFeature_RNA_Var = 300 ; nCount_RNA_Var = 10000 ; nfeatures = 2000 ; pcsCompute = 50 ; nDimsElbow = 50
  dimensions = c(1:15) ; resolution = 0.5 ; randomSeed = 0 ; nNeighbours = 30 ; spread = 0.5 ; seedUse = 1
}

if (GSM == 'GSM6021906_planula5d') {
  nFeature_RNA_Var = 250 ; nCount_RNA_Var = 20000 ; nfeatures = 2000 ; pcsCompute = 50 ; nDimsElbow = 50
  dimensions = c(1:20) ; resolution = 0.8 ; randomSeed = 0 ; nNeighbours = 30 ; spread = 0.5 ; seedUse = 1
}

if (GSM == 'GSM6021907_polyp8d') {
  nFeature_RNA_Var = 250 ; nCount_RNA_Var = 15000 ; nfeatures = 2000 ; pcsCompute = 50 ; nDimsElbow = 50
  dimensions = c(1:20) ; resolution = 0.5 ; randomSeed = 0 ; nNeighbours = 30 ; spread = 0.7 ; seedUse = 1
}

if (GSM == 'GSM6021908_polyp16d') {
  nFeature_RNA_Var = 250 ; nCount_RNA_Var = 20000 ; nfeatures = 2000 ; pcsCompute = 50 ; nDimsElbow = 50
  dimensions = c(1:15) ; resolution = 0.8 ; randomSeed = 0 ; nNeighbours = 30 ; spread = 0.5 ; seedUse = 1
}

if (GSM == 'GSM6021909_phbw') {
  nFeature_RNA_Var = 250 ; nCount_RNA_Var = 20000 ; nfeatures = 2000 ; pcsCompute = 50 ; nDimsElbow = 50
  dimensions = c(1:20) ; resolution = 0.8 ; randomSeed = 0 ; nNeighbours = 30 ; spread = 1 ; seedUse = 1
}

if (GSM == 'GSM6021910_mesF') {
  nFeature_RNA_Var = 200 ; nCount_RNA_Var = 20000 ; nfeatures = 2000 ; pcsCompute = 50 ; nDimsElbow = 50
  dimensions = c(1:23) ; resolution = 0.2 ; randomSeed = 0 ; nNeighbours = 10L ; spread = 1 ; seedUse = 1
}

if (GSM == 'GSM4663943_NvSubAdultTentacle') {
  nFeature_RNA_Var = 250 ; nCount_RNA_Var = 10000 ; nfeatures = 2000 ; pcsCompute = 50 ; nDimsElbow = 50
  dimensions = c(1:20) ; resolution = 1.1 ; randomSeed = 0 ; nNeighbours = 20 ; spread = 0.4 ; seedUse = 0
}

if (GSM == 'GGSM4663944_NvSubAdultMesentery') {
  nFeature_RNA_Var = 250 ; nCount_RNA_Var = 15000 ; nfeatures = 2000 ; pcsCompute = 50 ; nDimsElbow = 50
  dimensions = c(1:20) ; resolution = 0.9 ; randomSeed = 0 ; nNeighbours = 15 ; spread = 0.5 ; seedUse = 42
}

if (GSM == 'GSM4663945_NvSubAdultPharynx') {
  nFeature_RNA_Var = 250 ; nCount_RNA_Var = 10000 ; nfeatures = 2000 ; pcsCompute = 50 ; nDimsElbow = 50
  dimensions = c(1:20) ; resolution = 1 ; randomSeed = 0 ; nNeighbours = 25 ; spread = 0.5 ; seedUse = 1
}

if (GSM == 'GSM4663946_NvSubAdultBodywall') {
  nFeature_RNA_Var = 250 ; nCount_RNA_Var = 5000 ; nfeatures = 2000 ; pcsCompute = 50 ; nDimsElbow = 50
  dimensions = c(1:20) ; resolution = 0.6 ; randomSeed = 0 ; nNeighbours = 30 ; spread = 0.5 ; seedUse = 42
}

# Checking architecture and working environment 
scriptPath<-dirname(rstudioapi::getSourceEditorContext()$path)
scriptPath <- sub("/[^/]+$", "", scriptPath)
scriptPath <- sub("/[^/]+$", "", scriptPath)
scriptPath <- sub("/[^/]+$", "", scriptPath)
dataPath<-paste(scriptPath,'species',species,'analysis','STARmapping',paper,GEO,GSM,'Solo.out','Gene','filtered', sep = '/')
outputPath<-paste(scriptPath,'species',species,'analysis','Seurat',paper,GEO,GSM, sep='/')

while (file.exists(outputPath) == FALSE) {
  if (file.exists(paste(scriptPath,'species',sep='/'))){
    wdDirectory <- paste(scriptPath,'species', sep ='/')
    setwd(wdDirectory)
    if (file.exists(paste(wdDirectory,species,sep='/'))){
      wdDirectory <- paste(wdDirectory,species, sep ='/')
      setwd(wdDirectory)
      if (file.exists(paste(wdDirectory,'analysis',sep='/'))){
        wdDirectory <- paste(wdDirectory,'analysis', sep ='/')
        setwd(wdDirectory)
        if (file.exists(paste(wdDirectory,'Seurat',sep='/'))){
          wdDirectory <- paste(wdDirectory,'Seurat', sep ='/')
          setwd(wdDirectory)
          if (file.exists(paste(wdDirectory,paper,sep='/'))){
            wdDirectory <- paste(wdDirectory,paper, sep ='/')
            setwd(wdDirectory)
            if (file.exists(paste(wdDirectory,GEO,sep='/'))){
              wdDirectory <- paste(wdDirectory,GEO, sep ='/')
              setwd(wdDirectory)
              if (file.exists(paste(wdDirectory,GSM,sep='/'))){
                wdDirectory <- paste(wdDirectory,GSM, sep ='/')
                setwd(wdDirectory)
                print('Output path is setup')
              } else {
                dir.create(file.path(wdDirectory, GSM))
              }
            } else {
              dir.create(file.path(wdDirectory, GEO))
            }
          } else {
            dir.create(file.path(wdDirectory, paper))
          }
        } else {
          dir.create(file.path(wdDirectory, 'Seurat'))
        }
      } else {
        dir.create(file.path(wdDirectory, 'analysis'))
      }
    } else {
      dir.create(file.path(wdDirectory, species))
    }
  } else {
    dir.create(file.path(scriptPath, 'species'))
  }
}

setwd(dataPath)

# Data importation 
files <- list.files(path = ".")
for (i in files) {
  gzip(i,remove=FALSE)
}
projectName <- paste(species,paper,GEO,GSM,'seurat',sep='_')
project.data <- Read10X(data.dir = dataPath)
project <- CreateSeuratObject(counts=project.data,project=species,min.cells=3,min.features=200)
filesGz <- list.files(path = ".",pattern = ".gz")
for (i in filesGz) {
  file.remove(i)
}

# Pre-processing
project[["percent.mt"]] <- PercentageFeatureSet(project, pattern = "^MT-") # Mitochondrial gene detection
VlnSave=paste(outputPath,paste(projectName,'QC.png',sep='_'),sep='/')
png(file=VlnSave, width = 465, height = 225, units='mm', res = 300)
plot(VlnPlot(project, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)) # Visualize QC metrics
graphics.off()
VlnPlot(project, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(project, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(project, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
featureScatterSave=paste(outputPath,paste(projectName,'featureScatter.png',sep='_'),sep='/')
png(file=featureScatterSave, width = 465, height = 225, units='mm', res = 300)
plot(plot1 + plot2)
graphics.off()
plot1 + plot2
project <- subset(project, subset = nFeature_RNA > nFeature_RNA_Var & nCount_RNA < nCount_RNA_Var & percent.mt < 10)

# Normalization
project <- NormalizeData(project, normalization.method = "LogNormalize", scale.factor = 10000)

# Feature selection
project <- FindVariableFeatures(project, selection.method = "vst", nfeatures = nfeatures)
top10 <- head(VariableFeatures(project), 10) # Identify the 10 most highly variable genes
plot1 <- VariableFeaturePlot(project)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
variableFeatureSave=paste(outputPath,paste(projectName,'variableFeatures.png',sep='_'),sep='/')
png(file=variableFeatureSave, width = 465, height = 225, units='mm', res = 300)
plot(plot1 + plot2) # plot variable features with and without labels
graphics.off()
plot1 + plot2

# Data scaling 
all.genes <- rownames(project)
project <- ScaleData(project, features = all.genes)

# Linear dimension reduction
project <- RunPCA(project, features = VariableFeatures(object = project), pcs.compute = pcsCompute)
print(project[["pca"]], dims = 1:5, nfeatures = 5)
vizDimSave=paste(outputPath,paste(projectName,'dimensionVisualisation.png',sep='_'),sep='/')
png(file=vizDimSave, width = 465, height = 225, units='mm', res = 300)
plot(VizDimLoadings(project, dims = 1:2, reduction = "pca"))
graphics.off()
VizDimLoadings(project, dims = 1:2, reduction = "pca")
dimPlotSave=paste(outputPath,paste(projectName,'dimensionPlot.png',sep='_'),sep='/')
png(file=dimPlotSave, width = 465, height = 225, units='mm', res = 300)
plot(DimPlot(project, reduction = "pca"))
graphics.off()
DimPlot(project, reduction = "pca")
#dimheatmapSave=paste(outputPath,paste(projectName,'dimensionHeatmap.png',sep='_'),sep='/')
#png(file=dimheatmapSave, width = 465, height = 225, units='mm', res = 300)
#plot(DimHeatmap(project, dims = 1:15, cells = 500, balanced = TRUE))
#graphics.off()
DimHeatmap(project, dims = 1:15, cells = 500, balanced = TRUE)

# Determine dimensionality 
project <- JackStraw(project, num.replicate = 100)
project <- ScoreJackStraw(project, dims = 1:20)
jackStrawPlotSave=paste(outputPath,paste(projectName,'jackStrawPlot.png',sep='_'),sep='/')
png(file=jackStrawPlotSave, width = 465, height = 225, units='mm', res = 300)
plot(JackStrawPlot(project, dims = 1:20))
graphics.off()
JackStrawPlot(project, dims = 1:20)
elbowPlotSave=paste(outputPath,paste(projectName,'elbowPlot.png',sep='_'),sep='/')
png(file=elbowPlotSave, width = 465, height = 225, units='mm', res = 300)
plot(ElbowPlot(project), ndims=nDimsElbow)
graphics.off()
ElbowPlot(project, ndims=nDimsElbow)

# Cell clusterization
project <- FindNeighbors(project, dims = dimensions, nn.method = 'annoy', annoy.metric = 'cosine')
project <- FindClusters(project, resolution = resolution, random.seed = randomSeed)
head(Idents(project), 5) # cluster IDs of the first 5 cells

# Non-linear dimensional reduction (UMAP/tSNE)
project <- RunUMAP(project, n.neighbors = nNeighbours,dims = dimensions, spread = spread, seed.use = seedUse)
UMAPsave=paste(outputPath,paste(projectName,'UMAP.png',sep='_'),sep='/')
png(file=UMAPsave, width = 465, height = 225, units='mm', res = 300)
plot(DimPlot(project, reduction = 'umap',label=T))
graphics.off()
DimPlot(project, reduction = 'umap',label=T)

# Cluster biomarkers finding
#cluster2.markers <- FindMarkers(project, ident.1 = "2", min.pct = 0.25) # find all markers of cluster 2
#head(cluster2.markers, n = 5) 
#cluster5.markers <- FindMarkers(project, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25) # find all markers distinguishing cluster 5 from clusters 0 and 3
#head(cluster5.markers, n = 5)
project.markers <- FindAllMarkers(project, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # find markers for every cluster compared to all remaining cells, report only the positive ones
project.markers <- project.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)
write.csv(project.markers, paste(outputPath,paste(projectName,'biomarkers.csv',sep='_'),sep='/'))
#cluster11.markers <- FindMarkers(project, ident.1 = 11, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
#head(cluster5.markers)
#VlnPlot(project, features = c("NV2g013117000.1", "NV2g022411000.1"))
#VlnPlot(project, features = c("NV2g013117000.1", "NV2g022411000.1"), slot = "counts", log = TRUE) # Raw counts
#FeaturePlot(project, features = c("NV2g014941000.1","NV2g004530000.1","NV2g004680000.1","NV2g000160000.1","NV2g024159000.1",
#                                  "NV2g016733000.1","NV2g001189000.1","NV2g023847000.1")) #NV2g023847000.1 : GABA
#cluster5.markers <- FindMarkers(project, ident.1 = 5, logfc.threshold = 2, test.use = "roc", only.pos = TRUE)
#FeaturePlot(project, features = c("NV2g018311000.1","NV2g011815000.1","NV2g009746000.1","NV2g012621000.1","NV2g020491000.1")) #NV2g020491000.1 : Nematocyst expressed protein
#cluster3.markers <- FindMarkers(project, ident.1 = 3, logfc.threshold = 2, test.use = "roc", only.pos = TRUE)
#FeaturePlot(project, features = c("NV2g006264000.1","NV2g023230000.1","NV2g019682000.1","NV2g005353000.1","NV2g007911000.1")) #Astrocyte marker / Neural ? 
#cluster4.markers <- FindMarkers(project, ident.1 = 4, logfc.threshold = 2, test.use = "roc", only.pos = TRUE)
#FeaturePlot(project, features = c("NV2g002652000.1","NV2g010686000.1","NV2g018480000.1","NV2g004235000.1")) #Agrin marker / Muscle marker ?
#cluster9.markers <- FindMarkers(project, ident.1 = 9, logfc.threshold = 2, test.use = "roc", only.pos = TRUE)
#FeaturePlot(project, features = c("NV2g005216000.1","NV2g001200000.1","NV2g017740000.1","NV2g012348000.1")) #Tetraspanin, pathogenesis-related / Immune cluster ?
#cluster8.markers <- FindMarkers(project, ident.1 = 8, logfc.threshold = 2, test.use = "roc", only.pos = TRUE)
#FeaturePlot(project, features = c("NV2g016402000.1","NV2g013863000.1","NV2g008421000.1")) #TF Sox11 - Neurogenesis cluster ?
#cluster0.markers <- FindMarkers(project, ident.1 = 0, logfc.threshold = 1, test.use = "roc", only.pos = TRUE)


# Fulldataset
FeaturePlot(project, features = c("NV2g024886000.1")) # Embryonic endomesoderm : ENSEMBL_ID lacking => NvSnailB
FeaturePlot(project, features = c("NV2g009345000.1")) # Gastrodermis => FRIS-like4
FeaturePlot(project, features = c("NV2g019134000.1")) # Retractor muscle => Nve-Tpm2
FeaturePlot(project, features = c("NV2g019682000.1")) # Embryonic ectoderm => TXT51-like
FeaturePlot(project, features = c("NV2g009154000.1")) # Embryonic epidermis : ENSEMBL_ID lacking => NVE4448
FeaturePlot(project, features = c("NV2g011441000.1")) # Pharyngeal ectoderm => NvFoxA
FeaturePlot(project, features = c("NV2g004477000.1")) # NPC => NvSoxB.2
FeaturePlot(project, features = c("NV2g009665000.1")) # Neuronal : ENSEMBL_ID lacking => NvAshA
FeaturePlot(project, features = c("NV2g008165000.1")) # Secretory : NVE lacking => VKT52-like
FeaturePlot(project, features = c("NV2g012902000.1")) # Mucous gland
FeaturePlot(project, features = c("NV2g005200000.1")) # Cnidocyte => EVA1C-like3
FeaturePlot(project, features = c("NV2g019749000.1")) # Cnidocyte mature => FOS-like

#FeaturePlot(project, features = c("NV2g011657000.1")) # Cnidocyte gastrula
                                 
                                  
#project.markers %>%
#  group_by(cluster) %>%
#  top_n(n = 10, wt = avg_log2FC) -> top10
#DoHeatmap(project, features = top10$gene) + NoLegend()

# Assigning cell type identity to clusters
#new.cluster.ids <- c("")
#names(new.cluster.ids) <- levels(project)
#project <- RenameIdents(project, new.cluster.ids)
#DimPlot(project, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



