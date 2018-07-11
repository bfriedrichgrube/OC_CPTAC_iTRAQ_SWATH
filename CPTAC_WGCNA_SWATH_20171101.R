
library(MSnbase)
library(WGCNA)
library(RColorBrewer)
library(org.Hs.eg.db)
library(knitr)
library(dplyr)
library(ReactomePA)
library(clusterProfiler)

setwd('Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_Classification/')
> 
> options(stringsAsFactors=F)
> enableWGCNAThreads()
Allowing parallel execution with up to 3 working processes.
> load('CPTAC_SWATH_final_Set_20171101.Rdata.Rdata')
Error in readChar(con, 5L, useBytes = TRUE) : cannot open the connection
In addition: Warning message:
In readChar(con, 5L, useBytes = TRUE) :
  cannot open compressed file 'CPTAC_SWATH_final_Set_20171101.Rdata.Rdata', probable reason 'No such file or directory'
> load('CPTAC_SWATH_final_Set_20171101.Rdata')
> powers =1:20
> sft <- pickSoftThreshold(t(exprs(CPTAC_SWATH_Set)), powerVector=powers, networkType='signed', verbose=0)
   Power SFT.R.sq  slope truncated.R.sq mean.k. median.k.  max.k.
1      1   0.6080  7.110          0.796 946.000   955.000 1070.00
2      2   0.4140  2.670          0.805 571.000   578.000  735.00
3      3   0.1660  1.090          0.865 352.000   355.000  517.00
4      4   0.0162  0.257          0.910 221.000   221.000  372.00
5      5   0.0377 -0.353          0.957 141.000   139.000  272.00
6      6   0.1690 -0.715          0.988  91.100    87.900  202.00
7      7   0.3590 -1.120          0.974  59.900    56.400  153.00
8      8   0.4850 -1.410          0.969  39.900    36.600  116.00
9      9   0.6040 -1.640          0.971  27.000    23.800   90.00
10    10   0.6810 -1.830          0.973  18.500    15.600   70.40
11    11   0.7520 -1.930          0.984  12.800    10.400   55.50
12    12   0.8110 -2.000          0.988   8.980     6.980   44.20
13    13   0.8250 -2.140          0.983   6.370     4.690   35.60
14    14   0.8560 -2.180          0.985   4.560     3.200   28.80
15    15   0.8730 -2.210          0.984   3.310     2.210   23.50
16    16   0.8820 -2.200          0.978   2.420     1.540   19.30
17    17   0.8770 -2.230          0.962   1.790     1.090   16.00
18    18   0.8910 -2.200          0.966   1.340     0.775   13.30
19    19   0.9040 -2.160          0.975   1.010     0.550   11.20
20    20   0.9160 -2.150          0.978   0.768     0.390    9.41
> 
> pdf(file='softthreshold_SWATH_20171101.pdf',width=20, height=10)
> par(mfrow=c(1,2))
>  cex1=0.9
> 
>  plot(sft$fitIndices[,1],
+  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
+  xlab='Soft Threshold (power)',
+  ylab='Scale Free Topology Model Fit, signed R^2', type='n', 
+  main=paste('Scale indipendence'))
> 
>  text(sft$fitIndices[,1], 
+  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
+  labels=powers, cex=cex1, col='red')
> 
>  abline(h=0.9, col='red')
> 
>  plot(sft$fitIndices[,1],
+  sft$fitIndices[,5],
+  xlab='Soft Threshold (power)',
+  ylab='Mean connectivity',
+  type='n', main=paste('Mean connectivity'))
> 
>  text(sft$fitIndices[,1], 
+  sft$fitIndices[,5],
+  labels=powers,
+  cex=cex1,
+  col='red')
> 
> dev.off()
null device 
          1 
> -sign(sft$fitIndices[19,3])*sft$fitIndices[19,2],
Error: unexpected ',' in "-sign(sft$fitIndices[19,3])*sft$fitIndices[19,2],"
> -sign(sft$fitIndices[19,3])*sft$fitIndices[19,2]
[1] 0.9038705
> -sign(sft$fitIndices[19,3])*sft$fitIndices[20,2]
[1] 0.9157039
> -sign(sft$fitIndices[20,3])*sft$fitIndices[20,2]
[1] 0.9157039
> -sign(sft$fitIndices[18,3])*sft$fitIndices[18,2]
[1] 0.8911925
> softPower=19
> 
>  adjacency.mat <- adjacency(t(exprs(CPTAC_SWATH_Set)), power=softPower, type='signed')
> dim(adjacency.mat)
[1] 1599 1599
> [1] 1916 1916
Error: unexpected '[' in "["
> TOM <- TOMsimilarity(adjacency.mat, TOMType='unsigned')
..connectivity..
..matrix multiplication..
..normalization..
..done.
> 
> dissTOM <- 1-TOM
> 
> geneTree <- hclust(as.dist(dissTOM), method='average')
> 
> pdf(file='genetree_SWATH_20171101.pdf', width=20, height=10)
> 
> plot(geneTree, xlab='', sub='', 
+ main='Gene clustering on TOM-based dissimilarity',
+ labels=F, hang=0.04)
> 
> dev.off()
null device 
          1 
> 
> 
> minModuleSize <- 20
> 
> dynamicMods <- cutreeDynamic(dendro=geneTree, 
+  distM=dissTOM, 
+  minClusterSize=minModuleSize)
 ..cutHeight not given, setting it to 0.999  ===>  99% of the (truncated) height range in dendro.
 ..done.
> 
> dynamicColors <- labels2colors(dynamicMods, colorSeq=standardColors())
> 
> pdf(file='genetree_colors_SWATH_20171101.pdf', width=20, height=10)
> 
> plotDendroAndColors(geneTree, dynamicColors, 'Dynamic Tree Cut',
+  dendroLabels=F, hang=0.03, 
+  addGuide=T, guideHang=0.05,
+  main='Gene dendrogram and module colors')
> 
> dev.off()
null device 
          1 
> 
> 
> MEList <- moduleEigengenes(t(exprs(CPTAC_SWATH_Set)), colors=dynamicColors)
> MEs <- MEList$eigengenes
> 
> MEDiss <- 1- cor(MEs)
> 
> METree <- hclust(as.dist(MEDiss), method='average')
> 
> 
> pdf(file='module_eigengenes_SWATH_20171101.pdf', width=10, height=10)
> 
> plot(METree, main='Clustering of module eigengenes', 
+  xlab='', sub='')
> 
> MEDissThres <- 0.4
> 
> abline(h=MEDissThres, col='red')
> 
> dev.off()
null device 
          1 
> 
> merge <- mergeCloseModules(t(exprs(CPTAC_SWATH_Set)),
+ dynamicColors, 
+ cutHeight=MEDissThres, 
+ verbose=3)
 mergeCloseModules: Merging modules whose distance is less than 0.4
   multiSetMEs: Calculating module MEs.
     Working on set 1 ...
     moduleEigengenes: Calculating 8 module eigengenes in given set.
^[[D   multiSetMEs: Calculating module MEs.
     Working on set 1 ...
     moduleEigengenes: Calculating 5 module eigengenes in given set.
   Calculating new MEs...
   multiSetMEs: Calculating module MEs.
     Working on set 1 ...
     moduleEigengenes: Calculating 5 module eigengenes in given set.
> 
> moduleColors <- merge$colors
> 
> pdf(file='genetree_colors_merged_SWATH_20171101.pdf', height=10, width=20)
> 
> plotDendroAndColors(geneTree, cbind(dynamicColors, moduleColors), 
+  c('Dynamic Tree Cut', 'Merged dynamic'),
+  dendroLabels=F, hang=0.03, 
+  addGuide=T, guideHang=0.05, 
+  marAll=c(1,7,3,1))
> 
> dev.off()
null device 
          1 
> 
> 
> MEs0 <- moduleEigengenes(t(exprs(CPTAC_SWATH_Set)),moduleColors)$eigengenes
> MEs <- orderMEs(MEs0)
> datTraits <- pData(CPTAC_SWATH_Set)[,grep('binary', varLabels(CPTAC_SWATH_Set))]
> 
> colnames(datTraits)<- sub('.binary', '', colnames(datTraits))
> 
> head(datTraits)
            1     2     3
HRD_1    TRUE FALSE FALSE
HRD_10  FALSE  TRUE FALSE
HRD_101  TRUE FALSE FALSE
HRD_102  TRUE FALSE FALSE
HRD_103  TRUE FALSE FALSE
HRD_11  FALSE  TRUE FALSE
> moduleTraitCor <- cor(MEs, datTraits, use='p')
> moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, ncol(CPTAC_SWATH_Set))
> textMatrix <- paste(signif(moduleTraitCor, 2),
+ '\n(',
+ signif(moduleTraitPvalue, 1),
+ ')', sep='')
> 
> textMatrix<- matrix(textMatrix, ncol=ncol(moduleTraitCor))
> 
> pdf(file='module_trait_association_SWATH_20171101.pdf', height=10, width=10)
> 
> par(mar=c(6,8,3,2))
> 
> labeledHeatmap(Matrix=moduleTraitCor,
+  xLabels=colnames(datTraits),
+  yLabels=names(MEs),
+  ySymbols=names(MEs), 
+  colorLabels=F,
+  colors=blueWhiteRed(50),
+  textMatrix=textMatrix, 
+  setStdMargins=F,
+  cex.text=0.75,
+  zlim=c(-1,1),
+  main=paste('Module-Trait Association (p-value)'))
> 
> dev.off()
null device 
          1 
> 
> 
> features <- featureNames(CPTAC_SWATH_Set)
> 
> clustList <- tapply(features, moduleColors, c)
> 
> clustList <- lapply(clustList, intersect, mappedRkeys(org.Hs.egALIAS2EG))
> 
> clustEG <- lapply(clustList, function(x) sapply(as.list(org.Hs.egALIAS2EG[x]), '[', 1))
> 
> summary(clustEG)
          Length Class  Mode     
black       86   -none- character
green       42   -none- character
grey       360   -none- character
red         41   -none- character
turquoise 1063   -none- character
> xx <- compareCluster(clustEG, fun='enrichPathway',
+  organism='human', pvalueCutoff=0.05, 
+  qvalueCutoff=0.05, 
+  pAdjustMethod='none',
+  universe=unlist(clustEG))
^[[D> 
> pdf(file='cluster_enrichment_SWATH_20171101.pdf', height=10, width=10)
> plot(xx, colorBy='qvalue')
> dev.off()
null device 
          1 
> 
> yy <- compareCluster(clustEG, fun='enrichKEGG', 
+  organism='human', pvalueCutoff=0.05, 
+  qvalueCutoff=0.05, 
+  pAdjustMethod='none', 
+  universe=unlist(clustEG))
> 
> pdf(file='cluster_enrichmentKEGG_SWATH_20171101.pdf', height=10, width=10)
> plot(yy, colorBy='qvalue')
> dev.off()
null device 
          1 
> 
> zz <- compareCluster(clustEG, fun='enrichGO', 
+  ont='BP', organism='human', pvalueCutoff=0.05, 
+  qvalueCutoff=0.05, 
+  pAdjustMethod='none', 
+  universe=unlist(clustEG))
Loading required package: GO.db
> 
> pdf(file='cluster_enrichmentGO_SWATH_20171101.pdf', height=10, width=10)
> plot(zz, colorBy='qvalue')
> dev.off()
null device 
          1 
> 
> save(clustEG, file='CPTAC_SWATH_WGCNA_20171101.Rdata')
> 

