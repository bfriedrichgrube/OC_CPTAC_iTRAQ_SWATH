## Make heatmap for Figure 3A and S2A

## Step 1: Prepare R session
setwd('')
load('CPTAC_SWATH_final_Set_20171101.Rdata')

library(RColorBrewer)
library(pheatmap)
library(MSnbase)

## Step 2: Load R data

load('CPTAC_DDA_Set.Rdata')
load('CPTAC_SWATH_Set.Rdata')
load('CPTAC_SWATH_WGCNA.Rdata')

## Step 3: Prepare heatmap input

heatmap_black <- exprs(CPTAC_SWATH_Set)[names(clustEG$black),]
heatmap_green <- exprs(CPTAC_SWATH_Set)[names(clustEG$green),]
heatmap_grey <- exprs(CPTAC_SWATH_Set)[names(clustEG$grey),]
heatmap_grey <- exprs(CPTAC_SWATH_Set)[names(clustEG$grey),]
heatmap_red <- exprs(CPTAC_SWATH_Set)[names(clustEG$red),]
heatmap_turquoise <- exprs(CPTAC_SWATH_Set)[names(clustEG$turquoise),]
heatmap_input <- rbind(heatmap_black, heatmap_green, heatmap_grey, heatmap_red, heatmap_turquoise)

SWATH_class <- pData(CPTAC_SWATH_Set)[,18]
names(SWATH_class)<- rownames(pData(CPTAC_SWATH_Set))
SWATH_class <- sort(SWATH_class)
heatmap_input <- as.matrix(heatmap_input[,names(SWATH_class)])
DDA_class <- pData(CPTAC_DDA_Set)[,18]

anno_cols_all <- data.frame(
	SWATH_class  = pData(CPTAC_SWATH_Set)[,18],
	DDA_class = tmp[,18],
	mRNA_class = tmp[,17])
rownames(anno_cols_all) <- rownames(tmp)

gene_modules <- c(rep.int('ECM Organization', times=86), rep.int('Immune response', times=42), rep.int('Metabolism',times=360), rep.int('Complement cascade', times=41), rep.int('Translation',times=1063))

anno_row <- data.frame(Gene_module = gene_modules)

rownames(anno_row) <- rownames(heatmap_input)

color_anno_all <- list(
	SWATH_class = c('1'='palegreen', '2'='lightpink', '3'='lightblue1'),
	DDA_class = c('1' = 'lightgoldenrod', '2'='lightpink', '3'='lightblue1'),
	mRNA_class = c('Differentiated'='red', 'Immunoreactive'='chartreuse', 'Mesenchymal'='deepskyblue3', 'Proliferative'='darkorchid4', 'Data not available'='white'),
	Original_Class = c('Differentiated'='indianred1', 'Immunoreactive'='darkolivegreen1', 'Mesenchymal'='skyblue2', 'Proliferative'='mediumorchid1', 'Stromal'='gold'),
	Gene_module=c('Complement cascade'='red', 'Immune response'='green', 'ECM Organization'='black', 'Translation'='lightblue', 'Metabolism'='grey'))

colors <- colorRampPalette(c('blue', 'white', 'red'))(n=10)
breaks <- c(seq(-8, -0.5, length.out=5), 0, seq(0.5, 8, length.out=5))

## Step 4: Make heatmap

pheatmap(heatmap_input, color=colors, breaks=breaks, show_rownames=F, show_colnames=F,
filename='heatmap_SWATH_overlap_WGCNA.pdf', border_color=NA, cluster_cols=F, cluster_rows=F,
annotation_col=anno_cols_all, annotation_colors=color_anno_all, annotation_row=anno_row)


