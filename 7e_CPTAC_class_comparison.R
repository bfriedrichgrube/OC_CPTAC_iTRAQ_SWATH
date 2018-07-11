## Make group comparisons based on original  CPTAC classification ##
## Mesechymal to others should show significantly higher overlap in proteins regulated than other groups## 


## Step 1: Prepare R session 

#common proteins from mesenchymal class will be loaded ->  a and d
load('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_common_proteins_mesenchymal_20180312.Rdata')

#load DDA data -> DDA
load('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_Data/CPTAC_processed_data/CPTAC_DDA_protein_level_CDAP.Rdata')

#load SWATH data -> SWATH_imputed
load('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_Classification/SWATH_classification_workflow_20171101.Rdata')


## Step 2: Prepare data

dim(SWATH_imputed)
[1] 1659  103

dim(DDA)
[1] 8600  103
missing <- NULL
for(i in 1:nrow(DDA)){tmp <-length(which(is.na(DDA[i,])))
missing <- c(missing, tmp)
}

DDA_missing0 <- DDA[(missing==0),] 
dim(DDA_missing0)
[1] 4366  103

gene.ids <- intersect(rownames(SWATH_imputed), rownames(DDA_missing0))

SWATH_overlap <- SWATH_imputed[gene.ids, ]
DDA_overlap <- as.matrix(DDA_missing0[gene.ids,])

## Step 3: Make group comparison for Mesenchymal group and check overlap

z <- c(80, 23, 27, 92, 97, 51, 60, 8:9, 71:72, 11, 15, 18, 86, 28, 30, 32:33, 37, 40, 93, 42, 99)
length(z)
[1] 24

Mes_SWATH_pvalue <- NULL
Mes_SWATH_log2FC <- NULL

for(i in 1:nrow(SWATH_overlap)){
	tmp <- t.test(SWATH_overlap[i, z], SWATH_overlap[i,-z])
	tmp2 <- mean(SWATH_overlap[i, z])-mean(SWATH_overlap[i, -z])
	Mes_SWATH_pvalue <- c(Mes_SWATH_pvalue, tmp$p.value)
	Mes_SWATH_log2FC <- c(Mes_SWATH_log2FC, tmp2)
}

Mes_SWATH_adjust <- p.adjust(Mes_SWATH_pvalue, method='fdr')
Mes_SWATH_output <- data.frame(Protein=rownames(SWATH_overlap), log2FC=Mes_SWATH_log2FC, pvalue=Mes_SWATH_pvalue, padjust=Mes_SWATH_adjust)

upreg <- Mes_SWATH_output[which(Mes_SWATH_output$padjust<0.05 & Mes_SWATH_output$log2FC>log2(1.3)),]
downreg <- Mes_SWATH_output[which(Mes_SWATH_output$padjust<0.05 & Mes_SWATH_output$log2FC<(-log2(1.3))),]

dim(upreg)
[1] 113   4
dim(downreg)
[1] 46  4

p <- c(4:9, 13:14, 23, 26:27, 29, 36, 43, 51, 60, 62, 64:65, 69, 72, 74:75, 81)

Mes_DDA_pvalue <- NULL
Mes_DDA_log2FC <- NULL
for(i in 1:nrow(DDA_overlap)){
tmp <- t.test(DDA_overlap[i, p], DDA_overlap[i,-p])
tmp2 <- mean(DDA_overlap[i, p])-mean(DDA_overlap[i, -p])
Mes_DDA_pvalue <- c(Mes_DDA_pvalue, tmp$p.value)
Mes_DDA_log2FC <- c(Mes_DDA_log2FC, tmp2)
}
Mes_DDA_adjust <- p.adjust(Mes_DDA_pvalue, method='fdr')
Mes_DDA_output <- data.frame(Protein=rownames(DDA_overlap), log2FC=Mes_DDA_log2FC, pvalue=Mes_DDA_pvalue, padjust=Mes_DDA_adjust)

upreg_DDA <- Mes_DDA_output[which(Mes_DDA_output$padjust<0.05 & Mes_DDA_output$log2FC>log2(1.3)),]
downreg_DDA <- Mes_DDA_output[which(Mes_DDA_output$padjust<0.05 & Mes_DDA_output$log2FC<(-log2(1.3))),]

dim(upreg_DDA)
[1] 80  4
dim(downreg_DDA)
[1] 9 4

#overlap of upregulated proteins from mesenchymal 
length(intersect(upreg$Protein, upreg_DDA$Protein))
[1] 73
#overlap of downregulated proteins from mesenchymal
length(intersect(downreg$Protein, downreg_DDA$Protein))
[1] 6

#the mesenchymal protein module was positively correlated with the patient group, therefore the upregulated proteins will be chekced for overlap to the protein modules from iTRAQ and SWATH

#overlap of common module proteins and SWATH upreg
az <- intersect(names(d), upreg$Protein)
length(az)
[1] 69

#overlap of common module proteins and iTRAQ upreg
ax <- intersect(names(d), upreg_DDA$Protein)
length(ax)
[1] 66

#overlap of all 3 groups (common module proteins, and both upregulated vectors)
length(intersect(ax, az))
[1] 62

length(d)
[1] 84

#plot

pdf('Mesenchymal_SWATH_Volcano_20180314.pdf', width=10, height=10)
with(Mes_SWATH_output, plot(log2FC, -log10(padjust), pch=20, cex=2, cex.lab=0.1, yaxt='n', cex.axis=2.2, xlab='', ylab='', xlim=c(-2, 2), ylim=c(0,22)))
with(subset(Mes_SWATH_output, padjust<0.05 & log2FC>log2(1.3)), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='red'))
with(subset(Mes_SWATH_output, padjust<0.05 & log2FC<(-log2(1.3))), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='dodgerblue'))
axis(2, las=2, cex.axis=2.2)
abline(v=c(log2(1.3),-log2(1.3)), col=c("darkgray","darkgray"), lty=c(2,2))
abline(h=-log10(0.05), col="darkgray", lty=2)
dev.off()


pdf('Mesenchymal_DDA_Volcano_20180314.pdf', width=10, height=10)
with(Mes_DDA_output, plot(log2FC, -log10(padjust), pch=20, cex=2, cex.lab=0.1, yaxt='n', cex.axis=2.2, xlab='', ylab='', xlim=c(-2, 2), ylim=c(0,22)))
with(subset(Mes_DDA_output, padjust<0.05 & log2FC>log2(1.3)), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='red'))
with(subset(Mes_DDA_output, padjust<0.05 & log2FC<(-log2(1.3))), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='dodgerblue'))
axis(2, las=2, cex.axis=2.2)
abline(v=c(log2(1.3),-log2(1.3)), col=c("darkgray","darkgray"), lty=c(2,2))
abline(h=-log10(0.05), col="darkgray", lty=2)
dev.off()



pdf('venn_mesenchymal_Volcano_20180312.pdf', width=10, height=10)
 draw.triple.venn(84, 80, 113, 66, 73, 69, 62, lty=rep('blank', 3), fill=c('lightpink', 'lightblue', 'green'), cex=3)
dev.off()


common <- intersect(ax,az)
head(common)
[1] "SEPT11"   "COL3A1"   "CNN2"     "C1S"      "MYL9"     "SERPINF1"
write.table(common, 'common_proteins_20180313.txt' , quote=F, row.names=F)


## Step 4: Make group comparison for Proliferative group and check overlap


z <- c(74, 6, 68:69, 76:77, 17, 20, 29, 31, 35, 95, 50, 52, 54, 56)

Pro_SWATH_pvalue <- NULL
Pro_SWATH_log2FC <- NULL

for(i in 1:nrow(SWATH_overlap)){
	tmp <- t.test(SWATH_overlap[i, z], SWATH_overlap[i,-z])
	tmp2 <- mean(SWATH_overlap[i, z])-mean(SWATH_overlap[i, -z])
	Pro_SWATH_pvalue <- c(Pro_SWATH_pvalue, tmp$p.value)
	Pro_SWATH_log2FC <- c(Pro_SWATH_log2FC, tmp2)
}

Pro_SWATH_adjust <- p.adjust(Pro_SWATH_pvalue, method='fdr')
Pro_SWATH_output <- data.frame(Protein=rownames(SWATH_overlap), log2FC=Pro_SWATH_log2FC, pvalue=Pro_SWATH_pvalue, padjust=Pro_SWATH_adjust)

upreg <- Pro_SWATH_output[which(Pro_SWATH_output$padjust<0.05 & Pro_SWATH_output$log2FC>log2(1.3)),]
downreg <- Pro_SWATH_output[which(Pro_SWATH_output$padjust<0.05 & Pro_SWATH_output$log2FC<(-log2(1.3))),]

dim(upreg)
[1] 210   4
dim(downreg)
[1] 80  4


p <- c(3, 11, 21, 22, 31, 35, 40, 45, 61, 63, 67, 78, 88, 91, 94, 96)

Pro_DDA_pvalue <- NULL
Pro_DDA_log2FC <- NULL
for(i in 1:nrow(DDA_overlap)){
tmp <- t.test(DDA_overlap[i, p], DDA_overlap[i,-p])
tmp2 <- mean(DDA_overlap[i, p])-mean(DDA_overlap[i, -p])
Pro_DDA_pvalue <- c(Pro_DDA_pvalue, tmp$p.value)
Pro_DDA_log2FC <- c(Pro_DDA_log2FC, tmp2)
}
Pro_DDA_adjust <- p.adjust(Pro_DDA_pvalue, method='fdr')
Pro_DDA_output <- data.frame(Protein=rownames(DDA_overlap), log2FC=Pro_DDA_log2FC, pvalue=Pro_DDA_pvalue, padjust=Pro_DDA_adjust)

upreg_DDA <- Pro_DDA_output[which(Pro_DDA_output$padjust<0.05 & Pro_DDA_output$log2FC>log2(1.3)),]
downreg_DDA <- Pro_DDA_output[which(Pro_DDA_output$padjust<0.05 & Pro_DDA_output$log2FC<(-log2(1.3))),]

dim(upreg_DDA)
[1] 30  4
dim(downreg_DDA)
[1] 131   4



#overlap of upregulated proteins from proliferative 
length(intersect(upreg$Protein, upreg_DDA$Protein))
[1] 24
#overlap of downregulated proteins from proliferative
length(intersect(downreg$Protein, downreg_DDA$Protein))
[1] 76


#plot

pdf('Proliferative_SWATH_Volcano_20180314.pdf', width=10, height=10)
with(Pro_SWATH_output, plot(log2FC, -log10(padjust), pch=20, cex=2, cex.lab=0.1, yaxt='n', cex.axis=2.2, xlab='', ylab='', xlim=c(-2, 2), ylim=c(0,22)))
with(subset(Pro_SWATH_output, padjust<0.05 & log2FC>log2(1.3)), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='red'))
with(subset(Pro_SWATH_output, padjust<0.05 & log2FC<(-log2(1.3))), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='dodgerblue'))
axis(2, las=2, cex.axis=2.2)
abline(v=c(log2(1.3),-log2(1.3)), col=c("darkgray","darkgray"), lty=c(2,2))
abline(h=-log10(0.05), col="darkgray", lty=2)
dev.off()


pdf('Proliferative_DDA_Volcano_20180314.pdf', width=10, height=10)
with(Pro_DDA_output, plot(log2FC, -log10(padjust), pch=20, cex=2, cex.lab=0.1, yaxt='n', cex.axis=2.2, xlab='', ylab='', xlim=c(-2, 2), ylim=c(0,22)))
with(subset(Pro_DDA_output, padjust<0.05 & log2FC>log2(1.3)), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='red'))
with(subset(Pro_DDA_output, padjust<0.05 & log2FC<(-log2(1.3))), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='dodgerblue'))
axis(2, las=2, cex.axis=2.2)
abline(v=c(log2(1.3),-log2(1.3)), col=c("darkgray","darkgray"), lty=c(2,2))
abline(h=-log10(0.05), col="darkgray", lty=2)
dev.off()





## Step 5: Make group comparison for Stromal group and check overlap

z <- c(84, 88, 24, 96, 5)

Stro_SWATH_pvalue <- NULL
Stro_SWATH_log2FC <- NULL

for(i in 1:nrow(SWATH_overlap)){
	tmp <- t.test(SWATH_overlap[i, z], SWATH_overlap[i,-z])
	tmp2 <- mean(SWATH_overlap[i, z])-mean(SWATH_overlap[i, -z])
	Stro_SWATH_pvalue <- c(Stro_SWATH_pvalue, tmp$p.value)
	Stro_SWATH_log2FC <- c(Stro_SWATH_log2FC, tmp2)
}

Stro_SWATH_adjust <- p.adjust(Stro_SWATH_pvalue, method='fdr')
Stro_SWATH_output <- data.frame(Protein=rownames(SWATH_overlap), log2FC=Stro_SWATH_log2FC, pvalue=Stro_SWATH_pvalue, padjust=Stro_SWATH_adjust)

upreg <- Stro_SWATH_output[which(Stro_SWATH_output$padjust<0.05 & Stro_SWATH_output$log2FC>log2(1.3)),]
downreg <- Stro_SWATH_output[which(Stro_SWATH_output$padjust<0.05 & Stro_SWATH_output$log2FC<(-log2(1.3))),]

dim(upreg)
[1] 4 4
dim(downreg)
[1] 40  4


p <- c(49, 53, 56, 79, 103)

Stro_DDA_pvalue <- NULL
Stro_DDA_log2FC <- NULL
for(i in 1:nrow(DDA_overlap)){
tmp <- t.test(DDA_overlap[i, p], DDA_overlap[i,-p])
tmp2 <- mean(DDA_overlap[i, p])-mean(DDA_overlap[i, -p])
Stro_DDA_pvalue <- c(Stro_DDA_pvalue, tmp$p.value)
Stro_DDA_log2FC <- c(Stro_DDA_log2FC, tmp2)
}
Stro_DDA_adjust <- p.adjust(Stro_DDA_pvalue, method='fdr')
Stro_DDA_output <- data.frame(Protein=rownames(DDA_overlap), log2FC=Stro_DDA_log2FC, pvalue=Stro_DDA_pvalue, padjust=Stro_DDA_adjust)

upreg_DDA <- Stro_DDA_output[which(Stro_DDA_output$padjust<0.05 & Stro_DDA_output$log2FC>log2(1.3)),]
downreg_DDA <- Stro_DDA_output[which(Stro_DDA_output$padjust<0.05 & Stro_DDA_output$log2FC<(-log2(1.3))),]

dim(upreg_DDA)
[1] 3 4
dim(downreg_DDA)
[1] 1 4


#overlap of upregulated proteins from stromal 
length(intersect(upreg$Protein, upreg_DDA$Protein))
[1] 1
#overlap of downregulated proteins from stromal
length(intersect(downreg$Protein, downreg_DDA$Protein))
[1] 0


#plot

pdf('Stromal_SWATH_Volcano_20180314.pdf', width=10, height=10)
with(Stro_SWATH_output, plot(log2FC, -log10(padjust), pch=20, cex=2, cex.lab=0.1, yaxt='n', cex.axis=2.2, xlab='', ylab='', xlim=c(-2, 2), ylim=c(0,22)))
with(subset(Stro_SWATH_output, padjust<0.05 & log2FC>log2(1.3)), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='red'))
with(subset(Stro_SWATH_output, padjust<0.05 & log2FC<(-log2(1.3))), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='dodgerblue'))
axis(2, las=2, cex.axis=2.2)
abline(v=c(log2(1.3),-log2(1.3)), col=c("darkgray","darkgray"), lty=c(2,2))
abline(h=-log10(0.05), col="darkgray", lty=2)
dev.off()


pdf('Stromal_DDA_Volcano_20180314.pdf', width=10, height=10)
with(Stro_DDA_output, plot(log2FC, -log10(padjust), pch=20, cex=2, cex.lab=0.1, yaxt='n', cex.axis=2.2, xlab='', ylab='', xlim=c(-2, 2), ylim=c(0,22)))
with(subset(Stro_DDA_output, padjust<0.05 & log2FC>log2(1.3)), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='red'))
with(subset(Stro_DDA_output, padjust<0.05 & log2FC<(-log2(1.3))), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='dodgerblue'))
axis(2, las=2, cex.axis=2.2)
abline(v=c(log2(1.3),-log2(1.3)), col=c("darkgray","darkgray"), lty=c(2,2))
abline(h=-log10(0.05), col="darkgray", lty=2)
dev.off()


## Step 6: Make group comparison for immunoreactive group and check overlap

z <- c(7, 61, 64:67, 75, 78:79, 82, 19, 90, 25, 38, 41, 98, 44:45, 49, 53, 55, 58:59, 4)

Immu_SWATH_pvalue <- NULL
Immu_SWATH_log2FC <- NULL

for(i in 1:nrow(SWATH_overlap)){
	tmp <- t.test(SWATH_overlap[i, z], SWATH_overlap[i,-z])
	tmp2 <- mean(SWATH_overlap[i, z])-mean(SWATH_overlap[i, -z])
	Immu_SWATH_pvalue <- c(Immu_SWATH_pvalue, tmp$p.value)
	Immu_SWATH_log2FC <- c(Immu_SWATH_log2FC, tmp2)
}

Immu_SWATH_adjust <- p.adjust(Immu_SWATH_pvalue, method='fdr')
Immu_SWATH_output <- data.frame(Protein=rownames(SWATH_overlap), log2FC=Immu_SWATH_log2FC, pvalue=Immu_SWATH_pvalue, padjust=Immu_SWATH_adjust)

upreg <- Immu_SWATH_output[which(Immu_SWATH_output$padjust<0.05 & Immu_SWATH_output$log2FC>log2(1.3)),]
downreg <- Immu_SWATH_output[which(Immu_SWATH_output$padjust<0.05 & Immu_SWATH_output$log2FC<(-log2(1.3))),]

dim(upreg)
[1] 54  4
dim(downreg)
[1] 27  4


p <- c(2, 12, 15, 18:20, 30, 37:38, 42, 44, 55, 58, 70, 73, 80, 82:83, 87, 93, 95, 98, 100, 102)

Immu_DDA_pvalue <- NULL
Immu_DDA_log2FC <- NULL
for(i in 1:nrow(DDA_overlap)){
tmp <- t.test(DDA_overlap[i, p], DDA_overlap[i,-p])
tmp2 <- mean(DDA_overlap[i, p])-mean(DDA_overlap[i, -p])
Immu_DDA_pvalue <- c(Immu_DDA_pvalue, tmp$p.value)
Immu_DDA_log2FC <- c(Immu_DDA_log2FC, tmp2)
}
Immu_DDA_adjust <- p.adjust(Immu_DDA_pvalue, method='fdr')
Immu_DDA_output <- data.frame(Protein=rownames(DDA_overlap), log2FC=Immu_DDA_log2FC, pvalue=Immu_DDA_pvalue, padjust=Immu_DDA_adjust)

upreg_DDA <- Immu_DDA_output[which(Immu_DDA_output$padjust<0.05 & Immu_DDA_output$log2FC>log2(1.3)),]
downreg_DDA <- Immu_DDA_output[which(Immu_DDA_output$padjust<0.05 & Immu_DDA_output$log2FC<(-log2(1.3))),]

dim(upreg_DDA)
[1] 31  4
dim(downreg_DDA)
[1] 54  4


#overlap of upregulated proteins from immunoreactive 
length(intersect(upreg$Protein, upreg_DDA$Protein))
[1] 24

#overlap of downregulated proteins from immunoreactive
length(intersect(downreg$Protein, downreg_DDA$Protein))
[1] 23


#plot

pdf('Immunoreactive_SWATH_Volcano_20180314.pdf', width=10, height=10)
with(Immu_SWATH_output, plot(log2FC, -log10(padjust), pch=20, cex=2, cex.lab=0.1, yaxt='n', cex.axis=2.2, xlab='', ylab='', xlim=c(-2, 2), ylim=c(0,22)))
with(subset(Immu_SWATH_output, padjust<0.05 & log2FC>log2(1.3)), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='red'))
with(subset(Immu_SWATH_output, padjust<0.05 & log2FC<(-log2(1.3))), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='dodgerblue'))
axis(2, las=2, cex.axis=2.2)
abline(v=c(log2(1.3),-log2(1.3)), col=c("darkgray","darkgray"), lty=c(2,2))
abline(h=-log10(0.05), col="darkgray", lty=2)
dev.off()


pdf('Immunoreactive_DDA_Volcano_20180314.pdf', width=10, height=10)
with(Immu_DDA_output, plot(log2FC, -log10(padjust), pch=20, cex=2, cex.lab=0.1, yaxt='n', cex.axis=2.2, xlab='', ylab='', xlim=c(-2, 2), ylim=c(0,22)))
with(subset(Immu_DDA_output, padjust<0.05 & log2FC>log2(1.3)), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='red'))
with(subset(Immu_DDA_output, padjust<0.05 & log2FC<(-log2(1.3))), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='dodgerblue'))
axis(2, las=2, cex.axis=2.2)
abline(v=c(log2(1.3),-log2(1.3)), col=c("darkgray","darkgray"), lty=c(2,2))
abline(h=-log10(0.05), col="darkgray", lty=2)
dev.off()



## Step 7: Make group comparison for Differentiated group and check overlap

z <- c(1:3, 62:63, 10, 70, 73, 12:14, 16, 81, 83, 21:22, 85, 87, 89, 91, 26, 34, 36, 39, 94, 43, 46:48, 100:103, 57)

Diff_SWATH_pvalue <- NULL
Diff_SWATH_log2FC <- NULL

for(i in 1:nrow(SWATH_overlap)){
	tmp <- t.test(SWATH_overlap[i, z], SWATH_overlap[i,-z])
	tmp2 <- mean(SWATH_overlap[i, z])-mean(SWATH_overlap[i, -z])
	Diff_SWATH_pvalue <- c(Diff_SWATH_pvalue, tmp$p.value)
	Diff_SWATH_log2FC <- c(Diff_SWATH_log2FC, tmp2)
}

Diff_SWATH_adjust <- p.adjust(Diff_SWATH_pvalue, method='fdr')
Diff_SWATH_output <- data.frame(Protein=rownames(SWATH_overlap), log2FC=Diff_SWATH_log2FC, pvalue=Diff_SWATH_pvalue, padjust=Diff_SWATH_adjust)

upreg <- Diff_SWATH_output[which(Diff_SWATH_output$padjust<0.05 & Diff_SWATH_output$log2FC>log2(1.3)),]
downreg <- Diff_SWATH_output[which(Diff_SWATH_output$padjust<0.05 & Diff_SWATH_output$log2FC<(-log2(1.3))),]

dim(upreg)
[1] 1 4
dim(downreg)
[1] 0 4



p <- c(1, 10, 16:17, 24:25, 28, 32:34, 39, 41, 46:48, 50, 52, 54, 57, 59, 66, 68, 71, 76:77, 84:86, 89:90, 92, 97, 99, 101)

Diff_DDA_pvalue <- NULL
Diff_DDA_log2FC <- NULL
for(i in 1:nrow(DDA_overlap)){
tmp <- t.test(DDA_overlap[i, p], DDA_overlap[i,-p])
tmp2 <- mean(DDA_overlap[i, p])-mean(DDA_overlap[i, -p])
Diff_DDA_pvalue <- c(Diff_DDA_pvalue, tmp$p.value)
Diff_DDA_log2FC <- c(Diff_DDA_log2FC, tmp2)
}
Diff_DDA_adjust <- p.adjust(Diff_DDA_pvalue, method='fdr')
Diff_DDA_output <- data.frame(Protein=rownames(DDA_overlap), log2FC=Diff_DDA_log2FC, pvalue=Diff_DDA_pvalue, padjust=Diff_DDA_adjust)

upreg_DDA <- Diff_DDA_output[which(Diff_DDA_output$padjust<0.05 & Diff_DDA_output$log2FC>log2(1.3)),]
downreg_DDA <- Diff_DDA_output[which(Diff_DDA_output$padjust<0.05 & Diff_DDA_output$log2FC<(-log2(1.3))),]

dim(upreg_DDA)
[1] 12  4
dim(downreg_DDA)
[1] 1 4

#overlap of upregulated proteins from differentiated
length(intersect(upreg$Protein, upreg_DDA$Protein))
[1] 0
#overlap of downregulated proteins from differentiated
length(intersect(downreg$Protein, downreg_DDA$Protein))
[1] 0

#plot

pdf('Differentiated_SWATH_Volcano_20180314.pdf', width=10, height=10)
with(Diff_SWATH_output, plot(log2FC, -log10(padjust), pch=20, cex=2, cex.lab=0.1, yaxt='n', cex.axis=2.2, xlab='', ylab='', xlim=c(-2, 2), ylim=c(0,22)))
with(subset(Diff_SWATH_output, padjust<0.05 & log2FC>log2(1.3)), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='red'))
with(subset(Diff_SWATH_output, padjust<0.05 & log2FC<(-log2(1.3))), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='dodgerblue'))
axis(2, las=2, cex.axis=2.2)
abline(v=c(log2(1.3),-log2(1.3)), col=c("darkgray","darkgray"), lty=c(2,2))
abline(h=-log10(0.05), col="darkgray", lty=2)
dev.off()


pdf('Differentiated_DDA_Volcano_20180314.pdf', width=10, height=10)
with(Diff_DDA_output, plot(log2FC, -log10(padjust), pch=20, cex=2, cex.lab=0.1, yaxt='n', cex.axis=2.2, xlab='', ylab='', xlim=c(-2, 2), ylim=c(0,22)))
with(subset(Diff_DDA_output, padjust<0.05 & log2FC>log2(1.3)), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='red'))
with(subset(Diff_DDA_output, padjust<0.05 & log2FC<(-log2(1.3))), points(log2FC, -log10(padjust), pch=20, cex=2.2, col='dodgerblue'))
axis(2, las=2, cex.axis=2.2)
abline(v=c(log2(1.3),-log2(1.3)), col=c("darkgray","darkgray"), lty=c(2,2))
abline(h=-log10(0.05), col="darkgray", lty=2)
dev.off()







