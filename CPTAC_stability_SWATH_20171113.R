## Classification stability ## 
## Samples stability will be performed on same set of proteins as in the standard analysis (DDA_z, SWATH_z)
## Protein stability will start with all possible proteins (DDA_missing0, SWATH_imputed)




p-value <0.05 -> *

p-value <0.01 -> **

p-value <0.001 -> ***

p-value <0.0001 -> ****


## Step 1: Prepare R session ## 

setwd('Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_Classification/CPTAC_stability/')
library(mclust)
library(ggplot2)

## Step 2: Load data ##

load('../SWATH_classification_workflow_20171101.Rdata')
load('../classification_result.Rdata')

SWATH_class <- classification_result[,7]
names(SWATH_class) <- classification_result[,2]
SWATH_class <- SWATH_class[order(names(SWATH_class))]


## Step 3: Stability analysis on samples

#90%
ARI_90 <- NULL

for(i in 1:100){
set.seed(i)
tmp_take <- sample(1:103, 93, replace=F)
tmp_sample <- SWATH_z[, tmp_take]
tmp_original <- SWATH_class[tmp_take]
set.seed(0)
tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
tmp_class <- tmp_cluster$classification
tmp_ARI <- adjustedRandIndex(tmp_class, tmp_original)
ARI_90 <- c(ARI_90, tmp_ARI)
}

#80%
ARI_80 <- NULL

for(i in 1:100){
set.seed(i)
tmp_take <- sample(1:103, 83, replace=F)
tmp_sample <- SWATH_z[, tmp_take]
tmp_original <- SWATH_class[tmp_take]
set.seed(0)
tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
tmp_class <- tmp_cluster$classification
tmp_ARI <- adjustedRandIndex(tmp_class, tmp_original)
ARI_80 <- c(ARI_80, tmp_ARI)
}

#70%
ARI_70 <- NULL

for(i in 1:100){
set.seed(i)
tmp_take <- sample(1:103, 73, replace=F)
tmp_sample <- SWATH_z[, tmp_take]
tmp_original <- SWATH_class[tmp_take]
set.seed(0)
tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
tmp_class <- tmp_cluster$classification
tmp_ARI <- adjustedRandIndex(tmp_class, tmp_original)
ARI_70 <- c(ARI_70, tmp_ARI)
}

#60%
ARI_60 <- NULL

for(i in 1:100){
set.seed(i)
tmp_take <- sample(1:103, 62, replace=F)
tmp_sample <- SWATH_z[, tmp_take]
tmp_original <- SWATH_class[tmp_take]
set.seed(0)
tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
tmp_class <- tmp_cluster$classification
tmp_ARI <- adjustedRandIndex(tmp_class, tmp_original)
ARI_60 <- c(ARI_60, tmp_ARI)
}

#50%
ARI_50 <- NULL

for(i in 1:100){
set.seed(i)
tmp_take <- sample(1:103, 52, replace=F)
tmp_sample <- SWATH_z[, tmp_take]
tmp_original <- SWATH_class[tmp_take]
set.seed(0)
tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
tmp_class <- tmp_cluster$classification
tmp_ARI <- adjustedRandIndex(tmp_class, tmp_original)
ARI_50 <- c(ARI_50, tmp_ARI)
}



## Step 4: Stability analysis on proteins

# Prepare z scores of proteins without missing values

SWATH_z_all <- SWATH_imputed
SWATH_z_all <- sweep(SWATH_z_all, 1, apply(SWATH_z_all, 1, mean, na.rm=T), "-")
SWATH_z_all <- sweep(SWATH_z_all, 1, apply(SWATH_z_all, 1, sd, na.rm=T), "/")
dim(SWATH_z_all)


#90%
ARI_90protein <- NULL

for (i in 1:100){
set.seed(i)
tmp_take <- sample(1:1659, 1494, replace=F)
tmp_sample <- SWATH_z_all[tmp_take,]
set.seed(0)
tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
tmp_class <- tmp_cluster$classification
tmp_ARI <- adjustedRandIndex(tmp_class, SWATH_class)
ARI_90protein <- c(ARI_90protein, tmp_ARI)
}

#80%
ARI_80protein <- NULL

for (i in 1:100){
set.seed(i)
tmp_take <- sample(1:1659, 1328, replace=F)
tmp_sample <- SWATH_z_all[tmp_take,]
set.seed(0)
tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
tmp_class <- tmp_cluster$classification
tmp_ARI <- adjustedRandIndex(tmp_class, SWATH_class)
ARI_80protein <- c(ARI_80protein, tmp_ARI)
}

#70%
ARI_70protein <- NULL

for (i in 1:100){
set.seed(i)
tmp_take <- sample(1:1659, 1162, replace=F)
tmp_sample <- SWATH_z_all[tmp_take,]
set.seed(0)
tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
tmp_class <- tmp_cluster$classification
tmp_ARI <- adjustedRandIndex(tmp_class, SWATH_class)
ARI_70protein <- c(ARI_70protein, tmp_ARI)
}

#60%
ARI_60protein <- NULL

for (i in 1:100){
set.seed(i)
tmp_take <- sample(1:1659, 996, replace=F)
tmp_sample <- SWATH_z_all[tmp_take,]
set.seed(0)
tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
tmp_class <- tmp_cluster$classification
tmp_ARI <- adjustedRandIndex(tmp_class, SWATH_class)
ARI_60protein <- c(ARI_60protein, tmp_ARI)
}

#50%
ARI_50protein <- NULL

for (i in 1:100){
set.seed(i)
tmp_take <- sample(1:1659, 830, replace=F)
tmp_sample <- SWATH_z_all[tmp_take,]
set.seed(0)
tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
tmp_class <- tmp_cluster$classification
tmp_ARI <- adjustedRandIndex(tmp_class, SWATH_class)
ARI_50protein <- c(ARI_50protein, tmp_ARI)
}


## Step 5: Plot 
plot_SWATH_samples <- data.frame(ARI=c(ARI_90, ARI_80, ARI_70, ARI_60, ARI_50), Percentage=c(rep("90%",100), rep("80%", 100), rep("70%", 100), rep("60%", 100), rep("50%",100)))
plot_SWATH_samples$Percentage <- factor(plot_SWATH_samples$Percentage, levels=rev(levels(plot_SWATH_samples$Percentage)), ordered=T)

t.test(ARI_90, ARI_50)$p.value
[1] 1.260265e-19
t.test(ARI_90, ARI_80)$p.value
[1] 0.003726007
t.test(ARI_80, ARI_70)$p.value
[1] 0.009751242
t.test(ARI_70, ARI_60)$p.value
[1] 0.009068426
t.test(ARI_60, ARI_50)$p.value
[1] 0.08994563


d <- ggplot(plot_SWATH_samples, aes(Percentage, ARI))
sig1 <- data.frame(a=c(1,1:5,5), b=c(1.19, 1.21, 1.21, 1.21, 1.21,1.21,1.19))
sig2 <- data.frame(a=c(1,1:2,2), b=c(1.13, 1.15, 1.15, 1.13))
sig3 <- data.frame(a=c(2,2:3,3), b=c(1.09, 1.11, 1.11, 1.09))
sig4 <- data.frame(a=c(3,3:4,4), b=c(1.05, 1.07, 1.07, 1.05))
sig5 <- data.frame(a=c(4,4:5,5), b=c(1.01, 1.03, 1.03, 1.01))

pdf('boxplot_SWATH_stability_samples_20180207.pdf', height=10, width=10)
d+geom_boxplot()+theme_bw()+scale_y_continuous(limits=c(0,1.25))+
labs(title="Sub-sampling of 103 samples of the SWATH analysis")+
geom_line(data=sig1, aes(x=a, y=b))+annotate("text", x=3, y=1.22, label="****", size=20)+
geom_line(data=sig2, aes(x=a, y=b))+annotate("text", x=1.5, y=1.16, label="**", size=20)+
geom_line(data=sig3, aes(x=a, y=b))+annotate("text", x=2.5, y=1.12, label="**", size=20)+
geom_line(data=sig4, aes(x=a, y=b))+annotate("text", x=3.5, y=1.08, label="**", size=20)+
geom_line(data=sig5, aes(x=a, y=b))+annotate("text", x=4.5, y=1.10, label="n.s.", size=17)+
theme(axis.text.x=element_text(colour='black', size=35), axis.text.y=element_text(colour='black', size=35), axis.title=element_text(size=0), plot.title=element_text(size=0))
dev.off()


plot_SWATH_proteins <- data.frame(ARI=c(ARI_90protein, ARI_80protein, ARI_70protein, ARI_60protein, ARI_50protein), Percentage=c(rep("90%",100), rep("80%", 100), rep("70%", 100), rep("60%", 100), rep("50%",100)))
plot_SWATH_proteins$Percentage <- factor(plot_SWATH_proteins$Percentage, levels=rev(levels(plot_SWATH_proteins$Percentage)), ordered=T)

t.test(ARI_90protein, ARI_50protein)$p.value
[1] 0.3030113
t.test(ARI_90protein, ARI_80protein)$p.value
[1] 0.5972161
t.test(ARI_80protein, ARI_70protein)$p.value
[1] 0.5887578
t.test(ARI_70protein, ARI_60protein)$p.value
[1] 0.3167277
t.test(ARI_60protein, ARI_50protein)$p.value
[1] 0.8821407

e <- ggplot(plot_SWATH_proteins, aes(Percentage, ARI))
sig6 <- data.frame(a=c(1,1:2,2,2:3,3,3:4,4,4:5,5), b=c(0.97, 0.99, 0.99, 0.97, 0.99, 0.99, 0.97, 0.99, 0.99, 0.97, 0.99, 0.99, 0.97))

pdf('boxplot_SWATH_stability_proteins_20180207.pdf', height=10, width=10) 
e+geom_boxplot()+theme_bw()+scale_y_continuous(limits=c(0,1.25))+
labs(title="Sub-sampling of 1659 proteins of the SWATH analysis")+
geom_line(data=sig6, aes(x=a, y=b))+annotate("text", x=3, y=1.06, label="n.s.", size=17)+
theme(axis.text.x=element_text(colour='black', size=35), axis.text.y=element_text(colour='black', size=35), axis.title=element_text(size=0), plot.title=element_text(size=0))
dev.off()


save(ARI_90, ARI_80, ARI_70, ARI_60, ARI_50, plot_SWATH_samples, file='CPTAC_SWATH_stability_samples_20171113.Rdata')
save(ARI_50protein, ARI_60protein, ARI_70protein, ARI_80protein, ARI_90protein, plot_SWATH_proteins, SWATH_z_all, SWATH_class, file='CPTAC_SWATH_stability_proteins_20171113.Rdata')


