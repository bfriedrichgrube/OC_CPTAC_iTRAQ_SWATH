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


## Step 3: Stability analysis on samples of SWATH

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


## Step 4: Stability analysis on proteins of SWATH

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

## Step 5: Plot for SWATH
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

pdf('boxplot_SWATH_stability_samples.pdf', height=10, width=10)
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

pdf('boxplot_SWATH_stability_proteins.pdf', height=10, width=10) 
e+geom_boxplot()+theme_bw()+scale_y_continuous(limits=c(0,1.25))+
labs(title="Sub-sampling of 1659 proteins of the SWATH analysis")+
geom_line(data=sig6, aes(x=a, y=b))+annotate("text", x=3, y=1.06, label="n.s.", size=17)+
theme(axis.text.x=element_text(colour='black', size=35), axis.text.y=element_text(colour='black', size=35), axis.title=element_text(size=0), plot.title=element_text(size=0))
dev.off()


save(ARI_90, ARI_80, ARI_70, ARI_60, ARI_50, plot_SWATH_samples, file='CPTAC_SWATH_stability_samples_20171113.Rdata')
save(ARI_50protein, ARI_60protein, ARI_70protein, ARI_80protein, ARI_90protein, plot_SWATH_proteins, SWATH_z_all, SWATH_class, file='CPTAC_SWATH_stability_proteins_20171113.Rdata')

## Step 2: Load data ##

load('../DDA_classification_workflow_20171101.Rdata')
# load('../SWATH_classification_workflow_20171101.Rdata')
load('../classification_result.Rdata')

DDA_class <- classification_result[,6]
names(DDA_class) <- rownames(classification_result)
colnames(DDA_z) <- rownames(classification_result)

## Step 3: Stability analysis on samples

#90% 
ARI_90 <- NULL

for(i in 1:100){
set.seed(i)
tmp_take <- sample(1:103, 93, replace=F)
tmp_sample <- DDA_z[, tmp_take]
tmp_original <- DDA_class[tmp_take]
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
tmp_sample <- DDA_z[, tmp_take]
tmp_original <- DDA_class[tmp_take]
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
tmp_sample <- DDA_z[, tmp_take]
tmp_original <- DDA_class[tmp_take]
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
tmp_sample <- DDA_z[, tmp_take]
tmp_original <- DDA_class[tmp_take]
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
tmp_sample <- DDA_z[, tmp_take]
tmp_original <- DDA_class[tmp_take]
set.seed(0)
tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
tmp_class <- tmp_cluster$classification
tmp_ARI <- adjustedRandIndex(tmp_class, tmp_original)
ARI_50 <- c(ARI_50, tmp_ARI)
}

## Step 4: Stability analysis on proteins

# Prepare z scores of proteins without missing values
DDA_z_all <- as.matrix(DDA_missing0[4:4366,])
DDA_z_all <- sweep(DDA_z_all, 1, apply(DDA_z_all, 1, mean, na.rm=T), "-")
DDA_z_all <- sweep(DDA_z_all, 1, apply(DDA_z_all, 1, sd, na.rm=T), "/")

#90%
ARI_90protein <- NULL

for (i in 1:100){
set.seed(i)
tmp_take <- sample(1:4363, 3927, replace=F)
tmp_sample <- DDA_z_all[tmp_take,]
set.seed(0)
tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
tmp_class <- tmp_cluster$classification
tmp_ARI <- adjustedRandIndex(tmp_class, DDA_class)
ARI_90protein <- c(ARI_90protein, tmp_ARI)
}

#80%
ARI_80protein <- NULL

for (i in 1:100){
set.seed(i)
tmp_take <- sample(1:4363, 3491, replace=F)
tmp_sample <- DDA_z_all[tmp_take,]
set.seed(0)
tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
tmp_class <- tmp_cluster$classification
tmp_ARI <- adjustedRandIndex(tmp_class, DDA_class)
ARI_80protein <- c(ARI_80protein, tmp_ARI)
}

#70%
ARI_70protein <- NULL

for (i in 1:100){
set.seed(i)
tmp_take <- sample(1:4363, 3055, replace=F)
tmp_sample <- DDA_z_all[tmp_take,]
set.seed(0)
tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
tmp_class <- tmp_cluster$classification
tmp_ARI <- adjustedRandIndex(tmp_class, DDA_class)
ARI_70protein <- c(ARI_70protein, tmp_ARI)
}

#60%
ARI_60protein <- NULL

for (i in 1:100){
set.seed(i)
tmp_take <- sample(1:4363, 2618, replace=F)
tmp_sample <- DDA_z_all[tmp_take,]
set.seed(0)
tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
tmp_class <- tmp_cluster$classification
tmp_ARI <- adjustedRandIndex(tmp_class, DDA_class)
ARI_60protein <- c(ARI_60protein, tmp_ARI)
}

#50%
ARI_50protein <- NULL

for (i in 1:100){
set.seed(i)
tmp_take <- sample(1:4363, 2182, replace=F)
tmp_sample <- DDA_z_all[tmp_take,]
set.seed(0)
tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
tmp_class <- tmp_cluster$classification
tmp_ARI <- adjustedRandIndex(tmp_class, DDA_class)
ARI_50protein <- c(ARI_50protein, tmp_ARI)
}

#40%
ARI_40protein <- NULL

for (i in 1:100){
set.seed(i)
tmp_take <- sample(1:4363, 1746, replace=F)
tmp_sample <- DDA_z_all[tmp_take,]
set.seed(0)
tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
tmp_class <- tmp_cluster$classification
tmp_ARI <- adjustedRandIndex(tmp_class, DDA_class)
ARI_40protein <- c(ARI_40protein, tmp_ARI)
}

#30%
ARI_30protein <- NULL

for (i in 1:100){
set.seed(i)
tmp_take <- sample(1:4363, 1309, replace=F)
tmp_sample <- DDA_z_all[tmp_take,]
set.seed(0)
tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
tmp_class <- tmp_cluster$classification
tmp_ARI <- adjustedRandIndex(tmp_class, DDA_class)
ARI_30protein <- c(ARI_30protein, tmp_ARI)
}

## Step 5: Plot

plot_DDA_samples <- data.frame(ARI=c(ARI_90, ARI_80, ARI_70, ARI_60, ARI_50), Percentage=c(rep("90%",100), rep("80%", 100), rep("70%", 100), rep("60%", 100), rep("50%",100)))
plot_DDA_samples$Percentage<- factor(plot_DDA_samples$Percentage, levels=rev(levels(plot_DDA_samples$Percentage)), ordered=T)

t.test(ARI_90, ARI_50)$p.value
[1] 3.025899e-21
t.test(ARI_90, ARI_80)$p.value
[1] 0.02863022
t.test(ARI_80, ARI_70)$p.value
[1] 0.001287279
t.test(ARI_70, ARI_60)$p.value
[1] 0.001064435
t.test(ARI_60, ARI_50)$p.value
[1] 0.01155905

d <- ggplot(plot_DDA_samples, aes(Percentage, ARI))
sig1 <- data.frame(a=c(1,1:5,5), b=c(1.17, 1.19, 1.19, 1.19, 1.19,1.19,1.17))
sig2 <- data.frame(a=c(1,1:2,2), b=c(1.11, 1.13, 1.13, 1.11))
sig3 <- data.frame(a=c(2,2:3,3), b=c(1.05, 1.07, 1.07, 1.05))
sig4 <- data.frame(a=c(3,3:4,4), b=c(0.99, 1.01, 1.01, 0.99))
sig5 <- data.frame(a=c(4,4:5,5), b=c(0.93, 0.95, 0.95, 0.93))

pdf('boxplot_stability_samples.pdf', height=10, width=10)
d+geom_boxplot()+theme_bw()+scale_y_continuous(limits=c(0,1.25))+
labs(title="Sub-sampling of 103 samples of the iTRAQ DDA analysis")+
geom_line(data=sig1, aes(x=a, y=b))+annotate("text", x=3, y=1.21, label="****", size=10)+
geom_line(data=sig2, aes(x=a, y=b))+annotate("text", x=1.5, y=1.15, label="*", size=10)+
geom_line(data=sig3, aes(x=a, y=b))+annotate("text", x=2.5, y=1.09, label="**", size=10)+
geom_line(data=sig4, aes(x=a, y=b))+annotate("text", x=3.5, y=1.03, label="**", size=10)+
geom_line(data=sig5, aes(x=a, y=b))+annotate("text", x=4.5, y=0.97, label="*", size=10)+
theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20), axis.title=element_text(size=22), plot.title=element_text(size=24))
dev.off()

pdf('boxplot_stability_samples_20180206.pdf', height=10, width=10)
d+geom_boxplot()+theme_bw()+scale_y_continuous(limits=c(0,1.25))+
labs(title="Sub-sampling of 103 samples of the iTRAQ DDA analysis")+
geom_line(data=sig1, aes(x=a, y=b))+annotate("text", x=3, y=1.2, label="****", size=20)+
geom_line(data=sig2, aes(x=a, y=b))+annotate("text", x=1.5, y=1.14, label="*", size=20)+
geom_line(data=sig3, aes(x=a, y=b))+annotate("text", x=2.5, y=1.08, label="**", size=20)+
geom_line(data=sig4, aes(x=a, y=b))+annotate("text", x=3.5, y=1.02, label="**", size=20)+
geom_line(data=sig5, aes(x=a, y=b))+annotate("text", x=4.5, y=0.96, label="*", size=20)+
theme(axis.text.x=element_text(colour='black', size=35), axis.text.y=element_text(colour='black', size=35), axis.title=element_text(size=0), plot.title=element_text(size=0))
dev.off()





plot_DDA_proteins <- data.frame(ARI=c(ARI_90protein, ARI_80protein, ARI_70protein, ARI_60protein, ARI_50protein, ARI_40protein, ARI_30protein),
Percentage=c(rep("90%",100), rep("80%", 100), rep("70%", 100), rep("60%", 100), rep("50%",100), rep("40%", 100), rep("30%", 100)))
plot_DDA_proteins$Percentage <- factor(plot_DDA_proteins$Percentage, levels=rev(levels(plot_DDA_proteins$Percentage)), ordered=T)

t.test(ARI_90protein, ARI_30protein)$p.value
[1] 0.711517
t.test(ARI_90protein, ARI_80protein)$p.value
[1] 0.6098863
t.test(ARI_80protein, ARI_70protein)$p.value
[1] 0.8650532
t.test(ARI_70protein, ARI_60protein)$p.value
[1] 0.5708462
t.test(ARI_60protein, ARI_50protein)$p.value
[1] 0.3360966
t.test(ARI_50protein, ARI_40protein)$p.value
[1] 0.2741676
t.test(ARI_40protein, ARI_30protein)$p.value
[1] 0.1474785


e <- ggplot(plot_DDA_proteins, aes(Percentage, ARI))
sig6 <- data.frame(a=c(1,1:2,2,2:3,3,3:4,4,4:5,5,5:6,6,6:7,7), b=c(1.01, 1.03, 1.03, 1.01, 1.03, 1.03, 1.01, 1.03, 1.03, 1.01, 1.03, 1.03, 1.01, 1.03, 1.03, 1.01, 1.03, 1.03, 1.01))

pdf('boxplot_stability_proteins_20180206.pdf', height=10, width=10) 
e+geom_boxplot()+theme_bw()+scale_y_continuous(limits=c(0,1.25))+
labs(title="Sub-sampling of 4363 proteins of the iTRAQ DDA analysis")+
geom_line(data=sig6, aes(x=a, y=b))+annotate("text", x=4, y=1.10, label="n.s.", size=17)+
theme(axis.text.x=element_text(colour='black', size=35), axis.text.y=element_text(colour='black', size=35), axis.title=element_text(size=0), plot.title=element_text(size=0))
dev.off()

save(ARI_90, ARI_80, ARI_70, ARI_60, ARI_50, plot_DDA_samples, file='CPTAC_DDA_stability_samples_20171113.Rdata')
save(ARI_30protein, ARI_40protein, ARI_50protein, ARI_60protein, ARI_70protein, ARI_80protein, ARI_90protein, plot_DDA_proteins, DDA_z_all, DDA_class,
+ file='CPTAC_DDA_stability_proteins_20171113.Rdata')


## Step 2: Load data ##

load('../SWATH_classification_workflow_20171101.Rdata')
load('../classification_result.Rdata')

SWATH_class <- classification_result[,7]
names(SWATH_class) <- classification_result[,2]
SWATH_class <- SWATH_class[order(names(SWATH_class))]


## Step 3: Stability analysis on samples

Fisher_90_mes <- NULL
Fisher_90_other <- NULL

for (i in 1:100){
	set.seed(i)
		tmp_take <- sample(1:103, 93, replace=F)
		tmp_sample <- SWATH_z[,tmp_take]
		tmp_original <- SWATH_class[tmp_take]
	set.seed(0)
		tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
		tmp_class <- tmp_cluster$classification
		tmp_table <- table(tmp_class, tmp_original)
		tmp_fish_all <- NULL
	for(j in 1:nrow(tmp_table)){
		tmp_con <- matrix(data=c(sum(tmp_table[j,1:2]), tmp_table[j,3], sum(tmp_table[-j,1:2]), sum(tmp_table[-j,3]) ), ncol=2, byrow=T)
		tmp_fish <- fisher.test(tmp_con)$p.value
		tmp_fish_all <- c(tmp_fish_all, tmp_fish)		
	}
	Fisher_90_mes <- c(Fisher_90_mes, min(tmp_fish_all))
	Fisher_90_other <- c(Fisher_90_other, tmp_fish_all[-which(tmp_fish_all==min(tmp_fish_all))])
}


Fisher_80_mes <- NULL
Fisher_80_other <- NULL

for (i in 1:100){
	set.seed(i)
		tmp_take <- sample(1:103, 83, replace=F)
		tmp_sample <- SWATH_z[,tmp_take]
		tmp_original <- SWATH_class[tmp_take]
	set.seed(0)
		tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
		tmp_class <- tmp_cluster$classification
		tmp_table <- table(tmp_class, tmp_original)
		tmp_fish_all <- NULL
	for(j in 1:nrow(tmp_table)){
		tmp_con <- matrix(data=c(sum(tmp_table[j,1:2]), tmp_table[j,3], sum(tmp_table[-j,1:2]), sum(tmp_table[-j,3]) ), ncol=2, byrow=T)
		tmp_fish <- fisher.test(tmp_con)$p.value
		tmp_fish_all <- c(tmp_fish_all, tmp_fish)		
	}
	Fisher_80_mes <- c(Fisher_80_mes, min(tmp_fish_all))
	Fisher_80_other <- c(Fisher_80_other, tmp_fish_all[-which(tmp_fish_all==min(tmp_fish_all))])
}



Fisher_70_mes <- NULL
Fisher_70_other <- NULL

for (i in 1:100){
	set.seed(i)
		tmp_take <- sample(1:103, 73, replace=F)
		tmp_sample <- SWATH_z[,tmp_take]
		tmp_original <- SWATH_class[tmp_take]
	set.seed(0)
		tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
		tmp_class <- tmp_cluster$classification
		tmp_table <- table(tmp_class, tmp_original)
		tmp_fish_all <- NULL
	for(j in 1:nrow(tmp_table)){
		tmp_con <- matrix(data=c(sum(tmp_table[j,1:2]), tmp_table[j,3], sum(tmp_table[-j,1:2]), sum(tmp_table[-j,3]) ), ncol=2, byrow=T)
		tmp_fish <- fisher.test(tmp_con)$p.value
		tmp_fish_all <- c(tmp_fish_all, tmp_fish)		
	}
	Fisher_70_mes <- c(Fisher_70_mes, min(tmp_fish_all))
	Fisher_70_other <- c(Fisher_70_other, tmp_fish_all[-which(tmp_fish_all==min(tmp_fish_all))])
}


Fisher_60_mes <- NULL
Fisher_60_other <- NULL

for (i in 1:100){
	set.seed(i)
		tmp_take <- sample(1:103, 62, replace=F)
		tmp_sample <- SWATH_z[,tmp_take]
		tmp_original <- SWATH_class[tmp_take]
	set.seed(0)
		tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
		tmp_class <- tmp_cluster$classification
		tmp_table <- table(tmp_class, tmp_original)
		tmp_fish_all <- NULL
	for(j in 1:nrow(tmp_table)){
		tmp_con <- matrix(data=c(sum(tmp_table[j,1:2]), tmp_table[j,3], sum(tmp_table[-j,1:2]), sum(tmp_table[-j,3]) ), ncol=2, byrow=T)
		tmp_fish <- fisher.test(tmp_con)$p.value
		tmp_fish_all <- c(tmp_fish_all, tmp_fish)		
	}
	Fisher_60_mes <- c(Fisher_60_mes, min(tmp_fish_all))
	Fisher_60_other <- c(Fisher_60_other, tmp_fish_all[-which(tmp_fish_all==min(tmp_fish_all))])
}


Fisher_50_mes <- NULL
Fisher_50_other <- NULL

for (i in 1:100){
	set.seed(i)
		tmp_take <- sample(1:103, 52, replace=F)
		tmp_sample <- SWATH_z[,tmp_take]
		tmp_original <- SWATH_class[tmp_take]
	set.seed(0)
		tmp_cluster <- Mclust(t(tmp_sample), modelNames='VII')
		tmp_class <- tmp_cluster$classification
		tmp_table <- table(tmp_class, tmp_original)
		tmp_fish_all <- NULL
	for(j in 1:nrow(tmp_table)){
		tmp_con <- matrix(data=c(sum(tmp_table[j,1:2]), tmp_table[j,3], sum(tmp_table[-j,1:2]), sum(tmp_table[-j,3]) ), ncol=2, byrow=T)
		tmp_fish <- fisher.test(tmp_con)$p.value
		tmp_fish_all <- c(tmp_fish_all, tmp_fish)		
	}
	Fisher_50_mes <- c(Fisher_50_mes, min(tmp_fish_all))
	Fisher_50_other <- c(Fisher_50_other, tmp_fish_all[-which(tmp_fish_all==min(tmp_fish_all))])
}


plot_SWATH_other <- data.frame(p.value=c(Fisher_90_other_adj, Fisher_80_other_adj, Fisher_70_other_adj, Fisher_60_other_adj, Fisher_50_other_adj),Percentage=c(rep("90%",197), rep("80%", 190), rep("70%", 166), rep("60%", 147), rep("50%",113)))
plot_SWATH_other$Percentage <- factor(plot_SWATH_other$Percentage, levels=rev(levels(plot_SWATH_other$Percentage)), ordered=T)
d <- ggplot(plot_SWATH_other, aes(Percentage, p.value))



pdf('boxplot_SWATH_otherclusters.pdf', height=10, width=10)
d+geom_boxplot()+theme_bw()+scale_y_reverse()+geom_hline(aes(yintercept=0.01), colour="red", linetype="dashed")+labs(title="")
dev.off()


plot_SWATH_mes <- data.frame(p.value=c(Fisher_90_mes_adj, Fisher_80_mes_adj, Fisher_70_mes_adj, Fisher_60_mes_adj, Fisher_50_mes_adj),Percentage=c(rep("90%",100), rep("80%", 100), rep("70%", 100), rep("60%", 100), rep("50%",100)))
plot_SWATH_mes$Percentage <- factor(plot_SWATH_mes$Percentage, levels=rev(levels(plot_SWATH_mes$Percentage)), ordered=T)

e <- ggplot(plot_SWATH_mes, aes(Percentage, p.value))

pdf('boxplot_SWATH_mesenchymalcluster.pdf', height=10, width=10)
e+geom_boxplot()+theme_bw()+scale_y_reverse()+geom_hline(aes(yintercept=0.01), colour="red", linetype="dashed")+labs(title="")
dev.off()




values<-c(median(Fisher_90_mes_adj), median(Fisher_80_mes_adj), median(Fisher_70_mes_adj), median(Fisher_60_mes_adj), median(Fisher_50_mes_adj), 
median(Fisher_90_other_adj), median(Fisher_80_other_adj), median(Fisher_70_other_adj), median(Fisher_60_other_adj), median(Fisher_50_other_adj))
> length(values)
[1] 10
> Percentage <-c('90%', '80%', '70%', '60%', '50%', '90%', '80%', '70%', '60%', '50%')
> length(Percentage)
[1] 10
> Group <- c(rep('Mesenchymal cluster', 5), rep('Other clusters', 5))
> length(Group)
[1] 10
> plot_cluster <- data.frame(values, Percentage, Group)
> head(plot_cluster)
        values Percentage               Group
1 7.044070e-18        90% Mesenchymal cluster
2 6.581833e-14        80% Mesenchymal cluster
3 4.274643e-11        70% Mesenchymal cluster
4 3.267498e-07        60% Mesenchymal cluster
5 1.871104e-02        50% Mesenchymal cluster
6 6.382241e-04        90%      Other clusters
f <- ggplot(plot_cluster, aes(Percentage, values))
f+geom_point(aes(colour=factor(Group)))
f+geom_point(aes(colour=factor(Group)))+scale_y_reverse()
plot_cluster$Percentage <- factor(plot_cluster$Percentage, levels=rev(levels(plot_cluster$Percentage)), ordered=T)
f <- ggplot(plot_cluster, aes(Percentage, values))
f+geom_point(aes(colour=factor(Group)), size=4)+scale_y_reverse()+geom_hline(aes(yintercept=0.01), colour="red", linetype="dashed", size=2)+labs(title="")

pdf('cluster_stability_20180221.pdf', height=10, width=10) > f+geom_point(aes(shape=factor(Group)), size=5)+scale_y_reverse()+geom_hline(aes(yintercept=0.01), colour="red", linetype="dashed", size=2)+labs(title="")+theme_bw()+theme(axis.text.x=element_text(colour='black', size=35), axis.text.y=element_text(colour='black', size=35), axis.title=element_text(size=0), plot.title=element_text(size=0), legend.position='none')
dev.off()


