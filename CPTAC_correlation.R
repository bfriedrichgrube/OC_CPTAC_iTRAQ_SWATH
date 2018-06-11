## calculate correlation of iTRAQ vs. SWATH data and show distribution of values ###

## Step 1: Prepare R session 

setwd('Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/)
setwd('Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/')
load('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_Data/CPTAC_processed_data/CPTAC_protein_level_filtered_20171101.Rdata') 
load('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_Data/CPTAC_processed_data/CPTAC_DDA_protein_level_CDAP.Rdata')
load('CPTAC_Classification/SWATH_classification_workflow_20171101.Rdata')
ls()
[1] "CPTAC_SWATH_cluster_final" "DDA"                      
[3] "SWATH_filtered"            "SWATH_imputed"            
[5] "SWATH_overlap"             "SWATH_z"                  

missing <- NULL
for(i in 1:nrow(DDA)){tmp <-length(which(is.na(DDA[i,])))
missing <- c(missing, tmp)
}

DDA_missing0 <- DDA[(missing==0),] 
dim(DDA_missing0)
[1] 4366  103

gene.ids <- intersect(rownames(SWATH_imputed), rownames(DDA_missing0))
length(gene.ids)
[1] 1599

dim(SWATH_imputed)
[1] 1659  103
SWATH_overlap <- SWATH_imputed[gene.ids, ]
DDA_overlap <- DDA_missing0[gene.ids,]
dim(SWATH_overlap)
[1] 1599  103
dim(DDA_overlap)
[1] 1599  103
SWATH_overlap <- SWATH_overlap[,c(1,66,74,80,23,27,92,97,51,2,6:7,60,8,61:65,67:69,9:10,70:73,11,75:76,12:14,77,15,78:79,16:17,81:82,18:20,83,21,
22,84:90,24,91,25:26,28:41,93,42,94,43,95:96,98:99,44:50,100:101,52,102,53:58,103,59,3:5)]
class(SWATH_overlap)
[1] "matrix"
class(DDA_overlap)
[1] "data.frame"
DDA_overlap <- as.matrix(DDA_overlap)

swath.dda.overlap.corr <- sapply(seq_len(nrow(DDA_overlap)), function(i)
cor.test(DDA_overlap[i,],SWATH_overlap[i,], use="pairwise.complete.obs", method="spearman"))


padjust.prot <- p.adjust(swath.dda.overlap.corr[3,], method='fdr')
swath.dda.overlap.corr<- rbind(swath.dda.overlap.corr, padjust.prot)
mean(as.numeric(swath.dda.overlap.corr[4,]))
[1] 0.5926679
positive.corr.overlap.prot <- swath.dda.overlap.corr[,which(as.numeric(swath.dda.overlap.corr[4,])>0)]
dim(positive.corr.overlap.prot)
[1]    9 1594
1594/1599
[1] 0.996873

length(which(positive.corr.overlap.prot[9,]<0.01))
[1] 1541
1541/1599
[1] 0.9637273

a <- unlist(swath.dda.overlap.corr[4,])
head(a)
PDS5B.rho SRP72.rho EIF3M.rho  MTOR.rho FARP1.rho LRRC1.rho 
0.6896169 0.3394434 0.5231736 0.2793569 0.4783091 0.6038308 

b <- a[order(a)]
head(b)
  OSBPL8.rho    HERC2.rho  RABGGTA.rho     MDN1.rho HLA-DRB1.rho  SPECC1L.rho 
-0.036627422 -0.030389228 -0.029236041 -0.019801871 -0.019648113  0.000582085 
plot(b)
mean(b)
[1] 0.5926679
median(b)
[1] 0.6125621

which(b==median(b))
ATP2B4.rho 
       800 
length(b)
[1] 1599

max(b)
[1] 0.9542569
which(b==max(b))
FN1.rho 
   1599 


pdf('SWATH_iTRAQ_correlation_20180226.pdf', width=10, height=10)
par(cex.axis=2.5, las=1)
plot(b, xlab='', ylab='', xaxt='n')
points(x=(which(b==median(b))), y=median(b), col='red', pch=16, cex=3)
text(x=400, y=0.7, 'median=0.61', cex=3, col='red')
axis(1, padj=0.3)
dev.off()

##adjust colors 
pdf('Protein_abundances_distribution.pdf', height=10, width=10)
par(cex.axis=2.5, las=1)
vioplot(as.numeric(DDA_overlap), as.numeric(SWATH_overlap), col='grey', names=c('', ''))
dev.off()


