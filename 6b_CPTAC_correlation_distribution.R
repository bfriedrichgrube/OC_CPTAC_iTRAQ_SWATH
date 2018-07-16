## Calculate correlation of iTRAQ vs. SWATH data and show distribution of values ###

## Step 1: Prepare R session 

setwd('')
load('CPTAC_DDA_protein_level_CDAP.Rdata')
load('CPTAC_SWATH_imputed.Rdata')

## Step 2: Get only common proteins
gene.ids <- intersect(rownames(SWATH_imputed), rownames(DDA_missing0))
length(gene.ids)

dim(SWATH_imputed)
SWATH_overlap <- SWATH_imputed[gene.ids, ]
DDA_overlap <- DDA_missing0[gene.ids,]
dim(SWATH_overlap)
dim(DDA_overlap)

# sort SWATH in the same order as iTRAQ DDA
SWATH_overlap <- SWATH_overlap[,c(1,66,74,80,23,27,92,97,51,2,6:7,60,8,61:65,67:69,9:10,70:73,11,75:76,12:14,77,15,78:79,16:17,81:82,18:20,83,21,
22,84:90,24,91,25:26,28:41,93,42,94,43,95:96,98:99,44:50,100:101,52,102,53:58,103,59,3:5)]
class(SWATH_overlap)
class(DDA_overlap)
DDA_overlap <- as.matrix(DDA_overlap)

save(SWATH_overlap, file='CPTAC_SWATH_overlap.Rdata')
save(DDA_overlap, file='CPTAC_DDA_overlap.Rdata')

## Step 3: Calculate corrleation

swath.dda.overlap.corr <- sapply(seq_len(nrow(DDA_overlap)), function(i)
cor.test(DDA_overlap[i,],SWATH_overlap[i,], use="pairwise.complete.obs", method="spearman"))
padjust.prot <- p.adjust(swath.dda.overlap.corr[3,], method='fdr')
swath.dda.overlap.corr<- rbind(swath.dda.overlap.corr, padjust.prot)
mean(as.numeric(swath.dda.overlap.corr[4,]))
[1] 0.5926679
median(as.numeric(swath.dda.overlap.corr[4,]))
[1] 0.6125621
                                 
#check for positively correlated proteins
positive.corr.overlap.prot <- swath.dda.overlap.corr[,which(as.numeric(swath.dda.overlap.corr[4,])>0)]
dim(positive.corr.overlap.prot)
[1]    9 1594
1594/1599
[1] 0.996873

#check how many positive correlated ones are also significant
length(which(positive.corr.overlap.prot[9,]<0.01))
[1] 1541
1541/1599
[1] 0.9637273


## Step 4: Plot (Figure 2B)
a <- unlist(swath.dda.overlap.corr[4,])
head(a)
b <- a[order(a)]
head(b)
plot(b)
which(b==median(b))
length(b)
[1] 1599
max(b)
which(b==max(b))
                         
save(b, file="SWATH_iTRAQ_correlations.Rdata")

pdf('SWATH_iTRAQ_correlation.pdf', width=10, height=10)
par(cex.axis=2.5, las=1)
plot(b, xlab='', ylab='', xaxt='n')
points(x=(which(b==median(b))), y=median(b), col='red', pch=16, cex=3)
text(x=400, y=0.7, 'median=0.61', cex=3, col='red')
axis(1, padj=0.3)
dev.off()

## Step 5: Plot protein values distibution (Figure 2A)

pdf('Protein_abundances_distribution.pdf', height=10, width=10)
par(cex.axis=2.5, las=1)
vioplot(as.numeric(DDA_overlap), as.numeric(SWATH_overlap), col='grey', names=c('', ''))
dev.off()


