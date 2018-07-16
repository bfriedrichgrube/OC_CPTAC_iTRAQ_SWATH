## Network smoothing algorithm 

## Step 1: Prepare R session 

library(BioNetSmooth)
setwd('Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_HRD_networksmoothing/')


## Step 2: Load and prepare data 

network <- read.table('Downloads/9606.protein.links.detailed.v10.5.txt.gz', quote='\"', comment.char="", stringsAsFactors=F, header=T)
network[,1] <- as.character(substr(network[,1],6,20))
network[,2] <- as.character(substr(network[,2],6,20))
head(network[,1:2])
         protein1        protein2
1 ENSP00000000233 ENSP00000263431
2 ENSP00000000233 ENSP00000353863
3 ENSP00000000233 ENSP00000342026
4 ENSP00000000233 ENSP00000240874
5 ENSP00000000233 ENSP00000379847
6 ENSP00000000233 ENSP00000400088

network <- network[network$experimental>=800,]
dim(network)
[1] 31350    10

expr_mat <- read.table('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/CPTAC_mapDIA_analysis/CPTAC_mapDIA302/CPTAC_103HRD_newlib_jumbo_min2max3_quantnorm_batchcorr/analysis_output.txt', sep='\t', dec='.', header=T, as.is=T)
dim(expr_mat)
[1] 2724   10

mapping_table <- read.delim('M20180308AAFB7E4D2F1D05654627429E83DA5CCEF20450Z.tab', sep='\t', dec='.', header=T, as.is=T)
head(mapping_table)
    From              To
1 A0AV96 ENSP00000295971
2 A0AV96 ENSP00000371212
3 A0AV96 ENSP00000371214
4 A0AVT1 ENSP00000313454
5 A0AVT1 ENSP00000399234
6 A0FGR8 ENSP00000251527

mapping_table <- mapping_table[!duplicated(mapping_table$From),]
dim(mapping_table)
[1] 2717    2

#remove unmapped proteins
expr_mat <- expr_mat[grep('A2VCL2', expr_mat[,1], invert=T),]
expr_mat <- expr_mat[grep('P01834', expr_mat[,1], invert=T),]
expr_mat <- expr_mat[grep('P01860', expr_mat[,1], invert=T),]
expr_mat <- expr_mat[grep('P58107', expr_mat[,1], invert=T),]
expr_mat <- expr_mat[grep('P80303', expr_mat[,1], invert=T),]
expr_mat <- expr_mat[grep('Q58FF7', expr_mat[,1], invert=T),]
expr_mat <- expr_mat[grep('Q9UN81', expr_mat[,1], invert=T),]

dim(expr_mat)
[1] 2717   10

expr_mat <- data.frame(as.character(mapping_table[,2]), as.numeric(expr_mat[,6]), as.numeric(expr_mat[,6]))
colnames(expr_mat)<- c('gene_id', 'log2FC_1', 'log2FC_2')

## Step 3: Network mapping

netmap <- network_mapping(network, expr_mat, type = "laplacian", merge.by = "gene_id", global = TRUE)
Network with: 4424 nodes

## Step 4: Network smoothing


netmap$smoothmat <- network_smoothing(net = netmap$G, mat_intensities = netmap$mat_intensities, conditions = colnames (netmap$mat_intensities), iter = 20, alpha = 0.5, network_type = "laplacian")
length(intersect(network[,1],expr_mat[,1]))
[1] 1000

summary(netmap)
                Length   Class     Mode     
G               19571776 dgCMatrix S4       
mat_intensities     8848 -none-    numeric  
gene_names          4424 -none-    character
smoothmat           8848 -none-    numeric  

##Step 5: Extract top scored proteins
#5% top positive and 5% top negative proteins
# 5% of 4424 is 221

a <- netmap$smoothmat[,1]
names(a)<- netmap$gene_names

b<- a[order(a)]
x <- head(b, 221)
b <- a[order(a, decreasing=T)]
y <- head(b, 221)

write.table(names(y), file='positive_prots.txt', quote=F, row.names=F)
write.table(names(x), file='negative_prots.txt', quote=F, row.names=F)

# insert positive correlated in STRING (14.03.2018) 
## GO enrichment DNA repair and chromosome organization among top 50 GOTerms
## 

save(netmap, file='smoothed_network_20180314.Rdata')


## trial not promising
#top 10% absolute scores =442 proteins 

d <- abs(a)
e <- d[order(d, decreasing=T)]
z <- head(e, 442)
write.table(names(z), file='absolute_topprots.txt', quote=F, row.names=F)
