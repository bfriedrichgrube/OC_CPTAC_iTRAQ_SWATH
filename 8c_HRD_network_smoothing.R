## Network smoothing algorithm (Figure 7B)

## Step 1: Prepare R session 

library(BioNetSmooth)
setwd('')

## Step 2: Load and prepare data 

network <- read.table('9606.protein.links.detailed.v10.5.txt.gz', quote='\"', comment.char="", stringsAsFactors=F, header=T)
network[,1] <- as.character(substr(network[,1],6,20))
network[,2] <- as.character(substr(network[,2],6,20))

network <- network[network$experimental>=800,]

expr_mat <- read.table('analysis_output.txt', sep='\t', dec='.', header=T, as.is=T)

mapping_table <- read.delim('M20180308AAFB7E4D2F1D05654627429E83DA5CCEF20450Z.tab', sep='\t', dec='.', header=T, as.is=T)
mapping_table <- mapping_table[!duplicated(mapping_table$From),]

#remove unmapped proteins
expr_mat <- expr_mat[grep('A2VCL2', expr_mat[,1], invert=T),]
expr_mat <- expr_mat[grep('P01834', expr_mat[,1], invert=T),]
expr_mat <- expr_mat[grep('P01860', expr_mat[,1], invert=T),]
expr_mat <- expr_mat[grep('P58107', expr_mat[,1], invert=T),]
expr_mat <- expr_mat[grep('P80303', expr_mat[,1], invert=T),]
expr_mat <- expr_mat[grep('Q58FF7', expr_mat[,1], invert=T),]
expr_mat <- expr_mat[grep('Q9UN81', expr_mat[,1], invert=T),]

expr_mat <- data.frame(as.character(mapping_table[,2]), as.numeric(expr_mat[,6]), as.numeric(expr_mat[,6]))
colnames(expr_mat)<- c('gene_id', 'log2FC_1', 'log2FC_2')

## Step 3: Network mapping

netmap <- network_mapping(network, expr_mat, type = "laplacian", merge.by = "gene_id", global = TRUE)
Network with: 4424 nodes

## Step 4: Network smoothing

netmap$smoothmat <- network_smoothing(net = netmap$G, mat_intensities = netmap$mat_intensities, conditions = colnames (netmap$mat_intensities), iter = 20, alpha = 0.5, network_type = "laplacian")
length(intersect(network[,1],expr_mat[,1]))

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

save(netmap, file='smoothed_network.Rdata')
