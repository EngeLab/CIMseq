source("/home/martin/pancreas/RNAseq_2/RNAseq_functs.R") # Substitute path to the right one...

# Don't run the following lines, skip until line 20
counts <- load.count.data("/home/martin/correlative/GP5d/CIRM_master_expression_count_table.tsv", mincount=1e5)[[1]]
counts.ercc <- load.count.data("/home/martin/correlative/GP5d/CIRM_master_expression_count_table.tsv", mincount=1e5)[[2]]

o <- grep('10001029', colnames(counts)) # Experiment with doublet and singlet sorts.
counts <- counts[,o]
counts.ercc <- counts.ercc[,o]
counts.log <- norm.log.counts(counts)

good.cells <- counts.log['ACTB',] > get.cutoff.lognorm(counts.log)
table(good.cells)


counts.log <- counts.log[,good.cells]
counts <- counts[,good.cells]
counts.ercc <- counts.ercc[,good.cells]

# Skip until here

# Create annotations from plate names
dbl <- rep("Singlet", length=length(colnames(counts)))
dbl[grepl('1000102901', colnames(counts))] <- "Doublet" # Plate name 1000102901 contains the doublet sorted cells
dbl <- grepl('1000102901', colnames(counts))
frac.ercc <- colSums(counts.ercc) / (colSums(counts.ercc)+colSums(counts))


# Fraction reads that map to ERCC control RNA is lower in doublet cells as expected
sapply(split(frac.ercc, dbl), mean) # Mean is around 30% -> mostly triplets?
ps <- boxplot(split(frac.ercc, dbl), ylab="Fraction ERCC reads")
points(frac.ercc ~ jitter(as.numeric(as.factor(dbl)), factor=0.4), main="Fraction ERCC in singlets/doublets", pch=19, cex=0.6, col="brown")

# Alternative plot to the above, without boxplot
plot(frac.ercc ~ jitter(as.numeric(as.factor(dbl)), factor=1), main="Fraction ERCC in singlets/doublets", pch=19, cex=1, col="black")
dev2bitmap("ERCC_fraction.pdf", type="pdfwrite")


# Check co-expression of some cell type markers. Color doublet cells red. Might be a good angle to investigate these combinations further. For example what subtype of Mesenchymal cell sits next to an endothelial cell.
plot(counts.log['INS',] ~ counts.log['GCG',], col=as.factor(dbl == "Doublet")) # No INS/GCG doublets
plot(counts.log['INS',] ~ counts.log['THY1',], col=as.factor(dbl == "Doublet")) # Only doublets, should be bona fide... These might be an interesting subpopulation to check out
plot(counts.log['GCG',] ~ counts.log['THY1',], col=as.factor(dbl == "Doublet")) # Also interesting
plot(counts.log['INS',] ~ counts.log['NEUROG3',], col=as.factor(dbl == "Doublet")) # Frequently coexpressed in single cells
plot(counts.log['NEUROG3',] ~ counts.log['THY1',], col=as.factor(dbl == "Doublet")) # One doublet neurogenin+ / mesenchymal cells. 
plot(counts.log['FLT1',] ~ counts.log['THY1',], col=as.factor(dbl == "Doublet")) # Endothelial / Mesenchymal doublets
plot(counts.log['EPCAM',] ~ counts.log['THY1',], col=as.factor(dbl == "Doublet")) # Epithelial / Mesenchymal doublets
plot(counts.log['PROM1',] ~ counts.log['THY1',], col=as.factor(dbl == "Doublet")) # Ductal / Mesenchymal doublets
plot(counts.log['PROM1',] ~ counts.log['FLT1',], col=as.factor(dbl == "Doublet")) # Ductal / Endothelial doublets

# Same but with fancy 3d-plots.
library(rgl)
plot3d(counts.log['THY1',], counts.log['FLT1',],counts.log['EPCAM',],  col=c("black", "red")[(dbl == "Doublet")+1], size=10) # Mesenchymal / epithelial / endothelial multicells!
plot3d(counts.log['THY1',], counts.log['FLT1',],counts.log['INS',],  col=c("black", "red")[(dbl == "Doublet")+1], size=10) # Very low INS signal, unclear if there are any Mesenchymal / endothelial /endocrine multicells
col <- colorRampPalette(c("red", "white"))(102)[rank(frac.ercc)/length(frac.ercc)*100+1] 
plot3d(counts.log['THY1',], counts.log['FLT1',],counts.log['INS',],  col=col, size=10) # Colored by fraction ERCC: less ERCC=more red



library(RColorBrewer)
celltype.genes <- c("THY1", "INS", "GCG", "NEUROG3", "PROM1", "FLT1", "SST")
celltype.genes <- c("EPCAM", "THY1")
dev2bitmap("tSNE_doublets_EPCAM_THY1.pdf", type="pdfwrite")


# Plan: use tSNE and hierarchical clustering to determine cell type clusters (several clusters per cell type, actually). Then determine the most likely composition of (sub) cell types for each doublet/multicell. If there are sub-celltypes that preferentially sit next to a certain other cell type, we should be able to see it this way.

# Use the most highly expressed genes for distance calculations. Max instead of mean since we don't want to miss a cell type that is present in only small numbers.
maxs <- order(apply(counts.log, 1, max), decreasing=T)

# First, do a 2d projection of all cells (singlets + doublets). Just for visualization
good.tsne <- run.tsne(as.dist(1-cor(2^counts.log[maxs[1:2000],],method="p")), max_iter=6000, k=2, plot.callback=function(x, confusion = 30) {
    plot.nice(x, counts.log, celltype.genes)
    points(x, col=dbl=="Doublet", cex=1.5)
})

# Restricted to singlets
good.tsne.nodoublets <- run.tsne(as.dist(1-cor(2^counts.log[maxs[1:2000],dbl!="Doublet"],method="p")), max_iter=10000, k=2, perplexity = 10, plot.callback=function(x) {
    plot.nice(x, counts.log[,dbl != "Doublet"], celltype.genes)
#    points(x, col=dbl=="Doublet", cex=1.5)
})

# Use mclust to find (sub) celltypes.
library(mclust)
mod1 <- Mclust(d=dist(good.tsne.nodoublets), G=1:12)
plot(good.tsne.nodoublets, col=jet.colors(max(mod1$classification))[mod1$classification], pch=19)
#identify(good.tsne.nodoublets, labels=mod1$classification)


celltype.genes <- c("EPCAM", "THY1", "FLT1")
plot.nice(good.tsne.nodoublets, counts.log[,dbl!="Doublet"], celltype.genes, cex=1.6) # Most cells are either epithelial, mesenchymal or endothelial.
dev2bitmap("tSNE_singlets_classification_gene_expression.png", type="png16m", res=300)

counts.log.nodoubles <- counts.log[,dbl != "Doublet"]

# Make a matrix of average gene expression by cell types.
cell.types <- as.matrix(as.data.frame(lapply(unique(mod1$classification), function(x) {
    ingroup <- mod1$classification == x
    log2(rowMeans(2^counts.log.nodoubles[,ingroup]))
})))
heatmap(cell.types[maxs[1:2000],], scale = 'none')

# Use optimx to find the combination of cell types that best would approximate the doublet/mutliplet expression data.

fractions <- rep(1.0/(dim(cell.types)[2]), (dim(cell.types)[2]))

# Support function for the optimization
make.synthetic.slice <- function(cell.types, fractions) {
#    fractions <- fractions-min(fractions)
    fractions <- fractions/sum(fractions)
    res <- apply(cell.types, 1, function(x) {sum(x*fractions)})
#    res[is.na(res)] <- 0
}

# Distance function for the optimization - could possibly be better.
dist.to.slice <- function(fractions, cell.types, slice) {
    a <- make.synthetic.slice(cell.types, fractions)
    if(any(is.na(a))) {
        cat("NA in make synthetic slice!\n")
        cat(paste(fractions, sep="\t"), "\n")
    }
    cost <- sum(abs(a - slice))
    cost
}

# Expression matrix with only doublets (so we don't have to subset all the time)
counts.log.onlydoubles <- counts.log[,dbl =="Doublet"]
library(optimx)

cell.types.top2000 <- 2^cell.types[maxs[1:2000],]
counts.top2000 <- 2^counts.log.onlydoubles[maxs[1:2000],]

optim.res.all.2 <- lapply(1:(dim(counts.log.onlydoubles)[2]), function(i) {
    optimx(par=fractions, fn=dist.to.slice, gr=NULL, cell.types=cell.types.top2000, slice=counts.top2000[,i], method=c("L-BFGS-B"), lower=0.0, upper=1.0)
})


groups.mat <- as.matrix(data.frame(lapply(optim.res.all.2, function(x) {unlist(x[1:12])})))
heatmap(t(groups.mat), col=colorRampPalette(c("white", "red"))(100), ColSideColors=jet.colors(max(mod1$classification))[unique(mod1$classification)])
#dev2bitmap("heatmap_classification_by_fractions.pdf", type="pdfwrite")
#dev2bitmap("heatmap_classification_by_fractions.png", type="png16m", res=300)

colnames(groups.mat) <- colnames(counts.log.onlydoubles)
heatmap(t(counts.log.onlydoubles[maxs[1:2000],]), col=colorRampPalette(c("white", "red"))(100), scale='none') # Compare with clustering by raw expression levels.

# This bit is -NOT- working well. Might be a bad route to take.
majordiffs <- function(doublet.vec, cluster.vecs.mat, weights) {
    synslice <- make.synthetic.slice(2^cell.types, weights)
    (synslice - 2^doublet.vec)*2/(synslice + 2^doublet.vec)
}

my.resid1 <- majordiffs(counts.log.onlydoubles[,'X1000102901.G12'],  cell.types, groups.mat[,'X1000102901.G12'])
my.resid2 <- majordiffs(counts.log.onlydoubles[,'X1000102901.D9'],  cell.types, groups.mat[,'X1000102901.D9'])

my.resid3 <- majordiffs(counts.log.onlydoubles[,'X1000102901.E11'],  cell.types, groups.mat[,'X1000102901.E11'])
plot(my.resid1 ~ my.resid2)
points(log2(my.resid1) ~ log2(my.resid3), col="red")


          
# Find best 4 cells by brute force. Not good, noisy singl-cell expression leads to overfitting

dist.to.cell <- function(mixture, cell.types) {
    synth.cell <- rowSums(cell.types)
    synth.cell <- synth.cell/sum(synth.cell)*1022673
    cost <- sum((synth.cell - mixture)^2)
    cost
}

doublet.counts <- 2^counts.log[,dbl=="Doublet"]
singlet.counts <- 2^counts.log[,dbl=="Singlet"]


o <- c(1, 40, 35, 50)
dist.to.cell(doublet.counts[,1], singlet.counts[,o])
dsts <- lapply(1:100000, function(x) {
#    rcells <- sample(1:(dim(singlet.counts)[2]), sample(2:4, 1))
    rcells <- sample(1:(dim(singlet.counts)[2]), 2)
    return(list(rcells, dist.to.cell(doublet.counts[,1], singlet.counts[,rcells])))
})
o <- order(sapply(dsts, function(x) {x[[2]]}))
# Correct answer: THY1 + FLT1:
plot(doublet.counts[celltype.genes,1])
barplot(pp)
barplot(rowSums(singlet.counts[celltype.genes,dsts[[o[[3]]]][[1]]]))
