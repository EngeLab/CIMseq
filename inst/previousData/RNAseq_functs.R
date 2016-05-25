# Example run:
# Load data from counts table
#counts <- load.count.data("CIRM_data/CIRM_master_expression_count_table.tsv")[[1]]
#dim(counts)

# Calculate pairwise distances
#mode.fail.d <- calc.mode.fail.dist(counts)

# Do dimensionality reduction using tSNE. 
##my.tsne <- run.tsne(mode.fail.d, plot.callback=NULL)
#my.tsne <- run.tsne(mode.fail.d, plot.callback=my.plot.callback) # If you have rgl

# Normalize counts for diplay (eg. log2(counts per million+1)
#counts.log <- norm.log.counts(counts)

# Make a snazzy html 3d graph
#make.html(counts.log, c('VIM', 'CD74'), my.tsne, file.suffix="_expression_in_CIRM.html")

#dmat <- as.matrix(mode.fail.d)
#a.o <- match(colnames(dmat), annotation$rmangled.key)
#tissues <- annotation$inferred_tissue[a.o][o]
#tissue.source <- annotation$tissue[a.o][o]
#tissue.type <- annotation$tissue_type[a.o][o]
#experimenter <- annotation$experimenter[a.o]
#chip <- annotation$c1_chip_id[a.o]
#experiment <- annotation$experiment_name[a.o]
#week <- annotation$gestational_week[a.o][o]
#week.cols <- as.integer(week)
##week.cols <- week.cols-min(week.cols))/(max(week.cols)-min(week.cols)
#week.cols <- week.cols-min(week.cols[!is.na(week.cols)])+1
#week.cols[is.na(week.cols)] <- 0
#week.cols <- brewer.pal("Blues", n=max(week.cols+1))[week.cols+1]
#week.cols <- c("black", brewer.pal("PuBuGn", n=max(week.cols)))[week.cols+1]
#week.cols <- c("black", colorRampPalette(c("red", "yellow"))(max(week.cols)))[week.cols+1]
#plot3d(tsne.goodcov, col=week.cols, size=10)

#library(RColorBrewer)
#plot.sequential('VIM', counts.log, my.tsne)

#source.and.type <- paste(tissue.type, tissue.source, sep=".")
#source.and.type[source.and.type=='NA.pancreas'] <- 'fetal.pancreas'
#source.and.type[source.and.type=='NA.NA'] <- 'fetal.liver'
#source.and.type.i <- as.integer(as.factor(source.and.type))
#cols <- colors()[as.integer(source.and.type.i)+1]
#plot3d(tsne.goodcov, col=cols, size=10)

#o <- !is.na(annotation$gene_body_coverage[a.o]) &  annotation$gene_body_coverage[a.o] > 0.6

#tsne.goodcov <- run.tsne(as.dist(as.matrix(mode.fail.d)[o,o]))

gm.mean <- function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / sum(x > 0))
  }
}

jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

norm.log.counts <- function(counts) {
    norm.fact <- colSums(counts)
    counts.norm <- t( apply(counts, 1, function(x) {x/norm.fact*1000000+1}))
    counts.log <- log2(counts.norm)
    counts.log
}

load.count.data <- function(fname, mincount = 4e5, omit.bad.genes=TRUE, omit.bad.cells=TRUE) {
    counts <- read.table(fname, header=T, sep="\t")

    gene.names <- counts[[1]]
    counts <- as.matrix(counts[,-1])
    rownames(counts) <- gene.names

    ercc.counts <- counts[grepl("^ERCC", gene.names),]
    counts <- counts[!grepl("^ERCC", gene.names),]
    counts[is.na(counts)] <- 0
    last3.counts <- counts[(dim(counts)[1]-5):dim(counts)[1],]
    counts <- counts[1:(dim(counts)[1]-5),]


    # omit genes with no counts
    if(omit.bad.genes) {
        counts <- counts[rowSums(counts)>0,];
        
    }
    # omit cells with very poor coverage
    hist(colSums(counts), breaks=100)
    abline(v=mincount, col="red")
    ercc.counts <- ercc.counts[,colSums(counts)>mincount];
    last3.counts <- last3.counts[,colSums(counts)>mincount];
    counts <- counts[,colSums(counts)>mincount];
    
    list(counts=counts, ercc.counts=ercc.counts, last3.counts=last3.counts)
}


calc.mode.fail.dist <- function(gene.counts, n.cores = 8, ...) {
    require("scde")
    require("boot")
    require("parallel")
    o.ifm <- scde.error.models(counts=gene.counts,n.cores=n.cores,threshold.segmentation=T,save.crossfit.plots=F,save.model.plots=F,verbose=1, ...)
    valid.cells <- o.ifm$corr.a >0
    o.prior <- scde.expression.prior(models=o.ifm,counts=gene.counts,length.out=400,show.plot=T)
    o.fpm <- scde.expression.magnitude(o.ifm,counts=gene.counts);

    o.fail.curves <- scde.failure.probability(o.ifm,magnitudes=log((10^o.prior$x)-1))
    par(mfrow=c(1,1),mar = c(3.5,3.5,0.5,0.5), mgp = c(2.0,0.65,0), cex = 1);
    plot(c(),c(),xlim=range(o.prior$x),ylim=c(0,1),xlab="expression magnitude (log10)",ylab="drop-out probability")
    invisible(apply(o.fail.curves,2,function(y) lines(x=o.prior$x,y=y,col="orange")))

    p.self.fail <- scde.failure.probability(models=o.ifm,counts=gene.counts)
    cell.names <- colnames(gene.counts); names(cell.names) <- cell.names;

    # reclculate posteriors with the individual posterior modes 
    jp <- scde.posteriors(models=o.ifm,gene.counts,o.prior,return.individual.posterior.modes=T,n.cores=n.cores)
    # find joint posterior modes for each gene - a measure of MLE of group-average expression
    jp$jp.modes <- log(as.numeric(colnames(jp$jp)))[max.col(jp$jp)]
    p.mode.fail <- scde.failure.probability(models=o.ifm,magnitudes=jp$jp.modes)
    # weight matrix
    matw <- 1-sqrt(p.self.fail*sqrt(p.self.fail*p.mode.fail))
    # magnitude matrix (using individual posterior modes here)
    mat <- log10(exp(jp$modes)+1);
    # weighted distance
    mode.fail.dist <- as.dist(1-do.call(rbind,mclapply(cell.names,function(nam1) {
        unlist(lapply(cell.names,function(nam2) {
            corr(cbind(mat[,nam1],mat[,nam2]),w=sqrt(sqrt(matw[,nam1]*matw[,nam2])))
        }))
    },mc.cores=n.cores)),upper=F);
#    rownames(mode.fail.dist) <- colnames(gene.counts)
#    colnames(mode.fail.dist) <- colnames(gene.counts)    
    mode.fail.dist
}

my.plot.callback <- function(x) {
    require(rgl)
    plot3d(x, size=10, col=cols)
}

run.tsne <- function(my.dist, plot.callback=my.plot.callback, k=3, max_iter=5000, ...) {
    require(tsne)
    my.tsne <- tsne(my.dist, k=k, epoch_callback=plot.callback, initial_dims=50, max_iter=max_iter, ...)
    rownames(my.tsne) <- rownames(my.dist)
    my.tsne
}

plot.sequential <- function(geneset, counts.log, my.tsne, cols=NULL) {
for(gene in geneset) {
    invisible(readline(prompt="Press [enter] to continue"))
    if(!is.na(match(gene, rownames(counts.log)))) {
#        panc.genes.count <- counts[gene,]/colSums(counts)*1000000
#        panc.genes.max <- panc.genes.count
        panc.g.col <- counts.log[gene,]

#        sum(panc.genes.max <= 1)
#        panc.genes.max[panc.genes.max <= 1] <- 1
#        panc.genes.max[panc.genes.max == 0] <- 0.0001
#        panc.g.col <- as.numeric(log2(panc.genes.max))
        panc.g.col <- (panc.g.col-min(panc.g.col))/(max(panc.g.col)-min(panc.g.col))
        if(is.null(cols)) {
            cols <- colorRampPalette(c("Black", "Blue", "Yellow"))(100)
        }
        panc.g.col <- cols[panc.g.col*(length(cols)-2)+1]
#        plot3d(my.rtsne$Y[,1], my.rtsne$Y[,2], my.rtsne$Y[,3], col=panc.g.col, size=10, main=gene)
        plot3d(my.tsne, col=panc.g.col, size=10, main=gene)
    }
    else {
        cat(gene, " not found\n");
    }
}
}

geneset.colors <- function(gene, counts.log) {
    if(!is.na(match(gene[1], rownames(counts.log)))) {        
#        panc.genes.count <- counts[gene,]/colSums(counts)*1000000
#        panc.genes.max <- panc.genes.count
        panc.g.col1 <- counts.log[gene[1],]
        panc.g.col1 <- (panc.g.col1-min(panc.g.col1))/(max(panc.g.col1)-min(panc.g.col1))
        panc.g.col2 <- rep(0, dim(counts.log)[2])
        panc.g.col3 <- rep(0, dim(counts.log)[2])
        if(length(gene) > 1) {
            cat("two colors!")
            panc.g.col2 <- counts.log[gene[2],]
            panc.g.col2 <- (panc.g.col2-min(panc.g.col2))/(max(panc.g.col2)-min(panc.g.col2))
        }
        if(length(gene) > 2) {
            cat("three colors!")
            panc.g.col3 <- counts.log[gene[3],]
            panc.g.col3 <- (panc.g.col3-min(panc.g.col3))/(max(panc.g.col3)-min(panc.g.col3))
        }
#        sum(panc.genes.max <= 1)
#        panc.genes.max[panc.genes.max <= 1] <- 1
#        panc.genes.max[panc.genes.max == 0] <- 0.0001
#        panc.g.col <- as.numeric(log2(panc.genes.max))
        col <- rgb(panc.g.col1, panc.g.col2, panc.g.col3)
#        cat(col)
        return(col)
    }
    return(NA)
}

vecs2rgb <- function(r, g=NULL, b=NULL, NA.ret = NA) {
    r[is.na(r)] <- NA.ret
    r <- (r-min(r))/(max(r)-min(r))
    G <- rep(0, length(r))
    B <- rep(0, length(r))
    if(!is.null(g)) {
        G[is.na(G)] <- NA.ret
        G <- (g-min(g))/(max(g)-min(g))
    }
    if(!is.null(b)) {
        B[is.na(B)] <- NA.ret
        B <- (b-min(b))/(max(b)-min(b))
    }
    rgb(r, G, B)
}


plot.coreg <- function(geneset, counts.log, my.tsne) {
    require(rgl)
for(gene in geneset) {
#    invisible(readline(prompt="Press [enter] to continue"))
    if(!is.na(match(gene[1], rownames(counts.log))) & !is.na(match(gene[2], rownames(counts.log)))) {
#        panc.genes.count <- counts[gene,]/colSums(counts)*1000000
#        panc.genes.max <- panc.genes.count
        panc.g.col1 <- counts.log[gene[1],]
        panc.g.col1 <- (panc.g.col1-min(panc.g.col1))/(max(panc.g.col1)-min(panc.g.col1))
        panc.g.col2 <- rep(0, dim(counts.log)[1])
        panc.g.col3 <- rep(0, dim(counts.log)[1])
        if(length(gene) > 1) {
            cat("two colors!")
            panc.g.col2 <- counts.log[gene[2],]
            panc.g.col2 <- (panc.g.col2-min(panc.g.col2))/(max(panc.g.col2)-min(panc.g.col2))
        }
        if(length(gene) > 2) {
            cat("three colors!")
            panc.g.col3 <- counts.log[gene[3],]
            panc.g.col3 <- (panc.g.col3-min(panc.g.col3))/(max(panc.g.col3)-min(panc.g.col3))
        }
#        sum(panc.genes.max <= 1)
#        panc.genes.max[panc.genes.max <= 1] <- 1
#        panc.genes.max[panc.genes.max == 0] <- 0.0001
#        panc.g.col <- as.numeric(log2(panc.genes.max))
        col <- rgb(panc.g.col1, panc.g.col2, panc.g.col3)
#        plot3d(my.rtsne$Y[,1], my.rtsne$Y[,2], my.rtsne$Y[,3], col=panc.g.col, size=10, main=gene)
        plot3d(my.tsne, col=col, size=10, main=gene)
    }
    else {
        cat(gene, " not found\n");
    }
}
}

make.html <- function(counts.log, geneset, my.tsne, file.suffix="_expression_in_CIRM.html") {
    for(gene in geneset) {
#    invisible(readline(prompt="Press [enter] to continue"))
        if(!is.na(match(gene, rownames(counts.log)))) {
#        panc.genes.count <- counts[gene,]/colSums(counts)*1000000
#        panc.genes.max <- panc.genes.count
            panc.g.col <- counts.log[gene,]

#        sum(panc.genes.max <= 1)
#        panc.genes.max[panc.genes.max <= 1] <- 1
#        panc.genes.max[panc.genes.max == 0] <- 0.0001
#        panc.g.col <- as.numeric(log2(panc.genes.max))
            panc.g.col <- (panc.g.col-min(panc.g.col))/(max(panc.g.col)-min(panc.g.col))
            panc.g.col <- colorRampPalette(c("Black", "Blue", "Yellow"))(100)[panc.g.col*98+1]
#        plot3d(my.rtsne$Y[,1], my.rtsne$Y[,2], my.rtsne$Y[,3], col=panc.g.col, size=10, main=gene)
#        plot3d(my.tsne, col=panc.g.col, size=10, main=gene)
            make.3dplot(my.tsne, panc.g.col, paste(gene, file.suffix, sep=""))

        }
        else {
            cat(gene, " not found\n");
        }
    }
}

#clust.grps <- cutree(my.clust, k=25) # Gives nice clusters!
# Interesting group 17!
#plot.sequential(names(sort(all.groups.mean[[17]], decreasing=T)[1:20]), counts.log[,o], tsne.goodcov) # Mesenchymal?
# SULT1E1 sulfotransferase catalyzes sulfate conj. to hormones, etc.
# SEMA5A axonal guidance, neuronal dev.
# TFPI Coagulation inhibitor. High in placenta/leukaemia, not in blood.
# ACTA2 alpha-actin. Smooth muscle actin. myofibroblast formation
# CALD1 gene - calmodulin, actin binding protein. Smooth muscle/non-muscle contraction
# CDH6 cadherin. Hepatocellular/renal carcinoma
# CAV1 caveolin. Keep signalling molecules in vesicles...? Something like that.
# PCDH18 Cadherin related. Specific cell-cell connections in the brain.

#plot.sequential(names(sort(all.groups.mean[[17]], decreasing=T)[1:20]), counts.log[,o], tsne.goodcov) #Vasculature? Endothelial
# MYCT1 Myc target. Up in cancer...
# CD93 Anticoagulant
# EBF1 Maintain B-cell identity(!) and prevention of alternative fates in committed cells.


make.duplicate.predictor <- function(counts) {
    require(randomForest)
    n.samples <- dim(counts)[2]
    cat("Number of samples: ", n.samples)
    n.genes <- dim(counts)[1]
    cat("Number of genes: ", n.genes)
    mixed.df <- as.data.frame(lapply(1:n.samples, function(x) {
        cat("x: ", x)
        r <- sample(1:n.samples, 2, replace=FALSE)
        merge.cells.wercc(counts[,r[1]], counts[,r[2]])
    }))
    both.df <- cbind(counts, mixed.df)
    both.df.norm <- norm.log.counts(both.df)
    rownames(both.df.norm) <- rownames(counts)
    y <- as.factor(c(rep(1,n.samples ), rep(2, n.samples)))
    my.rf1 <- randomForest(x=t(both.df.norm), y=y, ntree=2000, importance=TRUE, do.trace=T)
    cat(summary(my.rf1))
    my.rf1
}

merge.cells <- function(x, y) {
    # pick equal number of random counts from x and y
    names <- 1:length(x)
    ox <- x[names] > 0
    lx <- unlist(lapply(names[ox], function(i) {rep(i, x[i])}))

    oy <- y[names] > 0
    ly <- unlist(lapply(names[oy], function(i) {rep(i, y[i])}))

    num.to.pick <- mean(length(lx), length(ly))

    lnew <- c(sample(lx, as.integer(num.to.pick/2)), sample(lx, as.integer(num.to.pick/2)))
    zt <- table(lnew)
    z <- x-x
    z[as.integer(names(zt))] <- as.integer(zt)
    z
}


merge.cells.wercc <- function(xorig, yorig) {
    # pick equal number of random counts from x and y
    ercc.o <- grepl("^ERCC", names(xorig))
    x.ercc <- xorig[ercc.o]
    z <- xorig-xorig
    z[!ercc.o] <- xorig[!ercc.o]+yorig[!ercc.o]
    z[ercc.o] <- xorig[ercc.o]
    z
}

# Select a cutoff for bad cells based on actin expression
get.cutoff.lognorm <- function(my.counts.log, quantile.cut=0.001, gene.name='ACTB') {
    cl.act <- my.counts.log[gene.name,]
    cl.act.m <- median(cl.act)
    cl.act.sd <- sqrt(sum((cl.act[cl.act > cl.act.m] - cl.act.m)^2)/(sum(cl.act  > cl.act.m)-1))
    my.cut <- qnorm(p=quantile.cut, mean=cl.act.m, sd=cl.act.sd)
    my.cut
}

pick.cluster <- function(layout, classification) {
    plot(layout, col=classification, pch=19, cex=0.8)
#    invisible(readline(prompt="Press [enter] to pick clusters"))
    selected.clusts <- unique(classification[identify(layout, labels=classification)])
#    invisible(readline(prompt="Press [enter] to show selected cells"))
    selected.cells <- !is.na(match(classification, selected.clusts))
    plot(layout, col=selected.cells+1, pch=19, cex=0.8)
    selected.cells
}

pick.clusters <- function(layout, classification) {
    plot(layout, col=classification, pch=19, cex=0.8)
#    invisible(readline(prompt="Press [enter] to pick clusters"))
    selected.clusts <- unique(classification[identify(layout, labels=classification)])
    selected.cells <- rep("undef", length=length(classification))
    a <- invisible(readline(prompt="Write name of cluster, [Enter] to stop"))
    while(a != "") {
        selected.cells[!is.na(match(classification, selected.clusts))] <- a
        plot(layout, col=as.integer(as.factor(selected.cells))+1, pch=19, cex=0.8)
        invisible(readline(prompt="Press [enter] to pick clusters"))
        plot(layout, col=classification, pch=19, cex=0.8)
        selected.clusts <- unique(classification[identify(layout, labels=classification)])
        a <- invisible(readline(prompt="Write name of cluster, [Enter] to stop"))
    }
    plot(layout, col=as.integer(as.factor(selected.cells))+1, pch=19, cex=0.8)
    selected.cells
}


sel.by.feature <- function(counts, annot, feature, featureVals) {
    o <- !is.na(match(annot[,feature], featureVals))
    o.n <- annot$rmangled.name[o]
    o <- !is.na(match(colnames(counts), o.n))
}

revcomp <- function(dna) {
    sapply(strsplit(chartr("ATGC", "TACG", dna), NULL), function(x) {paste(rev(x), collapse='')})
}

vioplot.list <- function(l, ...) {
    require(vioplot)
    l$names <- names(l)
    names(l)[1] <- 'x'
 #   l <- c(l, as.list(match.call(expand.dots = TRUE)[-1:-2]))
#    return(l)
    do.call(vioplot, l)
}

scatterhist <- function(x, y, xlab="", ylab="", onlyNulls = FALSE){
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, plot=FALSE)
  yhist = hist(y, plot=FALSE)
  if(onlyNulls) {
      xhist = hist(x[y==0], plot=FALSE)
      yhist = hist(y[x==0], plot=FALSE)
  }
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(3,3,1,1))
  plot(x,y)
  par(mar=c(0,3,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
  par(mar=c(3,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0, 
    at=.8 * (mean(x) - min(x))/(max(x)-min(x)))
  mtext(ylab, side=2, line=1, outer=TRUE, adj=0, 
    at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
}


# targets: vector, ngenes long
# values: matrix, ngenes x ncells long
# result: char vector: ncells long color list
col.from.targets <- function(targets, values) {
    targets <- targets[1:(dim(values)[1])]
    if(is.matrix(values)) {
        v <- t(apply(values, 1, function(x) {(x-min(x))/(max(x)-min(x))}))
        fractions <- apply(v, 2, function(x) {x/sum(x)})
        fractions[is.nan(fractions)] <- 1.0/(dim(fractions)[1])
        targets.rgb <- col2rgb(targets)
        res <- vector("character", length=length(targets))
        for(i in 1:(dim(values)[2])) {
            mytarget.rgb <- 255-t(apply(targets.rgb, 1, function(x) {(255-x) * v[,i]}))
            mytarget.rgb <- rowSums(t(apply(mytarget.rgb, 1, function(x) {x * fractions[,i]})))
#        cat(fractions)
#        cat(i)
            res[i] <- rgb(red=mytarget.rgb['red']/256, green=mytarget.rgb['green']/256, blue=mytarget.rgb['blue']/256)
        }
        return(res)
    } else {
        v <- (values-min(values))/(max(values)-min(values))
#        fractions <- apply(v, 2, function(values) {values/sum(values)})
#        fractions[is.nan(fractions)] <- 1.0/(dim(fractions)[1])
        targets.rgb <- col2rgb(targets)
        res <- vector("character", length=length(targets))
        for(i in 1:length(values)) {
            mytarget.rgb <- 255-t(apply(targets.rgb, 1, function(values) {(255-values) * v[i]}))
#            mytarget.rgb <- rowSums(t(apply(mytarget.rgb, 1, function(values) {values * fractions[,i]})))
            res[i] <- rgb(red=mytarget.rgb[1]/256, green=mytarget.rgb[2]/256, blue=mytarget.rgb[3]/256)
        }
    return(res)
    }
}

plot.nice <- function(layout, gene.vals, genes, pal=NULL, maxVal=NULL, cex=1, rim.modifier=0.85) {
    if(is.null(pal)) {
        require(RColorBrewer)
        pal <- brewer.pal(9, "Set1")
    }
    targets <- pal[1:length(genes)]
    values <- gene.vals[genes,]
    if(!is.null(maxVal)) {
        if(maxVal < 1) {
            values <- t(apply(values, 1, function(x) {
                x[x > quantile(x[x != min(x)], maxVal)] <- quantile(x[x != min(x)],maxVal)
                x
            }))
        } else {
            values <- apply(values, 2, function(x) {x[x > maxVal] <- maxVal})
        }
#        values[values > max] <- max
    }
    cols <- col.from.targets(targets, values)
    plot(layout, col="black", pch=19, cex=cex)
    points(layout, col=cols, pch=19, cex=cex*rim.modifier)
    legend("bottomright", legend=genes, fill=pal)
}

sel.by.cv <- function(counts.nodups) {
    require(statmod)# library(pcaMethods); library(fastICA)
    ed <- counts.nodups*1000000/colSums(counts.nodups) # Second pass, no duplicates
    means <- rowMeans(ed)
    vars <- apply(ed,1,var)
    cv2 <- vars/means^2
    winsorize <- function (x, fraction=0.05) {
        if(length(fraction) != 1 || fraction < 0 ||
           fraction > 0.5) {
            stop("bad value for 'fraction'")
        }
        lim <- quantile(x, probs=c(fraction, 1-fraction))
        x[ x < lim[1] ] <- lim[1]
        x[ x > lim[2] ] <- lim[2]
        x
    }
    wed <- t(apply(ed, 1, winsorize, fraction=2/ncol(ed))) 
    means = rowMeans(wed); vars = apply(wed,1,var); cv2 <- vars/means^2
    useForFit <- means >= unname( quantile( means[ which( cv2 > .3 ) ], .95 ) ) 
    fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
    xg <- exp(seq( min(log(means[means>0])), max(log(means), na.rm=T), length.out=1000 ))
    afit <- fit$coef["a1tilde"]/means+fit$coef["a0"]
    vfit <- fit$coef["a1tilde"]/xg+fit$coef["a0"]
    varFitRatio <- vars/(afit*means^2)
    varorder <- order(varFitRatio,decreasing=T)
    return(varorder)
}

random.sel.graph <- function(x, startcell, k=5, frac=0.5) {
    require(igraph)
    numcells <- dim(x)[1]
    o <- c(startcell, sample((1:numcells)[-startcell], size=(numcells*frac)-1))
    adjc.mat<-apply(as.matrix(x[o,o]), 1, function(x) {
        x[x > (sort(x, decreasing=F)[k+1])] <- 0
        x
    })
#    adjc.mat[lower.tri(adjc.mat, diag=TRUE)]<-0
    g_islets <- graph.adjacency(adjc.mat > 0,weighted=TRUE,diag=FALSE,mode="undirected")
#    E(g_islets)$weight<-rep(round(t(adjc.mat)[(t(adjc.mat))>0],2), 2)
    res <- sapply(2:(numcells*frac), function(i) {
        length(V(g_islets)[get.shortest.paths(g_islets, from=1, to=i)[[1]][[1]]])
    })
    names(res) <- o[-1]
    res[res == 0] <- max(res)+1
    res
}


error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

plot.nice.whist <- function(coords, gene.vals, genes, pal=NULL) {
    if(is.null(pal)) {
        require(RColorBrewer)
        pal <- brewer.pal(9, "Set1")
    }
    ngenes <- length(genes)
    ll <- rep(1, length.out=ngenes*2)
    ll[seq(1:ngenes)*2] <- (1:ngenes)+1
    zones=matrix(ll, ncol=2, byrow=TRUE)
    layout(zones, widths=c(4/5,1/5))
    targets <- pal[1:length(genes) ]
    values <- gene.vals[genes,]
    cols <- col.from.targets(targets, values)
    plot(coords, col="black", pch=19, cex=1)
    points(coords, col=cols, pch=19, cex=0.85)
    legend("bottomright", legend=genes, fill=pal)
    sapply(genes, function(gene) {hist(gene.vals[gene,], main=gene)})
}


spatial.median <- function(matx) {
    dsums <- apply(matx, 2, function(x) {
        sum(apply(matx, 2, function(y) {
            sqrt(sum((x-y)^2))
        }))
    })
    which.min(dsums)
}

