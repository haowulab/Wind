#########################################################
## functions for weighted normalized mutual information
#########################################################

## This function creates reference structure,
## given an expression matrix and true cell classes
createRef <- function(Y, classes) {
    ## normalize by total
    ss = colMeans(Y)
    ss = ss / median(ss)
    Ynorm = sweep(Y, 2, ss, FUN="/")
    ## get cluster means
    allClasses = unique(classes)
    mY = matrix(0, nrow=nrow(Y), ncol=length(allClasses))
    for(i in 1:length(allClasses)) {
        ix = classes == allClasses[i]
        mY[,i] = rowMeans(Ynorm[,ix])
    }
    colnames(mY) = allClasses

    ## select some marker genes to hclust
    ix = selectMarker(mY)
    mY2 = mY[ix,]
    hc = hclust(dist(t(log(mY2+1))))

    ## create output
    weights = rev(hc$height) / max(hc$height)
    Ref = matrix(0, nrow=length(allClasses), ncol=length(allClasses)-1)
    for(i in 2:length(allClasses))
        Ref[,i-1] = cutree(hc,k=i)
    rownames(Ref) = allClasses

    list(Ref=Ref, weights=weights, hc=hc)

}


## function to select marker genes for creating reference structure
## I'm using the ones with the largest variance (of log counts), and the mean greater than 15
selectMarker <- function(YY) {
    rr = matrixStats::rowVars(log(YY+1))
    ## rr = rowMeans(mY)
    flag = rep(FALSE, length(rr))
    ix = order(rr, decreasing=TRUE)[1:1000] ## select 1000 or so
    flag[ix] = TRUE
    flag[rowMeans(YY)<15] = FALSE
    ix = which(flag)
    ix
}


## Weighte normalized mutual information
## Given the celltype structure, true classes and cluster result
## Note trueclass don't have to be the groud truth.
## It can just be a clustering of cells
wNMI <- function(ctStruct, trueclass, cluster, use.weight=TRUE) {
    Ref = ctStruct$Ref
    if(use.weight)
        weights = ctStruct$weights
    else
        weights = rep(1, length(ctStruct$weights))

    ## create extended reference data with the same number of rows as Y
    X = Ref[trueclass,]
    X1 = cbind(1, X)
    HX = rep(NA,ncol(X))
    for(i in 2:ncol(X1)) {
        HX[i-1] = condentropy(as.factor(X1[,i]), as.factor(X1[,i-1]) )
    }

    HX.Y = apply( X, 2, function(xxx) condentropy(as.factor(xxx),cluster) )
    HX.Y = diff(c(0,HX.Y))
    HX0 = sum(HX)
    HY0 = entropy(cluster)
    MI1 = (sum(HX*weights)-sum(HX.Y*weights)) / sum(HX*weights) #fraction of HX* explained by Y
    NMI1 = MI1*HX0*2 / (HX0+HY0)
    NMI1
}

## adjusted entropy
aH <- function(X,d=rep(1,ncol(X))) {
    X1 = cbind(rep(1,nrow(X)),X)
    HX = rep(NA,ncol(X))

    for(i in 2:ncol(X1)){
        HX[i-1] = condentropy(X1[,i],X1[,i-1])
    }
    sum(HX*d)
}

