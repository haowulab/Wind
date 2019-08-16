############################################
## functions for weighted Rand Index (wRI)
############################################

## This function computes the weights used for weighted Rand Index,
## given an expression matrix and true cell classes
createWeights <- function(Y, classes) {
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

    ## select some marker genes - this step is tricky
    ix = selectMarker(mY)
    mY2 = mY[ix,]

    ## W1 is the correlation of cluster centers
    W1 = cor(mY2) ## note it can be negative

    ## W0 is the average within group variances
    W0 = matrix(1, nrow=nrow(W1), ncol=ncol(W1))
    rownames(W0) = rownames(W1)
    colnames(W0) = colnames(W1)
    for( i in 1:nrow(W1) ) {
        ix = classes==allClasses[i]
        tmp = cor(Y[,ix])
        diag(tmp) = NA
        W0[i,i] = 1- mean(tmp, na.rm=TRUE)
    }

    list(W0=W0, W1=W1)
}


## The main Weighted Rand Index function
wRI <- function(trueclass, cluster, w0=NULL, w1=NULL) {
    J = length(unique(trueclass)) #number of true classes

    if (is.null(w0)) {
        w0 = 1-diag(J)
        colnames(w0) = names(table(trueclass))
    } # credit for keeping celltype j1 and j2  separate
    if (is.null(w1)) {
        w1 = diag(J)
        colnames(w1) = names(table(trueclass))
    }

    T1 = table(cluster, trueclass)
    T1 = T1[,colnames(w0)]
    T0 = table(trueclass, trueclass)
    T0 = T0[colnames(w0), colnames(w0)]

    S1 = GetScore(T1, w0=w0, w1=w1)
    S0 = GetScore(T0, w0=w0, w1=w1)
    wRI1 = sum(S1[1:4])/sum(S0[1:4])  ## weighted Rand Index

    NI1 = (S1["N11"]+S1["N01"])/(S1["N11"]+S1["d0"]) #adjusted positive predicative value
    NI2 = (S1["N00"]+S1["N10"])/(S1["b0"]+S1["c0"]) #adjusted negative predicative value
    n = length(trueclass) 
    output = c("wRI"=wRI1,"NI1"=NI1,"NI2"=NI2,"p1"=(S1["N11"]+S1["d0"])/choose(n,2),
               "p0"=(S1["b0"]+S1["c0"])/choose(n,2))
    names(output) = c("wRI","NI1","NI2","p1", "p0")
    output

}

GetScore <- function(T1, w0, w1) {
    n = sum(T1)
    J = ncol(T1)
    I = nrow(T1)

    a = sum(choose(T1,2)) # friends remain friends (N11)
    c = d = 0 # c:friends separated (N10); d: foes called friends (N01)
    b0 = c0 = d0 = 0

    for(i in 1:I){
        for( j in 1:J) {
            c = c + sum(T1[i,j]*T1[-i,j]*w0[j,j])/2 #partial score in N10
            d = d + sum(T1[i,j]*T1[i,-j]*w1[j,-j])/2 #partial score in N01

            b0 = b0 + sum(T1[i,j]*colSums(matrix(T1[-i,-j],ncol=J-1)))/2
            c0 = c0 + sum(T1[i,j]*T1[-i,j])/2
            d0 = d0 + sum(T1[i,j]*T1[i,-j])/2
        }
    }
    c(N11=a,N00=b0,N10=c,N01=d,b0=b0,c0=c0,d0=d0)
}
