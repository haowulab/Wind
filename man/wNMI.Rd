\name{wNMI}
\alias{wNMI}
\title{
  Weighted normalized mutual information (wNMI)
}
\description{
  This function computes the weighted normalized mutual information
  (wNMI),  which computes the agreement of two partitions (often a gold
  standard reference and a clustering result) based on a modified
  normalized mutual information definition to account for cell type
  hierarchy. 
}

\usage{
wNMI(ctStruct, trueclass, cluster, use.weight = TRUE)
}

\arguments{
  \item{ctStruct}{The reference hierarchy, returned from `createRef' function.}
  \item{trueclass}{A vector for the true classes of cells. Often from an outside experiement.}
  \item{cluster}{A vector for the clustering results.}
  \item{use.weight}{A logical for using the weight or not. If TRUE,
	weighted normalized mutual information will be computed. If FALSE,
	the weights are not used and the result will just be the traditional
	normalized mutual information.} 
}


\value{
  A numeric for the value of weighted normalized mutual information. 
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
  Zhijin Wu <zhijin_wu@brown.edu>, Hao Wu <hao.wu@emory.edu>
}

\seealso{
  createRef, wRI
}

\examples{
data(Zhengmix8eq)
ctStruct = createRef(Y, trueclass)

methods = names(clusterRes)
allNMI = matrix(0, nrow=length(methods), ncol=2)
rownames(allNMI) = methods
colnames(allNMI) = c("NMI", "wNMI")
for(i in 1:length(clusterRes)) {
    allNMI[i,1] = wNMI(ctStruct, trueclass, clusterRes[[i]], FALSE)
    allNMI[i,2] = wNMI(ctStruct, trueclass, clusterRes[[i]])
}
barplot(t(allNMI), beside=TRUE, ylim=c(0.4,1.05), legend.text=TRUE, xpd=FALSE)
}
