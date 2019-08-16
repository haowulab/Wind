\name{wRI}
\alias{wRI}

\title{
Weighted Rand Index (wRI)
}

\description{
  This function computes the weighted Rand index (wRI), which computes
  the agreement of two partitions (often a gold standard reference and a
  clustering result) based on a modified Rand index definition to
  account for cell type hierarchy.  
}

\usage{
wRI(trueclass, cluster, w0 = NULL, w1 = NULL)
}

\arguments{
  \item{trueclass}{A vector for the true classes of cells. Often from an
	outside experiement.}
  \item{cluster}{A vector for the clustering results.}
  \item{w0}{An optional square matrix of dimension JxJ, where J is the
	number of classes provided by 'trueclass'. The (i,j)-th entry is the
	score of separating two cells of types i and j in different
	clusters.}
  \item{w1}{An optional square matrix of dimension JxJ, where J is the
	number of classes provided by 'trueclass'. The(i,j)-th entry is the
	score for putting two cells of different types in the same cluster.}
}

\details{
  The w0 and w1 matrics can be specified by user based on based on prior
  biological knowledge. If not given, they will be estimated using
  'createWeights' function.
}

\value{
  A vector of numerics for values related to the Weighted Rand Index (wRI),
  including:
  \item{wRI}{The wRI value.}
  \item{NI1}{The index with the positive prediction value interpretation: given two subjects are put in the same cluster, the probability that they truely belong to the same type.}
  \item{NI2}{The index with the negative prediction value interpretation: given two subjects are separated in different clusters, the probability they belong to different types.}
  \item{p1}{The proportaion of "in-cluster" relationships in the reference, among all pair-wise relationships. }
  \item{p0}{The proportaion of "seprated in different cluster" relationships in the reference, among all pair-wise relationships.}
}

\author{
  Zhijin Wu <zhijin_wu@brown.edu>, Hao Wu <hao.wu@emory.edu>
}

\seealso{
  createWeights, wNMI
}

\examples{
data(Zhengmix8eq)
weights = createWeights(Y, trueclass)
methods = names(clusterRes)
allRI = matrix(0, nrow=length(methods), ncol=6)
rownames(allRI) = methods
colnames(allRI) = c("RI", "NI1","NI2","wRI","wNI1","wNI2")
for(i in 1:length(clusterRes)) {
    allRI[i,1:3] = wRI(trueclass, clusterRes[[i]]) [1:3]
    allRI[i,4:6] = wRI(trueclass, clusterRes[[i]], weights$W0, weights$W1)[1:3]
}
barplot(t(allRI[,c(1,4)]), beside=TRUE, ylim=c(0.7,1.05), legend.text=TRUE, xpd=FALSE)
}

