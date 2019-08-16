\name{createWeights}
\alias{createWeights}

\title{
  Create weights used for computing weighted Rand Index (wRI)
  
}
\description{
  This function takes an expression matrix and true cell classes, and
  returns two square (JxJ) matrices, where J is the
  number of classes provided by true cell class. 
}

\usage{
createWeights(Y, classes)
}

\arguments{
  \item{Y}{The expression matrix. Rows are genes, columns are cells.}
  \item{classes}{True cell classes.}
}

\details{
  There are two weight matrices W1 (score for putting two cells in the
  same cluster) and W0 (score for separating two cells in different
  clusters). We set W1 has diagonal values 1 and off-diagonal less than
  1. To obtain off-diagonal values entry (i,j) in W1, we compute the
  mean expression profiles for all cell types from single-cell data, and
  then use the Pearson’s correlation of mean expression from cell types
  i and j as W1[i,j]. Note that in general, correlations among
  expression from different cell types are high because a majority of
  genes are not differentially expressed among cell types. To make W1
  scores more distinctive, we take a marker selection step to pick the
  top 1000 genes with the largest variances of log expressions. The
  Pearson’s correlation of mean expression from these genes are then
  taken as the off-diagonal values for W1[i,j].

  We make W0 has off diagonal values 1 and diagonal between 0 and 1,
  reflecting that keeping the
  separation existing in the reference receives full credit, but
  breaking a (weak) tie may not reduce the score completely to 0. We
  compute W0[i,i] based on the inter-cellular expression variances within
  cell type i. To be specific, we take expressions for all cells in cell
  type i and compute their Pearson’s correlations. W0[i,i] is defined as
  1 minus the average inter-cellular Pearson’s correlations from all
  cells. Thus, for a tight cluster where the within cell type
  correlation is high, the score will be closer to 0, indicating
  stronger penalty for separating cells for such a cluster. On the other
  hand, for a loose cluster where the within cell type correlation is
  low, the score will be larger, indicating weaker penalty.
}

\value{
  A list with following fields:
  \item{W0}{A square matrix of dimension JxJ, where J is the
	number of classes provided by true cell class. The (i,j)-th entry is
	the	score of separating two cells of types i and j in different
	clusters.}
  \item{W1}{A square matrix of dimension JxJ, where J is the
	number of classes provided by by true cell class The(i,j)-th entry is
	the score for putting two cells of different types in the same
	cluster.}
}

\author{
  Zhijin Wu <zhijin_wu@brown.edu>, Hao Wu <hao.wu@emory.edu>
}

\seealso{
wRI
}

\examples{
data(Zhengmix8eq)
weights = createWeights(Y, trueclass)
}
