\name{createRef}
\alias{createRef}
\title{
  Create reference hierarchy for computing weighted normalized mutual
  information (wNMI)
}

\description{
  This function takes an expression matrix and true cell classes, and
  returns two square (JxJ) matrices, where J is the
  number of classes provided by true cell class.
  
}
\usage{
createRef(Y, classes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{The expression matrix. Rows are genes, columns are cells.}
  \item{classes}{True cell classes.}
}

\details{
The weights used in wNMI are computed based on the heights of the branches
of the hierarchical tree of different cell types. We first compute mean
expression profiles for all cell types, and then select 1000 marker
genes using the same procedure as described in 'createWeights' function. A hierarchical tree
is then constructed based on the mean expression from marker genes. We
then obtain the tree height (from the bottom) at each branching
point. The ratios of these heights to the maximum tree height are used
as the weight in computing wNMI.
}

\value{
  A list with following fields:
  \item{Ref}{A matrix for the reference cell type clustering when
	cutting the full hierarchical tree at different branching point. Each column
	represents a cut. For example, if a column is [1,1,1,2,2,3], it
	means the first 3 cell types are deemed in the same cluster, the
	next two are in another cluster, and the last cell type is in a
	cluster by itself.
  }
  \item{weights}{A vector for the relative tree heights at different
	cut. This will be used as weights in computing wNMI.}
  \item{hc}{An object of class 'hclust'. This is the hierarchical tree
	of the reference cell types.} 
}

\author{
Zhijin Wu <zhijin_wu@brown.edu>, Hao Wu <hao.wu@emory.edu>
}


\seealso{

}
\examples{
data(Zhengmix8eq)
ctStruct = createRef(Y, trueclass)
plot(ctStruct$hc)
}
