\name{fastclime}
\alias{fastclime}

\title{
The main solver for fastclime package
}

\description{
A fast parametric simplex solver for constrainted l1 minimization approach to sparse precision matrix estimation
}

\usage{
fastclime(x, lambda.min.ratio = NULL)
}

\arguments{
  \item{x}{
There are 2 options: (1) \code{x} is an \code{n} by \code{d} data matrix (2) a \code{d} by \code{d} sample covariance matrix. The program automatically identifies the input matrix by checking the symmetry. (\code{n} is the sample size and \code{d} is the dimension)
}

  \item{lambda.min.ratio}{
It is the largest sparsity of each column for estimated graphs. This program makes the sparsity level of the graph path increase from \code{0} to \code{lambda.min.ratio}.The default value is \code{0.1}. When you have a large data set, make sure the sparsity level is small enough, otherwise the program might take a long time.
}

}

\details{
This program uses a parametric simplex linear programming method to solve a CLIME (a constrained l1 minimization sparse precision matrix estimation) problem. The solution path of the problem is corresponding to the parameter in the simplex method
}

\note{
The maximum number of nlambda is 50, so it can only store a solution path of length 50. When the data set is large, the solution path is set to be less than 15.
}

\value{
An object with S3 class \code{"fastclime"} is returned:  
  \item{data}{
The \code{n} by \code{d} data matrix or \code{d} by \code{d} sample covariance matrix from the input
}
  \item{cov.input}{
An indicator of the sample covariance. 
}
  \item{sigmahat}{
The empirical covariance of the data. If cov.inpu is TRUE, sigmahat = data
}
  \item{nlambda}{
The length of the path. The maximum value is 50 when d is small and it is 15 when d is large. The program will stop when it reaches this value
}
  \item{lambda}{
The sequence of regularization parameters for each column, it is a \code{nlambda} by \code{d} matrix. It will
be filled with 0 when the program finds the required sparsity level in that column.
}
  \item{path}{
A \code{nlambda} list of \code{d} by \code{d} adjacency matrices of estimated graphs as a graph path corresponding to \code{lambda}.
}
  \item{sparsity}{
The sparsity levels of the graph path. It is a vector of length \code{nlambda}. The program takes the average weight of the two direction, so it is possible the final sparsity level is a bit higher than the required sparsity level.
}
  \item{icov}{
A \code{nlambda} list of \code{d} by \code{d} precision matrices as an alternative graph path (numerical path) corresponding to \code{lambda}. 
}
  \item{df}{
It is a \code{d} by \code{nlambda} matrix. Each row contains the number of nonzero coefficients along the solution path. 
}
}

\author{
Haotian Pang, Han Liu and Robert Vanderbei \cr
Maintainers: Haotan Pang<hpang@princeton.edu>
}

\seealso{
\code{\link{fastclime.generator}}, \code{\link{fastclime.npn}}, \code{\link{fastclime.plot}}, \code{\link{fastclime.roc}}, and \code{\link{fastclime-package}}.
}

\examples{
#generate data
L = fastclime.generator(n = 100, d = 20)

#graph path estimation
out1 = fastclime(L$data,0.4)
fastclime.plot(out1$path[[4]])

#graph path estimation using the sample covariance matrix as the input.
out1 = fastclime(cor(L$data),0.4)
fastclime.plot(out1$path[[4]])
}
