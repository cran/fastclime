\name{fastclime}
\alias{fastclime}

\title{
The main solver for fastclime package
}

\description{
A fast parametric simplex solver for constrainted l1 minimization approach to sparse precision matrix estimation
}

\usage{
fastclime(x, lambda.min = 0.1, nlambda = 50)
}

\arguments{
  \item{x}{
There are 2 options: (1) \code{x} is an \code{n} by \code{d} data matrix (2) a \code{d} by \code{d} sample covariance matrix. The program automatically identifies the input matrix by checking the symmetry. (\code{n} is the sample size and \code{d} is the dimension)
}

  \item{lambda.min}{
This is the smallest value of lambda you would like the solver to explorer. The default value is \code{0.1}. If \code{nlambda} is large enough, the precision matrix selector function \code{\link{fastclime.lambda}} will be able to find all precision matrix corresponding to all lambda values ranging from \code{1} to \code{lambda.min}.
}

  \item{nlambda}{
It is the number of the path length one would like to achieve. The default length is 50. Note if d is large and nlambda is also large, it is possible that the program will fail to allocate memory for the path.
}
}

\details{
This program uses parametric simplex linear programming method to solve CLIME (Constrained l1 Minimization Sparse Precision Matrix Estimation) problem. The solution path of the problem corresponds to the parameter in the parametric simplex method.
}

\note{
The program will stop when either the maximum number of iteration for each column \code{nlambda} is achieved or when the required \code{lambda.min} is achieved for each column. When the dimension is huge, make sure \code{nlambda} is small so that there are enough memory to allocate the solution path. \code{lambdamtx} and \code{icovlist} will be used in \code{\link{fastclime.lambda}}.
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
  \item{maxnlambda}{
The length of the path. If the program finds \code{lambda.min} in less than \code{nlambda} iterations for all columns, then the acutal maximum lenth for all columns will be returned. Otherwise it equals \code{nlambda}.
}
  \item{lambdamtx}{
The sequence of regularization parameters for each column, it is a \code{nlambda} by \code{d} matrix. It will be filled with 0 when the program finds the required \code{lambda.min} value for that column. This parameter is required for \code{\link{fastclime.lambda}}.
}
  \item{icovlist}{
A \code{nlambda} list of \code{d} by \code{d} precision matrices as an alternative graph path (numerical path) corresponding to \code{lambdamtx}. This parameter is also required for \code{\link{fastclime.lambda}}.
}

}

\author{
Haotian Pang, Han Liu and Robert Vanderbei \cr
Maintainer: Haotan Pang<hpang@princeton.edu>
}

\seealso{
\code{\link{fastclime.generator}}, \code{\link{fastclime.plot}}, \code{\link{fastclime.lambda}} and \code{\link{fastclime-package}}.
}

\examples{
#generate data
L = fastclime.generator(n = 100, d = 20)

#graph path estimation
out1 = fastclime(L$data,0.1)
O = fastclime.lambda(out1$lambdamtx, out1$icovlist,0.2)
fastclime.plot(O$path)

#graph path estimation using the sample covariance matrix as the input.
out1 = fastclime(cor(L$data),0.1)
O = fastclime.lambda(out1$lambdamtx, out1$icovlist,0.2)
fastclime.plot(O$path)
}
