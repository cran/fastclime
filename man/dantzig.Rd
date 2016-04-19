\name{dantzig}
\alias{dantzig}

\title{
A solver for the Dantzig selector estimator
}

\description{
Implementation of the Primal Dual (i.e. Self Dual) Simplex Method on Dantzig selector
}

\usage{
dantzig(X, y, lambda = 0.01, nlambda = 50)
}

\arguments{
 \item{X}{
\code{x} is an \code{n} by \code{d} data matrix
}
 \item{y}{
\code{y} is a length \code{n} response vector
}
  \item{lambda}{
The parametric simplex method will stop when the calculated parameter is smaller than lambda. The default value is \code{0.01}.
}
  \item{nlambda}{
 This is the number of the maximum path length one would like to achieve. The default length is 50.
}
}

\details{
This program applies the parametric simplex linear programming method to the Dantzig selector to solve for the regression coefficient vector. The solution path of the problem corresponds to the parameter in the parametric simplex method.
}

\note{
The program will stop when either the maximum number of iterations for each column \code{nlambda} is achieved or when the required \code{lambda} is achieved for each column. Note if d is large and nlambda is also large, it is possible that the program will fail to allocate memory for the path.
}



\value{
An object with S3 class \code{"dantzig"} is returned:  
  \item{X}{
\code{X} is the \code{n} by \code{d} data matrix.
}
  \item{y}{
\code{y} is a length \code{n} response vector.
}
  \item{BETA0}{
\code{BETA0} is a \code{d} by \code{validn} matrix where each column has an estimated regression coefficient vector given a lambda interval.
}

  \item{n0}{
\code{n0} is the number of rows in the \code{n} by \code{d} data matrix.
}
  \item{d0}{
\code{d0} is the number of columns in the \code{n} by \code{d} data matrix.
}
  \item{validn}{
\code{validn} is the number of solutions along the solution path. The maximum is \code{nlambda}.
}
  \item{lambdalist}{
\code{lambdalist} is the decrementing path of the lambda solution values.
}

}

\author{
Haotian Pang, Han Liu, Robert Vanderbei and Di Qi \cr
Maintainer: Haotian Pang<hpang@princeton.edu>
}

\seealso{
\code{\link{dantzig.selector}}
}

\examples{

#generate data
a = dantzig.generator(n = 200, d = 100, sparsity = 0.1)

#regression coefficient estimation
b = dantzig(a$X0, a$y, lambda = 0.1, nlambda = 100)

}
