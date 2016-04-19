\name{dantzig.selector}
\alias{dantzig.selector}

\title{
Dantzig selector
}

\description{
Function used to select the solution path for a given lambda
}

\usage{
dantzig.selector(lambdalist, BETA0, lambda)
}

\arguments{
  \item{lambdalist}{
\code{lambdalist} is the decrementing path of the lambda solution values. It is obtained from the \code{dantzig} function.
}
  \item{BETA0}{
\code{BETA0} is a \code{d} by \code{validn} matrix where each column has an estimated regression coefficient vector given a given lambda interval. It is obtained from the \code{dantzig} function.
}
  \item{lambda}{
\code{lambda} is the lambda solution value the user wishes to estimate a regression coefficient vector with.
}
}

\details{
Finds the estimated regression coefficient vector associated with a given \code{lambda}
}

\value{
\item{beta0}{
\code{beta0} is the estimated regression coefficient vector for the given \code{lambda}.
}
}

\author{
Di Qi, Haotian Pang \cr
Maintainer: Di Qi <dqi@princeton.edu>
}

\seealso{
\code{\link{dantzig}}
}

\examples{

# generate data
a = dantzig.generator(n = 200, d = 100, sparsity = 0.1)

# regression coefficient estimation
b = dantzig(a$X0, a$y, lambda = 0.1, nlambda = 100)

# estimated regression coefficient vector
c = dantzig.selector(b$lambdalist, b$BETA0, 15)

}
