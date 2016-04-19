\name{dantzig.generator}
\alias{dantzig.generator}

\title{
Dantzig data generator
}

\description{
Generates sparse linear regression model for testing \code{dantzig} function
}

\usage{
dantzig.generator(n = 50, d = 100, sparsity = 0.1, sigma0=1)
}

\arguments{
  \item{n}{
The number of observations (sample size). The default value is \code{50}.
}
  \item{d}{
The number of variables (dimension). The default value is \code{100}.
}
  \item{sparsity}{
\code{d} is either the number of nonzero entries out of \code{d} or the proportion of nonzero entries in \code{BETA}
}
  \item{sigma0}{
\code{sigma0} is the standard deviation of the noise vector
}
}

\details{
Generates sparse linear regression model for testing \code{dantzig} function
}

\value{
An object with S3 class "dantzig.generator" is returned:
\item{X0}{
\code{X0} is the \code{n} by \code{d} matrix for the generated data
}
\item{y}{
A \code{y} is a \code{n} response vector for the generated data
}
  \item{BETA}{
\code{BETA} is a length \code{d} regression coefficient vector
}
  \item{s}{
\code{s} is the number of nonzero entries out of \code{d}
}
  \item{pos}{
A vector containing the indices of the nonzero entries (may contain repeats)
}
}

\author{
Haotian Pang, Han Liu, Robert Vanderbei and Di Qi \cr
Maintainer: Di Qi <dqi@princeton.edu>
}

\seealso{
\code{\link{dantzig}}
}

\examples{
	
##
L = dantzig.generator(n = 50, d = 100, sparsity = 0.1)	

}
