\name{fastclime.lambda}
\alias{fastclime.lambda}

\title{
A precision matrix and path selector function for fastclime
}

\description{
Select the precision matrix and solution path for a given parameter lambda
}

\usage{
fastclime.lambda(lambdamtx, icovlist, lambda)
}

\arguments{
  \item{lambdamtx}{
The sequence of regularization parameters for each column, it is a \code{nlambda} by \code{d} matrix.}

  \item{icovlist}{
A \code{nlambda} list of \code{d} by \code{d} precision matrices as an alternative graph path (numerical path) corresponding to \code{lambdamtx}. 
}
  \item{lambda}{
The user specified parameter lambda. The function will return the solution path corresponding to this value. Note \code{lambda} has to be larger than or equal to \code{lambda.min} input in fastclime.
}
}

\details{
The output from \code{\link{fastclime}} stores a list of precision matrice and a matrix of parameters. This program will select the required solution path and precision matrix for a given parameter \code{lambda}.
}

\note{
The function is able to estimate the precision matrices corresponding to all \code{lambda} values ranging from \code{1} to \code{lambda.min}, provided a large enough \code{nlambda} is used in \code{\link{fastclime}}. The function will give a message if the program could not find the solution path corresponding to the required \code{lamba}. The user may want to incease \code{nlambda} in \code{\link{fastclime}} in order to find the required solution path.
}

\value{
An object with S3 class "fastclime.lambda" is returned:
 \item{icov}{
The estimated precision matrix corresponding to \code{lambda}.
}
  \item{path}{
The estimated graph path corresponding to \code{lambda}.
}
  \item{sparsity}{
The sparsity level of this estimated graph path.
}

}

\author{
Haotian Pang, Han Liu and Robert Vanderbei \cr
Maintainer: Haotan Pang<hpang@princeton.edu>
}

\seealso{
\code{\link{fastclime}} and \code{\link{fastclime-package}}
}

\examples{
#generate data
L = fastclime.generator(n = 100, d = 20)

#graph path estimation
out1 = fastclime(L$data,0.1)
O = fastclime.lambda(out1$lambdamtx, out1$icovlist,0.2)
fastclime.plot(O$path)
}
