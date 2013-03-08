\name{fastclime.npn}
\alias{fastclime.npn}

\title{
Nonparanormal(npn) transformation
}

\description{
Implements the Gausianization to help relax the assumption of normality.
}

\usage{
fastclime.npn(x, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)
}

\arguments{
  \item{x}{
The \code{n} by \code{d} data matrix representing \code{n} observations in \code{d} dimensions
}
  \item{npn.func}{
The transformation function used in the npn transformation. If \code{npn.func = "truncation"}, the truncated ECDF is applied. If \code{npn.func = "shrinkage"}, the shrunken ECDF is applied. The default is \code{"shrinkage"}. If \code{npn.func = "skeptic"}, the nonparanormal skeptic is applied.
}
  \item{npn.thresh}{
The truncation threshold used in nonparanormal transformation, ONLY applicable when \code{npn.func = "truncation"}. The default value is \code{1/(4*(n^0.25)*} \code{sqrt(pi*log(n)))}.
}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}

\details{
The nonparanormal extends Gaussian graphical models to semiparametric Gaussian copula models.Motivated by sparse additive models, the nonparanormal method estimates the Gasussian copula by marginally transforming the variables using smooth functions.Computationally, the estimation of a nonparanormal transformation is very efficient and only requires one pass of the data matrix.
}

\value{
\item{data}{
A \code{d} by \code{d} nonparanormal correlation matrix if \code{npn.func = "skeptic"}, and A \code{n} by \code{d} data matrix representing \code{n} observations in \code{d} transformed dimensions other wise.
}
}

\author{
Haotian Pang, Han Liu and Robert Vanderbei \cr
Maintainers: Haotan Pang<hpang@princeton.edu>
}

\seealso{
  \code{\link{fastclime}} and \code{\link{fastclime-package}}.
}

\examples{
# generate nonparanormal data
L = fastclime.generator(graph = "cluster", g = 5)
L$data = L$data^5

# transform the data using the shrunken ECDF
Q = fastclime.npn(L$data)

# transform the non-Gaussian data using the truncated ECDF
Q = fastclime.npn(L$data, npn.func = "truncation")

# transform the non-Gaussian data using the truncated ECDF
Q = fastclime.npn(L$data, npn.func = "skeptic")

}
