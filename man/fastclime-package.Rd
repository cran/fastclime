\name{fastclime-package}
\alias{fastclime-package}
\docType{package}
\title{
Fast Parametric Simplex Solver for CLIME
}
\description{
A package for constrainted l1 minimization approach to sparse precision matrix estimation
}
\details{
\tabular{ll}{
Package: \tab fastclime\cr
Type: \tab Package\cr
Version: \tab 0.0.7\cr
Date: \tab 2013-03-15\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The package "fastclime" provides 3 main functions:\cr
(1) the data generator creates random samples from multivariate normal distributions with different graph structures. Please refer to \code{\link{fastclime.generator}}.\cr
(2) the nonparanormal (npn) transformation helps relax the normality assumption. Please refer to \code{\link{fastclime.npn}}.\cr
(3) The parametric simplex solver for constrainted l1 minimization approach to sparse precision matrix estimation. Please refer to \code{\link{fastclime}}.\cr
}
\author{
Haotian Pang, Han Liu and Robert Vanderbei \cr
Maintainers: Haotan Pang<hpang@princeton.edu>;
}
\seealso{
\code{\link{fastclime.generator}}, \code{\link{fastclime}}, \code{\link{fastclime.plot}} and \code{\link{fastclime.roc}}
}
