\name{fastclime-package}
\alias{fastclime-package}
\docType{package}
\title{
Fast Parametric Simplex Solver for CLIME and Linear Programming
}
\description{
A package for generic linear programming, parameterized linear programming and constrainted l1 minimization approach to sparse precision matrix estimation
}
\details{
\tabular{ll}{
Package: \tab fastclime\cr
Type: \tab Package\cr
Version: \tab 1.2.4\cr
Date: \tab 2014-04-25\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The package "fastclime" provides 4 main functions:\cr
(1) the data generator creates random samples from multivariate normal distributions with different graph structures. Please refer to \code{\link{fastclime.generator}}.\cr
(2) The parametric simplex solver for constrainted l1 minimization approach to sparse precision matrix estimation. Please refer to \code{\link{fastclime}}.\cr
(3) The path selector function gives the path and precision matrix for a given parameter in CLIME. Please refer to \code{\link{fastclime.lambda}}.\cr
(4) A generic linear programming solver and a parameterized linear programming solver. Please refer to \code{\link{fastlp}} and \code{\link{paralp}}.\cr
}
\author{
Haotian Pang, Han Liu and Robert Vanderbei \cr
Maintainer: Haotan Pang<hpang@princeton.edu>;
}

\seealso{
\code{\link{fastclime.generator}}, \code{\link{fastclime}}, \code{\link{fastclime.plot}}, \code{\link{fastclime.lambda}}, \code{\link{fastlp}} and\code{\link{paralp}}
}
