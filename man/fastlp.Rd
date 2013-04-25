\name{fastlp}
\alias{fastlp}

\title{
A generic LP solver
}

\description{
A generic linear programming solver using parametric simplex method
}

\usage{
fastlp(obj,mat,rhs,lambda=0)
}

\arguments{
  \item{obj}{
The objective vector of the coefficient with length n.
}
  \item{mat}{
The constraint matrix of the linear programming with dimension m*n. Note this argument must be in matrix form even it is a vector. 
}
  \item{rhs}{
The right hand side vector of the constraint with length m.
}
  \item{lambda}{
The parametric simplex method will stop when the calculated paramter is smaller than lambda. The default value is zero and it corresponds to the optimal value.
}
}

\details{
This function is used to solve a general linear programming in standard inequality form: "maximize obj*x, subject to: mat*x<=rhs, x>=0"
}

\note{
The linear programming should be in the form "maximize obj*x, subject to: mat*x<=rhs, x>=0". If the original problem is not in this form. The user has to convert it into this form. For example, the equality constrants can be separated into two inequality constraints. 
}

\value{
The optimal value will be returned if it exists. Otherwise the function will indicate the problem is infeasible or unbounded. 
}

\author{
Haotian Pang, Han Liu and Robert Vanderbei \cr
Maintainer: Haotan Pang<hpang@princeton.edu>
}

\seealso{
\code{\link{fastclime}} and \code{\link{fastclime-package}}
}

\examples{
#generate an LP problem and solve it
A=matrix(c(-1,-1,0,1,-2,1),nrow=3)
b=c(-1,-2,1)
c=c(-2,3)
fastlp(c,A,b)
}
