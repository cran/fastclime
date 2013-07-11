\name{paralp}
\alias{paralp}

\title{
A solver for parameterized LP problems
}

\description{
A parameterized linear programming solver using parametric simplex method
}

\usage{
paralp(obj,mat,rhs,obj_bar,rhs_bar,lambda=0)
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

  \item{obj_bar}{
The vector used to time the parameter and added to the objective vector, with length n. This pertubation vector must be nonnegative.
}

  \item{rhs_bar}{
The vector used to time the parameter and added to the right hand side vector, with length m. This pertubation vector must be nonnegative.
}

  \item{lambda}{
The parametric simplex method will stop when the calculated paramter is smaller than lambda. The default value is zero and it corresponds to the optimal value.
}
}

\details{
This function is used to solve a general linear programming in standard inequality form: "maximize obj*x+obj_bar*lambda, subject to: mat*x<=rhs+rhs_bar*lambda, x>=0"
}


\value{
The optimal value will be returned if it exists with a proper value of chosen lambda. Otherwise the function will indicate the problem is infeasible or unbounded. 
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
b_bar=c(1,1,1)
c_bar=c(1,1)
paralp(c,A,b,c_bar,b_bar)
}
