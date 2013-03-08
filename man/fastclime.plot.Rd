\name{fastclime.plot}
\alias{fastclime.plot}

\title{
Graph visualization
}

\description{
Implements the graph visualization using adjacency matrix. It can automatic organize 2D embedding layout.
}

\usage{
fastclime.plot(G, epsflag = FALSE, graph.name = "default", cur.num = 1, 
location)
}

\arguments{
  \item{G}{
The adjaceny matrix corresponding to the graph.
}
  \item{epsflag}{
If \code{epsflag = TRUE}, save the plot as an eps file in the target directory. The default value is \code{FALSE}.
}
  \item{graph.name}{
The name of the output eps files. The default value is "default".
}
  \item{cur.num}{
The number of plots saved as eps files. Only applicale when \code{epsflag = TRUE}. The default value is 1.
}
  \item{location}{
Target directory. The default value is the current working directory.
}
}

\details{
The user can change \code{cur.num} to plot several figures and select the best one. The implementation is based on the popular package "igraph".
}

\author{
Haotian Pang, Han Liu and Robert Vanderbei \cr
Maintainers: Haotan Pang<hpang@princeton.edu>
}

\seealso{
\code{\link{fastclime}} and \code{\link{fastclime-package}}
}

\examples{
## visualize the hub graph
L = fastclime.generator(graph = "hub")
fastclime.plot(L$theta)

## visualize the band graph
L = fastclime.generator(graph = "band",g=5)
fastclime.plot(L$theta)

## visualize the cluster graph
L = fastclime.generator(graph = "cluster")
fastclime.plot(L$theta)

#show working directory
getwd()
#plot 5 graphs and save the plots as eps files in the working directory  
fastclime.plot(L$theta, epsflag = TRUE, cur.num = 5)
}
