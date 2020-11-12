\name{generateX0}
\Rdversion{1.1}
\alias{generateX0}
\title{
Generate some initial points for the hybrid explore phase
}
\description{
A function to randomly generate \code{n0}-many points
within requested bounds, \code{lb}, \code{ub}.
The points selected are to have the largest, minimum
distance between any two points.
}
\usage{
generateX0 (lb, ub, n0 = 10, npool = 100)

}
\arguments{
  \item{lb}{Lower (finite) bounds for points generated.
}
  \item{ub}{Upper (finite) bounds for points generated.
}
  \item{n0}{The number of points to be generated.
}
  \item{npool}{The size of the pool of sets of randomly
generated points in determining the best set of points.
}
}
\value{
A matrix is returned which can then be supplied to
\code{\link{hybrid.explore}}.
}
\references{
"Efficient MCMC Schemes for Computationally Expensive Posterior Distributions",
Fielding, Nott and Liong (2011).
}
\author{
Mark James Fielding <mark.fielding@gmx.com>
}
\seealso{
\code{\link{hybrid.explore}},
}
\examples{
lb <- c(-3,-3)
ub <- c(3,3)
X0 <- generateX0(lb, ub)
}

