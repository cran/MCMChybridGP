\name{generateX0}
\Rdversion{1.1}
\alias{generateX0}
\title{
Generate some initial points for the hybrid Exploratory phase
}
\description{
A function to randomly generate \code{n0}-many points
within requested bounds, \code{lb}, \code{ub}.
The points selected are to have the largest, minimum
distance between any two points.
}
\usage{
generateX0 (lb, ub, n0 = 10, npool = 100, ranX = NULL) 

}
\arguments{
  \item{lb}{Lower (finite) bounds for points generated.
}
  \item{ub}{Upper (finite) bounds for points generated.
}
  \item{n0}{The number of points to be generated.
}
  \item{npool}{The number of sets of randomly generated points
to be considered in determining the best set.
}
  \item{ranX}{An optional function to generate a random set
of points to be considered.
}
}
\value{
A matrix is returned which can then be supplied to \code{\link{hybrid.explore}}.
}
\references{
"Efficient MCMC schemes for Bayesian calibration of computer models",
Fielding, Mark, Nott, David J. and Liong Shie-Yui (2009), in preparation.
}
\author{
Mark James Fielding <stafmj@nus.edu.sg>
}
\seealso{
\code{\link{hybrid.explore}},
}
\examples{
lb <- c(-3,-3)
ub <- c(3,3)
X0 <- generateX0(lb, ub)

## Or with 'ranX' supplied:
ranX <- function() matrix(rnorm(40),ncol=2)
X1 <- generateX0(lb, ub, n0=20, ranX=ranX)
}
