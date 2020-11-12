\name{hybrid.sample}
\Rdversion{1.1}
\alias{hybrid.sample}
\alias{.hybrid.template}
\title{
Sampling phase applying results from Exploratory phase.
}
\description{
Sampling phase in hybrid MCMC, which takes the output from
\code{\link{hybrid.explore}} to samples from the target
distribution supplied.  The number of chains,
Leapfrog moves and Gaussian process parameters are the
same as used in \code{\link{hybrid.explore}},
or the values may be updated here.  For a target distribution
time consuming to evaluate, the target can be replaced
completely by the Gaussian process approximation in some
or all of the chains.  Bounds supplied act as reflecting
barriers.
}
\usage{
hybrid.sample(Explore, n = 1000, replace.target = c(0, 1, 2),
              lb = NULL, ub = NULL, L = NULL, delta = NULL,
              nchains = NULL, T.mult = NULL, maxleap = NULL,
              r = 5, nswaps = NULL, graph = FALSE)
}
\arguments{
  \item{Explore}{
Output from \code{\link{hybrid.explore}}.
A list object consisting of \code{X}, \code{y}, \code{f} and \code{GPfit}.
\code{X} is the set of points from which a Gaussian process
fit is gained. Then \code{y} is the corresponding values
of the target log-density function \code{f}, with further arguments
\code{...} as supplied to \code{\link{hybrid.explore}}.
If \code{GPfit} is \code{NULL} then a new Gaussian process fit is determined.
}
  \item{n}{
The number of sampling iterations.
}
  \item{L}{
An optional integer argument passed from \code{\link{hybrid.explore}}.
The number of steps used in Leapfrog moves.
}
  \item{delta}{
An optional numerical argument passed from \code{\link{hybrid.explore}}.
The size of steps used in Leapfrog moves.
}
  \item{lb}{
An optional numeric argument passed from \code{\link{hybrid.explore}}.
Lower (finite) bounds placed on X.
(If supplied, then \code{ub} must also be supplied and finite.)
}
  \item{ub}{
An optional numeric argument passed from \code{\link{hybrid.explore}}.
Upper (finite) bounds placed on X.
(If supplied, then \code{lb} must also be supplied and finite.)
}
  \item{nchains}{
An optional integer argument passed from \code{\link{hybrid.explore}}.
The number of MCMC parallel chains to be used.
}
  \item{T.mult}{
An optional integer argument passed from \code{\link{hybrid.explore}}.
The number of Parallel Tempering chains to be used.
}
  \item{maxleap}{
An optional numerical argument passed from \code{\link{hybrid.explore}}.
Gives the maximum standard deviation of \code{f(x)} that a point can vary
from points in Gaussian process fit.
}
  \item{r}{
An optional numerical argument (default 5).
A penalty factor on points straying from the region of Gaussian process fit,
when the target distribution is replaced.
}
  \item{nswaps}{
An optional integer argument passed from \code{\link{hybrid.explore}}.
The number of repititions of swaps proposed between MCMC chains.
}
  \item{replace.target}{
The sampling scheme to be used (0, 1 or 2) in acceptance of
MCMC proposals.  Where 0 represents using the true target
distribution in all chains. 1 (default) represents using the true target
distribution only in the primary chain (having temperature 1).
2 represents replacing the target distribution in all chains by the
Gaussian process approximation.
}
  \item{graph}{
An optional boolean argument (default is FALSE).
Request graphical progress display during the sample phase.
}
}
\details{
The method used in \code{\link{hybrid.sample}} 
is described in Fielding, Nott and Liong (2011).
}
\value{
A list is returned consisting of the following.
\item{SAMP }{A matrix with rows corresponding to sampled points generated from
the target distribution.}
\item{y }{A column of the corrresponding values of the log-density
of the target distribution.}
\item{acceptance }{A column of 0 (rejected) and 1 (accepted) giving a record
of sampling proposal acceptance.}
\item{function.calls }{The number of function calls to evaluate the true log-density.}
}
\references{
"Efficient MCMC Schemes for Computationally Expensive Posterior Distributions",
Fielding, Nott and Liong (2011).
}
\author{
Mark J. Fielding <mark.fielding@gmx.com>
}
\note{
A record is kept throughout a run of \code{\link{hybrid.sample}}
stored as a global variable list, \code{hybrid.sample.out}.
Useful for a run stopped prematurely.

The method used in \code{\link{hybrid.sample}} gives extensions
to the work of Rasmussen (2003) and is described in
Fielding, Nott and Liong (2011).
}
\seealso{
\code{\link{hybrid.explore}}
}
\examples{
    mu1 <- c(-1, -1)
    mu2 <- c(+1, +1)
    sigma.sq <- 0.16
    ub <- c(1.5, 3)
    X0 <- matrix(c(-2,-1, 0,-2, 0, 1, 0, 1, 1,
                   -2,-1,-2, 0, 0, 0, 2, 1, 2), ncol = 2)
    f <- function(x) {
        px <- 1/4/pi/sqrt(sigma.sq) * exp(-1/2/sigma.sq *
            sum((x - mu1)^2)) + 1/4/pi/sqrt(sigma.sq) *
            exp(-1/2/sigma.sq * sum((x - mu2)^2))
        return(log(px))
    }

\donttest{
    explore.out <- hybrid.explore(f, X0, ub=ub, n=150, graph=TRUE)
    sample.out <- hybrid.sample(explore.out, n=500, graph=TRUE)

    opar <- par(mfrow=c(2,1))
    plot(density(sample.out$SAMP[,1]), xlab="x1", ylab="f(x)")
    plot(density(sample.out$SAMP[,2]), xlab="x2", ylab="f(x)")
    par(opar)
}
}
