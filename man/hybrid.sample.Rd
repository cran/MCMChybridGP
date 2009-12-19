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
or all of the chains.
}
\usage{
hybrid.sample(Explore = NULL, n = 1000, L = NULL, delta = NULL,
nchains = NULL, T.mult = NULL, cor.max=0.9, r=5, nswaps = NULL,
replace.target = c(0, 1, 2), graph = FALSE, graph.template = NULL,
key.vars = NULL, delta.increase = NULL, nrefresh.graph = NULL)
}
\arguments{
  \item{Explore}{
The ouput from \code{\link{hybrid.explore}}.
}
  \item{n}{
The number of sampling iterations.
}
  \item{L}{
An optional integer argument.  The number of steps used in Leapfrog moves.
(Defaults to the same value used in \code{\link{hybrid.explore}}.)
}
  \item{delta}{
An optional numerical argument.  The size of steps used in Leapfrog moves.
(Defaults to the same value used in \code{\link{hybrid.explore}}.)
}
  \item{nchains}{
An optional integer argument.  The number of MCMC parallel chains to be used.
(Defaults to the same value used in \code{\link{hybrid.explore}}.)
}
  \item{T.mult}{
An optional integer argument.  The number of Parallel Tempering chains to be used.
(Defaults to the same value used in \code{\link{hybrid.explore}}.)
}
  \item{cor.max}{
An optional numerical argument between 0 and 1 (default 0.9).  Gives
the maximum amount a point can vary from points in Gaussian process fit.
Gauge this by \code{sd_y} output from \code{\link{hybrid.explore}}.
}
  \item{r}{
An optional numerical argument (default 5).
A penalty factor on points varying from points in Gaussian process fit,
when replacing the target distribution.
}
  \item{nswaps}{
An optional integer argument.  The number of repititions in swaps proposed
between MCMC chains (defaulting to value in \code{\link{hybrid.explore}}).
}
  \item{replace.target}{
The sampling scheme to be used (0, 1 or 2) in acceptance of
MCMC proposals.  Where 0 represents using the true target
distribution in all chains. 1 (default) represents using the true target
distribution only in the primary chain (with temperature 1).
2 represents replacing the target distribution in all chains by the
Gaussian process approximation.
}
  \item{graph}{
An optional boolean argument (default is FALSE).
To request a graphical progress display during the Sampling phase.
}
  \item{graph.template}{
An optional function argument.
A template on which to draw points generated during the Sampling phase.
By default a graph would be produced scaled appropriately for points sampled.
Its form can be seen by \code{\link{.hybrid.template}}.
}
  \item{key.vars}{
An optional integer argument.
A column denoting any variables of particular interest to be graphed.
}
  \item{delta.increase}{
An optional boolean argument for a proportionally larger Leapfrog 
step size with larger temperatures.
}
  \item{nrefresh.graph}{
An optional integer argument.
The number of cycles through the MCMC chains after which the graphical
display should be refreshed.
}
}
\details{
The method used in \code{\link{hybrid.sample}} 
is fully described in Fielding, Mark, Nott, David J. and Liong Shie-Yui (2009).
}
\value{
A list is returned consisting of the following.
\item{SAMP }{A matrix with rows corresponding to sampled points generated from
the target distribution.}
\item{y }{A column of the corrresponding values of the log-density
of the target distribution.}
\item{Ey }{Corresponding values of the (zero-mean) Gaussian process fit.}
\item{fcalls }{The number of function calls to evaluate the true log-density.}
\item{acceptance }{A column of 0's (for proposals not accepted) and
1's (for proposals accepted) for outcomes of Leapfrog proposals.}
\item{swap.acceptance }{A list of outcomes of swaps proposed,
consisting of \code{Primary}, swaps with the primary chain,
and \code{Others}, swaps amongst remaining chains.
Each gives a column of 0's (for proposals not accepted) and
1's (for proposals accepted).}
}
\references{
"Efficient MCMC schemes for Bayesian calibration of computer models",
Fielding, Mark, Nott, David J. and Liong Shie-Yui (2009), in preparation.
}
\author{
Mark J. Fielding <stafmj@nus.edu.sg>
}
\note{
The method used in \code{\link{hybrid.sample}} gives extensions
to the work of Rasmussen (2003), as described in
Fielding, Mark, Nott, David J. and Liong Shie-Yui (2009).
}
\seealso{
\code{\link{hybrid.explore}}
}
\examples{
\dontrun{

    mu1 <- c(-1, -1)
    mu2 <- c(+1, +1)
    sigma.sq <- 0.16
    X0 <- matrix(c(-2,-1,0,-2,0,2,0,1,2, -2,-1,-2,0,0,0,2,1,2), ncol = 2)
    f <- function(x) {
        px <- 1/4/pi/sqrt(sigma.sq) * exp(-1/2/sigma.sq *
            sum((x - mu1)^2)) + 1/4/pi/sqrt(sigma.sq) *
            exp(-1/2/sigma.sq * sum((x - mu2)^2))
        return(log(px))
    }

    explore.out <- hybrid.explore(f, X0, n=100, graph=TRUE)

    sample.out <- hybrid.sample(explore.out, n=500, graph=TRUE)

    opar <- par(mfrow=c(2,1))
    plot(density(sample.out$SAMP[,1]), xlab="x1", ylab="f(x)")
    plot(density(sample.out$SAMP[,2]), xlab="x2", ylab="f(x)")
    par(opar)

}
}
