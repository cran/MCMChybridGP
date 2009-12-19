\name{hybrid.explore}
\Rdversion{1.1}
\alias{hybrid.explore}
\title{
Exploratory phase to determine points used for a Gaussian process fit.
}
\description{
The Exploratory phase of hybrid MCMC using a
Gaussian process approximation of derivatives.
The user must provide the log-density of a target distribution
and a small set of (say 10) points forming the matrix X0.
Using a Gaussian process approximation, points are added
until a appropriate set of points are determined
for a final Gaussian process fit to the target distribution.
Results can then be passed to the Sampling phase
or the Gaussian process (zero mean)
approximation \code{Ef} can be used to approximate \code{f}.
}
\usage{
hybrid.explore(f, X0, n = 200, L = 1000, delta = 0.001,
nchains = 5, T.mult = 1.5, lb = NULL, ub = NULL, cor.max = 0.9,
delta.increase = FALSE, nswaps = choose(nchains, 2), spread = 1,
graph = FALSE, graph.template = NULL, key.vars = 1:dim(X)[2],
nrefresh.graph = ceiling(20/nchains))
}
\arguments{
  \item{f}{
A function giving the log-density of the target distribution.
}
  \item{X0}{
A matrix with rows giving say 10 points for an initial
Gaussian process fit.
}
  \item{n}{
An optional integer argument.
The number of points to be generated for the final
Gaussian process fit.
}
  \item{L}{
An optional integer argument.
The number of steps to be used in Leapfrog moves in MCMC proposals.
}
  \item{delta}{
An optional numerical argument.
The size of Leapfrog steps to be made, with a suitable
balance between \code{L} and \code{delta}.
}
  \item{nchains}{
An optional integer argument.
The number of Parallel Tempering chains to be used (default 5).
}
  \item{T.mult}{
An optional numerical argument.
The Temperature multiple to be used in Parallel Tempering (default 1.5).
}
  \item{lb}{
An optional numerical argument.  Lower (finite) bounds placed on X.
(If supplied, then \code{ub} must also be supplied.)
}
  \item{ub}{
An optional numerical argument.  Upper (finite) bounds placed on X.
(If supplied, then \code{lb} must also be supplied.)
}
  \item{cor.max}{
An optional numerical argument. Used in early stops for Leapfrog
moves (default 0.9).  The Leapfrog is stopped when
correlation with respect to f(x) exceeds this value.
If \code{cor.max} greater than 1 then the Leapfrog is stopped instead
when the standard devation exceeds the value supplied.
}
  \item{delta.increase}{
An optional boolean argument.
For chains with a larger Temperature to use proportionally larger
Leapfrog step sizes.
}
  \item{nswaps}{
An optional integer argument.
The number of repetitions of proposed swaps between the parallel chains.
}
  \item{spread}{
An optional numerical argument.
To encourage greater spread of values generated for Gaussian process,
increasingly with standard deviation with respect to \code{f}.
}
  \item{graph}{
An optional boolean argument.
To request graphical progress display during the Exploratory phase.
}
  \item{graph.template}{
A template on which to draw points generated during the Exploratory phase.
When not supplied a graph will be produced appropriately for points generated.
}
  \item{key.vars}{
An optional integer argument.
A column denoting any variables of particular interest to be graphed.
}
  \item{nrefresh.graph}{
An optional integer argument.
The number of cycles through the MCMC chains after which
the graphical display should be refreshed
}
}
\details{
The method used in \code{\link{hybrid.explore}} is fully described in
Fielding, Mark, Nott, David J. and Liong Shie-Yui (2009).
}
\value{
A list is returned to be used as input to \code{\link{hybrid.sample}}.
\item{X }{A matrix made up of the set of points to be
used in the final Gaussian process fit.}
\item{y }{A column of the corresponding values of \code{f}.}
\item{Ey }{A column of the corresponding (zero-mean) Gaussian process fit.}
\item{sd_y }{List of 2 columns of corrresponding standard deviations
\code{sd} with respect to \code{f} (and standardised values \code{cor}).
Useful for setting \code{cor.max} in \code{\link{hybrid.sample}}.}
\item{f }{The original log-density function as supplied.}
\item{fcalls }{The number of function calls used to evaluate \code{f}.}
\item{acceptance }{A column of 0's (for proposals not accepted) and
1's (for proposals accepted) for the outcomes of Leapfrog proposals.}
\item{swap.acceptance }{A column of 0's (for proposals not accepted) and
1's (for proposals accepted) for the outcomes of swaps proposed
between parallel chains.}
\item{details }{A list containing information particular to \code{\link{hybrid.sample}}.}
\item{GPfit }{Ouput from \code{\link{GProcess}} with the final
Gaussian process fit, and in particular the (zero mean) approximation \code{Ef}.}
}
\references{
"Efficient MCMC schemes for Bayesian calibration of computer models",
Fielding, Mark, Nott, David J. and Liong Shie-Yui (2009), in preparation.
}
\author{
Mark J. Fielding <stafmj@nus.edu.sg>
}
\note{
The methods used in \code{\link{hybrid.explore}} and \code{\link{hybrid.sample}}
give extensions to the work of Rasmussen (2003), as described in
Fielding, Mark, Nott, David J. and Liong Shie-Yui (2009).
}
\seealso{
\code{\link{hybrid.sample}}, \code{\link{GProcess}}, \code{\link{generateX0}}.
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
