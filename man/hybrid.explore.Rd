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
or the Gaussian process approximation \code{Ef} can be
used to approximate \code{f}.
}
\usage{
hybrid.explore(f, X0, ..., y0 = NULL, n = 200, L = 1000, delta = 0.003,
                 nchains = 5, T.mult = 1.5, lb = NULL, ub = NULL,
                 maxleap=0.95, finetune = FALSE, Tinclude = nchains,
                 nswaps = choose(nchains, 2), graph = FALSE)
}
\arguments{
  \item{f}{
A function giving the log-density of the target distribution.
}
  \item{X0}{
A matrix with rows giving say 10 points for an initial
Gaussian process fit.
}
  \item{...}{
Further arguments to be passed to \code{f}.
}
  \item{y0}{
Corresponding function values if already known.
Otherwise the corresponding function values will first be evaluated.
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
An optional numeric argument (default 0.003).
The size of Leapfrog steps to be made, with a suitable
balance between \code{L} and \code{delta}.
}
  \item{nchains}{
An optional integer argument.
The number of Parallel Tempering chains to be used (default 5).
}
  \item{T.mult}{
An optional numeric argument.
The Temperature multiple to be used in Parallel Tempering (default 1.5).
}
  \item{lb}{
An optional numeric argument.  Lower bounds placed on X.
Only points within bounds are included in Gaussian process fit.
}
  \item{ub}{
An optional numeric argument.  Upper bounds placed on X.
Only points within bounds are included in Gaussian process fit.
}
  \item{maxleap}{
An optional numeric 0-1 argument (default 0.95).
Used in early stops in Leapfrog
moves.  The Leapfrog is stopped when the
standard deviation of \code{f(x)} exceeds an appropriate
corresponding value.
}
  \item{finetune}{
An optional boolean argument (default FALSE).
Option for fine-tuning Gaussian process parameters.
This may be useful when a good fit is difficult to achieve.
}
  \item{Tinclude}{
An optional integer argument (default \code{nchains}).
Limit the number of temperatures to contribute points to
the Gaussian process fit.  With low acceptance
rates for \code{T=1}, \code{Tinclude=1} might be used.
Otherwise a larger value gives a faster explore phase.
}
  \item{nswaps}{
  An optional integer argument.
  The number of repetitions of proposed swaps between the parallel chains.
  }
  \item{graph}{
  An optional boolean argument (default FALSE).
  To request a graphical display of progress during the explore phase.
  }
}
\details{
The method used in \code{\link{hybrid.explore}} is described in
Fielding, Nott and Liong (2011).
}
\value{
A list is returned to be used as input to \code{\link{hybrid.sample}}.
\item{X }{A matrix made up of the set of points
used in the final Gaussian process fit.}
\item{y }{A column of the corresponding evaluations of \code{f}.}
\item{f }{The original log-density function supplied.}
\item{maxleap }{The value of \code{maxleap} as default or supplied.}
\item{function.calls }{The number of function calls used to evaluate \code{f}.}
\item{details }{A list containing information particular to \code{\link{hybrid.sample}}.}
\item{GPfit }{Ouput from \code{\link{GProcess}} with the final
Gaussian process fit and in particular the (zero mean) approximation \code{Ef}.}
}
\references{
"Efficient MCMC Schemes for Computationally Expensive Posterior Distributions",
Fielding, Nott and Liong (2011).
}
\author{
Mark J. Fielding <mark.fielding@gmx.com>
}
\note{
The methods used in \code{\link{hybrid.explore}} and
\code{\link{hybrid.sample}}
give extensions to the work of Rasmussen (2003), as described in
Fielding, Nott and Liong (2011).

For very low acceptance rates, points included at later stages are likely
to be more useful with a fit only deteriorated by the earlier points.
In such a case a second run of \code{\link{hybrid.explore}} might be useful,
taking values for \code{X0} and \code{y0} as those output for \code{X} and
\code{y} from the first run.
}
\seealso{
\code{\link{hybrid.sample}}, \code{\link{GProcess}}, \code{\link{generateX0}}.
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
