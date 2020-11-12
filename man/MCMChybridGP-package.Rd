\name{MCMChybridGP-package}
\Rdversion{1.1}
\alias{MCMChybridGP-package}
\alias{MCMChybridGP}
\docType{package}
\title{
Hybrid MCMC for a multimodal density with derivatives replaced by Gaussian process
}
\description{
Hybrid Markov chain Monte Carlo (MCMC) to simulate from
a multimodal target distribution with derivatives unknown.
A Gaussian process fit is used to approximate derivatives.
The Package consists of an Exploratory phase,
with \code{\link{hybrid.explore}}, followed by a Sampling
phase, with \code{\link{hybrid.sample}}.
The user is to supply the log-density \code{f} of the
target distribution along with a small number of (say 10)
points to get things started.
The Sampling phase allows replacement of the true
target in high temperature chains, or complete replacement
of the target.  A full description of the method is given in
Fielding, Nott and Liong (2011).

The authors gratefully acknowledge the support &
contributions of the Singapore-Delft Water Alliance.
The research presented in this work was carried
out as part of the Multi-Objective Multi-Reservoir
Management research programme (R-264-001-272).
}
\details{
\tabular{ll}{
Package: \tab MCMChybridGP\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-09-15\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Mark James Fielding <mark.fielding@gmx.com>

Maintainer: Mark James Fielding <mark.fielding@gmx.com>
}
\references{
"Efficient MCMC Schemes for Computationally Expensive Posterior Distributions",
Fielding, Nott and Liong (2011).
}
\keyword{ package }
\seealso{
\code{\link{hybrid.explore}}, \code{\link{hybrid.sample}}
}
\examples{
    mu1 <- c(-1, -1)
    mu2 <- c(+1, +1)
    sigma.sq <- 0.1225
    ub <- c(1.5, 3)
    X0 <- generateX0(lb=c(-2,-2), ub=ub)
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
