\name{GProcess}
\Rdversion{1.1}
\alias{GProcess}
\title{
Determine a Gaussian process fit to a multivariate log-density function.
}
\description{
A function to determine a Gaussian process fit to a set of
points forming a matrix \code{X}, given a column
of corresponding values of the log-density
of a target distribution.  Returned is a (zero mean)
approximation \code{Ef} of the log-density and various
components of the Gaussian process fit as used by
\code{\link{hybrid.explore}} and \code{\link{hybrid.sample}}.
}
\usage{
GProcess(X, y, dsq_eta = NULL, request.functions = TRUE)
}
\arguments{
  \item{X}{
A matrix of at least 2 columns with rows representing the
points (nodes) for the Gaussian process.
}
  \item{y}{
A column of corresponding values of the log-density.
Each entry corresponds to the log-density evaluated
at the respective row in \code{X}.
}
  \item{dsq_eta}{
Gaussian process parameters as used in \code{\link{hybrid.explore}}.
(Optimal values will be determined when not supplied.)
}
  \item{request.functions}{
Boolean argument (with default TRUE) to request the return of function
approximations \code{Ef} (and \code{sigmaf}) of the log-density.
}
}
\value{
Returned is a list as requested consisting of:
\item{Ef }{Upon request, the Gaussian process (zero mean) approximation
of the log-density function.}
\item{sigmaf }{Upon request, a function giving the Gaussian process
approximation of the standard deviation with respect to the
log-density for a given point.}
\item{Sigma }{Covariance matrix used in generating the Gaussian process fit.}
\item{Sigma.inv }{The inverse of the Covariance matrix.}
\item{inverseOK }{Boolean flag to indicate successful calculation
of \code{Sigma.inv}.}
\item{X }{The original \code{X} matrix as supplied.}
\item{ystar }{The zero mean form of \code{y} as supplied.}
\item{dsq_eta }{Parameter values determined for the Gaussian process.}
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
    mu1 <- c(-1, -1)
    mu2 <- c(+1, +1)
    sigma.sq <- 0.1225
    X <- matrix(c(-2,-1,0,-2,0,2,0,1,2, -2,-1,-2,0,0,0,2,1,2), ncol = 2)
    f <- function(x) {
        px <- 1/4/pi/sqrt(sigma.sq) * exp(-1/2/sigma.sq *
            sum((x - mu1)^2)) + 1/4/pi/sqrt(sigma.sq) *
            exp(-1/2/sigma.sq * sum((x - mu2)^2))
        return(log(px))
    }

    y <- rep(NA, 9)
    for(i in 1:9) y[i] <- f(X[i,])

    Ef <- GProcess(X, y, request.functions = TRUE)$Ef
    Ey <- NA*y
    for(i in 1:9) Ey[i] <- Ef(X[i,])
    data.frame(X, ystar=y-mean(y), Ey)
    ## Gaussian process close to exact at points supplied.

}
