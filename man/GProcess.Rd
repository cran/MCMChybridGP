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
GProcess(X, y, params = NULL, request.functions = TRUE, finetune = FALSE)
}
\arguments{
  \item{X}{
A matrix of at least 2 columns with rows representing the
points (nodes) for a Gaussian process fit.
}
  \item{y}{
A column of corresponding values of the log-density.
Each entry corresponds to the log-density evaluated
at the respective row in \code{X}.
}
  \item{params}{
Gaussian process parameters as used in \code{\link{hybrid.explore}}.
(Optimal values are determined (starting from a supplied value).)
}
  \item{request.functions}{
Optional boolean argument (default TRUE) to request the return of function
approximations \code{Ef} (and \code{sigmaf}) of the log-density.
}
  \item{finetune}{
Optional boolean argument (default FALSE) to determine
fine-tuned optimal values in \code{params}.
}
}
\value{
Returned is a list as requested consisting of:
\item{Ef }{The Gaussian process approximation of the log-density function.}
\item{sigmaf }{Upon request, a function giving the Gaussian process
approximation of the standard deviation with respect to \code{y=f(x)}
for a given point.}
\item{Sigma }{Covariance matrix used in the Gaussian process fit.}
\item{Sigma.inv }{The inverse of the Covariance matrix.}
\item{inverseOK }{Boolean flag to indicate successful calculation
of \code{Sigma.inv}.}
\item{X }{The original \code{X} matrix as supplied.}
\item{y }{\code{y} as supplied. (Note calculations use \code{y-mean(y)})}
\item{params }{Parameter values for the Gaussian process fit.}
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
    data.frame(X, y, Ey)
    ## Gaussian process close to exact at points supplied.
}

