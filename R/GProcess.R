GProcess <-
function (X, y, dsq_eta = NULL, request.functions = TRUE) 
{
    X <- as.matrix(X)
    n <- dim(X)[1]
    d <- dim(X)[2]
    if(d < 2) stop("Points with dimension of at least 2 assumed for X")
    ystar <- y - mean(y)
    if (is.null(dsq_eta)) 
        dsq_eta <- rep(1, 1 + d)
    Lreject <- -1e6
    if (exists(".MCMChybridGP.Lreject")) Lreject <- .MCMChybridGP.Lreject
    .MCMChybridGP.Lreject <- Lreject
    assign(".MCMChybridGP.Lreject", Lreject, envir=.GlobalEnv)
    Lok <- TRUE
    L <- function(dsq_eta) {
        dsq <- dsq_eta[1]
        eta <- dsq_eta[-1]
        c.scalar <- function(xi, xj) {
            h <- abs(xi - xj)
            return(dsq * prod((1 + eta * h) * exp(-eta * h)))
        }
        c.vec <- function(x) {
            vec <- rep(NA, n)
            for (j in 1:n) vec[j] <- c.scalar(x, X[j, ])
            return(vec)
        }
        Sigma <- .Call("calcSigma", X, dsq_eta, PACKAGE = "MCMChybridGP")
        cholOK <- FALSE
        try({
            R <- chol(Sigma)
            cholOK <- TRUE
        })
        if (!cholOK) {
            val <- .MCMChybridGP.Lreject * runif(1, 1, 10)
            cat("GProcess: L(dsq, eta): chol(Sigma) failed. Returning", 
                val, "\n")
            Lok <- FALSE
            return(val)
        }
        ldet.Sigma <- 2 * sum(log(diag(R)))
        Sigma.inv <- chol2inv(R)
        val <- -1/2 * ldet.Sigma - 1/2 * t(ystar) %*% Sigma.inv %*% 
            ystar
        if (is.na(val)) {
            val <- .MCMChybridGP.Lreject * runif(1, 1, 10)
            Lok <- FALSE
            cat("L(dsq, eta) NaN set to", val, "\n")
            return(val)
        }
        val <- as.double(val)
        if (val < .MCMChybridGP.Lreject) 
            assign(".MCMChybridGP.Lreject", val, envir = .GlobalEnv)
        return(val)
    }
    cat("GProcess: Maximum likelihood for delta.sq and eta..\n")
    cat("initial L(", dsq_eta, ") =", L(dsq_eta), "\n")
    try({
        solveL <- optim(dsq_eta, L, control = list(factr = 1e+15, 
            fnscale = -1), method = "L-BFGS-B", lower = 1e-06)
        if (Lok) 
            dsq_eta <- solveL$par
        cat("optim => L(", round(1e+05 * dsq_eta)/1e+05, ") =", 
            as.double(solveL$value), "\n")
    })
    dsq <- dsq_eta[1]
    eta <- dsq_eta[-1]
    Sigma <- .Call("calcSigma", X, dsq_eta, PACKAGE = "MCMChybridGP")
    inverseOK <- FALSE
    try({
        Sigma.inv <- ginv(Sigma)
        inverseOK <- TRUE
    })
    if (!request.functions) 
        return(list(Sigma = Sigma, Sigma.inv = Sigma.inv, inverseOK = inverseOK, 
            X = X, ystar = ystar, dsq_eta = dsq_eta))
    c.scalar <- function(xi, xj) {
        h <- abs(xi - xj)
        val <- dsq * prod((1 + eta * h) * exp(-eta * h))
        return(val)
    }
    c.vec <- function(x) {
        vec <- rep(NA, dim(X)[1])
        for (j in 1:dim(X)[1]) vec[j] <- c.scalar(x, X[j, ])
        return(vec)
    }
    if (!inverseOK) 
        Ef <- function(x) NA
    else Ef <- function(x) {
        Ey_x <- as.double(t(c.vec(x)) %*% Sigma.inv %*% ystar)
        return(Ey_x)
    }
    if (!inverseOK) 
        sigmaf <- function(x) NA
    else sigmaf <- function(x) {
        c.vec_x <- c.vec(x)
        val <- dsq - as.double(t(c.vec_x) %*% Sigma.inv %*% c.vec_x)
        if (is.na(val)) 
            return(Inf)
        if (val < 0) {
            if (val > -0.001) 
                val <- -val
            else {
                rprt <- paste("sigmaf(")
                for (i in 1:length(x)) rprt <- paste(rprt, x[i])
                rprt <- paste(rprt, ") =", val)
                stop(rprt)
            }
        }
        return(sqrt(val))
    }
    return(list(Sigma = Sigma, Sigma.inv = Sigma.inv, inverseOK = inverseOK, 
        X = X, ystar = ystar, dsq_eta = dsq_eta, Ef = Ef, sigmaf = sigmaf))
}

