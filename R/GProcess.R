GProcess <-
function (X, y, params = NULL, request.functions = TRUE, finetune=FALSE)
{
    X <- as.matrix(X)
    n <- dim(X)[1]
    d <- dim(X)[2]
    if(d < 2) stop("Points with dimension of at least 2 assumed for X")
    mY <- mean(y)
    ystar <- y - mY
    if (is.null(params)) {
        params <- rep(1, 1 + d)
        finetune <- TRUE
    } else factr <- 1e10
    if(finetune) factr <- 1e15
    Lreject <- -1e13
    Lok <- TRUE
    Sigma <- matrix(rep(0, n*n), ncol=n)
    L <- function(params) {
        lsq <- params[1]
        eta <- params[-1]
        c.scalar <- function(xi, xj) {
            h <- abs(xi - xj)
            return(lsq * prod((1 + eta * h) * exp(-eta * h)))
        }
        c.vec <- function(x) {
            vec <- rep(NA, n)
            for (j in 1:n) vec[j] <- c.scalar(x, X[j, ])
            return(vec)
        }
        Sigma <- .C("GPcovar", as.double(X), as.integer(n), as.integer(d),
                    as.double(params), Sigma = as.double(Sigma))$Sigma
        dim(Sigma) <- c(n, n)
        cholOK <- FALSE
        try({
            R <- chol(Sigma)
            cholOK <- TRUE
        }, silent=TRUE)
        if (!cholOK) {
            val <- Lreject * runif(1, 1, 10)
            Lok <- FALSE
            return(val)
        }
        ldet.Sigma <- 2 * sum(log(diag(R)))
        Sigma.inv <- chol2inv(R)
        val <- -1/2 * ldet.Sigma - 1/2 * t(ystar) %*% Sigma.inv %*% 
            ystar
        if (is.na(val)) {
            val <- Lreject * runif(1, 1, 10)
            Lok <- FALSE
            message(paste("L(lsq, eta) NaN set to", val))
            return(val)
        }
        val <- as.double(val)
        return(val)
    }
    try({
        solveL <- optim(params, L, control = list( fnscale = -1, factr=factr ),
          method = "L-BFGS-B", lower = 1e-06)
        if (Lok) 
            params <- solveL$par
    })
    lsq <- params[1]
    eta <- params[-1]
    Sigma <- .C("GPcovar", as.double(X), as.integer(n), as.integer(d),
                as.double(params), Sigma = as.double(Sigma))$Sigma
    dim(Sigma) <- c(n, n)
    inverseOK <- FALSE
    try({
        Sigma.inv <- ginv(Sigma)
        inverseOK <- TRUE
    })
    if (!request.functions) 
        return(list(Sigma = Sigma, Sigma.inv = Sigma.inv, inverseOK = inverseOK, 
            X = X, y = y, params = params))
    c.scalar <- function(xi, xj) {
        h <- abs(xi - xj)
        val <- lsq * prod((1 + eta * h) * exp(-eta * h))
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
        return(mY + Ey_x)
    }
    if (!inverseOK) 
        sigmaf <- function(x) NA
    else sigmaf <- function(x) {
        c.vec_x <- c.vec(x)
        val <- lsq - as.double(t(c.vec_x) %*% Sigma.inv %*% c.vec_x)
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
        X = X, y = y, params = params, Ef = Ef, sigmaf = sigmaf))
}

