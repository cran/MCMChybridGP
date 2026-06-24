GProcess <-
function (X, y, params = NULL, request.functions = TRUE, finetune = FALSE)
{
    X <- as.matrix(X)
    n <- nrow(X)
    d <- ncol(X)

    if (d < 2)
        stop("Points with dimension of at least 2 assumed for X")

    mY    <- mean(y)
    ystar <- y - mY

    if (is.null(params)) {
        params   <- rep(1, 1 + d)
        finetune <- TRUE
    } else {
        factr <- 1e10
    }

    if (finetune)
        factr <- 1e15

    Lreject <- -1e13
    Lok     <- TRUE

    # --------------------------------------------------
    # Log-likelihood for given params
    # --------------------------------------------------
    L <- function(params) {

        lsq <- params[1]
        eta <- params[-1]

        # covariance matrix via Rcpp
        Sigma <- GPcovar(X, params)

        cholOK <- FALSE
        val    <- NA_real_

        # Cholesky for log-det and inverse
        try({
            R      <- chol(Sigma)
            cholOK <- TRUE
        }, silent = TRUE)

        if (!cholOK) {
            val <- Lreject * runif(1, 1, 10)
            Lok <- FALSE
            return(val)
        }

        ldet.Sigma <- 2 * sum(log(diag(R)))
        Sigma.inv  <- chol2inv(R)

        val <- -0.5 * ldet.Sigma - 0.5 * drop(t(ystar) %*% Sigma.inv %*% ystar)

        if (is.na(val)) {
            val <- Lreject * runif(1, 1, 10)
            Lok <- FALSE
            message(paste("L(lsq, eta) NaN set to", val))
            return(val)
        }

        as.double(val)
    }

    # --------------------------------------------------
    # Optimise hyperparameters
    # --------------------------------------------------
    try({
        solveL <- optim(
            params, L,
            control = list(fnscale = -1, factr = factr),
            method  = "L-BFGS-B",
            lower   = 1e-06
        )
        if (Lok)
            params <- solveL$par
    })

    lsq <- params[1]
    eta <- params[-1]

    # Final covariance and inverse
    Sigma <- GPcovar(X, params)

    inverseOK <- FALSE
    Sigma.inv <- NULL

    try({
        Sigma.inv <- MASS::ginv(Sigma)
        inverseOK <- TRUE
    })

    if (!request.functions)
        return(list(
            Sigma     = Sigma,
            Sigma.inv = Sigma.inv,
            inverseOK = inverseOK,
            X         = X,
            y         = y,
            params    = params
        ))

    # --------------------------------------------------
    # Covariance helpers
    # --------------------------------------------------
    c.scalar <- function(xi, xj) {
        h   <- abs(xi - xj)
        val <- lsq * prod((1 + eta * h) * exp(-eta * h))
        val
    }

    c.vec <- function(x) {
        vec <- numeric(nrow(X))
        for (j in seq_len(nrow(X)))
            vec[j] <- c.scalar(x, X[j, ])
        vec
    }

    # --------------------------------------------------
    # Predictive mean
    # --------------------------------------------------
    if (!inverseOK) {
        Ef <- function(x) NA
    } else {
        Ef <- function(x) {
            Ey_x <- as.double(t(c.vec(x)) %*% Sigma.inv %*% ystar)
            mY + Ey_x
        }
    }

    # --------------------------------------------------
    # Predictive standard deviation
    # --------------------------------------------------
    if (!inverseOK) {
        sigmaf <- function(x) NA
    } else {
        sigmaf <- function(x) {
            c.vec_x <- c.vec(x)
            val     <- lsq - as.double(t(c.vec_x) %*% Sigma.inv %*% c.vec_x)

            if (is.na(val))
                return(Inf)

            if (val < 0) {
                if (val > -0.001) {
                    val <- -val
                } else {
                    rprt <- paste0(
                        "sigmaf(",
                        paste(x, collapse = " "),
                        ") = ", val
                    )
                    stop(rprt)
                }
            }

            sqrt(val)
        }
    }

    list(
        Sigma     = Sigma,
        Sigma.inv = Sigma.inv,
        inverseOK = inverseOK,
        X         = X,
        y         = y,
        params    = params,
        Ef        = Ef,
        sigmaf    = sigmaf
    )
}

