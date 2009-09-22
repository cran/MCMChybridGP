generateX0 <-
function (lb, ub, n0 = 10, npool = 100, ranX = NULL) 
{
    if(!is.null(dim(lb))) stop("lb in error")
    if(!is.null(dim(ub))) stop("ub in error")
    d <- length(lb)
    if(length(ub) != d) stop("Dimension mismatch for lb, ub")
    if(sum(abs(lb)+abs(ub)) == Inf) stop("Finite lb, ub required")
    if (is.null(ranX)) 
        ranX <- function() {
            X <- matrix(rep(NA, n0 * d), ncol = d)
            for (i in 1:d) {
                Delta <- (ub[i] - lb[i])/2/n0
                X[, i] <- sample(seq(lb[i] + Delta, ub[i] - Delta, 
                  length.out = n0))
            }
            return(X)
        }
    X <- Xbest <- as.matrix(ranX())
    if (n0 < 2) {
        return(Xbest)
    }
    dmax <- -Inf
    DIST <- function(x, y) sum(((x - y)/(ub - lb))^2)
    for (pool in 1:npool) {
        dmin <- Inf
        for (i in 1:(n0 - 1)) for (j in (i + 1):n0) if (DIST(X[i, 
            ], X[j, ]) < dmin) 
            dmin <- DIST(X[i, ], X[j, ])
        if (dmin > dmax) {
            dmax <- dmin
            Xbest <- X
        }
        if (pool < npool) 
            X <- as.matrix(ranX())
    }
    return(Xbest)
}

