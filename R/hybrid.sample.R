hybrid.sample <-
function (Explore, n = 1000, replace.target = c(0, 1, 2),
    lb=NULL, ub=NULL, L = NULL, delta = NULL, nchains = NULL,
    T.mult = NULL, maxleap=NULL, r = 5, nswaps = NULL, graph = FALSE)
{
    negInf <- function(X)
      as.double(apply(X, 2, function(x) 101*min(x)-100*max(x)))
    posInf <- function(X)
      as.double(apply(X, 2, function(x) 101*max(x)-100*min(x)))
    X <- Explore$X
    y <- Explore$y
    f <- Explore$f
    d <- dim(X)[2]
    if (d < 2)
        graph <- FALSE
    nexplore <- dim(X)[1]
    if (is.null(L)) L <- Explore$details$L
    if (is.null(L)) L <- 1000
    if (is.null(delta)) Explore$details$delta
    if (is.null(delta)) delta <- 0.003
    if (is.null(nchains)) Explore$details$nchains
    if (is.null(nchains)) nchains <- 5
    if (is.null(nswaps)) nswaps <- Explore$details$nswaps
    if (is.null(nswaps)) nswaps <- choose(nchains, 2)
    if (is.null(T.mult)) T.mult <- Explore$details$nchains
    if (is.null(T.mult)) T.mult <- 1.5
    if (T.mult == 1) 
        if (replace.target[1] == 1) 
            replace.target <- 0
    Temp <- T.mult^(1:nchains - 1)
    if(is.null(lb)) lb <- Explore$details$lb
    if(is.null(ub)) ub <- Explore$details$ub
    BOUNDS <- FALSE
    if ((!is.null(lb)) || (!is.null(ub))) BOUNDS <- TRUE
    if(is.null(lb)) lb <- rep(-Inf, d)
    low <- lb
    low[!(low>-Inf)] <- negInf(X)[!(low>-Inf)]
    if(is.null(ub)) ub <- rep(Inf, d)
    upp <- ub
    upp[!(upp<Inf)] <- posInf(X)[!(upp<Inf)]
    if(is.null(lb)) lb <- rep(-Inf,d)
    if(is.null(ub)) ub <- rep(Inf,d)
    if (BOUNDS) {
      if (!is.null(lb)) 
          if (length(lb) == 1) 
              lb <- rep(lb, d)
      if (!is.null(ub)) 
          if (length(ub) == 1) 
              ub <- rep(ub, d)
      if ((length(lb) != d) || (length(lb) != d)) 
          stop("lb, ub incorrect dimension")
      if (is.na(sum(abs(lb)+abs(ub))))
          stop("NaN for lb, ub")
    }
    W <- function(p) return(sum(p^2)/2)
    if (!replace.target[1] %in% 0:2) 
        stop(paste("hybrid.sample: Choose replace.target from:\n", 
            "1: Use true target distribution in primary chain only\n", 
            "0: Use true target distribution in all chains\n", 
            "2: Replace target by Gaussian process target approximation\n"))
    select <- 1:nexplore
    for(i in nexplore:1) for(j in 1:d) {
      if ((X[i,j] < lb[j]) || (X[i,j] > ub[j]) )
        select <- select[-i]
    }
    nselect <- length(select)
    if(nselect < nchains) stop("Not enough points supplied within bounds")
    y_select <- y[select]
    X_select <- X[select,]
    indx_best <- sort(y_select, decreasing=TRUE, index.return=TRUE)$ix[1:nchains]
    x <- X_select[indx_best, ]
    fx <- y_select[indx_best]
    if(is.null(dim(x))) x <- t(x)
    fcalls <- 0
    x1 <- x[1, ]
    fx1 <- fx[1]
    digits <- rep(NA, d)
    for (j in 1:d) {
        rangeX <- IQR(X[, j])
        if (rangeX == 0) 
            digits[j] <- 4
        else digits[j] <- max(1, round(3 - log10(rangeX)))
    }
    rnd <- function(x) return(round(10^digits * x)/10^digits)
    KEY.VARS <- NULL
    key.vars <- 1:d
    for (i in 1:(length(key.vars) - 1)) for (j in (i + 1):length(key.vars)) {
        pair <- c(key.vars[i], key.vars[j])
        KEY.VARS <- rbind(KEY.VARS, pair)
    }
    if (graph) {
        if(dim(X)[2] <= 2) nrefresh.graph <- Inf else
          nrefresh.graph <- ceiling(10/nchains)
        graph.template <- .hybrid.template
        count <- key.indx <- 1
        if (is.null(graph.template)) 
            graph.template <- Explore$details$graph.template
        keyvars <- KEY.VARS[key.indx, ]
        graph.template(X, key.vars = keyvars, lb, ub)
        title(c("", "", "", "Sampling phase"), adj = 0, cex.main = 0.8)
    }
    clrs <- c("blue", rep(c("green", "turquoise", "orange", "yellow", 
        "red", "grey", "pink"), nchains))
    reds <- c("blue", rep(c("green", "turquoise", "orange", "yellow", 
        "red", "grey", "pink"), nchains))
    if (T.mult == 1) 
        Rank <- 1:nchains
    else Rank <- rank(Temp)
    if (graph) {
        keyvars <- KEY.VARS[key.indx, ]
        points(X[, keyvars[1]], X[, keyvars[2]], col = "grey", 
            pch = 20)
    }
    if (T.mult == 1) 
        1:nchains
    else Rank <- rank(Temp)
    for (chain in 1:nchains) {
        path <- Rank[chain]
        if (graph) 
            points(x[chain, keyvars[1]], x[chain, keyvars[2]], 
                pch = 19, col = clrs[path])
    }
    if(is.null(Explore$GPfit))
      Explore$GPfit <- GProcess(X, y, params = NULL, request.functions = FALSE)
    GPfit <- Explore$GPfit
    if (!GPfit$inverseOK) 
        stop("hybrid.sample: Gaussian Process inverse covariance matrix failed")
    params <- GPfit$params
    Sigma.inv <- GPfit$Sigma.inv
    ystar <- y - mean(y)
    lsq <- params[1]
    eta <- params[-1]
    c.scalar <- function(xi, xj) {
        h <- abs(xi - xj)
        val <- lsq * prod((1 + eta * h) * exp(-eta * h))
        return(val)
    }
    c.vec <- function(x) {
        vec <- rep(NA, nexplore)
        for (j in 1:nexplore) vec[j] <- c.scalar(x, X[j, ])
        return(vec)
    }
    Ef <- function(x) {
        Ey_x <- as.double(t(c.vec(x)) %*% Sigma.inv %*% ystar)
        return(Ey_x)
    }
    sigmaf <- function(x) {
        c.vec_x <- c.vec(x)
        val <- lsq - as.double(t(c.vec_x) %*% Sigma.inv %*% c.vec_x)
        if (is.na(val)) {
            warning("hybrid.sample: NaN sigmaf(", rnd(x), ")\n")
            return(Inf)
        }
        if (val < 0) {
            if (abs(val) < 0.001) 
                val <- abs(val)
            else {
                rprt <- paste("sigmaf(")
                for (i in 1:length(x)) rprt <- paste(rprt, x[i])
                rprt <- paste(rprt, ") =", val)
                stop(rprt)
            }
        }
        return(sqrt(val))
    }
    Ey <- cory <- NULL
    for(i in 1:nexplore) {
      Ey <- c(Ey, Ef(X[i,]))
      cory <- c(cory, sigmaf(X[i,])/sqrt(lsq))
    }
    if (is.null(maxleap)) maxleap <- max(cory)
    SAMP <- ySAMP <- NULL
    its <- 1
    acceptance <- NULL
    swaps.record <- swaps.record1 <- NULL
    repeat {
        if (nchains >= 3) 
            for (reps in 1:nswaps) for (index1 in sample(2:(nchains - 
                1))) {
                index2 <- index1 + 1
                if (T.mult == 1) 
                  Tindex <- 1:nchains
                else Tindex <- sort(Temp, index.return = TRUE)$ix
                swap <- Tindex[c(index1, index2)]
                T1 <- Temp[swap[1]]
                T2 <- Temp[swap[2]]
                if (replace.target[1] == 0) {
                  Efswap1 <- fx[swap[1]]
                  Efswap2 <- fx[swap[2]]
                }
                else {
                  corf1 <- sigmaf(x[swap[1], ])/sqrt(lsq)
                  corf2 <- sigmaf(x[swap[1], ])/sqrt(lsq)
                  Efswap1 <- Ef(x[swap[1], ]) - r * (corf1 - maxleap)^2/(1 - corf1)^2
                  Efswap2 <- Ef(x[swap[2], ]) - r * (corf2 - maxleap)^2/(1 - corf2)^2
                }
                alpha <- (Efswap2/T1 - Efswap1/T1 + Efswap1/T2 - Efswap2/T2)
                if(is.na(alpha)) alpha <- -Inf
                if (log(runif(1)) < alpha) {
                  Temp[swap[1]] <- T2
                  Temp[swap[2]] <- T1
                  swaps.record <- c(swaps.record, 1)
                }
                else {
                  swaps.record <- c(swaps.record, 0)
                }
            }
        if (nchains >= 2) {
            if (T.mult == 1) 
                Tindex <- 1:nchains
            else Tindex <- sort(Temp, index.return = TRUE)$ix
            swap <- Tindex[c(1, 2)]
            T1 <- Temp[swap[1]]
            if (T1 != 1) 
                stop(paste("hybrid.sample swaps: T1 =", T1))
            T2 <- Temp[swap[2]]
            if (replace.target[1] == 0) {
                fswap1 <- Efswap1 <- fx[swap[1]]
                fswap2 <- Efswap2 <- fx[swap[2]]
            }
            else if (replace.target[1] == 2) {
                corf1 <- sigmaf(x[swap[1], ])/sqrt(lsq)
                corf2 <- sigmaf(x[swap[1], ])/sqrt(lsq)
                Efswap1 <- fswap1 <- Ef(x[swap[1], ]) - r * (corf1 - maxleap)^2/(1 - corf1)^2
                Efswap2 <- fswap2 <- Ef(x[swap[2], ]) - r * (corf2 - maxleap)^2/(1 - corf2)^2
            }
            else if (replace.target[1] == 1) {
                fswap1 <- fx1
                xswap2 <- x[swap[2], ]
                corf1 <- sigmaf(x[swap[1], ])/sqrt(lsq)
                corf2 <- sigmaf(x[swap[1], ])/sqrt(lsq)
                Efswap1 <- Ef(x[swap[1], ]) - r * (corf1 - maxleap)^2/(1 - corf1)^2
                Efswap2 <- Ef(x[swap[2], ]) - r * (corf2 - maxleap)^2/(1 - corf2)^2
                date1 <- date()
                fswap2 <- f(x[swap[2], ])
                date2 <- date()
                mean.acceptance <- 0
                if (!is.null(acceptance))
                  mean.acceptance <- mean(acceptance)
                if (is.na(fswap2)) 
                  stop("f(x) NaN")
                if (fswap2 == Inf) 
                  stop(paste("Swaps: Infinite f(x) not allowed:\n", 
                    "f(", deparse(xswap2), ") =", fswap2))
                fcalls <- fcalls + 1
                if (graph) 
                  .show.time(date1, date2)
            }
            alpha <- ((fswap2 - fswap1)/1 + (Efswap1 - Efswap2)/T2)
            if(is.na(alpha)) alpha <- -Inf
            if (log(runif(1)) < alpha) {
                Temp[swap[1]] <- T2
                Temp[swap[2]] <- 1
                swaps.record1 <- c(swaps.record1, 1)
                if (replace.target[1] == 1) {
                  x1 <- xswap2
                  fx1 <- fswap2
                }
            }
            else
                swaps.record1 <- c(swaps.record1, 0)
        }
        if (T.mult == 1) 
            Rank <- 1:nchains
        else Rank <- rank(Temp)
        for (chain in 1:nchains) {
            path <- Rank[chain]
            x.C <- x[chain, ]
            p.C <- NULL
            if (replace.target[1] == 0) 
                fx.C <- fx[chain]
            else if ((replace.target[1] == 1) && (path == 1)) 
                fx.C <- fx1
            else fx.C <- Ef(x[chain, ])
            xp_Early <- NULL
            repeat {
                p.C <- rnorm(d)
                xp_Early <- matrix(c(x.C, p.C, rep(0, d)), ncol = 3)
                xOK <- FALSE
                delta_adj <- delta
                try({
                  out <- .C("Leap", as.double(X), as.integer(dim(X)[1]),
                    as.integer(d), as.double(y), as.double(x.C),
                    as.double(p.C), as.integer(L), as.double(delta),
                    as.double(GPfit$Sigma.inv), as.double(params),
                    as.double(Temp[chain]), as.double(sqrt(lsq)*maxleap),
                    as.double(low), as.double(upp),
                    as.double(rep(0,n)), as.double(rep(0,n)),
                    as.double(rep(0,d)), as.double(rep(0,d)),
                    xp_Early = as.double(xp_Early) )
                  xp_Early <- out$xp_Early
                  dim(xp_Early) <- c(d,3)
                  xOK <- TRUE
                  xOK <- TRUE
                })
                xout <- xp_Early[, 1]
                if (is.na(sum(xp_Early))) 
                  xOK <- FALSE
                if (xOK) 
                  break
            }
            out <- list(x = xp_Early[, 1], p = xp_Early[, 2])
            H_xp.C <- -fx.C + W(p.C)
            x.P <- out$x
            p.P <- out$p
            if (graph) 
              lines(c(x.C[keyvars[1]], x.P[keyvars[1]]), 
                c(x.C[keyvars[2]], x.P[keyvars[2]]), lwd = 0.5, 
                lty = 2, col = clrs[path])
            if (replace.target[1] == 2) 
              fx.P <- Ef(x.P)
            else if ((replace.target[1] == 1) && (path > 1)) 
              fx.P <- Ef(x.P)
            else {
              date1 <- date()
              fx.P <- f(x.P)
              date2 <- date()
              if (is.na(fx.P)) 
                stop("f(x) NaN")
              if (fx.P == Inf) 
                stop(paste("Leapfrog: Infinite f(x) not allowed:\n", 
                  "f(", deparse(x.P), ") =", fx.P))
              fcalls <- fcalls + 1
              if (graph) 
                if (path == 1) 
                  .show.time(date1, date2)
            }
            H_xp.P <- -fx.P + W(p.P)
            if ((replace.target[1] == 2) || ((replace.target == 1) && (path > 1))) {
              corf_x.C <- sigmaf(x.C)/sqrt(lsq)
              corf_x.P <- sigmaf(x.P)/sqrt(lsq)
              H_xp.C <- H_xp.C + r * (corf_x.C - maxleap)^2/(1 - corf_x.C)^2
              H_xp.P <- H_xp.P + r * (corf_x.P - maxleap)^2/(1 - corf_x.P)^2
            }
            alpha <- min(((-H_xp.P + H_xp.C)/Temp[chain]), 0)
            if(is.na(alpha)) alpha <- -Inf
            if (log(runif(1)) < alpha) {
              x[chain, ] <- x.P
              if (replace.target[1] == 0) 
                fx[chain] <- fx.P
              else if (replace.target[1] == 1) 
                if (path == 1) 
                  fx1 <- fx.P
              if (path == 1) {
                SAMP <- rbind(SAMP, x.P)
                ySAMP <- c(ySAMP, fx.P)
                acceptance <- c(acceptance, 1)
                if(dim(SAMP)[1] %% floor(n/10) == 0)
                  message(paste(sep="","hybrid.sample: ",
                    round(100*dim(SAMP)[1]/n),"% Complete with ",
                    round(100*mean(acceptance),1),"% acceptance of proposals"))
              }
              if (graph) 
                points(x.P[keyvars[1]], x.P[keyvars[2]], 
                  pch = 20, col = clrs[path])
            }
            else {
              if (graph)
                points(x.P[keyvars[1]], x.P[keyvars[2]], 
                  pch = "X", col = reds[path])
              if (path == 1) {
                SAMP <- rbind(SAMP, x.C)
                ySAMP <- c(ySAMP, fx.C)
                acceptance <- c(acceptance, 0)
              }
            }
            mean.acceptance <- 0
            if (!is.null(acceptance))
              mean.acceptance <- mean(acceptance)
        }
        if (graph) {
            if (count >= nrefresh.graph) {
                count <- 1
                key.indx <- key.indx + 1
                if (key.indx > dim(KEY.VARS)[1]) 
                  key.indx <- 1
                keyvars <- KEY.VARS[key.indx, ]
                graph.template(X, key.vars = keyvars, lb, ub)
                points(X[, keyvars[1]], X[, keyvars[2]], col = "grey", pch = 20)
                points(x.P[keyvars[1]], x.P[keyvars[2]], pch = 20, col = "grey")
                title(c("", "", "", "Sampling phase"), adj = 0, cex.main = 0.8)
            }
            points(x[, keyvars[1]], x[, keyvars[2]], pch = 20, col = "grey")
            points(SAMP[acceptance>0, keyvars[1]], SAMP[acceptance>0, keyvars[2]], pch = 20, col = "black")
            count <- count + 1
        }
        its <- its + 1
        if (length(acceptance) > n) 
            break
    }
    return(list(SAMP = SAMP, y = ySAMP, acceptance = acceptance, function.calls = fcalls))
}

