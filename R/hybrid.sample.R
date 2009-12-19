hybrid.sample <-
function (Explore = NULL, n = 1000, L = NULL, delta = NULL, nchains = NULL, 
    T.mult = NULL, cor.max=0.9, r=5, nswaps = NULL, replace.target = c(0, 1, 2),
    graph = FALSE, graph.template = NULL, key.vars = NULL, delta.increase = NULL,
    nrefresh.graph=NULL) 
{
    f <- Explore$f
    y <- Explore$y
    X <- Explore$X
    bsf <- Explore$details$bsf
    f.bsf <- Explore$details$f.bsf
    if (is.null(dim(bsf))) 
        bsf <- t(bsf)
    dsq_eta <- Explore$details$dsq_eta
    if (is.null(L)) 
        L <- Explore$details$L
    if (is.null(delta)) 
        delta <- Explore$details$delta
    if (is.null(nchains)) 
        nchains <- Explore$details$nchains
    if (is.null(nswaps)) 
        nswaps <- Explore$details$nswaps
    if (is.null(T.mult)) 
        T.mult <- Explore$details$T.mult
    if (T.mult == 1) 
        if (replace.target[1] == 1) 
            replace.target <- 0
    Temp <- T.mult^(1:nchains - 1)
    lb <- Explore$details$lb
    ub <- Explore$details$ub
    BOUNDS <- Explore$details$BOUNDS
    KEY.VARS <- Explore$details$KEY.VARS
    if (is.null(delta.increase)) 
        delta.increase <- Explore$details$delta.increase
    if (!is.logical(delta.increase)) 
        stop("delta.increase not TRUE/FALSE")
    if(is.null(nrefresh.graph))
      nrefresh.graph <- Explore$details$nrefresh.graph
    nrefresh.graph <- as.integer(nrefresh.graph)
    if(is.na(nrefresh.graph)) stop("'nrefresh.graph' gave NaN")
    if(nrefresh.graph < 0) stop("'nrefresh.graph' negative not asllowed")
    if(nrefresh.graph == 0) nrefresh.graph <- Inf
    W <- function(p) return(sum(p^2)/2)
    if (!replace.target[1] %in% 0:2) 
        stop(paste("hybrid.sample: Choose replace.target from:\n", 
            "1: Use true target distribution in primary chain only\n", 
            "0: Use true target distribution in all chains\n", 
            "2: Replace target by Gaussian process target approximation\n"))
    nbsf <- dim(bsf)[1]
    indx.bsf <- sort(f.bsf, index.return = TRUE, decreasing = TRUE)$ix
    bsf <- bsf[indx.bsf, ]
    f.bsf <- f.bsf[indx.bsf]
    if (nchains <= nbsf) {
        select <- indx.bsf[1:nchains]
        x <- bsf[select, ]
        if (is.null(dim(x))) 
            x <- t(x)
        fx <- f.bsf[select]
    }
    else {
        x <- bsf
        fx <- f.bsf
        nremain <- nchains - nbsf
        while (nremain > 0) {
            select <- sample(1:nbsf, 1)
            x <- rbind(x, bsf[select, ])
            fx <- c(fx, f.bsf[select])
            nremain <- nremain - 1
        }
    }
    cat("\n")
    cat("hybrid.sample: L =", L, "and delta =", delta, "\n")
    cat("hybrid.sample: Iterations n =", n, "with", nchains, 
        "chains: Tj =", Temp, "\n")
    d <- dim(X)[2]
    if (d < 2) 
        graph <- FALSE
    fcalls <- 0
    T1.loc <- 0
    if (dim(x)[1] == 1) 
        T1.loc <- 1
    else repeat {
        T1.loc <- T1.loc + 1
        if (Temp[T1.loc] == 1) 
            break
    }
    x1 <- x[T1.loc, ]
    fx1 <- fx[T1.loc]
    digits <- rep(NA, d)
    for (j in 1:d) {
        rangeX <- IQR(X[, j])
        if (rangeX == 0) 
            digits[j] <- 4
        else digits[j] <- max(1, round(3 - log10(rangeX)))
    }
    rnd <- function(x) return(round(10^digits * x)/10^digits)
    if (graph) {
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
        cat(sep="", "x(", chain, ", T =", Temp[chain], ") =\n")
        cat(rnd(x[chain, ]), "\n")
        if (graph) 
            points(x[chain, keyvars[1]], x[chain, keyvars[2]], 
                pch = 19, col = clrs[path])
    }
    GPfit <- GProcess(X, y, dsq_eta = dsq_eta, request.functions = FALSE)
    if (!GPfit$inverseOK) 
        stop("Gaussian Process failed")
    dsq_eta <- GPfit$dsq_eta
    Sigma.inv <- GPfit$Sigma.inv
    ystar <- y - mean(y)
    dsq <- dsq_eta[1]
    eta <- dsq_eta[-1]
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
    Ef <- function(x) {
        Ey_x <- as.double(t(c.vec(x)) %*% Sigma.inv %*% ystar)
        return(Ey_x)
    }
    sigmaf <- function(x) {
        c.vec_x <- c.vec(x)
        val <- dsq - as.double(t(c.vec_x) %*% Sigma.inv %*% c.vec_x)
        if (is.na(val)) {
            cat("hybrid.sample: NaN sigmaf(", rnd(x), ")\n")
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
    acceptance <- NULL
    SAMP <- yReport <- NULL
    its <- 1
    swaps.record <- swaps.record1 <- NULL
    repeat {
        cat("\n")
        if (nchains == 1) {
            cat("hybrid.sample: Single chain used ")
            if (replace.target[1] == 0) 
                cat("using true target\n")
            else cat("with target replaced\n")
        }
        else {
            cat("Parallel Tempering with", nchains, "chains and T =", 
                Temp, "\n")
            if (replace.target[1] == 1) 
                cat("hybrid.sample(replace.target=1): Replace target in chains with T>1\n")
            else if (replace.target[1] == 2) 
                cat("hybrid.sample(replace.target=2): Replace target in all chains\n")
            else cat("hybrid.sample: Using true target\n")
        }

        if (nchains >= 3) 
            for (reps in 1:nswaps) for (index1 in sample(2:(nchains - 1))) {
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
                  cor1 <- sigmaf(x[swap[1], ])/sqrt(dsq)
                  cor2 <- sigmaf(x[swap[1], ])/sqrt(dsq)
                  Efswap1 <- Ef(x[swap[1], ]) - r*(cor1-cor.max)^2/(1-cor1)^2
                  Efswap2 <- Ef(x[swap[2], ]) - r*(cor2-cor.max)^2/(1-cor2)^2
                }
                alpha <- (Efswap2/T1 - Efswap1/T1 + Efswap1/T2 - Efswap2/T2)
                if (log(runif(1)) < alpha) {
                  Temp[swap[1]] <- T2
                  Temp[swap[2]] <- T1
                  cat("hybrid.sample: Swap T =", T1, "with", T2, "\n")
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
                cor1 <- sigmaf(x[swap[1], ])/sqrt(dsq)
                cor2 <- sigmaf(x[swap[1], ])/sqrt(dsq)
                Efswap1 <- fswap1 <- Ef(x[swap[1], ]) - r*(cor1-cor.max)^2/(1-cor1)^2
                Efswap2 <- fswap2 <- Ef(x[swap[2], ]) - r*(cor2-cor.max)^2/(1-cor2)^2
            }
            else if (replace.target[1] == 1) {
                fswap1 <- fx1
                xswap2 <- x[swap[2], ]
                cor1 <- sigmaf(x[swap[1], ])/sqrt(dsq)
                cor2 <- sigmaf(x[swap[1], ])/sqrt(dsq)
                Efswap1 <- Ef(x[swap[1], ]) - r*(cor1-cor.max)^2/(1-cor1)^2
                Efswap2 <- Ef(x[swap[2], ]) - r*(cor2-cor.max)^2/(1-cor2)^2
                date1 <- date()
                fswap2 <- f(x[swap[2], ])
                date2 <- date()
                if (is.null(acceptance)) 
                  mean.acceptance <- NA
                else mean.acceptance <- mean(acceptance)
                cat(sep = "", "Sampling x(T=1) of ", nchains, " chains: ", its, " / ",
                  n, " iterations with ", round(10000 * mean.acceptance)/100, "% acceptance\n")
                if (is.na(fswap2)) 
                  stop("f(x) NaN")
                if (fswap2 == Inf) 
                  stop(paste("Swaps: Infinite f(x) not allowed:\n", 
                    "f(", deparse(xswap2), ") =", fswap2))
                cat("Swaps: f(x) =", fswap2, "\n")
                fcalls <- fcalls + 1
                if (graph) 
                  .show.time(date1, date2, graph)
            }
            alpha <- ((fswap2 - fswap1)/1 + (Efswap1 - Efswap2)/T2)
            if (log(runif(1)) < alpha) {
                cat("hybrid.sample: Swap T = 1 with", T2, "\n")
                Temp[swap[1]] <- T2
                Temp[swap[2]] <- 1
                swaps.record1 <- c(swaps.record1, 1)
                if (replace.target[1] == 1) {
                  x1 <- xswap2
                  fx1 <- fswap2
                }
            }
            else {
                cat("hybrid.sample: Reject swap between T = 1 &", T2, "\n")
                swaps.record1 <- c(swaps.record1, 0)
            }
        }
        if (T.mult == 1) 
            Rank <- 1:nchains
        else Rank <- rank(Temp)
        for (chain in 1:nchains) {
            path <- Rank[chain]
            x.C <- x[chain, ]
            p.C <- NULL
            cat("\n")
            cat(sep="", "Leapfrog from x(", chain, ", T = ", Temp[chain], "):\n")
            cat("", rnd(x.C), "\n")
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
                if (!delta.increase) 
                  delta_adj <- delta
                else delta_adj <- delta * Temp[chain]
                try({
                  xp_Early <- .Call("Leapfrog", X, y, x.C, p.C, 
                    L, delta_adj, GPfit$Sigma.inv, dsq_eta, Temp[chain], 
                    Inf, PACKAGE = "MCMChybridGP")
                  xOK <- TRUE
                })
                xout <- xp_Early[, 1]
                if (is.na(sum(xp_Early))) 
                  xOK <- FALSE
                if (xOK) 
                  break
                cat("hybrid.sample: Proposed x failed:\n", xout, "\n")
            }
            out <- list(x = xp_Early[, 1], p = xp_Early[, 2])
            H_xp.C <- -fx.C + W(p.C)
            x.P <- out$x
            p.P <- out$p
            OUTSIDE_BOUNDS <- FALSE
            if (BOUNDS) 
                for (i in 1:d) {
                  if (x.P[i] < lb[i]) 
                    OUTSIDE_BOUNDS <- TRUE
                  if (x.P[i] > ub[i]) 
                    OUTSIDE_BOUNDS <- TRUE
                }
            if (OUTSIDE_BOUNDS) {
                if (path == 1) {
                  cat(sep="", "Proposed x(", chain, ", T =", Temp[chain], ") outside bounds with:\n")
                  cat("", rnd(x.P), "\n")
                  SAMP <- rbind(SAMP, x.C)
                  yReport <- c(yReport, fx.C)
                  acceptance <- c(acceptance, 0)
                  if (graph) 
                    lines(c(x.C[keyvars[1]], x.P[keyvars[1]]), 
                      c(x.C[keyvars[2]], x.P[keyvars[2]]), lwd = 0.5, 
                      lty = 6, col = clrs[path])
                }
            }
            else {
                cat(sep="", "Proposed x(", chain, ", T =", Temp[chain], "):\n")
                cat("", rnd(x.P), "\n")
                if (graph) 
                  lines(c(x.C[keyvars[1]], x.P[keyvars[1]]), 
                    c(x.C[keyvars[2]], x.P[keyvars[2]]), lwd = 0.5, 
                    lty = 2, col = clrs[path])
                if (replace.target[1] == 2) 
                  fx.P <- Ef(x.P)
                else if ((replace.target[1] == 1) && (path > 1)) 
                  fx.P <- Ef(x.P)
                else {
                  cat("hybrid.sample: Calculating true f(x) ..\n")
                  date1 <- date()
                  fx.P <- f(x.P)
                  date2 <- date()
                  if (is.na(fx.P)) 
                    stop("f(x) NaN")
                  if (fx.P == Inf) 
                    stop(paste("Leapfrog: Infinite f(x) not allowed:\n", 
                      "f(", deparse(x.P), ") =", fx.P))
                  fcalls <- fcalls + 1
                  cat("f(x) =", fx.P, "with leap from f(x) =", 
                    fx.C, "\n")
                  if (graph) 
                    if (path == 1) 
                      .show.time(date1, date2, graph)
                }
                H_xp.P <- -fx.P + W(p.P)
                if ((replace.target[1] == 2) || ((replace.target == 1) && (path > 1))) {
                  cor_x.C <- sigmaf(x.C)/sqrt(dsq)
                  cor_x.P <- sigmaf(x.P)/sqrt(dsq)
                  H_xp.C <- H_xp.C + r*(cor_x.C-cor.max)^2/(1-cor_x.C)^2
                  H_xp.P <- H_xp.P + r*(cor_x.P-cor.max)^2/(1-cor_x.P)^2
                }
                alpha <- min(((-H_xp.P + H_xp.C)/Temp[chain]), 
                  0)
                if (log(runif(1)) < alpha) {
                  x[chain, ] <- x.P
                  if (replace.target[1] == 0) 
                    fx[chain] <- fx.P
                  else if (replace.target[1] == 1) 
                    if (path == 1) 
                      fx1 <- fx.P
                  cat("ACCEPT x.P(", chain, ")\n")
                  if (path == 1) {
                    acceptance <- c(acceptance, 1)
                    SAMP <- rbind(SAMP, x.P)
                    yReport <- c(yReport, fx.P)
                  }
                  if (graph) 
                    points(x.P[keyvars[1]], x.P[keyvars[2]], 
                      pch = 20, col = clrs[path])
                }
                else {
                  if (graph) 
                    points(x.P[keyvars[1]], x.P[keyvars[2]], 
                      pch = "X", col = reds[path])
                  cat("REJECT x.P(", chain, ")\n")
                  if (path == 1) {
                    SAMP <- rbind(SAMP, x.C)
                    yReport <- c(yReport, fx.C)
                    acceptance <- c(acceptance, 0)
                  }
                }
                if (is.null(acceptance)) 
                  mean.acceptance <- NA
                else mean.acceptance <- mean(acceptance)
                cat(sep = "", "Sampling x(T=1) of ", nchains, " chains: ", its, " / ",
                  n, " iterations with ", round(10000 * mean.acceptance)/100, "% acceptance\n")
                cat("Leapfrog used: L =", L, ", delta =", delta_adj, "\n")
                if (nchains >= 2) 
                  cat(sep = "", "Swaps (T=1 with T>1) ", round(10000 * 
                    mean(swaps.record1))/100, "%\n")
                if (nchains >= 3) 
                  cat(sep = "", "Swaps (T>1) ", round(10000 * 
                    mean(swaps.record))/100, "%\n")
            }
        }
        if (graph) {
            if (count >= nrefresh.graph) {
                count <- 1
                key.indx <- key.indx + 1
                if (key.indx > dim(KEY.VARS)[1]) 
                  key.indx <- 1
                keyvars <- KEY.VARS[key.indx, ]
                graph.template(X, key.vars = keyvars, lb, ub)
                points(x.P[keyvars[1]], x.P[keyvars[2]], pch = 20, 
                  col = "grey")
                title(c("", "", "", "Sampling phase"), adj = 0, 
                  cex.main = 0.8)
            }
            points(x[, keyvars[1]], x[, keyvars[2]], pch = 20, 
                col = "grey")
            points(SAMP[, keyvars[1]], SAMP[, keyvars[2]], pch = 20, 
                col = "black")
            count <- count + 1
        }
        its <- its + 1
        assign("hybrid.sample.out", list(SAMP = SAMP, y = yReport, 
            fcalls = fcalls, acceptance = acceptance, swap.acceptance = list(Primary = swaps.record1, 
                Other = swaps.record)), envir = .GlobalEnv)
        if (its > n) 
            break
    }
    return(list(SAMP = SAMP, y = yReport, fcalls = fcalls, acceptance = acceptance, 
        swap.acceptance = list(Primary = swaps.record1, Other = swaps.record)))
}

