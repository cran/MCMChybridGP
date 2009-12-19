hybrid.explore <-
function (f, X0, n = 200, L = 1000, delta = 0.001, nchains = 5, T.mult = 1.5,
    lb = NULL, ub = NULL, cor.max = 0.9, delta.increase = FALSE,
    nswaps = choose(nchains,2), spread = 1, graph = FALSE,
    graph.template = NULL, key.vars = 1:dim(X)[2],
    nrefresh.graph = ceiling(20/nchains)) 
{
    X <- as.matrix(X0)
    n0 <- dim(X)[1]
    d <- dim(X)[2]
    if(d < 2) stop("Points with dimension of at least 2 assumed for X")
    W <- function(p) return(sum(p^2)/2)
    sd_y <- data.frame(sd = NA, cor = NA)
    if (d < 2) 
        graph <- FALSE
    if (n0 > 1) 
        for (i in 2:n0) for (j in 1:(i - 1)) if (sum(abs(X[i, ] - X[j, ])) == 0) 
            stop(paste(sep = "", "Duplicate rows (", i, " & ", 
                j, ") in X not allowed"))
    if (!is.logical(delta.increase)) 
        stop("delta.increase not TRUE/FALSE")
    BOUNDS <- TRUE
    if (is.null(lb)) 
        if (is.null(ub)) 
            BOUNDS <- FALSE
    if (BOUNDS) 
        if (is.null(ub) || is.null(ub)) 
            stop("Both or neither bounds must be supplied")
    if (BOUNDS) {
        if (!is.null(lb)) 
            if (length(lb) == 1) 
                lb <- rep(lb, d)
        if (!is.null(ub)) 
            if (length(ub) == 1) 
                ub <- rep(ub, d)
        if ((length(lb) != d) || (length(lb) != d)) 
            stop("lb, ub incorrect dimension")
        if (is.na(sum(lb + ub))) 
            stop("hybrid.explore: NaN lb, ub")
        xOK <- TRUE
        for (i in 1:n0) {
            for (j in 1:d) if ((X[i, j] <= lb[j]) || (X[i, j] >= ub[j])) 
                xOK <- FALSE
        }
        if (!xOK) 
            stop("hybrid.explore: X not strictly within bounds")
        cat("hybrid.explore: Placing bounds lb =\n")
        cat(lb, "\n")
        cat("hybrid.explore: Placing bounds ub =\n")
        cat(ub, "\n")
    }
    digits <- rep(NA, d)
    for (j in 1:d) {
        rangeX <- IQR(X[, j])
        if (rangeX == 0) 
            digits[j] <- 4
        else digits[j] <- max(1, round(3 - log10(rangeX)))
    }
    rnd <- function(x) return(round(10^digits * x)/10^digits)
    KEY.VARS <- NULL
    if (graph) {
        key.vars <- as.integer(key.vars)
        if (is.na(sum(key.vars))) 
            stop("key.vars NaN")
        for (i in 1:(length(key.vars) - 1)) for (j in (i + 1):length(key.vars)) {
            if (key.vars[i] == key.vars[j]) 
                stop("Duplicate key.vars not allowed")
            pair <- c(key.vars[i], key.vars[j])
            KEY.VARS <- rbind(KEY.VARS, pair)
        }
        if (is.null(graph.template)) 
            graph.template <- .hybrid.template
        key.indx <- 1
        keyVars <- KEY.VARS[key.indx, ]
    }
    cat(sep = "", "hybrid.explore: ", n0, " initial ", d, "D points supplied for ", 
        n, " Exploratory points\n")
    cat("Leapfrogs: L =", L, "and delta =", delta, "\n")
    if (T.mult < 1) 
        stop("T.mult < 1 not allowed")
    swaps.record <- NULL
    y <- rep(NA, n0)
    Ey <- rep(0, n0)
    Temp <- T.mult^(1:nchains - 1)
    cat("Parallel Tempering:", nchains, "chains with T:", Temp, 
        "\n")
    date1 <- date2 <- date()
    if (graph) {
        plot.new()
        title(c("", "", "", "", "", "", "", "", "", "", "", 
            "", "", "Evaluate initial points .."))
    }
    for (j in 1:n0) {
        functionOK <- FALSE
        try({
            date1 <- date()
            y[j] <- f(X[j, ])
            date2 <- date()
            cat("\n")
            cat(sep="", "f(x", j, ") = ", y[j], " for x", j, ":\n")
            cat("", X[j, ], "\n")
            .show.time(date1, date2, graph=FALSE)
            functionOK <- TRUE
        })
        if (is.na(y[j])) 
            stop("f(x) NaN")
        if (y[j] == Inf) 
            stop("Infinite f(x) not allowed")
    }
    fcalls <- n0
    x <- matrix(rep(NA, nchains * d), ncol = d)
    fx <- rep(NA, nchains)
    target.bsf <- f.bsf <- rep(NA, nchains)
    select <- sample(1:n0, nchains)
    x <- X[select, ]
    if (is.null(dim(x))) 
        x <- t(x)
    bsf <- x
    fx <- f.bsf <- y[select]
    if (graph) {
        graph.template(X, key.vars = keyVars, lb, ub)
        points(X[, keyVars[1]], X[, keyVars[2]], pch = 20, col = "black")
        title(c("", "", "", "Exploratory phase"), adj = 0, cex.main = 0.8)
        clrs <- c("blue", rep(c("green", "turquoise", "orange", 
            "yellow", "red", "grey", "pink"), nchains))
        reds <- c("blue", rep(c("green", "turquoise", "orange", 
            "yellow", "red", "grey", "pink"), nchains))
    }
    acceptance <- NULL
    its <- 1
    nacc <- 0
    dsq_eta <- rep(1, 1 + d)
    refresh_count <- 1
    repeat {
        cat("\n")
        GPfit <- GProcess(X, y, dsq_eta = dsq_eta, request.functions = FALSE)
        Sigma.inv <- GPfit$Sigma.inv
        dsq_eta <- GPfit$dsq_eta
        dsq <- dsq_eta[1]
        eta <- dsq_eta[-1]
        if (cor.max < 1) {
            maxsig <- cor.max * sqrt(dsq)
            sfactor <- sqrt(dsq)
        }
        else {
            maxsig <- cor.max
            sfactor <- 1
        }
        c.scalar <- function(xi, xj) {
            h <- abs(xi - xj)
            val <- dsq * prod((1 + eta * h) * exp(-eta * h))
            return(val)
        }
        c.vec <- function(x) {
            vec <- rep(NA, dim(X)[1])
            for (j in 1:dim(X)[1]) vec[j] <- c.scalar(x, X[j, 
                ])
            return(vec)
        }
        sigmaf <- function(x) {
            c.vec_x <- c.vec(x)
            val <- dsq - as.double(t(c.vec_x) %*% Sigma.inv %*% 
                c.vec_x)
            if (is.na(val)) 
                return(Inf)
            if (val < 0) {
                if (val > -0.001) 
                  val <- -val
                else {
                  xrep <- (x)
                  rprt <- paste("sigmaf(")
                  for (i in 1:length(xrep)) rprt <- paste(rprt, 
                    xrep[i])
                  rprt <- paste(rprt, ") =", val)
                  stop(rprt)
                }
            }
            return(sqrt(val))
        }
        if (its == 1) 
            for (i in 1:nchains) target.bsf[i] <- fx[i] + spread * 
                sigmaf(x[i, ])/sfactor
        if (nchains >= 2) {
            for (reps in 1:nswaps) for (index1 in sample(1:nchains)) for (index2 in sample((1:nchains)[-index1])) {
                if (T.mult == 1) 
                  Tindex <- 1:nchains
                else Tindex <- sort(Temp, index.return = TRUE)$ix
                swap <- Tindex[c(index1, index2)]
                T1 <- Temp[swap[1]]
                T2 <- Temp[swap[2]]
                fswap1 <- fx[swap[1]] + spread * sigmaf(x[swap[1]])/sfactor
                fswap2 <- fx[swap[2]] + spread * sigmaf(x[swap[2]])/sfactor
                alpha <- (fswap2/T1 - fswap1/T1 + fswap1/T2 - 
                  fswap2/T2)
                if (log(runif(1)) < alpha) {
                  Temp[swap[1]] <- T2
                  Temp[swap[2]] <- T1
                  swaps.record <- c(swaps.record, 1)
                  cat("hybrid.explore: Swap T =", T1, 
                    "with", T2, "\n")
                }
                else {
                  swaps.record <- c(swaps.record, 0)
                }
            }
        }
        Xnext <- X
        ynext <- y
        if (T.mult == 1) 
            Rank <- 1:nchains
        else Rank <- rank(Temp)
        for (chain in 1:nchains) {
          path <- Rank[chain]
          x.C <- as.double(x[chain, ])
          p.C <- NULL
          sigmaf_x.C <- sigmaf(x.C)
          cat("\n")
          cat(sep="", "Leapfrog from x(", chain, ", T = ", Temp[chain], "):\n")
          cat("", rnd(x.C), "\n")
          fx.C <- fx[chain]
          xp_Early <- NULL
          count <- 0
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
                  maxsig, PACKAGE = "MCMChybridGP")
                xOK <- TRUE
              })
              if (is.na(sum(xp_Early))) 
                xOK <- FALSE
              if (xOK) 
                break
              cat("hybrid.explore: Proposed x failed\n")
              count <- count + 1
              if (count == 10) {
                for (remv in 1:nchains) {
                  X <- X[-length(y), ]
                  y <- y[-length(y)]
                }
                GPfit <- GProcess(X, y, dsq_eta = dsq_eta, 
                  request.functions = TRUE)
                Sigma.inv <- GPfit$Sigma.inv
                dsq_eta <- GPfit$dsq_eta
                dsq <- dsq_eta[1]
                eta <- dsq_eta[-1]
                if (cor.max < 1) {
                  maxsig <- cor.max * sqrt(dsq)
                  sfactor <- sqrt(dsq)
                }
                sigmaf <- GPfit$sigmaf
                if (its == 1) 
                  for (i in 1:nchains) target.bsf[i] <- fx[i] + 
                    spread * sigmaf(x[i, ])/sfactor
                Xnext <- X
                ynext <- y
                count <- 0
              }
            }
            out <- list(x = xp_Early[, 1], p = xp_Early[, 2], 
                Early = xp_Early[1, 3])
            x.P <- out$x
            OUTSIDE_BOUNDS <- FALSE
            if (BOUNDS) 
                for (i in 1:d) {
                  if (x.P[i] < lb[i]) 
                    OUTSIDE_BOUNDS <- TRUE
                  if (x.P[i] > ub[i]) 
                    OUTSIDE_BOUNDS <- TRUE
                }
            if (OUTSIDE_BOUNDS) {
                cat(sep="", "Proposed x(", chain, ", T =", Temp[chain], ") outside bounds with:\n")
                cat("", rnd(x.P), "\n")
                if (graph) 
                  lines(c(x.C[keyVars[1]], x.P[keyVars[1]]), 
                    c(x.C[keyVars[2]], x.P[keyVars[2]]), lwd = 0.2, 
                    lty = 6, col = clrs[path])
            }
            else {
                p.P <- out$p
                cat(sep="", "Proposed x(", chain, ", T =", Temp[chain], "):\n")
                cat("", rnd(x.P), "\n")
                sigmaf_x.P <- sigmaf(x.P)
                if (out$Early < 0) {
                  cat("Early stop: cor(f(x)) =", sigmaf_x.P/sfactor, ">",
                    cor.max, "( with T =", Temp[chain], ")\n")
                  if (sigmaf_x.P < maxsig) 
                    stop(paste("Leapfrog.C gave sigmaf =", -out$Early, 
                      "but sigmaf(x) =", sigmaf_x.P))
                }
                else {
                  cat("Full Leapfrog: cor(f(x)) =", sigmaf_x.P/sfactor, "<", cor.max, "\n")
                  if (sigmaf_x.P > maxsig)
                    stop(paste("Caution: Leapfrog gave sigmaf =", out$Early, 
                      "when sigmaf(x) =", sigmaf_x.P))
                }
                if (graph) 
                  lines(c(x.C[keyVars[1]], x.P[keyVars[1]]), 
                    c(x.C[keyVars[2]], x.P[keyVars[2]]), lwd = 0.2, 
                    lty = 2, col = clrs[path])
                cat("Leapfrog: Calculate true f(x) for proposal.\n")
                date1 <- date()
                fx.P <- f(x.P)
                date2 <- date()
                if (is.na(fx.P)) 
                  stop("f(x) NaN")
                if (fx.P == Inf) 
                  stop("Infinite f(x) not allowed")
                fcalls <- fcalls + 1
                cat("f(x) =", fx.P, "leap from f(x)", fx.C, "\n")
                .show.time(date1, date2, graph)
                H_xp.C <- -fx.C - spread * sigmaf_x.C/sfactor + 
                  W(p.C)
                H_xp.P <- -fx.P - spread * sigmaf_x.P/sfactor + 
                  W(p.P)
                if (graph) 
                  if (sigmaf_x.P > maxsig) 
                    points(x.P[keyVars[1]], x.P[keyVars[2]], 
                      pch = 1, cex = 1.5, col = clrs[path])
                alpha <- min(((-H_xp.P + H_xp.C)/Temp[chain]), 
                  0)
                if (log(runif(1)) < alpha) {
                  sd_y <- rbind(sd_y, c(sigmaf_x.P, sigmaf_x.P/sqrt(dsq)))
                  x[chain, ] <- x.P
                  fx[chain] <- fx.P
                  if (fx.P + sigmaf_x.P/sfactor > target.bsf[path]) {
                    bsf[path, ] <- x.P
                    target.bsf[path] <- fx.P + spread * sigmaf_x.P/sfactor
                    f.bsf[path] <- fx.P
                  }
                  Xnext <- rbind(Xnext, x.P)
                  ynext <- c(ynext, fx.P)
                  acceptance <- c(acceptance, 1)
                  Ey <- c(Ey, sigmaf_x.P)
                  cat(sep="", "ACCEPT x(", chain, ", T = ", Temp[chain], "):\n")
                  cat("", rnd(x.P), "\n")
                  if (graph) 
                    points(x.P[keyVars[1]], x.P[keyVars[2]], 
                      pch = 20, col = clrs[path])
                }
                else {
                  cat(sep="", "REJECT x.P (", chain, "):\n")
                  cat("", rnd(x.P), "\n")
                  acceptance <- c(acceptance, 0)
                  if (graph) 
                    points(x.P[keyVars[1]], x.P[keyVars[2]], 
                      pch = "X", col = reds[path])
                }
                if (is.null(acceptance)) 
                  mean.acceptance <- NA
                else mean.acceptance <- mean(acceptance)
                cat(sep = "", "hybrid.explore: ", dim(Xnext)[1], 
                  " / ", n, " complete with ", round(10000 * mean.acceptance)/100, 
                  "% acceptance\n")
                cat("Leapfrog used: L =", L, ", delta =", delta_adj, "\n")
                if (nchains >= 2) 
                  cat("Swaps:", round(10000 * mean(swaps.record))/100,
                    "% acceptance\n")
                if (dim(Xnext)[1] >= n) 
                  break
            }
        }
        X <- Xnext
        y <- ynext
        if (graph)
          points(X[, keyVars[1]], X[, keyVars[2]], pch = 20)
        assign("hybrid.explore.out", list(X = X, X = X, y = y, 
            acceptance = acceptance, swap.acceptance = swaps.record, 
            X0 = X0, bsf = bsf, f.bsf = f.bsf), envir = .GlobalEnv)
        if (is.null(acceptance)) 
            mean.acceptance <- NA
        else mean.acceptance <- mean(acceptance)
        if (dim(X)[1] >= n) 
            break
        if (graph) {
            if (refresh_count < nrefresh.graph) 
                refresh_count <- refresh_count + 1
            else {
                refresh_count <- 1
                key.indx <- key.indx + 1
                if (key.indx > dim(KEY.VARS)[1]) 
                  key.indx <- 1
                keyVars <- KEY.VARS[key.indx, ]
                graph.template(X, key.vars = keyVars, lb, ub)
                title(c("", "", "", "Exploratory phase"), adj = 0, 
                  cex.main = 0.8)
                points(X[, keyVars[1]], X[, keyVars[2]], pch = 20)
            }
        }
        its <- 2
    }
    sd_y <- sd_y[-1, ]
    GPfit <- GProcess(X, y, dsq_eta, request.functions = TRUE)
    details <- list(bsf = bsf, f.bsf = f.bsf, nchains = nchains, 
        T.mult = T.mult, nswaps = nswaps, L = L, delta = delta, 
        BOUNDS = BOUNDS, lb = lb, ub = ub, KEY.VARS = KEY.VARS, 
        graph.template = graph.template, delta.increase = delta.increase, 
        nrefresh.graph = nrefresh.graph, cor.max = cor.max, dsq_eta = dsq_eta)
    return(list(X = X, y = y, Ey = Ey, sd_y = sd_y, f = f, fcalls = fcalls, 
        acceptance = acceptance, swap.acceptance = swaps.record, 
        details = details, GPfit=GPfit))
}

