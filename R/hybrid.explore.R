hybrid.explore <-
function (f, X0, y0=NULL, n = 200, L = 1000, delta = 0.003, nchains = 5, 
    T.mult = 1.5, lb = NULL, ub = NULL, maxleap=0.95, finetune = FALSE,
    Tinclude = nchains, nswaps = choose(nchains, 2), graph = FALSE) 
{
    X <- as.matrix(X0)
    yfinal <- Xfinal <- NULL
    n0 <- dim(X)[1]
    Taccept <- rep(Tinclude+1, n0)
    acceptance <- NULL
    d <- dim(X)[2]
    if (d < 2) 
        stop("Points with dimension of at least 2 assumed for X")
    W <- function(p) return(sum(p^2)/2)
    if (d < 2) 
        graph <- FALSE
    if (n0 > 1) 
        for (i in 2:n0) for (j in 1:(i - 1)) if (sum(abs(X[i, 
            ] - X[j, ])) == 0) 
            stop(paste(sep = "", "Duplicate rows (", i, " & ", 
                j, ") in X not allowed"))
    BOUNDS <- FALSE
    if ((!is.null(lb))  || (!is.null(ub)))
            BOUNDS <- TRUE
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
      if (is.na(sum(lb + ub))) 
          stop("hybrid.explore: NaN lb, ub")
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
    if (graph) {
        if(!exists(".Rpackage_DEMO")) {
          Rpackage_DEMO <- .Rpackage_DEMO <- FALSE
        } else {
          Rpackage_DEMO <- .Rpackage_DEMO
          assign(".Rpackage_DEMO", FALSE, envir=.GlobalEnv)
        }
        KEY.VARS <- NULL
        nrefresh.graph <- ceiling(10/nchains)
        if(is.logical(Rpackage_DEMO))
          if(Rpackage_DEMO) nrefresh.graph <- Inf
        key.vars <- 1:dim(X)[2]
        for (i in 1:(length(key.vars) - 1)) for (j in (i + 1):length(key.vars)) {
            if (key.vars[i] == key.vars[j]) 
                stop("Duplicate key.vars not allowed")
            pair <- c(key.vars[i], key.vars[j])
            KEY.VARS <- rbind(KEY.VARS, pair)
        }
        graph.template <- .hybrid.template
        key.indx <- 1
        keyVars <- KEY.VARS[key.indx, ]
    }
    Temp <- T.mult^(1:nchains - 1)
    swaps.record <- NULL
    fcalls <- 0
    if(!is.null(y0)) {
      y <- y0 
    } else {
      cat(sep = "", "hybrid.explore: ", n0, " initial ", d, "-D points supplied for ", 
          n, " Exploratory points\n")
      cat("Leapfrogs: L =", L, "and delta =", delta, "\n")
      if (T.mult < 1) 
          stop("T.mult < 1 not allowed")
      y <- rep(NA, n0)
      cat("Parallel Tempering:", nchains, "chains with T:", Temp, 
          "\n")
      date1 <- date2 <- date()
      for (j in 1:n0) {
          functionOK <- FALSE
          try({
            Xj <- X[j,]
            OUTSIDE_BOUNDS <- FALSE
            if (BOUNDS) {
              for (i in 1:d) {
                if (Xj[i] < lb[i]) {
                  OUTSIDE_BOUNDS <- TRUE
                  Xj[i] <- lb[i]
                }
                if (Xj[i] > ub[i]) {
                  OUTSIDE_BOUNDS <- TRUE
                  Xj[i] <- ub[i]
                }
              }
            }
            if(OUTSIDE_BOUNDS) stop("X0 violates bounds supplied")
            date1 <- date()
            cat("\n")
            y[j] <- f(X[j,])
            date2 <- date()
            cat(sep = "", "f(x", j, ") = ", y[j], " for x", j, 
                ":\n")
            cat("", X[j, ], "\n")
            .show.time(date1, date2, graph = FALSE)
            functionOK <- TRUE
          })
          if(!functionOK) stop("f(x) failed")
          if (is.na(y[j])) stop("f(x) NaN")
          if (y[j] == Inf) stop("Infinite f(x) not allowed")
      }
      fcalls <- n0
    }
    select <- sample(1:n0, nchains)
    x <- X[select, ]
    if (is.null(dim(x))) 
        x <- t(x)
    fx <- y[select]
    if (graph) {
        graph.template(X, key.vars = keyVars, lb, ub)
        points(X[, keyVars[1]], X[, keyVars[2]], pch = 20, col = "grey")
        points(Xfinal[, keyVars[1]], Xfinal[, keyVars[2]], pch = 20, col = "black")
        title(c("", "", "", "Exploratory phase"), adj = 0, cex.main = 0.8)
        clrs <- c("blue", rep(c("green", "turquoise", "orange", 
            "yellow", "red", "grey", "pink"), nchains))
        reds <- c("blue", rep(c("green", "turquoise", "orange", 
            "yellow", "red", "grey", "pink"), nchains))
    }
    count10 <- 0
    nacc <- 0
    refresh_count <- 1
    params <- rep(1, d+1)
    repeat {
        cat("\n")
        GPfit <- GProcess(X, y, params = params, request.functions = FALSE, finetune = finetune)
        if(!GPfit$inverseOK) stop("hybrid.explore: GProcess inverse covariance matrix failed")
        Sigma.inv <- GPfit$Sigma.inv
        params <- GPfit$params
        lsq <- params[1]
        eta <- params[-1]
        c.scalar <- function(xi, xj) {
            h <- abs(xi - xj)
            val <- lsq * prod((1 + eta * h) * exp(-eta * h))
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
            val <- lsq - as.double(t(c.vec_x) %*% Sigma.inv %*% 
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
        if (nchains >= 2) {
            for (reps in 1:nswaps) for (index1 in sample(1:nchains)) for (index2 in sample((1:nchains)[-index1])) {
                if (T.mult == 1) 
                  Tindex <- 1:nchains
                else Tindex <- sort(Temp, index.return = TRUE)$ix
                swap <- Tindex[c(index1, index2)]
                T1 <- Temp[swap[1]]
                T2 <- Temp[swap[2]]
                sdf1 <- sigmaf(x[swap[1]])
                sdf2 <- sigmaf(x[swap[2]])
                corf1 <- sdf1/sqrt(lsq)
                corf2 <- sdf2/sqrt(lsq)
                fswap1 <- fx[swap[1]] + sdf1
                fswap2 <- fx[swap[2]] + sdf2
                alpha <- (fswap2/T1 - fswap1/T1 + fswap1/T2 - fswap2/T2)
                if(is.na(alpha)) alpha <- -Inf
                if (log(runif(1)) < alpha) {
                  Temp[swap[1]] <- T2
                  Temp[swap[2]] <- T1
                  swaps.record <- c(swaps.record, 1)
                  cat("hybrid.explore: Swap T =", T1, "with", 
                    T2, "\n")
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
            sdf_x.C <- sigmaf(x.C)
            corf_x.C <- sdf_x.C/sqrt(lsq)
            cat("\n")
            cat(sep = "", "Leapfrog from x(", chain, ", T = ", 
                Temp[chain], "):\n")
            cat("", rnd(x.C), "\n")
            fx.C <- fx[chain]
            xp_Early <- NULL
            count <- 0
            repeat {
                p.C <- rnorm(d)
                xp_Early <- matrix(c(x.C, p.C, rep(0, d)), ncol = 3)
                xOK <- FALSE
                low <- lb
                upp <- ub
                try({
                  xp_Early <- .Call("Leapfrog", X, y, x.C, p.C, 
                    L, delta, GPfit$Sigma.inv, params, Temp[chain], sqrt(lsq)*maxleap, low, upp)
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
                  GPfit <- GProcess(X, y, params = params, request.functions = TRUE)
                  if(!GPfit$inverseOK)
                    stop("hybrid.explore: GProcess inverse covariance matrix failed")
                  Sigma.inv <- GPfit$Sigma.inv
                  params <- GPfit$params
                  lsq <- params[1]
                  eta <- params[-1]
                  sigmaf <- GPfit$sigmaf
                  Xnext <- X
                  ynext <- y
                  count <- 0
                }
            }
            out <- list(x = xp_Early[, 1], p = xp_Early[, 2], 
                Early = xp_Early[1, 3])
            x.P <- out$x
            OUTSIDE_BOUNDS <- FALSE
            if (BOUNDS) {
                for (i in 1:d) {
                  if (x.P[i] < lb[i]) {
                    OUTSIDE_BOUNDS <- TRUE
                  }
                  if (x.P[i] > ub[i]) {
                    OUTSIDE_BOUNDS <- TRUE
                  }
                }
            }
            p.P <- out$p
            cat(sep = "", "Proposed x(", chain, ", T =", 
              Temp[chain], "):\n")
            cat("", rnd(x.P), "\n")
            sdf_x.P <- sigmaf(x.P)
            corf_x.P <- sdf_x.P/sqrt(lsq)
            if(as.integer(count10/10) == count10/10) MAX10 <- 0
            if(corf_x.P > MAX10) MAX10 <- corf_x.P
            if (out$Early < 0) {
              cat("Early stop: sd(f(x))/lambda =", corf_x.P, "where maxleap =", maxleap, "\n")
            }
            else {
              cat("Full Leapfrog: sd(f(x))/lambda =", corf_x.P, "where maxleap =", maxleap, "\b)\n")
              if (corf_x.P > maxleap) 
                stop(paste("Caution: Leapfrog gave sigmaf =", 
                  out$Early/sqrt(lsq), "when sigmaf(x)/lambda =", corf_x.P))
            }
            if (graph) {
              clr <- clrs[path]
              clr_red <- reds[path]
              if(!OUTSIDE_BOUNDS)
                lines(c(x.C[keyVars[1]], x.P[keyVars[1]]), 
                c(x.C[keyVars[2]], x.P[keyVars[2]]), lwd = 0.2, 
                lty = 2, col = clr)
            }
            if(OUTSIDE_BOUNDS) {
              cat("T =", Temp[chain], "\b: Ignore point out of bounds\n")
              fx.P <- alpha <- -Inf
            } else {
              cat("Leapfrog: Calculate true f(x) for proposal.\n")
              date1 <- date()
              cat("\n")
              fx.P <- f(x.P)
              date2 <- date()
              if (is.na(fx.P)) 
                stop("f(x) NaN")
              if (fx.P == Inf) 
                stop("Infinite f(x) not allowed")
              fcalls <- fcalls + 1
              cat("f(x) =", fx.P, "leap from f(x)", fx.C, "\n")
              .show.time(date1, date2, graph)
              H_xp.C <- -fx.C - sdf_x.C + W(p.C)
              H_xp.P <- -fx.P - sdf_x.P + W(p.P)
              alpha <- min(((-H_xp.P + H_xp.C)/Temp[chain]), 0)
            }
            if (log(runif(1)) < alpha) {
              x[chain, ] <- x.P
              fx[chain] <- fx.P
              Xnext <- rbind(Xnext, x.P)
              ynext <- c(ynext, fx.P)
              if(path <= Tinclude) {
                Xfinal <- rbind(Xfinal, x.P)
                yfinal <- c(yfinal, fx.P)
                Taccept <- c(Taccept, path)
                acceptance <- c(acceptance, 1)
              }
              cat(sep = "", "ACCEPT x(", chain, ", T = ", 
                Temp[chain], "):\n")
              cat("", rnd(x.P), "\n")
              if(path > Tinclude) cat("Point marked for later removal",
                          "(Tinclude =", Tinclude, "\b)\n")
              if (graph) 
                points(x.P[keyVars[1]], x.P[keyVars[2]], 
                  pch = 20, col = clr)
            }
            else {
              cat(sep = "", "REJECT x.P (", chain, ", T=", Temp[chain], "):\n")
              if(path == 1) if(!OUTSIDE_BOUNDS) {
                Xnext <- rbind(Xnext, x.P)
                ynext <- c(ynext, fx.P)
                Taccept <- c(Taccept, Tinclude+1)
                cat("Primary chain point marked for later removal\n")
              }
              cat("", rnd(x.P), "\n")
              if(path <= Tinclude) acceptance <- c(acceptance, 0)
              if (graph) if(!OUTSIDE_BOUNDS)
                points(x.P[keyVars[1]], x.P[keyVars[2]], 
                  pch = "X", col = clr_red)
            }
            if(dim(Xnext)[1] > n) {
              i <- 0
              repeat {
                i <- i+1
                if(Taccept[i] > Tinclude) {
                  Xnext <- Xnext[-i,]
                  ynext <- ynext[-i]
                  i <- i-1
                }
                if(length(ynext) <= n) break
              }
            }
            if (nchains >= 2) 
              cat(sep="", "Swaps: ", round(10000 * mean(swaps.record))/100, "% acceptance\n")
            if (is.null(acceptance)) 
              mean.acceptance <- 0 else
              mean.acceptance <- mean(acceptance[ceiling(length(acceptance)/2):length(acceptance)])
            cat(sep = "", "hybrid.explore: ",
              round(10000 * mean.acceptance)/100, "% acceptance of required points\n")
            cat("Leapfrog used: L =", L, ", delta =", delta, "\n")
            cat("Present Gaussian process consists of", dim(Xnext)[1], "points\n")
            cat("So far", sum(Taccept<=Tinclude), "/", n, "required points in",
              Tinclude, "out of", nchains, "chains\n")
            if(!is.null(Xfinal)) if(dim(Xfinal)[1] >= n) break
            count10 <- count10+1
            cat("sd(f(x))/lambda =", corf_x.P, "where maxleap =", maxleap, "\n")
            cat("Maximum of 10 observed sd(f(x)) =", MAX10, "\n")
        }
        X <- Xnext
        y <- ynext
        assign("hybrid.explore.out", list(X = Xfinal, y = yfinal, Xfull=X, yfull=y, acceptance = acceptance,
            X0 = X0), envir = .GlobalEnv)
        if (is.null(acceptance)) 
            mean.acceptance <- NA
        else mean.acceptance <- mean(acceptance)
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
                points(X[, keyVars[1]], X[, keyVars[2]], pch = 20, col="grey")
                points(Xfinal[, keyVars[1]], Xfinal[, keyVars[2]], pch = 20)
            }
        }
        if(!is.null(Xfinal)) if (dim(Xfinal)[1] >= n) break
    }
    X <- Xfinal
    y <- yfinal
    GPfit <- GProcess(X, y, params, request.functions = TRUE, finetune=TRUE)
    if(!GPfit$inverseOK)
      cat("hybrid.explore: GProcess inverse covariance matrix failed\n")
    details <- list(nchains = nchains, T.mult = T.mult, nswaps = nswaps, L = L, delta = delta,
        lb = lb, ub = ub, KEY.VARS = KEY.VARS)
    return(list(X = X, y = y, f = f, maxleap=maxleap, function.calls = fcalls, details = details, GPfit = GPfit))
}

