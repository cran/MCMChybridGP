.hybrid.template <-
function (X, key.vars, lb, ub) 
{
    if(is.null(lb)) lb <- rep(-Inf, length(x1))
    if(is.null(ub)) ub <- rep(Inf, length(x1))
    x1 <- X[, key.vars[1]]
    x2 <- X[, key.vars[2]]
    low1 <- min(x1)
    upp1 <- max(x1)
    low2 <- min(x2)
    upp2 <- max(x2)
    extra1 <- extra2 <- 0.25
    if(sum(abs(lb)) < Inf) {
      low1 <- lb[key.vars[1]]
      low2 <- lb[key.vars[2]]
      extra1 <- 0.05
    }
    if(sum(abs(ub)) < Inf) {
      upp1 <- ub[key.vars[1]]
      upp2 <- ub[key.vars[2]]
      extra2 <- 0.05
    }
    l1 <- low1 - extra1 * (upp1 - low1)
    l2 <- low2 - extra1 * (upp2 - low2)
    u1 <- upp1 + extra2 * (upp1 - low1)
    u2 <- upp2 + extra2 * (upp2 - low2)
    xlim <- c(l1, u1)
    ylim <- c(l2, u2)
    plot(x1, x2, xlim=xlim, type="n", ylim=ylim,
         xlab=paste(sep="", "x", key.vars[1]),
         ylab=paste(sep="", "x", key.vars[2]))
    abline(v=lb[key.vars[1]], lty=2)
    abline(v=ub[key.vars[1]], lty=2)
    abline(h=lb[key.vars[2]], lty=2)
    abline(h=ub[key.vars[2]], lty=2)
    legend(xlim[2], ylim[1], c("Accepted points", "Rejected points"), 
        pch=c(20, 4), pt.cex=c(1.2, 1), 
        cex=0.7, col=c("black", "red"), xjust=1, 
        yjust=0)
}

.refresh.title <-
function (x, attempt=NULL) 
{
    for (ch in LETTERS) {
        str <- ch
        for (i in 1:80) str <- paste(sep="", str, ch)
        title(str, col.main="white")
        title(c(str, str), col.main="white")
        title(str, col.main="white", cex.main=1.5)
        title(c(str, str), col.main="white", cex.main=1.5)
    }
    if (!is.null(attempt)) 
        title(paste(sep="", attempt, ":"), cex.main=0.7, 
            adj=0)
    title(deparse(x), cex.main=0.7)
}

.runMCMCdemo <-
function () 
{
    d <- 2
    mu1 <- c(-1, -1)
    mu2 <- c(+1, +1)
    sigma.sq <- 0.16
    X <- matrix(c( -1, -1,
                   -1, -2,
                 -1.7, -1.7,
                   -2, -1,
                 -1.7, -0.3,
                   -1,  0,
                 -0.3, -0.3,
                    0, -1,
                 -0.3, -1.7,
                    1,  1,
                    1,  0,
                  0.3,  0.3,
                    0,  1,
                  0.3,  1.7,
                    1,  2,
                  1.2,  1.7,
                  1.4,  1,
                  1.2,  0.3,
                  NULL), byrow=TRUE, ncol=2)
    lb <- NULL
    ub <- c(1.5, 3)
    cat("\n")
    cat("mu1 <- c(-1, -1)\n")
    cat("mu2 <- c(+1, +1)\n")
    cat("sigma.sq <-", sigma.sq, "\n")
    lshow <- lb
    ushow <- ub
    nms <- paste(sep="", "x", 1:d)
    if(is.null(lshow)) lshow <- rep(-Inf, d)
    if(is.null(ushow)) ushow <- rep(Inf, d)
    showbounds <- data.frame(nms, lb=lshow, ub=ushow)
    colnames(showbounds)[1] <- ""
    print(showbounds)

    f <- function(x) {
        px <- 1/4/pi/sqrt(sigma.sq) * exp(-1/2/sigma.sq * sum((x - 
            mu1)^2)) + 1/4/pi/sqrt(sigma.sq) * exp(-1/2/sigma.sq * 
            sum((x - mu2)^2))
        return(log(px))
    }
    cat("\n")
    cat("f <- function(x) {\n")
    cat("  px <- 1/4/pi/sqrt(sigma.sq) * exp(-1/2/sigma.sq * sum((x-mu1)^2)) +\n")
    cat("        1/4/pi/sqrt(sigma.sq) * exp(-1/2/sigma.sq * sum((x-mu2)^2))\n")
    cat("  return(log(px))\n")
    cat("}\n")
    cat("\n")
    cat(sep="\t", "X0:", X[1,], "\n")
    cat(sep="\t", "   ", X[2,], "\n")
    cat(sep="\t", "   ", X[3,], "\n")
    cat(sep="\t", "   ", X[4,], "\n")
    cat(sep="\t", "   ", X[5,], "\n")
    cat("    ...\n")
    cat("\n")
    plot.density <- function(X, phase) {
        select1 <- list(1:131, 1:181)
        opar <- par(mfrow=c(2, 1), oma=c(0, 0, 3, 0))
        on.exit(par(opar))
        for (param in 1:2) {
            p1 <- function(x) return(1/4/pi/sqrt(sigma.sq) * 
                exp(-t(x - mu1[param]) * (x - mu1[param])/2/sigma.sq) + 
                1/4/pi/sqrt(sigma.sq) * exp(-t(x - mu2[param]) * 
                  (x - mu2[param])/2/sigma.sq))
            p2 <- function(x) return(1/4/pi/sqrt(sigma.sq) * 
                exp(-t(x - mu1[param]) * (x - mu1[param])/2/sigma.sq) + 
                1/4/pi/sqrt(sigma.sq) * exp(-t(x - mu2[param]) * 
                  (x - mu2[param])/2/sigma.sq))
            x <- seq(-5, 5, by=0.05)
            px <- 0.5 * p1(x) + 0.5 * p2(x)
            px <- px/sum(px * 0.05)
            densty <- density(X[, param])
            plot(x[select1[[param]]], px[select1[[param]]],
              type="l", main=paste(sep="", "Density x", param),
              lwd=0.5, col="blue", xlab="", ylab="",
              xlim=c(-3,3), ylim=c(0, max(c(densty$y, px))))
            abline(v=ub[param], lty=2)
            select2 <- (1:length(densty$x))[densty$x<ub[param]]
            densty$y <- densty$y[select2]
            densty$x <- densty$x[select2]
            lines(densty, lty=2, col="red")
        }
        title(phase, outer=TRUE)
    }
    x.Contour <- seq(-1 - 4 * sqrt(sigma.sq), 1 + 4 * sqrt(sigma.sq), 
        length.out=25)
    fx1.Contour <- fx2.Contour <- rep(NA, 25)
    fx.Contour <- matrix(rep(NA, 25^2), ncol=25)
    n <- 150
    nchains <- 5
    T.mult <- 2
    L <- 100
    delta <- 0.05
    cat("\n")
    cat(sep="", "> explore.out <- hybrid.explore(f, X0, n, L, delta,
        nchains, T.mult, ...)\n")
    input <- readline(prompt=paste(sep="",
      "Enter number of points for Gaussian process, n [", n, "]: "))
    if (input != "") 
        try({
            num <- as.integer(input)
            if (!is.na(num)) 
                n <- num
        })
    if(n > 400) {
      cat("n will be limited to 400\n")
      n <- 400
    }
    input <- readline(prompt=paste(sep="", "Enter number of parallel MCMC chains, nchains [", 
        nchains, "]: "))
    if (input != "") 
        try({
            num <- as.integer(input)
            if (!is.na(num)) 
                nchains <- num
        })
    if (nchains > 1) {
        input <- readline(prompt=paste(sep="", "Enter Temperature multiple, T.mult [", 
            T.mult, "]: "))
        if (input != "") 
            try({
                num <- as.double(input)
                if (!is.na(num)) 
                  T.mult <- num
            })
    }
    cat(sep="", "> explore.out <- hybrid.explore(f, X0, n, L=", L,
      ", delta=", delta, ",\n")
    cat(sep="", "+                  nchains=", nchains, ", T.mult=",
      T.mult, ", lb=lb, ub=ub)\n")
    explore.out <- hybrid.explore(f, X, L=L, delta=delta, 
        n=n, nchains=nchains, T.mult=T.mult, 
        lb=lb, ub=ub, graph=TRUE, nswaps=1)
    iterations <- 500
    replace.target <- 1
    cat("\n")
    cat("> sample.out <- hybrid.sample(explore.out, iterations, replace.target)\n")
    input <- readline(prompt=paste(sep="", "Enter number of sampling iterations [", 
        iterations, "]: "))
    if (input != "") {
        num <- as.integer(input)
        if (!is.na(num)) 
            iterations <- num
    }
    input <- readline(prompt=paste(sep="", "Enter replace.target (0, 1, 2) [", replace.target, "]: "))
    if (input != "") {
        num <- as.integer(input)
        if (input %in% 0:2) 
            replace.target <- num
    }
    cat(sep="", "> sample.out <- hybrid.sample(explore.out, n=", iterations,
      ", replace.target=", replace.target, ")\n")
    sample.out <- hybrid.sample(explore.out, n=iterations, replace.target = replace.target, graph=TRUE)
    cat("\n")
    cat("explore.out: ", names(explore.out), "\n")
    cat("sample.out: ", names(sample.out), "\n")
    plot.density(sample.out$SAMP, phase="sample.out$SAMP")
    return(list(explore.out=explore.out, sample.out=sample.out))
}

.show.time <-
function (date1, date2) 
{
    TIME <- 24 * 60 * 60 * as.integer(substr(date2, 9, 10)) + 
        60 * 60 * as.integer(substr(date2, 12, 13)) + 60 * as.integer(substr(date2, 
        15, 16)) + as.integer(substr(date2, 18, 19)) - 24 * 60 * 
        60 * as.integer(substr(date1, 9, 10)) - 60 * 60 * as.integer(substr(date1, 
        12, 13)) - 60 * as.integer(substr(date1, 15, 16)) - as.integer(substr(date1, 
        18, 19))
    for (ch in LETTERS) {
        str <- ch
        for (i in 1:15) str <- paste(sep="", str, ch)
        title(c("", "", "", str), adj=1, cex.main=0.7, 
            col.main="white")
        if (TIME > 0) 
            title(c("", "", "", paste("Call:", TIME, "seconds")), 
                adj=1, cex.main=0.7)
    }
}

