.hybrid.template <-
function (X, key.vars, lb, ub) 
{
    if (!is.null(X)) {
        x1 <- X[, key.vars[1]]
        x2 <- X[, key.vars[2]]
    }
    else x1 <- x2 <- NULL
    DRAW_BOUNDS <- FALSE
    if (!is.null(lb)) {
        l1 <- lb[key.vars[1]]
        u1 <- ub[key.vars[1]]
        l2 <- lb[key.vars[2]]
        u2 <- ub[key.vars[2]]
        DRAW_BOUNDS <- TRUE
    }
    else {
        low1 <- min(x1)
        upp1 <- max(x1)
        low2 <- min(x2)
        upp2 <- max(x2)
        l1 <- low1 - 0.5 * (upp1 - low1)
        u1 <- upp1 + 0.5 * (upp1 - low1)
        l2 <- low2 - 0.5 * (upp2 - low2)
        u2 <- upp2 + 0.5 * (upp2 - low2)
    }
    xlim <- c(l1, u1)
    ylim <- c(l2, u2)
    plot(x1, x2, xlim = xlim, type = "n", ylim = ylim, xlab = paste(sep = "", 
        "x", key.vars[1]), ylab = paste(sep = "", "x", key.vars[2]))
    if (DRAW_BOUNDS) {
        abline(v = l1, lty = 2)
        abline(v = u1, lty = 2)
        abline(h = l2, lty = 2)
        abline(h = u2, lty = 2)
    }
    legend(xlim[2], ylim[1], c("Accepted points", "Rejected points"), 
        pch = c(20, 4), pt.cex = c(1.2, 1), 
        cex = 0.7, col = c("black", "red"), xjust = 1, 
        yjust = 0)
}

.refresh.title <-
function (x, attempt = NULL) 
{
    for (ch in LETTERS) {
        str <- ch
        for (i in 1:80) str <- paste(sep = "", str, ch)
        title(str, col.main = "white")
        title(c(str, str), col.main = "white")
        title(str, col.main = "white", cex.main = 1.5)
        title(c(str, str), col.main = "white", cex.main = 1.5)
    }
    if (!is.null(attempt)) 
        title(paste(sep = "", attempt, ":"), cex.main = 0.7, 
            adj = 0)
    title(deparse(x), cex.main = 0.7)
}

.runMCMCdemo <-
function () 
{
    mu1 <- c(-1, -1)
    mu2 <- c(+1, +1)
    sigma.sq <- 0.16
    X <- matrix(c(-2, -2, -1, -1, 0, -2, -2, 0, 0, 0, 2, 0, 0, 
        2, 1, 1, 2, 2), byrow = TRUE, ncol = 2)
    cat("\n")
    cat("mu1 <- c(-1, -1)\n")
    cat("mu2 <- c(+1, +1)\n")
    cat("sigma.sq <-", sigma.sq, "\n")
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
    cat("t(X): -2   -1    0   -2    0    2    0    1    2 \n")
    cat("      -2   -1   -2    0    0    0    2    1    2 \n")
    plot.density <- function(X, phase) {
        oma <- par(mfrow = c(2, 1), oma = c(0, 0, 3, 0))
        for (param in 1:2) {
            p1 <- function(x) return(1/4/pi/sqrt(sigma.sq) * 
                exp(-t(x - mu1[param]) * (x - mu1[param])/2/sigma.sq) + 
                1/4/pi/sqrt(sigma.sq) * exp(-t(x - mu2[param]) * 
                  (x - mu2[param])/2/sigma.sq))
            p2 <- function(x) return(1/4/pi/sqrt(sigma.sq) * 
                exp(-t(x - mu1[param]) * (x - mu1[param])/2/sigma.sq) + 
                1/4/pi/sqrt(sigma.sq) * exp(-t(x - mu2[param]) * 
                  (x - mu2[param])/2/sigma.sq))
            x <- seq(-5, 5, by = 0.05)
            px <- 0.5 * p1(x) + 0.5 * p2(x)
            px <- px/sum(px * 0.05)
            densty <- density(X[, param])
            plot(x, px, type = "l", main = paste(sep = "", 
                "Density x", param), lwd = 0.5, col="blue",
                ylim = c(0, max(c(densty$y, px))))
            lines(densty, lty = 2, col = "red")
        }
        title(phase, outer = TRUE)
        par(oma)
    }
    x.Contour <- seq(-1 - 4 * sqrt(sigma.sq), 1 + 4 * sqrt(sigma.sq), 
        length.out = 25)
    fx1.Contour <- fx2.Contour <- rep(NA, 25)
    fx.Contour <- matrix(rep(NA, 25^2), ncol = 25)
    n <- 100
    nchains <- 5
    T.mult <- 1.5
    L <- 100
    delta = 0.005
    cat("\n")
    cat(sep="", "> explore.out <- hybrid.explore(f, X0, n, L, delta, nchains, T.mult, ...)\n")
    input <- readline(prompt = paste(sep = "", "Enter number of points for Gaussian process, n [", n, "]: "))
    if (input != "") 
        try({
            num = as.integer(input)
            if (!is.na(num)) 
                n <- num
        })
    if(n > 400) {
      cat("n will be limited to 400\n")
      n <- 400
    }
    input <- readline(prompt = paste(sep = "", "Enter number of parallel MCMC chains, nchains [", 
        nchains, "]: "))
    if (input != "") 
        try({
            num = as.integer(input)
            if (!is.na(num)) 
                nchains <- num
        })
    if (nchains > 1) {
        input <- readline(prompt = paste(sep = "", "Enter Temperature multiple, T.mult [", 
            T.mult, "]: "))
        if (input != "") 
            try({
                num = as.double(input)
                if (!is.na(num)) 
                  T.mult <- num
            })
    }
    cat(sep="", "> explore.out <- hybrid.explore(f, X0, n, L=", L,
      ", delta=", delta, ",\n")
    cat(sep="", "+                  nchains=", nchains, ", T.mult=",
      T.mult, ")\n")
    template.explore <- function(key.vars, ...) {
        if (nchains > 1) {
            cmd <- paste(sep = "", " hybrid.explore(f, X, n=", 
                n, ", nchains=", nchains, ", T.mult=", T.mult, ")")
        }
        else cmd <- paste(sep = "", " hybrid.explore(f, X, n=", 
            n, ", nchains=", nchains, ")")
        for (i in 1:25) for (j in 1:25) {
            xin <- rep(0, 2)
            xin[key.vars] <- c(x.Contour[i], x.Contour[j])
            fx.Contour[i, j] <- f(xin)
        }
        par(mfrow = c(1, 1), oma = c(0, 1, 0, 0))
        par(mfrow = c(1, 1), oma = c(0, 1, 0, 0))
        contour(x.Contour, x.Contour, exp(fx.Contour), xlab = paste(sep = "", 
            "x", key.vars[1]), ylab = paste(sep = "", "x", key.vars[2]), 
            cex.main = 1, xlim = c(-4, 4), ylim = c(-4, 4), levels = c(0.01, 
                1e-04), labels = rep("", 2))
        title(c("", "", "", "", "", "", cmd), cex.main = 0.7, 
            adj = 0)
        legend(4, -4, c("Included points", "Not included", "Early stops"), 
            pch = c(20, 4, 1), pt.cex = c(1.2, 1, 1.2), cex = 0.7, 
            col = c("black", "red", "blue"), xjust = 1, yjust = 0)
    }
    explore.out <- hybrid.explore(f, X, L = L, delta = delta, 
        n = n, nchains = nchains, cor.max = 0.7, T.mult = T.mult, 
        graph = TRUE, graph.template = template.explore, key.vars = 1:2,
        nswaps=5, spread=1.5, nrefresh.graph = 25)
    assign("explore.out", explore.out, envir = .GlobalEnv)
    iterations <- 200
    replace.target <- 0
    cat("\n")
    cat("> sample.out <- hybrid.sample(explore.out, iterations, replace.target)\n")
    input <- readline(prompt = paste(sep = "", "Enter number of sampling iterations [", 
        iterations, "]: "))
    if (input != "") {
        num = as.integer(input)
        if (!is.na(num)) 
            iterations <- num
    }
    input <- readline(prompt=paste(sep = "", "Enter replace.target (0, 1, 2) [", replace.target, "]: "))
    if (input != "") {
        num = as.integer(input)
        if (input %in% 0:2) 
            replace.target <- num
    }
    template.sample <- function(key.vars, ...) {
        cmd1 <- paste(sep = "", " hybrid.sample(explore.out, iterations=", 
            iterations, ",")
        cmd2 <- paste("                                  replace.target=",
            deparse(replace.target), ")")
        for (i in 1:25) for (j in 1:25) {
            xin <- rep(0, 2)
            xin[key.vars] <- c(x.Contour[i], x.Contour[j])
            fx.Contour[i, j] <- f(xin)
        }
        par(mfrow = c(1, 1), oma = c(0, 1, 0, 0))
        par(mfrow = c(1, 1), oma = c(0, 1, 0, 0))
        contour(x.Contour, x.Contour, exp(fx.Contour), xlab = paste(sep = "", 
            "x", key.vars[1]), ylab = paste(sep = "", "x", key.vars[2]), 
            cex.main = 1, xlim = c(-4, 4), ylim = c(-4, 4), levels = c(0.01, 
                1e-04), labels = rep("", 2))
        title(c("", "", "", "", "", "", cmd1, cmd2), cex.main = 0.79, 
            adj = 0)
        legend(4, -4, c("Sampled points", "Rejected points"), 
            pch = c(20, 4, 1), pt.cex = c(1.2, 1, 1.2), cex = 0.7, 
            col = c("black", "red", "blue"), xjust = 1, yjust = 0)
    }
    cat(sep="", "> sample.out <- hybrid.sample(explore.out, n=", iterations,
      ", replace.target=", replace.target, ")\n")
    sample.out <- hybrid.sample(explore.out, n = iterations,
      graph = TRUE, graph.template = template.sample,
      key.vars = 1:2, nrefresh.graph=50)
    assign("sample.out", sample.out, envir = .GlobalEnv)
    cat("sample.out: ", names(sample.out), "\n")
    plot.density(sample.out$SAMP, phase = "sample.out$SAMP")
    return(list(explore.out = explore.out, sample.out = sample.out))
}

.show.time <-
function (date1, date2, graph) 
{
    TIME <- 24 * 60 * 60 * as.integer(substr(date2, 9, 10)) + 
        60 * 60 * as.integer(substr(date2, 12, 13)) + 60 * as.integer(substr(date2, 
        15, 16)) + as.integer(substr(date2, 18, 19)) - 24 * 60 * 
        60 * as.integer(substr(date1, 9, 10)) - 60 * 60 * as.integer(substr(date1, 
        12, 13)) - 60 * as.integer(substr(date1, 15, 16)) - as.integer(substr(date1, 
        18, 19))
    if (TIME > 0) 
        cat("Function call:", TIME, "seconds\n")
    if (graph) {
        for (ch in LETTERS) {
            str <- ch
            for (i in 1:15) str <- paste(sep = "", str, ch)
            title(c("", "", "", str), adj = 1, cex.main = 0.7, 
                col.main = "white")
        }
        if (TIME > 0) 
            title(c("", "", "", paste("Call:", TIME, "seconds")), 
                adj = 1, cex.main = 0.7)
    }
}

