#Originally from David Clifford and Peter McCullagh (2011). 
#regress: Gaussian linear models with linear covariance 
#structure. R package version 1.3.
#https://CRAN.R-project.org/package=regress

regress <- function (formula, Vformula, identity = TRUE, kernel = NULL, 
    start = NULL, taper = NULL, pos, verbose = 0, gamVals = NULL, 
    maxcyc = 50, tol = 1e-04, data, fraction = NULL, print.level = NULL) 
{
    if (!is.null(print.level)) {
        cat("\nWarning: print.level has been replaced by verbose and has been deprecated.\nIt will be removed in the next version of regress\n\n")
        verbose <- print.level
    }
    if (!is.null(print.level)) {
        cat("\nWarning: fraction has been deprecated and replaced by taper - a vector of values from 0 to 1 giving a unique fraction for each step of the Newton Raphson algorithm.\n")
    }
    if (missing(data)) 
        data <- environment(formula)
    mf <- model.frame(formula, data = data, na.action = na.pass)
    mf <- eval(mf, parent.frame())
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
    model <- list()
    model <- c(model, mf)
    if (missing(Vformula)) 
        Vformula <- NULL
    if (!is.null(Vformula)) {
        V <- model.frame(Vformula, data = data, na.action = na.pass)
        V <- eval(V, parent.frame())
        Vcoef.names <- names(V)
        V <- as.list(V)
        k <- length(V)
    }
    else {
        V <- NULL
        k <- 0
        Vcoef.names = NULL
    }
    isNA <- is.na(y)
    y <- na.omit(y)
    n <- length(y)
    Xcolnames <- dimnames(X)[[2]]
    if (is.null(Xcolnames)) {
        Xcolnames <- paste("X.column", c(1:dim(as.matrix(X))[2]), 
            sep = "")
    }
    X <- X[isNA == F, ]
    X <- matrix(X, n, length(X)/n)
    qr <- qr(X)
    rankQ <- n - qr$rank
    if (qr$rank) {
        X <- matrix(X[, qr$pivot[1:qr$rank]], n, qr$rank)
        Xcolnames <- Xcolnames[qr$pivot[1:qr$rank]]
    }
    else {
        cat("\nERROR: X has rank 0\n\n")
    }
    if (missing(kernel)) {
        K <- X
        colnames(K) <- Xcolnames
        reml <- TRUE
        kernel <- NULL
    }
    else {
        if (length(kernel) == 1 && kernel > 0) {
            K <- matrix(rep(1, n), n, 1)
            colnames(K) <- c("1")
        }
        if (length(kernel) == 1 && kernel <= 0) {
            K <- Kcolnames <- NULL
            KX <- X
            rankQK <- n
        }
        if (length(kernel) > 1) {
            K <- kernel[isNA == F, ]
        }
        reml <- FALSE
    }
    if (!is.null(K)) {
        Kcolnames <- colnames(K)
        qr <- qr(K)
        rankQK <- n - qr$rank
        if (qr$rank == 0) 
            K <- NULL
        else {
            K <- matrix(K[, qr$pivot[1:qr$rank]], n, qr$rank)
            Kcolnames <- Kcolnames[qr$pivot[1:qr$rank]]
            KX <- cbind(K, X)
            qr <- qr(KX)
            KX <- matrix(KX[, qr$pivot[1:qr$rank]], n, qr$rank)
        }
    }
    if (missing(maxcyc)) 
        maxcyc <- 50
    if (missing(tol)) 
        tol <- 1e-04
    delta <- 1
    for (i in 1:k) {
        if (is.matrix(V[[i]])) {
            V[[i]] <- V[[i]][isNA == F, ]
            V[[i]] <- V[[i]][, isNA == F]
        }
        if (is.factor(V[[i]])) {
            V[[i]] <- V[[i]][isNA == F]
        }
    }
    In <- diag(rep(1, n), n, n)
    if (identity) {
        V[[k + 1]] <- as.factor(1:n)
        names(V)[k + 1] <- "In"
        k <- k + 1
        Vcoef.names <- c(Vcoef.names, "In")
        Vformula <- as.character(Vformula)
        Vformula[2] <- paste(Vformula[2], "+In")
        Vformula <- as.formula(Vformula)
    }
    model <- c(model, V)
    model$formula <- formula
    model$Vformula <- Vformula
    if (!missing(pos)) 
        pos <- as.logical(pos)
    if (missing(pos)) 
        pos <- rep(FALSE, k)
    pos <- c(pos, rep(FALSE, k))
    pos <- pos[1:k]
    SWsolveINDICATOR <- FALSE
    if (all(sapply(V, is.factor))) {
        SWsolveINDICATOR <- TRUE
        Z <- list()
        for (i in 1:length(V)) {
            if (is.factor(V[[i]])) {
                Vi <- model.matrix(~V[[i]] - 1)
                colnames(Vi) <- levels(V[[i]])
                Z[[i]] <- Vi
                V[[i]] <- tcrossprod(Vi)
            }
        }
        names(Z) <- names(V)
    }
    else {
        for (i in 1:length(V)) {
            if (is.factor(V[[i]])) {
                Vi <- model.matrix(~V[[i]] - 1)
                Vi <- tcrossprod(Vi)
                V[[i]] <- Vi
            }
        }
        Z <- V
    }
    A <- matrix(rep(0, k^2), k, k)
    entries <- expand.grid(1:k, 1:k)
    x <- rep(0, k)
    sigma <- c(1, rep(0, k - 1))
    stats <- rep(0, 0)
    if (missing(taper)) {
        taper <- rep(0.9, maxcyc)
        if (missing(start) && k > 1) 
            taper[1:2] <- c(0.5, 0.7)
    }
    else {
        taper <- pmin(abs(taper), 1)
        if ((l <- length(taper)) < maxcyc) 
            taper <- c(taper, rep(taper[l], maxcyc - l))
    }
    if (!is.null(start)) {
        start <- c(start, rep(1, k))
        start <- start[1:k]
    }
    if (k > 2 && is.null(start)) 
        start <- rep(var(y, na.rm = TRUE), k)
    if (k == 1 && is.null(start)) 
        start <- var(y, na.rm = TRUE)
    if (is.null(start) && k == 2) {
        if (missing(gamVals)) {
            gamVals <- seq(0.01, 0.02, length = 3)^2
            gamVals <- sort(c(gamVals, seq(0.1, 0.9, length = 3), 
                1 - gamVals))
            gamVals <- 0.5
        }
        if (length(gamVals) > 1) {
            if (verbose >= 1) 
                cat("Evaluating the llik at gamma = \n")
            if (verbose >= 1) 
                cat(gamVals)
            if (verbose >= 1) 
                cat("\n")
            reg.obj <- reml(gamVals, y, X, V[[1]], V[[2]], verbose = verbose)
            llik <- reg.obj$llik
            llik <- as.real(llik)
            if (verbose >= 2) 
                cat(llik, "\n")
            gam <- gamVals[llik == max(llik)]
            gam <- gam[1]
            if (verbose >= 2) 
                cat("MLE is near", gam, "and llik =", max(llik), 
                  "there\n")
        }
        if (length(gamVals) == 1) {
            gam <- gamVals[1]
            reg.obj <- list(rms = var(y))
        }
        start <- c(1 - gam, gam) * reg.obj$rms
        if (gam == 0.9999) {
            taper[1] <- taper[1]/100
            maxcyc <- maxcyc * 10
        }
        if (verbose >= 1) 
            cat(c("start algorithm at", round(start, 4), "\n"))
    }
    if (is.null(start) & k > 2) {
        LLvals <- NULL
        V2 <- V[[2]]
        for (ii in 3:k) V2 <- V2 + V[[ii]]
        LLvals <- c(LLvals, reml(0.5, y, X, V[[1]], V2)$llik)
        V2 <- V[[1]] + V2
        for (ii in 1:k) {
            V2 <- V2 - V[[ii]]
            LLvals <- c(LLvals, reml(0.75, y, X, V2, V[[ii]])$llik)
        }
        best <- which.max(LLvals)
        if (verbose) {
            cat("Checking starting points\n")
            cat("llik values of", LLvals, "\n")
        }
        if (best == 1) {
            start <- rep(var(y, na.rm = TRUE), k)
        }
        else {
            start <- rep(0.25, k)
            start[best] <- 0.75
        }
    }
    sigma <- coef <- start
    coef[pos] <- log(sigma[pos])
    coef[!pos] <- sigma[!pos]
    T <- vector("list", length = k)
    for (ii in 1:k) T[[ii]] <- matrix(NA, n, n)
    for (cycle in 1:maxcyc) {
        ind <- which(pos)
        if (length(ind)) {
            coef[ind] <- pmin(coef[ind], 20)
            coef[ind] <- pmax(coef[ind], -20)
            sigma[ind] <- exp(coef[ind])
        }
        if (verbose >= 1) {
            cat(cycle, " ")
        }
        if (!SWsolveINDICATOR) {
            Sigma <- 0
            for (i in 1:k) Sigma <- Sigma + V[[i]] * sigma[i]
            W <- solve(Sigma, In)
        }
        else {
            W <- SWsolve2(Z[1:(k - 1)], sigma)
        }
        if (is.null(K)) 
            WQK <- W
        else {
            WK <- W %*% K
            WQK <- W - WK %*% solve(t(K) %*% WK, t(WK))
        }
        if (reml) 
            WQX <- WQK
        else {
            WX <- W %*% KX
            WQX <- W - WX %*% solve(t(KX) %*% WX, t(WX))
        }
        rss <- as.numeric(t(y) %*% WQX %*% y)
        sigma <- sigma * rss/rankQK
        coef[!pos] <- sigma[!pos]
        coef[pos] <- log(sigma[pos])
        WQK <- WQK * rankQK/rss
        WQX <- WQX * rankQK/rss
        rss <- rankQK
        eig <- sort(eigen(WQK, symmetric = TRUE, only.values = TRUE)$values, 
            decreasing = TRUE)[1:rankQK]
        if (any(eig < 0)) {
            cat("error: Sigma is not pos def on contrasts: range(eig)=", 
                range(eig), "\n")
            WQK <- WQK + (tol - min(eig)) * diag(nobs)
            eig <- eig + tol - min(eig)
        }
        ldet <- sum(log(eig))
        llik <- ldet/2 - rss/2
        if (cycle == 1) 
            llik0 <- llik
        delta.llik <- llik - llik0
        llik0 <- llik
        if (verbose) 
            cat("sigma =", sigma, "(scale-adjusted)\n")
        if (verbose && reml) 
            cat("resid llik =", llik, "delta.llik =", delta.llik, 
                "\n")
        if (verbose && !reml) 
            cat("llik =", llik, "delta.llik =", delta.llik, "\n")
        x <- NULL
        var.components <- rep(1, k)
        ind <- which(pos)
        if (length(ind)) 
            var.components[ind] <- sigma[ind]
        if (!SWsolveINDICATOR) {
            if (identity) {
                T[[k]] <- WQK
                if (k > 1) {
                  for (ii in (k - 1):1) T[[ii]] <- WQK %*% V[[ii]]
                }
            }
            else {
                for (ii in 1:k) T[[ii]] <- WQK %*% V[[ii]]
            }
        }
        else {
            if (identity) {
                T[[k]] <- WQK
                if (k > 1) {
                  for (ii in (k - 1):1) T[[ii]] <- tcrossprod(WQK %*% 
                    Z[[ii]], Z[[ii]])
                }
            }
            else {
                for (ii in 1:k) T[[ii]] <- tcrossprod(WQK %*% 
                  Z[[ii]], Z[[ii]])
            }
        }
        x <- sapply(T, function(x) as.numeric(t(y) %*% x %*% 
            WQX %*% y - sum(diag(x))))
        x <- x * var.components
        ff <- function(x) sum(T[[x[1]]] * t(T[[x[2]]])) * var.components[x[1]] * 
            var.components[x[2]]
        aa <- apply(entries, 1, ff)
        A[as.matrix(entries)] <- aa
        stats <- c(stats, llik, sigma[1:k], x[1:k])
        if (verbose == -1) {
        }
        A.svd <- ginv(A)
        x <- A.svd %*% x
        if (qr(A)$rank < k) {
            if (cycle == 1) {
                if (verbose) {
                  cat("Warning: Non identifiable dispersion model\n")
                  cat(sigma)
                  cat("\n")
                }
            }
        }
        coef <- coef + taper[cycle] * x
        sigma[!pos] <- coef[!pos]
        sigma[pos] <- exp(coef[pos])
        if (cycle > 1 & abs(delta.llik) < tol * 10) 
            break
        if (max(abs(x)) < tol) 
            break
    }
    if (cycle == maxcyc) {
        if (verbose) 
            cat("WARNING:  maximum number of cycles reached before convergence\n")
    }
    stats <- as.numeric(stats)
    stats <- matrix(stats, cycle, 2 * k + 1, byrow = TRUE)
    colnames(stats) <- c("llik", paste("s^2_", Vcoef.names, sep = ""), 
        paste("der_", Vcoef.names, sep = ""))
    WX <- W %*% X
    XtWX <- crossprod(X, WX)
    cov <- XtWX
    cov <- solve(cov, cbind(t(WX), diag(1, dim(XtWX)[1])))
    beta.cov <- matrix(cov[, (dim(t(WX))[2] + 1):dim(cov)[2]], 
        dim(X)[2], dim(X)[2])
    cov <- matrix(cov[, 1:dim(t(WX))[2]], dim(X)[2], dim(X)[1])
    beta <- cov %*% y
    beta <- matrix(beta, length(beta), 1)
    row.names(beta) <- Xcolnames
    beta.se <- sqrt(abs(diag(beta.cov)))
    pos.cov <- (diag(beta.cov) < 0)
    beta.se[pos.cov] <- NA
    beta.se <- matrix(beta.se, length(beta.se), 1)
    row.names(beta.se) <- Xcolnames
    rms <- rss/rankQ
    fitted.values <- X %*% beta
    Q <- In - X %*% cov
    predicted <- NULL
    if (identity) {
        gam <- sigma[k]
        if (SWsolveINDICATOR) {
            Sigma <- 0
            for (i in 1:k) {
                Sigma <- Sigma + V[[i]] * sigma[i]
            }
        }
        predicted <- fitted.values + (Sigma - gam * In) %*% W %*% 
            (y - fitted.values)
    }
    sigma.cov <- (A.svd[1:k, 1:k] * 2)
    FI <- A/2
    FI.c <- matrix(0, dim(FI)[1], dim(FI)[2])
    FI.c <- FI/tcrossprod((sigma - 1) * pos + 1)
    sigma.cov <- ginv(FI.c)
    names(sigma) <- Vcoef.names
    rownames(sigma.cov) <- colnames(sigma.cov) <- Vcoef.names
    result <- list(trace = stats, llik = llik, cycle = cycle, 
        rdf = rankQ, beta = beta, beta.cov = beta.cov, beta.se = beta.se, 
        sigma = sigma[1:k], sigma.cov = sigma.cov[1:k, 1:k], 
        W = W, Q = Q, fitted = fitted.values, predicted = predicted, 
        pos = pos, Vnames = Vcoef.names, formula = formula, Vformula = Vformula, 
        Kcolnames = Kcolnames, model = model, Z = Z)
    class(result) <- "regress"
    result
}
