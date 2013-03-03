library(ensembleBMA)

ensembleBMAgamma

####################################

function (ensembleData, trainingDays, dates = NULL, control = controlBMAgamma(), 
    exchangeable = NULL) 
{
    if (!inherits(ensembleData, "ensembleData")) 
        stop("not an ensembleData object")
    if (missing(trainingDays)) 
        stop("trainingDays must be specified")
    call <- match.call()
    warmStart <- FALSE
    if (is.list(trainingDays)) 
        trainingsDays <- trainingDays[[1]]
    ensMemNames <- ensembleMembers(ensembleData)
    nForecasts <- length(ensMemNames)
    exchangeable <- getExchangeable(exchangeable, ensembleGroups(ensembleData), 
        nForecasts)
    M <- !dataNA(ensembleData)
    if (!all(M)) 
        ensembleData <- ensembleData[M, ]
    nObs <- nrow(ensembleData)
    if (!nObs) 
        stop("no data")
    Dates <- as.character(ensembleValidDates(ensembleData))
    DATES <- sort(unique(Dates))
    julianDATES <- ymdhTOjul(DATES)
    incr <- min(1, min(diff(julianDATES)))
    forecastHour <- ensembleFhour(ensembleData)
    lag <- ceiling(forecastHour/24)
    dates <- getDates(DATES, julianDATES, dates, trainingDays, 
        lag, incr)
    juliandates <- ymdhTOjul(dates)
    nDates <- length(dates)
    biasCoefs <- matrix(NA, nrow = 2, ncol = nDates, dimnames = list(NULL, 
        dates))
    varCoefs <- array(NA, c(2, nDates), dimnames = list(NULL, 
        dates))
    weights <- array(NA, c(nForecasts, nDates), dimnames = list(ensMemNames, 
        dates))
    trainTable <- rep(0, nDates)
    names(trainTable) <- dates
    nIter <- loglikelihood <- rep(0, nDates)
    names(nIter) <- names(loglikelihood) <- dates
    K <- 1:nForecasts
    L <- length(juliandates)
    twin <- 1:trainingDays
    cat("\n")
    l <- 0

    # compute training set 
    for (i in seq(along = juliandates)) {
        I <- (juliandates[i] - lag * incr) >= julianDATES
        if (!any(I)) 
            stop("insufficient training data")
        j <- which(I)[sum(I)]
        if (j != l) {
            twin <- (j + 1) - (1:trainingDays)
            D <- as.logical(match(Dates, DATES[twin], nomatch = 0))
            if (!any(D)) 
                stop("this should not happen")
            d <- ensembleValidDates(ensembleData[D, ])
            if (length(unique(d)) != trainingDays) 
                stop("wrong # of training days")
            cat("modeling for date", dates[i], "...")
            kNA <- apply(ensembleForecasts(ensembleData[D, ]), 
                2, function(x) all(is.na(x)))
            if (any(kNA)) {
                if (!is.null(x <- exchangeable)) 
                  x <- exchangeable[-K[kNA]]

                  ######################

                fit <- fitBMAgamma(ensembleData[D, -K[kNA]], 
                  control = control, exchangeable = x)
            }
            else {
                fit <- fitBMAgamma(ensembleData[D, ], control = control, 
                  exchangeable = exchangeable)

                ########################
            }
            l <- j
            trainTable[i] <- length(unique(Dates[D]))
            nIter[i] <- fit$nIter
            loglikelihood[i] <- fit$loglikelihood
            if (warmStart) 
                control$start$weights <- weights[, i]
            cat("\n")
        }
        else {
            trainTable[i] <- -abs(trainTable[i - 1])
            nIter[i] <- -abs(nIter[i - 1])
            loglikelihood[i] <- loglikelihood[i - 1]
        }
        biasCoefs[, i] <- fit$biasCoefs
        varCoefs[, i] <- fit$varCoefs
        weights[K[!kNA], i] <- fit$weights
    }
    structure(list(training = list(days = trainingDays, lag = lag, 
        table = trainTable), biasCoefs = biasCoefs, varCoefs = varCoefs, 
        weights = weights, nIter = nIter, exchangeable = exchangeable, 
        power = fit$power), forecastHour = forecastHour, initializationTime = ensembleItime(ensembleData), 
        call = match.call(), class = c("ensembleBMAgamma", "ensembleBMA"))
}
<environment: namespace:ensembleBMA>

####################################



fitBMAgamma
####################################
function (ensembleData, control = controlBMAgamma(), exchangeable = NULL) 
{
    ZERO <- 1e-100
    powfun <- function(x, power) x^power

    # check whether observations are exchangeable (in our dataset, NO)
    if (is.null(exchangeable)) 
        exchangeable <- ensembleGroups(ensembleData)
    if (length(unique(exchangeable)) == length(exchangeable)) 
        exchangeable <- NULL
    if (!(nullX <- is.null(exchangeable))) {
        namX <- as.character(exchangeable)
        uniqueX <- unique(namX)
        nX <- length(uniqueX)
    }

    # control object is result of controlBMAgamma -- see below
    maxIter <- control$maxIter
    tol <- eps <- control$tol
    ensembleData <- ensembleData[!dataNA(ensembleData, dates = FALSE), 
        ]
    nObs <- dataNobs(ensembleData)
    if (!nObs) 
        stop("no observations")
    obs <- dataVerifObs(ensembleData)
    if (is.null(startup <- dataStartupSpeed(ensembleData))) {
        if (is.null(control$startupSpeed)) 
            stop("default anemometer startup speed not specified")
        startup <- control$startupSpeed
    }
    if (length(startup) != nrow(ensembleData)) {
        startup <- rep(startup, length = nrow(ensembleData))
    }
    else if (length(startup) != 1) 
        stop("startup speed improperly specified")
    if (any(is.na(startup))) {
        if (is.null(control$startupSpeed)) 
            stop("default anemometer startup speed not specified")
        startup[is.na(startup)] <- control$startupSpeed
    }
    ensMemNames <- ensembleMembers(ensembleData)
    nForecasts <- length(ensMemNames)
    Y0 <- obs == 0
    if (sum(!Y0) < 2) 
        stop("less than 2 nonzero obs")
    ensembleData <- ensembleForecasts(ensembleData)
    ensembleData <- as.matrix(apply(ensembleData, 2, powfun, 
        power = control$power))

    # transform data with power function (we won't do this)
    obs <- as.vector(sapply(obs, powfun, power = control$power))

    # fit linear model
    lmFunc <- function(x, y) {
        beta0 <- min(y)
        x <- as.matrix(x)
        n <- ncol(x)
        x <- as.vector(x)
        nax <- is.na(x) # listwise deletion of missing obs
        x <- x[!nax]
        y <- rep(y, n)[!nax]
        if (all(!x)) {
            fit <- list(coefficients = c(mean(y), 0), fitted.values = rep(mean(y), 
                length(y)))
        }
        else {
            fit <- lm(y ~ x)
            coefs <- fit$coefficients
            if (coefs[1] <= 0) {
                coefs[1] <- beta0
                coefs[2] <- sum((y - beta0) * x)/sum(x * x)
                fit$coefficients <- coefs
                fit$fitted.values <- cbind(1, x) %*% coefs
            }
        }
        fit
    }

    # p. 3212, eqn 3
    meanFit <- lmFunc(ensembleData, obs)
    biasCoefs <- meanFit$coefficients
    meanVec <- as.vector(ensembleData)
    meanVec[!is.na(meanVec)] <- meanFit$fitted.values
    MEAN <- matrix(meanVec, nObs, nForecasts)
    miss <- is.na(ensembleData)
    Mzero <- miss[Y0, , drop = FALSE]
    Mnonz <- miss[!Y0, , drop = FALSE]
    completeDataLLmiss <- function(z, w, m, X, obs, startup) {
        objective <- function(par) {
            nObs <- length(obs)
            nFor <- ncol(X)
            miss <- is.na(X)
            Y0 <- obs == 0
            Mzero <- miss[Y0, , drop = FALSE]
            Mnonz <- miss[!Y0, , drop = FALSE]
            W <- matrix(w, nObs, nFor, byrow = TRUE)
            W[miss] <- 0
            W <- sweep(W, MARGIN = 1, FUN = "/", STATS = apply(W, 
                1, sum))
            v <- (par[1]^2 + (par[2]^2) * X)^2
            r <- m/v
            q <- array(NA, dim(z))
            q[Y0, ][!Mzero] <- log(pgamma(startup[Y0], shape = (r * 
                m)[Y0, , drop = FALSE][!Mzero], rate = r[Y0, 
                , drop = FALSE][!Mzero]))
            q[!Y0, ][!Mnonz] <- dgamma(matrix(obs, nObs, nFor)[!Y0, 
                , drop = FALSE][!Mnonz], shape = (r * m)[!Y0, 
                ][!Mnonz], rate = r[!Y0, , drop = FALSE][!Mnonz], 
                log = TRUE)
            Wzero <- W == 0
            include <- !miss & !Wzero
            -sum(z[include] * (q[include] + log(W[include])))
        }
        objective
    }

    # "we restricted the variance parameters to be constant across all ensemble members" (p. 3212, section (c) )
    varCoefs <- if (is.null(control$init$varCoefs)) 
        c(1, 1)
    else control$init$varCoefs
    varCoefs <- pmax(varCoefs, 1e-04)
    names(varCoefs) <- names(weights) <- NULL

    # weights are initially uniform unless otherwise specified
    weights <- if (is.null(control$init$weights)) 
        1
    else control$init$weights
    if (length(weights) == 1) 
        weights <- rep(weights, nForecasts)
    weights <- weights/sum(weights) # normalize
    weights <- pmax(weights, 1e-04) # make at least 0.0001
    weights <- weights/sum(weights) # normalize again
    if (!is.null(names(weights))) 
        weights <- weights[ensMemNames]
    if (!nullX) {
        for (labX in uniqueX) {
            I <- namX == labX
            weights[I] <- mean(weights[I]) # ? 
        }
    }
    nIter <- 0
    z <- matrix(1/nForecasts, ncol = nForecasts, nrow = nObs)
    objold <- 0
    newLL <- 0

    # "we estimate w_k; c_0; and c_1 by... maximum likelihood" - p. 3213
    while (TRUE) {
    		# "we maximize it numerically using the expectation-maximization algorithm"
    		# expectation step
        VAR = (varCoefs[1] + varCoefs[2] * ensembleData)^2
        RATE <- MEAN/VAR
        SHAPE <- RATE * MEAN
        # estimate z_j+1 given w_j
        z[Y0, ][!Mzero] <- pgamma(startup[Y0], shape = SHAPE[Y0, 
            ][!Mzero], rate = RATE[Y0, ][!Mzero])
        z[!Y0, ][!Mnonz] <- dgamma(matrix(obs, nObs, nForecasts)[!Y0, 
            ][!Mnonz], shape = SHAPE[!Y0, ][!Mnonz], rate = RATE[!Y0, 
            ][!Mnonz], log = TRUE)
        zmax = apply(z[!Y0, ], 1, max, na.rm = TRUE)
        z[!Y0, ] <- exp(sweep(z[!Y0, ], MARGIN = 1, FUN = "-", 
            STATS = zmax))
        z <- sweep(z, MARGIN = 2, FUN = "*", STATS = weights)
        oldLL <- newLL
        newLL <- sum(zmax + log(apply(z[!Y0, , drop = FALSE], 
            1, sum, na.rm = TRUE)))
        newLL <- newLL + sum(log(apply(z[Y0, , drop = FALSE], 
            1, sum, na.rm = TRUE)))
        z <- z/apply(z, 1, sum, na.rm = TRUE)
        z[z < ZERO] <- 0
        # estimate w_j+1 given z_j+1 
        wold <- weights
        zsum2 <- apply(z, 2, sum, na.rm = TRUE)
        weights <- zsum2/sum(zsum2) # weight the observations by how much they improve prediction (eqn. 5, p. 3213-4)
        weights[weights < ZERO] <- 0
        if (!nullX) {
            weights <- sapply(split(weights, namX), mean)[namX]
        }

        # maximization step 
        weps <- max(abs(wold - weights)/(1 + abs(weights)))
        fn <- completeDataLLmiss(z, weights, MEAN, ensembleData, 
            obs, startup) # specifying the loss function
        optimResult = optim(sqrt(varCoefs), fn = fn, method = "BFGS")
        if (optimResult$convergence) 
            warning("optim does not converge")
        varOld <- varCoefs
        varCoefs <- optimResult$par^2
        veps <- max(abs(varOld - varCoefs)/(1 + abs(varCoefs)))
        ERROR <- abs(objold - optimResult$value)/(1 + abs(optimResult$value))
        objold <- optimResult$value
        if (nIter > 0) {
            error <- abs(oldLL - newLL)/(1 + abs(newLL))
            if (error < eps) 
                break # iterate until convergence...
        }
        nIter <- nIter + 1
        if (nIter >= maxIter) 
            break # ... or until max iterations are reached
    }
    if (nIter >= maxIter && error > eps) 
        warning("iteration limit reached")
    names(biasCoefs) <- NULL
    names(weights) <- ensMemNames
    startup <- unique(startup)
    if (length(startup) > 1) 
        startup <- NA
    structure(list(biasCoefs = biasCoefs, varCoefs = varCoefs, 
        weights = weights, nIter = nIter, loglikelihood = newLL, 
        power = control$power, startupSpeed = startup), class = c("fitBMAgamma", 
        "fitBMA"))
}
<environment: namespace:ensembleBMA>

####################################



controlBMAgamma

####################################

# tol: smallest number the machine can represent
# maxIter: largest integer machine can represent (2147483647)

function (maxIter = Inf, tol = sqrt(.Machine$double.eps), power = 1, 
    startupSpeed = NULL, init = list(varCoefs = NULL, weights = NULL)) 
{
    if (is.infinite(maxIter) && maxIter > 0) 
        maxIter <- .Machine$integer.max
    list(maxIter = maxIter, tol = tol, power = power, startupSpeed = startupSpeed, 
        init = init)
}
<environment: namespace:ensembleBMA>

####################################