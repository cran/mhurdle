mhurdle.lnl <- function(param, X1, X2, X3, X4, y, gradient = FALSE,
                        fitted = FALSE, dist = NULL, corr = FALSE, robust = TRUE){
#    otime <- proc.time()
    if (robust){
        frho <- function(x) atan(x) * 2 / pi
        grho <- function(x) 2 / pi / (1 +  x ^ 2)
        fmu <- function(x) exp(x)
        gmu <- function(x) exp(x)
        fsd <- function(x) exp(x)
        gsd <- function(x) exp(x)
    }
    else{
        frho <- function(x) x
        grho <- function(x) rep(1, length(x))
        fmu <- function(x) x
        gmu <- function(x) 1
        fsd <- function(x) x
        gsd <- function(x) 1
    }

    N <- length(y)
    
    #  dummies for the existing equations hi and number of
    #  coefficients Ki
    h1 <- ! is.null(X1) ;  K1 <- ifelse(h1, ncol(X1), 0)
    h3 <- ! is.null(X3) ;  K3 <- ifelse(h3, ncol(X3), 0)
    h4 <- ! is.null(X4) ;  K4 <- ifelse(h4, ncol(X4), 0)
    # KR is either 1 (h1 or h3) or 3 (h1 and h3)
    KR <- ifelse(corr, h1 + h3 + h1 * h3, 0)
    K2 <- ncol(X2)
    # shape and scale parameters for box-cox and ihs transformations
    if (dist %in% c("ln2", "bc", "bc2", "ihs")){
        lambda <- param[K1 + K2 + K3 + 1 + K4 + KR + 1]
        if (dist %in% c("bc2", "ln2")){
            if (dist == "bc2") posmu <- K1 + K2 + K3 + 1 + K4 + KR + 2
            else posmu <- K1 + K2 + K3 + 1 + K4 + KR + 1
            mu <- fmu(param[posmu])
            gradientmu <- gmu(param[posmu])
        }
        if (dist == "bc") mu <- 0
    }
    if (dist %in% c("bc", "bc2")) sgn <- sign(lambda) else sgn <- + 1
    if (dist == "ln") mu <- 0
    # equation 1
    if (h1){
        beta1 <- param[1:K1]
        bX1 <- as.numeric(crossprod(t(X1), beta1))
        Phi1 <- pnorm(bX1) ; phi1 <- dnorm(bX1)
        Phi1 <- sanitize(Phi1, string = c("Phi1", "mhurdle.lnl"),
                         verbal = FALSE, replace = TRUE)
    }
    else{
        bX1 <- Inf ; beta1 <- NULL
        Phi1 <- 1 ; phi1 <- 0;
    }

    # equation 2
    beta2 <- param[(K1 + 1):(K1 + K2)]
    bX2 <- as.numeric(crossprod(t(X2), beta2))

    # equation 3
    if (h3){
        beta3 <- param[(K1 + K2 + 1):(K1 + K2 + K3)]
        bX3 <- as.numeric(crossprod(t(X3), beta3))
        Phi3 <- pnorm(bX3) ; phi3 <- dnorm(bX3)
        Phi3 <- sanitize(Phi3, string = c("Phi3", "mhurdle.lnl"),
                         verbal = FALSE, replace = TRUE)
    }
    else{
        bX3 <- Inf ; beta3 <- NULL
        Phi3 <- 1 ; phi3 <- 0
    }

    # standard deviation
    sd <- param[K1 + K2 + K3 + 1]
    gradientsd <- gsd(sd)
    sd <- fsd(sd)
    # equation 4
    if (h4){
        beta4 <- param[(K1 + K2 + K3 + 2):(K1 + K2 + K3 + K4 + 1)]
        bX4 <- as.numeric(crossprod(t(X4), beta4))
        pbX4 <- pnorm(bX4)
        sigma <- sd * pbX4
    }
    else{
        sigma <- sd
        pbX4 <- 1
    }

    # correlation coefficients
    if (corr){
        posrho <- (K1 + K2 + K3 + K4 + 2):(K1 + K2 + K3 + K4 + 1 + KR)
        rho <- param[posrho]
        # In case of only one correlation coefficient, use the whole
        # vector with only one non-null component
        if (h1 & ! h3) rho <- c(rho, 0, 0)
        if (! h1 & h3) rho <- c(0, 0, rho)
        gradientrho <- grho(rho)
        rho <- frho(rho)
        if (h1 & h3){
            # In case of a trivariate normal distribution, check the
            # joint relation of the three coefficients
            rho3 <- function(rho) rho[1] * rho[2] + c(-1, 1) *
                sqrt(1 + rho[1] ^ 2 * rho[2] ^ 2 - rho[1] ^ 2 - rho[2] ^ 2)
            if (rho[3] < rho3(rho)[1]) rho[3] <- rho3(rho)[1] + 1E-04 
            if (rho[3] > rho3(rho)[2]) rho[3] <- rho3(rho)[2] - 1E-04
        }
    }
    else rho <- rep(0, 3)
    # Transformation of the dependent variable
    Ty <- switch(dist,
                 "ln"  = log2(y) + log(Phi3),
                 "ln2" = log2(y * Phi3 + mu),
                 "bc"  = (exp(lambda * log(y * Phi3)) - 1) / lambda,
                 "bc2" = (exp(lambda * log(y * Phi3 + mu)) - 1) / lambda,
                 "ihs" = log(lambda * y * Phi3 + sqrt(1 + (lambda  * y * Phi3) ^ 2)) / lambda,
                 y * Phi3
                 )
    
    # logarithm of the jacobian
    lnJ <- switch(dist,
                  "ln"  = - log2(y),
                  "ln2" = log(Phi3) - log(mu + Phi3 * y),
                  "bc"  = (lambda - 1) * log2(y) + lambda * log(Phi3),
                  "bc2" = (lambda - 1) * log(Phi3 * y + mu) + log(Phi3),
                  "ihs" = - 0.5 * log(1 + (lambda * Phi3 * y) ^ 2) + log(Phi3),
                  log(Phi3)
                  )
    
    # derivative of lnJ respective with lambda
    lnJlb <- switch(dist,
                    "bc"  = log2(y) + log(Phi3),
                    "bc2" = log(Phi3 * y + mu),
                    "ihs" = - lambda * y ^ 2 * Phi3 ^ 2 / (1 + (lambda * y * Phi3) ^ 2)
                    )

    # derivative of lnJ respective with mu
    if (dist == "ln2") lnJmu <- - 1 / (mu + Phi3 * y)
    if (dist == "bc2") lnJmu <- (lambda - 1) / (Phi3 * y + mu)
    
    # the  residual of the consumption equation
    resid <- (Ty - bX2)
    # problem with bc and lambda < 0, for y = 0, resid = -inf and
    # log(dnorm(resid)) = -inf
    resid[y == 0] <- 0


    # The opposite of the correspondant z values are then computed

    mzn <- + Inf
    if (dist %in% c("bc", "bc2") && lambda > 0) mzn <- (1 / lambda + bX2) / sigma
    if (dist == "tn") mzn <- - bX2 / sigma
    # A verifier
    if (dist == "tn") mzn <- bX2 / sigma

    mzx <- - Inf
    if (dist %in% c("bc", "bc2") && lambda < 0) mzx <- (1 / lambda + bX2) / sigma
    
    mz0 <- + Inf
    if (dist == "bc2") mz0 <- ( - (mu ^ lambda - 1) / lambda + bX2) / sigma
    if (dist == "bc" && lambda > 0)  mz0 <- (1 / lambda + bX2) / sigma
    if (dist == "ln2") mz0 <- (bX2 - log(mu)) / sigma
    if (dist %in% c("n", "tn")) mz0 <- bX2 / sigma
    
    
    # compute the relevant bivariate and trivariate cumulative normals
    if (h1) arg1 <- (bX1 + rho[1] * resid / sigma) / sqrt(1 - rho[1] ^ 2) else arg1 <- Inf
    if (h3) arg3 <- (bX3 + rho[3] * resid / sigma) / sqrt(1 - rho[3] ^ 2) else arg3 <- Inf
    Pr123A <- PHI3(bX1, mz0, bX3, rho)
    Pr123B <- PHI3(bX1, mzx, bX3, rho)
    
    Pr13 <- PHI2(arg1, arg3,
                 (rho[2] - rho[1] * rho[3]) / sqrt(1 - rho[1] ^ 2) / sqrt(1 - rho[3] ^ 2)
                 )
    Pr13$f <- sanitize(Pr13$f, string = c("Pr13", "mhurdle.ln"),
                       verbal = FALSE, replace = TRUE)
    # PI is the correction of the truncature
    PI <- punorm(mzn) - punorm(mzx)
    Pplus <- (Pr123A$f - Pr123B$f) / PI
    Numerator <- PI - Pr123A$f + Pr123B$f
    lnL.null <- log(Numerator) - log(PI)
    lnL.one <- log(PI - Numerator) - log(PI)
    lnL.null[y != 0] <- 0
    lnL.pos <-
        - log(sigma) +
            dnorm(resid / sigma, log = TRUE) +
                log(Pr13$f) +
                    lnJ - log(PI)
    lnL.pos[y == 0] <- 0    
    lnL <- lnL.null * (y == 0) + lnL.pos * (y != 0)
    if (any(is.na(lnL) | any(is.infinite(lnL)))){
        lnL[is.na(lnL)] <- 0
        warnings("infinite or missing values of lnLi")
    }
    attr(lnL, "parts") <- c(lnLNull = sum(lnL.null[y == 0]),
                            lnLOne  = sum(lnL.one[y != 0]),
                            lnLPos  = sum(lnL.pos[y != 0]))

    if (gradient){
        gradi <- c()
        
        # derivatives respective to beta1
        if (h1){
            lnL.beta1 <- (y == 0) * ( (- Pr123A$d1 + Pr123B$d1) / Numerator) +
                (y != 0) * ( Pr13$d1 / sqrt(1 - rho[1] ^ 2) / Pr13$f)
            gradi <- cbind(gradi, lnL.beta1 * X1)
        }


        # derivatives respective to beta2
        PI2 <-  (dnorm(mzn) - dnorm(mzx) ) / sigma
        lnL.beta2 <- (y == 0) * ( (PI2 - Pr123A$d2 / sigma + Pr123B$d2 / sigma) / Numerator -
                                     PI2 / PI) +
            (y != 0) * ( resid / sigma ^ 2 -
                            ( Pr13$d1 * rho[1] / sigma / sqrt(1 - rho[1] ^ 2) +
                                 Pr13$d2 * rho[3] / sigma / sqrt(1 - rho[3] ^ 2) ) / Pr13$f -
                                     PI2 / PI)
        gradi <- cbind(gradi, lnL.beta2 * X2)
        
        # derivatives respective to beta3
        if (h3){
            # derivative of T(Phi3 y) with bX3
            Ty3 <- switch(dist,
                          "ln" = mills(bX3),
                          "ln2" = y * phi3 / (mu + y * Phi3),
                          "bc" = exp(lambda * log(y * Phi3)) * mills(bX3),
                          "bc2" = exp((lambda - 1) * log(y * Phi3 + mu)) * phi3 * y,
                          "ihs" = y * phi3 / sqrt( 1 + (lambda * y * Phi3) ^ 2),
                          y * phi3
                          )
            # derivative of lnJ with bX3
            lnJ3 <- switch(dist,
                           "ln" = 0,
                           "ln2" = mills(bX3) - y * phi3 / (mu + y * Phi3),
                           "bc" = lambda * mills(bX3),
                           "bc2" = (lambda - 1) * phi3 * y / (Phi3 * y + mu) + mills(bX3),
                           "ihs" = - phi3 * Phi3 * lambda ^ 2 * y ^ 2 /
                               (1 + (lambda * y * Phi3) ^ 2) + mills(bX3),
                           mills(bX3)
                          )
            
            lnL.beta3 <- (y == 0) * (- Pr123A$d3 + Pr123B$d3) / Numerator + 
                (y != 0) * (- resid / sigma ^ 2 * Ty3 +
                                ( Pr13$d1 * Ty3 * rho[1] / sigma / sqrt(1 - rho[1] ^ 2) +
                                     Pr13$d2 * (1 + Ty3 * rho[3] / sigma) /
                                         sqrt(1 - rho[3] ^ 2) ) / Pr13$f +
                                             lnJ3
                            )
            gradi <- cbind(gradi, lnL.beta3 * X3)
        }
        
        # derivatives respective to sigma
        # PI is only relevant for truncated normal and box-cox transformation
        PIs <- 0
        Pr123Bs <- 0
#        print(mz0);stop()
#        if (is.infinite(mz0)) Pr123As <- 0 else Pr123As <- Pr123A$d2 * (- mz0 / sigma)
        Pr123As <- ifelse(is.infinite(mz0), 0, Pr123A$d2 * (- mz0 / sigma))
        if (dist == "tn" | (dist %in% c("bc", "bc2") && lambda > 0)) PIs <- PIs + dnorm(mzn) * (- mzn / sigma)
        if (dist %in% c("bc", "bc2") && lambda < 0){
            PIs <- PIs - dnorm(mzx) * (- mzx / sigma)
            Pr123Bs <- Pr123B$d2 * (- mzx / sigma)
        }
        lnL.sigma <- (y == 0) * ( (PIs - Pr123As + Pr123Bs) / Numerator - PIs / PI) +
            (y != 0) * (- 1 / sigma + resid ^ 2 / sigma ^ 3 +
                            ( - Pr13$d1 * rho[1] * resid / sigma ^ 2 / sqrt(1 - rho[1] ^ 2) -
                                 Pr13$d2 * rho[3] * resid / sigma ^ 2 /
                                     sqrt(1 - rho[3] ^ 2))  / Pr13$f -
                                         PIs / PI)
        gradi <- cbind(gradi, sigma = lnL.sigma * pbX4  * gradientsd)
        
        # derivatives respective to beta4
        if (!is.null(X4)){
            gradi <- cbind(gradi,  lnL.sigma * sd * dnorm(bX4) * X4)
        }
        
        # derivatives respective to rho
        if (corr){
            if (h1 & ! h3){
                Pr123A$dr <- cbind(Pr123A$dr, 0, 0)
                Pr123B$dr <- cbind(Pr123B$dr, 0, 0)
            }
            if (! h1 & h3){
                Pr123A$dr <- cbind(0, 0, Pr123A$dr)
                Pr123B$dr <- cbind(0, 0, Pr123B$dr)
            }
            if ( (h1 & h3) & !is.matrix(Pr123A$dr) ) Pr123A$dr <- cbind(0, Pr123A$dr, 0)

            if (length(Pr123B$dr) == 1 && Pr123B$dr == 0) Pr123B$dr <- cbind(0, 0, 0)
            lnL.rho <- c()
            if (h1){
                Drho12 <- (rho[1] * rho[2] - rho[3]) /
                    (1 - rho[1] ^ 2) ^ 1.5  / sqrt(1 - rho[3] ^ 2)
                lnL.rho12 <- (y == 0) * (- Pr123A$dr[, 1] + Pr123B$dr[, 1]) / Numerator  +
                    (y != 0) * ( Pr13$d1 * (resid / sigma + rho[1] / (1 - rho[1] ^ 2) *
                                                (bX1 + rho[1] * resid / sigma) ) /
                                    (1 - rho[1] ^ 2) ^ 0.5 +
                                        Pr13$dr * Drho12) / Pr13$f
                lnL.rho <- cbind(lnL.rho, lnL.rho12 * gradientrho[1])
            }
            if (h1 & h3){
                Drho13 <- 1 / sqrt(1 - rho[1] ^ 2) / sqrt(1 - rho[3] ^ 2)
                lnL.rho13 <- (y == 0) * (- Pr123A$dr[, 2] + Pr123B$dr[, 2]) / Numerator  +
                    (y != 0) * ( Pr13$dr * Drho13) / Pr13$f
                lnL.rho <- cbind(lnL.rho, lnL.rho13 * gradientrho[2])
            }
            if (h3){
                Drho23 <- (rho[3] * rho[2] - rho[1]) /
                    (1 - rho[3] ^ 2) ^ 1.5  / sqrt(1 - rho[1] ^ 2)
                lnL.rho23 <- (y == 0) * (- Pr123A$dr[, 3] + Pr123B$dr[, 3]) / Numerator  +
                    (y != 0) * ( Pr13$d2 * (resid / sigma + rho[3] / (1 - rho[3] ^ 2) *
                                                (bX3 + rho[3] * resid / sigma) ) /
                                    (1 - rho[3] ^ 2) ^ 0.5 +
                                        Pr13$dr * Drho23) / Pr13$f
                lnL.rho <- cbind(lnL.rho, lnL.rho23 * gradientrho[3])
            }
            gradi <- cbind(gradi, lnL.rho)
        }
        # derivative respective to lambda
        if (dist %in% c("bc", "bc2")){
            Tylb <- (log(Phi3 * y + mu) * (Phi3 * y + mu) ^ lambda - Ty) / lambda
            if (mu == 0) T0lb <- ( 1 / lambda ^ 2 * (lambda > 0) + 0 * (lambda < 0))
            else T0lb <- (log(mu) * mu ^ lambda - (mu ^ lambda - 1) / lambda) / lambda
            Tymaxlb <- 0 * (lambda > 0) + (1 / lambda ^ 2) * (lambda < 0)
            PIl <- - sign(lambda) * dnorm( (bX2 + 1 / lambda) / sigma ) / (sigma * lambda ^ 2)
            lnL.lambda <- vector(mode = "numeric", length = length(y))
            
            lnL.lambda[y == 0] <- ( (PIl - Pr123A$d2 * (- T0lb /  sigma) + Pr123B$d2 *
                                         (- Tymaxlb / sigma)) / Numerator - PIl / PI)[y == 0]

            
            lnL.lambda[y > 0] <- ( (- resid / sigma ^ 2 +
                                        (Pr13$d1 * rho[1] / sqrt(1 - rho[1] ^ 2) / sigma +
                                             Pr13$d2 * rho[3] / sqrt(1 - rho[3] ^ 2) / sigma) / Pr13$f ) *
                                      Tylb + lnJlb - PIl / PI)[y > 0]
            gradi <- cbind(gradi, tr = lnL.lambda)
        }

        if (dist == "bc2"){
            Tymu <- (Phi3 * y + mu) ^ (lambda - 1)
            T0mu <- mu ^ (lambda - 1)
            lnL.mu <- vector(mode = "numeric", length = length(y))
            lnL.mu[y == 0] <- (- Pr123A$d2 * (- T0mu / sigma) / Numerator)[y == 0]
            lnL.mu[y > 0] <- (( - resid / sigma ^ 2 +
                                   (Pr13$d1 * rho[1] / sqrt(1 - rho[1] ^ 2) / sigma +
                                        Pr13$d2 * rho[3] / sqrt(1 - rho[3] ^ 2) / sigma) / Pr13$f
                               ) * Tymu + lnJmu)[y != 0]
            gradi <- cbind(gradi, mu = lnL.mu * gradientmu)
        }

        if (dist == "ihs"){
            Tylb <- (y * Phi3) / lambda / sqrt(1 + (lambda * y * Phi3) ^ 2) - Ty / lambda
            lnL.lambda <- vector(mode = "numeric", length = length(y))
            ## lnL.lambda[y != 0] <- (( -resid / sigma ^ 2 +  mills(z1$f) * rho1 / sigma / sqrt(1 - rho1 ^ 2) +
            ##                             mills(z3$f) * rho3 / sigma / sqrt(1 - rho3 ^ 2) ) * Tylb + lnJlb)[y != 0]
            lnL.lambda[y != 0]
            gradi <- cbind(gradi, lnL.lambda)
        }
        
        if (dist == "ln2"){
            Tymu <- 1 / (mu + Phi3 * y)
            T0mu <- 1 / mu
            lnL.mu <- vector(mode = "numeric", length = length(y))
            lnL.mu[y == 0] <- ( - Pr123A$d2 * (- T0mu / sigma ) / Numerator)[y == 0]
            lnL.mu[y != 0] <- (( - resid / sigma ^ 2 +
                                    (Pr13$d1 * rho[1] / sqrt(1 - rho[1] ^ 2) / sigma +
                                         Pr13$d2 * rho[3] / sqrt(1 - rho[3] ^ 2) / sigma) / Pr13$f
                                ) * Tymu + lnJmu)[y != 0]
            gradi <- cbind(gradi, mu = lnL.mu * gradientmu)
        }

        if (any(is.na(gradi))) gradi[is.na(gradi)] <- 0
        attr(lnL, "gradient") <- gradi
    }
    if (fitted){
        if (dist %in% c("bc", "bc2")){
            if (length(mzn) == 1) mzn <- rep(mzn, length(y))
            if (length(mzx) == 1) mzx <- rep(mzx, length(y))
            if (length(Phi3) == 1) Phi3 <- rep(Phi3, length(y))
            if (length(Phi1) == 1) Phi1 <- rep(Phi1, length(y))
            
            resid <- function(z, index){
                switch(dist,
                       "ln"  = log2(z) + log(Phi3[index]),
                       "ln2" = log2(z * Phi3[index] + mu),
                       "bc"  = (exp(lambda * log(z * Phi3[index])) - 1) / lambda,
                       "bc2" = (exp(lambda * log(z * Phi3[index] + mu)) - 1) / lambda,
                       "ihs" = log(lambda * z * Phi3[index] + sqrt(1 + (lambda  * z * Phi3[index]) ^ 2)) / lambda,
                       z * Phi3[index]
                       ) - bX2[index]
            }
            lnJ <- function(z, index){
                switch(dist,
                       "ln"  = - log2(z),
                       "ln2" = log(Phi3[index]) - log(mu + Phi3[index] * z),
                       "bc"  = (lambda - 1) * log2(z) + lambda * log(Phi3[index]),
                       "bc2" = (lambda - 1) * log(Phi3[index] * z + mu) + log(Phi3[index]),
                       "ihs" = - 0.5 * log(1 + (lambda * Phi3[index] * z) ^ 2) + log(Phi3[index]),
                       log(Phi3[index])
                       )
            }
            
            arg1 <- function(z, index) if (h1) (bX1[index] + rho[1] * resid(z, index) / sigma) / sqrt(1 - rho[1] ^ 2) else Inf
            arg3 <- function(z, index) if (h3) (bX3[index] + rho[3] * resid(z, index) / sigma) / sqrt(1 - rho[3] ^ 2) else Inf
            Pr13 <- function(z, index) PHI2(arg1(z, index), arg3(z, index),
                                            (rho[2] - rho[1] * rho[3]) / sqrt(1 - rho[1] ^ 2) / sqrt(1 - rho[3] ^ 2))$f
            # Pbnorm plante et pas PHI2        
            PI <- function(index) punorm(mzn[index]) - punorm(mzx[index])
            fplus <- function(z, index){
                x <- z * exp(- log(sigma) + dnorm(resid(z, index) / sigma, log = TRUE) + log(Pr13(z, index)) + lnJ(z, index) - log(PI(index)))
                x
            }        
            if (dist %in% c("bc", "bc2") && lambda < 0) maxint <- max(y) * 3
            else maxint <- +Inf
            E <- sapply(seq_len(length(y)), function(i) integrate(function(x) fplus(x, i), 0, maxint)$value)
        }

        if (dist %in% c("ln", "ln2")){
            if (h1) arg1 <- bX1 + rho[1] * sigma else arg1 <- + Inf
            if (h3) arg3 <- bX3 + rho[3] * sigma else arg3 <- + Inf
            E <- exp(bX2 + 0.5 * sigma ^ 2) / Phi3 * ptnorm(arg1, mz0 + sigma, arg3, rho)  -
                mu * ptnorm(bX1, mz0, bX3, rho)/ Phi3
        }
        if (dist %in% c("n", "tn")){
            phi2 <- dnorm(bX2 / sigma)
            phi1 <- dnorm(bX1)
            phi3 <- dnorm(bX3)
            Pr13 <- pbnorm((bX1 - rho[1] * bX2 / sigma) / sqrt(1 - rho[1] ^ 2),
                           (bX3 - rho[3] * bX2 / sigma) / sqrt(1 - rho[3] ^ 2),
                           (rho[2] - rho[1] * rho[3]) / sqrt( (1 - rho[1] ^ 2) * (1 - rho[3] ^ 2)))
            Pr23 <- pbnorm((bX2 / sigma - rho[1] * bX1) / sqrt(1 - rho[1] ^ 2),
                           (bX3         - rho[2] * bX1) / sqrt(1 - rho[2] ^ 2),
                           (rho[3] - rho[1] * rho[2]) / sqrt( (1 - rho[1] ^ 2) * (1 - rho[2] ^ 2)))
            
            Pr12 <- pbnorm((bX1         - rho[2] * bX3) / sqrt(1 - rho[2] ^ 2),
                           (bX2 / sigma - rho[3] * bX3) / sqrt(1 - rho[3] ^ 2),
                           (rho[1] - rho[2] * rho[3]) / sqrt(1 - rho[2] ^ 2) / sqrt(1 - rho[3] ^ 2)
                           )
            E <- bX2 / Phi3 + sigma / (Phi3 * Pplus) * (phi2 * Pr13 +
                                                                rho[1] * phi1 * Pr23 / sqrt( (1 - rho[1] ^ 2) * (1 - rho[3] ^ 2)) +
                                                                    rho[3] * phi3 * Pr12 / sqrt( (1 - rho[2] ^ 2) * (1 - rho[3] ^ 2)))

            E <- bX2 / Phi3 + sigma / (Phi3 * Pplus) * (phi2 * Pr13 + rho[1] * phi1 * Pr23 + rho[3] * phi3 * Pr12)

        }
        attr(lnL, "fitted") <- cbind(zero = 1 - Pplus, 
                                     pos = E / Pplus)
    }
    lnL
}
