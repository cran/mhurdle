inner_prod <- function(x, y = NULL){
    x <- drop(x)
    if (is.null(y)) y <- x else y <- drop(y)
    sum(x * y)
}

outer_prod <- function(x, y = NULL){
    x <- drop(x)
    if (is.null(y)) y <- x else y <- drop(y)
    outer(x, y)
}

cross_prod <- function(x, y = NULL){
    x <- drop(x)
    if (is.null(y)) y <- x else y <- drop(y)
    drop(crossprod(x, y))
}


UP <- function(B, dx, dg)
    B + (inner_prod(dx, dg) + inner_prod(dg, cross_prod(B, dg))) *
        outer_prod(dx) / inner_prod(dx, dg) ^ 2 -
            (outer_prod(cross_prod(B, dg), dx) +
             t(outer_prod(cross_prod(B, dg), dx))) / inner_prod(dx, dg)

get_gradient <- function(x) - apply(attr(x, "gradient"), 2, sum)
get_gradientObs <- function(x) - attr(x, "gradient")
get_hessian <- function(x) if(! is.null(attr(x, "hessian"))) - attr(x, "hessian") else NULL
get_value <- function(x) - sum(as.numeric(x))
get_opg <- function(x) crossprod(get_gradientObs(x))

optim.ml <- function(f, start, g = NULL, H = NULL,
                     direction = c("max", "min"),
                     method = c("bfgs", "bhhh", "nr"),
                     gradient_tol = 1E-05,#sqrt(.Machine$double.eps),
                     value_tol = 1E-05,#sqrt(.Machine$double.eps),
                     step_tol = 1E-05,#sqrt(.Machine$double.eps),
                     iter_lim = 100, trace = 1L, ...){
    
    init_time <- proc.time()
    direction <- match.arg(direction)
    method <- match.arg(method)
    iter <- 0
    gradient_crit <- 10
    compute_hessian <- (method == "nr")
    x <- unname(start)
    K <- length(x)
    lnL <- f(x, gradient = TRUE, hessian = compute_hessian, ...)
    if (method == "bfgs") B <- diag(K);solve(get_opg(lnL))#diag(K)
    value <- get_value(lnL)
    if (trace > 0L) cat(paste("Initial value of lnL:", round(value, 2), "\n"))
    while(gradient_crit > gradient_tol & iter <= iter_lim){
        gradient <- get_gradient(lnL)
        iter <- iter + 1
        if (method == "bfgs") d_x <- - cross_prod(B, gradient)
        if (method == "bhhh") d_x <- - solve(get_opg(lnL), gradient)
        if (method == "nr") d_x <- - solve(get_hessian(lnL), gradient)
        step <- 1
        previous_lnL <- lnL
        lnL <- f(x + d_x, hessian = compute_hessian, ...)
        increasing <- get_value(lnL) < get_value(previous_lnL)
        if (increasing){
            while(get_value(lnL) < get_value(previous_lnL)){
                step <- step * 2
                previous_lnL <- lnL
                lnL <- f(x + step * d_x, hessian = compute_hessian, ...)
            }
            step <- step / 2
            lnL <- previous_lnL
        }
        else{
            while(get_value(lnL) > get_value(previous_lnL)){
                step <- step / 2
                print(step)
                previous_lnL <- lnL
                lnL <- f(x + step * d_x, hessian = compute_hessian, ...)
            }
        }
        x <- x + step * d_x
        d_g <- get_gradient(lnL) - gradient
                
        if (method == "bfgs") B <- UP(B, step * d_x, d_g)
#        gradient_crit <- inner_prod(get_gradient(lnL), solve(get_opg(lnL), get_gradient(lnL)))
        gradient_crit <- sqrt(mean(get_gradient(lnL) ^ 2))
        print(get_value(previous_lnL))
        print(get_value(lnL))
        value_crit <- 1 - get_value(previous_lnL) / get_value(lnL)
        if (is.na(value_crit)) value_crit <- 10
        print(get_value(lnL))
        print(value_crit)
        
        if (value_crit < value_tol) code <- 2L
        if (step < step_tol) code <- 3L
        if (iter >= iter_lim) code <- 4L
        if (gradient_crit < gradient_tol) code <- 1L
        if (trace > 0L) cat(paste("iteration:", iter, "lnL = ", round(get_value(lnL), 3),
                                  "step: ", round(step, 4), "crit:", round(gradient_crit, 3), "\n"))        
        if (trace > 1L){
            cat("--------------------------------------------\n")
            resdet <- rbind(param = x, gradient = get_gradient(lnL))
            print(round(resdet,3))
            cat("--------------------------------------------\n")
        }

        est_stat <- structure(list(method = method,
                                   gradient_crit = gradient_crit,
                                   value_crit = value_crit,
                                   step_crit = step,
                                   iter_crit = iter,
                                   code = code,
                                   time = proc.time() - init_time),
                              class = "est_stat")
        
        result <- structure(list(lnL = get_value(lnL), coefficients = x, contributions = as.numeric(lnL),
                                 gradientObs = get_gradientObs(lnL), gradient = get_gradient(lnL),
                                 opg = get_opg(lnL), hessian = get_hessian(lnL),
                                 est_stat = est_stat),
                            class = "ml_max")
    }
    result
}

print.est_stat <- function(x, ...){

    s <- round(x$time[1],0)
    h <- s %/% 3600
    s <- s - 3600 * h
    m <- s %/% 60
    s <- s - 60 * m
    cat(paste(x$method, "method\n"))
    tstr <- paste(h, "h:", m, "m:", s, "s", sep="")
    cat(paste(x$iter_crit,"iterations,",tstr,"\n"))
    cat(paste("g'(-H)^-1g =", sprintf("%5.3G", as.numeric(x$gradient_crit)),"\n"))
    msg <- switch(x$code,
                  "1" = "gradient close to zero",
                  "2" = "successive function values within tolerance limits",
                  "3" = "last step couldn't find higher value",
                  "4" = "iteration limit exceeded"
                  )
    cat(paste(msg, "\n"))
}

print.ml_max <- function(x, ...) print(x$est_stat)


summary.ml_max <- function(x){
    se <- sqrt(diag(solve(x$opg)))
    round(cbind(x$coefficients, se, student = x$coefficients / se))
}
