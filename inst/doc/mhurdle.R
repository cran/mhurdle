## ----label = setup, include = FALSE----------------------------
knitr::opts_chunk$set(echo = TRUE, cache = FALSE, message = FALSE, warning = FALSE, widtht = 65)
options(width = 65)

## ----eval = TRUE, label = Data---------------------------------
library("mhurdle")
data("Interview", package = "mhurdle")

## --------------------------------------------------------------
round(c(mean = mean(Interview$shows),
        "% of 0" = mean(Interview$shows == 0), 
        "mean pos" = mean(Interview$shows) / mean(Interview$shows > 0)),
      2)

## ----label = tobit---------------------------------------------
N010I <- mhurdle(shows ~ 0 | linc + smsa + age + 
                     educ + size, data = Interview, 
                 h2 = TRUE, dist = "n", method = "bhhh")
L010I <- update(N010I, dist = "ln")

## ----label = selection-----------------------------------------
L100D <- mhurdle(shows ~ smsa + age + educ + size | 
                     linc, data = Interview, 
                 h2 = FALSE, dist = "ln", corr = TRUE, method = "bhhh", 
                 finalHessian = TRUE)
L100D2 <- update(L100D, start = coef(L100D), robust = FALSE)
L110D <- update(L100D, h2 = TRUE)
L110D2 <- update(L110D, start = coef(L110D), robust = FALSE)
L001D <- mhurdle(shows ~ 0 | linc | smsa + age +  
                       educ + size, data = Interview, 
                   h2 = FALSE, corr = TRUE, method = "bhhh", 
                 finalHessian = TRUE)
L001D2 <- update(L001D, start = coef(L001D), robust = FALSE)
L011D <- update(L001D, h2 = TRUE)
L011D2 <- update(L011D, start = coef(L011D), robust = FALSE)

## ----label = threeeq-------------------------------------------
L101D <- mhurdle(shows ~ educ + size | linc | 
                        smsa + age, data = Interview, 
                    h2 = FALSE, method = "bhhh", corr = TRUE, 
                 finalHessian = TRUE)
L101D2 <- update(L101D, start = coef(L101D), robust = FALSE)
L111D <- update(L101D, h2 = TRUE)
L111D2 <- update(L111D, start = coef(L111D), robust = FALSE)
L111I <- update(L111D, corr = FALSE)
L111I2 <- update(L111I, start = coef(L111I), robust = FALSE)

## ----label = estimations, echo = FALSE, eval = TRUE, results = 'asis'----
models <- list(L110D = L110D, L011D = L011D, 
               L101D = L101D, L111D = L111D, L111I = L111I)
coefs <- unique(names(Reduce("c", lapply(models, coef))))
o_h1 <- grep("h1", coefs, value = TRUE)
o_h2 <- grep("h2", coefs, value = TRUE)
o_h3 <- grep("h3", coefs, value = TRUE)
n_h1 <- substr(o_h1, 4, 100)
n_h2 <- paste(substr(o_h2, 4, 100), " ", sep = "")
n_h3 <- paste(substr(o_h3, 4, 100), "  ", sep = "")
mps <- c(n_h1, n_h2, n_h3,"$\\sigma$", "$\\alpha$", "$\\rho_{12}$", "$\\rho_{13}$", "$\\rho_{23}$")
names(mps) <- c(o_h1, o_h2, o_h3, "sd.sd", "pos", "corr12", "corr13", "corr23")

gm <- tibble::tribble(
    ~raw,        ~clean,          ~fmt,
    "nobs",      "$N$",             0,
    "dpar", "degree of par.",       1,
    "nobs.zero", "$N_o$",           0,
    "nobs.pos", "$N_+$",            0,
    "R2.zero", "$R^2_o$",           3,
    "R2.pos", "$R^2_+$",            3,
    "logLik", "log likelihood",     1)

v <- modelsummary::msummary(models,
                       statistic = NULL, estimate = "{estimate} ({std.error})", fmt = 2,
                       title = "Multiple hurdle regressions",
                       coef_map = mps,
                       gof_map = gm, escape = FALSE)
## v <- kableExtra::pack_rows(v, index = c(selection = 5, consumption = 2, expense = 5,
##                             "Miscellanous parameters" = 5, "Goodness of fit measures" = 5))
v

## ----label = loglik, results = 'hide'--------------------------
library("lmtest")
lrtest(L111D, L111I)

## ----label = vuongnested, results = 'hide'---------------------
vuongtest(L111D, L111I, type = "nested")

## ----label = hvuongnested, include = FALSE---------------------
pvlr <- round(lrtest(L111D, L111I)[["Pr(>Chisq)"]][2], 3)
pvvuong <- vuongtest(L111D, L111I, type = "nested")$p.value

## ----label = vuongdh, eval = FALSE-----------------------------
#  vuongtest(L110D, L011D)

## ----label hvuongunnested, include = FALSE---------------------
vt <- vuongtest(L110D, L011D)
stat <- as.numeric(round(vt$statistic, 3))
pvaltrad <- round(vt$p.value, 3)

