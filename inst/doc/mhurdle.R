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

## ----label = tabmodels, include = FALSE, eval = TRUE-----------
models <- list(L110D = L110D, L011D = L011D, 
               L101D = L101D, L111D = L111D, L111I = L111I)
coefs <- unique(Reduce("c", lapply(models, function(x) names(coef(x)))))
coefs1 <- grep("h1", coefs)
coefs2 <- grep("h2", coefs)
coefs3 <- grep("h3", coefs)
coefso <- (1:length(coefs))[-c (coefs1, coefs2, coefs3)]
custcoefs <- coefs
custcoefs[coefso] <- c("$\\sigma$", "$\\rho_{12}$", "$\\alpha$", "$\\rho_{23}$", "$\\rho_{13}$")
coefs.h <- grep("h.\\.", coefs)
#custcoefs[coefs.h] <- substring(custcoefs[coefs.h], 4, 1E4)
groups <- list("**Hurdle 1**" = 1:length(coefs1),
               "**Hurlde 2**" = (length(coefs1)+1):(length(coefs1) + length(coefs2)),
               "**Hurdle 3**"= (length(coefs1)+length(coefs2)+1):(length(coefs1) + length(coefs2) + length(coefs3)),
               "**Others**" = (length(coefs1) + length(coefs2) + length(coefs3) + 1) : (length(coefs1) + length(coefs2) + length(coefs3) + length(coefso)))
#coeflabels <- gsub("h[1-3]\\.(\\w*)", "\\1", coefs)

## ----estimations, echo = FALSE, results = 'asis', eval = TRUE----
texreg::texreg(models, reorder.coef = c(coefs1, coefs2, coefs3, coefso), 
       custom.coef.names = custcoefs,#[c(coefs1, coefs2, coefs3, coefso)],
       caption = "Estimation of cencored models for the fees and admissions good",
       label = "tab:estimations",
       groups = groups
       ) 

## ----label = loglik, results = 'hide'--------------------------
library("lmtest")
lrtest(L111D, L111I)

## ----label = vuongnested, results = 'hide'---------------------
vuongtest(L111D, L111I, type = "nested")

## ----label = hvuongnested, include = FALSE---------------------
pvlr <- round(lrtest(L111D, L111I)[["Pr(>Chisq)"]][2], 3)
pvvuong <- vuongtest(L111D, L111I, type = "nested")$p.value

## ----label = vuongdh, results = 'hide'-------------------------
ndvuongtest(L110D, L011D)

## ----label hvuongunnested, include = FALSE---------------------
unnested <- ndvuongtest(L110D, L011D)
stat <- as.numeric(round(unnested$trad$statistic, 3))
pvaltrad <- round(unnested$trad$p.value, 3)
pvalshi <- round(unnested$nd$p.value, 3)

