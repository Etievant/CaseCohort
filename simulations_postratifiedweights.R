## -----------------------------------------------------------------------------
## Filename: simulation.R
## -----------------------------------------------------------------------------
## Description: This script compares the estimation of relative hazard, 
##              cumulative baseline and pure risks when using the whole cohort, 
##              a stratified case-cohort with design weights, a stratified
##              case-cohort with weights poststratified on the number of 
##              non-cases, an un-stratified case cohort with design weights, and 
##              an un-stratified case-cohort with poststratified weights. When 
##              using the case-cohorts, we also compare the robust variance 
##              estimate and the estimate with two components, one for the 
##              variation due to sampling from the superpopulation, and one due 
##              to sampling from the cohort, as proposed in Etievant and Gail 
##              (2022).
##
##              This script gives the simulation results displayed in Web 
##              Appendix E.4 in Etievant and Gail (2022)
## -----------------------------------------------------------------------------
## Required Package: dplyr, ggplot2, grid, gridExtra, gtable, parallel, 
## survival, xtable
## -----------------------------------------------------------------------------
## Required Functions: use help.functions.R
## -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Read functions and packages --------------------------------------------------

library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(parallel)
library(survival)
library(xtable)

source("help.functions.R")

# ------------------------------------------------------------------------------
# Different scenarios to be investigated over simulations ----------------------

param <- scenarios(n = c(5000, 10000), prob.y = c(0.02, 0.05, 0.1), 
                   noncases.per.case = c(2, 4), part = 4)

# ------------------------------------------------------------------------------
# Fixed parameters -------------------------------------------------------------

time    <- 10 
alpha1  <- 0.05 
alpha2  <- -0.35
beta1   <- - 0.2 
beta2   <- 0.25
beta3   <- - 0.3
Prob.x2 <- matrix(c(0.7, 0.05, 0.25, 0.45, 0.2, 0.35, 0.40, 0.3, 0.3), nrow = 3,
                  byrow = TRUE)
x1.cat  <- 1:3
w       <- 0:3
Tau1    <- 0
Tau2    <- 8
x1      <- c(-1, 1, -0.6) # given covariate profiles for the pure risk
x2      <- c(1, -1, 0.6)
x3      <- c(1, 1, 0.6)
names.beta  <- paste0("beta", c(1,2,3))
names.Pi    <- paste0("Pi.x", c(1,2,3), ".Tau1Tau2")

# ------------------------------------------------------------------------------
# Quantities needed to compute the time-fixed baseline hazard ------------------

densityX1X3 <- Vectorize(densityX1X3)
E.exp       <- E.RH(alpha1, alpha2, beta1, beta2, beta3, Prob.x2) # E(exp(beta1 * X1 + beta2 * X2 + beta3 * X3))

E.exp.W   <- E.RH.w(alpha1, alpha2, beta1, beta2, beta3, Prob.x2)
E.exp.w0  <- E.exp.W$E.exp.w0 # E(exp(beta1 * X1 + beta2 * X2 + beta3 * X3) | W = w)
E.exp.w1  <- E.exp.W$E.exp.w1
E.exp.w2  <- E.exp.W$E.exp.w2
E.exp.w3  <- E.exp.W$E.exp.w3

Prob.w <- c(integrate(f = densityX1, lower = 0, upper = 1)$value * 
              sum((Prob.x2[2,])[1]) + 
              integrate(f = densityX1, lower = 1, upper = 10)$value * 
              sum((Prob.x2[3,])[1]),
            integrate(f = densityX1, lower = -10, upper = -2)$value * 
              sum((Prob.x2[1,])[1:2]) + 
              integrate(f = densityX1, lower = -2, upper = 0)$value * 
              sum((Prob.x2[2,])[1:2]),
            integrate(f = densityX1, lower = 0, upper = 1)$value * 
              sum((Prob.x2[2,])[2:3]) + 
              integrate(f = densityX1, lower = 1, upper = 10)$value * 
              sum((Prob.x2[3,])[2:3]),
            integrate(f = densityX1, lower = -10, upper = -2)$value * 
              sum((Prob.x2[1,])[3]) + 
              integrate(f = densityX1, lower = -2, upper = 0)$value * 
              sum((Prob.x2[2,])[3])) # P(W = w)

Nreplic <- 5000 / 4 # we break down the replications in 4 

set.seed(12345)
seed <- round(abs(rnorm(4) * 10^5))

Onerun <- function (p) {
  
  Prob.y            <- param[p, 1] 
  n                 <- param[p, 2] 
  K                 <- param[p, 3] # number of non-cases we wish to sample for each case
  part              <- param[p, 4]
  
  set.seed(seed[part])
  
  # Time-fixed baseline hazard -------------------------------------------------
  bh                <- Prob.y / E.exp / time # (time-fixed) baseline hazard 
  names(bh)         <- NULL
  
  # Fixed numbers of individuals that will be sampled to obtain an stratified 
  # case-cohort ----------------------------------------------------------------
  E.strata.n        <- n * Prob.w # expected strata sizes
  strata.m          <- round(E.strata.n * K * bh * time * 
    c(E.exp.w0, E.exp.w1, E.exp.w2, E.exp.w3) / (1 - bh * time * 
                                                   c(E.exp.w0, E.exp.w1, 
                                                     E.exp.w2, E.exp.w3))) 
  
  # Fixed number of individuals that will be sampled to obtain an unstratified 
  # case-cohort ----------------------------------------------------------------
  m                 <- round(n * K * bh * time * E.exp / 
                               (1 - bh * time * E.exp))
  
  # True values of the cumulative baseline hazard in [Tau1, Tau2] and of the 
  # pure risks in [Tau1, Tau2] and for covariate profiles x1, x2 and x3 --------
  Lambda0.Tau1Tau2  <- bh * (Tau2 - Tau1)
  Pi.x1.Tau1Tau2    <- c(1 - exp(-exp(x1 %*% c(beta1, beta2, beta3)) * bh * 
                                   (Tau2 - Tau1)))
  Pi.x2.Tau1Tau2    <- c(1 - exp(-exp(x2 %*% c(beta1, beta2, beta3)) * bh * 
                                   (Tau2 - Tau1)))
  Pi.x3.Tau1Tau2    <- c(1 - exp(-exp(x3 %*% c(beta1, beta2, beta3)) * bh * 
                                   (Tau2 - Tau1)))
  
  res               <- NULL
  
  for (nrep in 1:Nreplic) {
    
    # --------------------------------------------------------------------------
    # Generation of the cohort -------------------------------------------------
    
    gener <- cohort.generation(n, time = 10, bh, alpha1 = 0.05, alpha2 = -0.35, 
                                beta1 = - 0.2, beta2 = 0.25, beta3 = - 0.3,
                                Prob.x2 = matrix(c(0.7, 0.05, 0.25, 0.45, 0.2, 
                                                 0.35, 0.40, 0.3, 0.3), 
                                               nrow = 3, byrow = TRUE), 
                                strata.m = strata.m)
    
    cohort    <- gener$cohort
    cohort    <- cbind(cohort, n = rep(n, n), m = rep(m, n))

    # --------------------------------------------------------------------------
    # Sampling of the stratified case cohort -----------------------------------
  
    strata.n  <- gener$strata.n
    strata.m  <- gener$strata.m 
    strata.n.cases <- gener$strata.n.cases
    
    indiv.sampled <- NULL
    for (ww in w) {
      indiv.sampled <- c(indiv.sampled, sample(cohort$id[cohort$W == ww], 
                                               size = strata.m[ww + 1], 
                                               replace = FALSE))
    }
    subcohort <- rep(0, n)
    subcohort[indiv.sampled] <- 1
    subcohort[cohort$status == 1] <- 1 # we add all the remaining cases
    
    casecohort <- cohort
    casecohort <- casecohort[-which((subcohort == 0)), ]
    casecohort <- as.data.frame(casecohort)
    casecohort$weights <- 1
    casecohort$weights[casecohort$status == 0] <- casecohort$strata.n[casecohort$status == 0] / casecohort$strata.m[casecohort$status == 0] # 1 / stratum-specific sampling fraction
    
    # --------------------------------------------------------------------------
    # Sampling of the unstratified case cohort ---------------------------------
    
    indiv.sampled <- sample(cohort$id, size = m, replace = FALSE)
    unstrat.subcohort <- rep(0, n)
    unstrat.subcohort[indiv.sampled] <- 1
    unstrat.subcohort[cohort$status == 1] <- 1 # we add all the remaining cases
    
    unstrat.casecohort <- cohort
    unstrat.casecohort <- unstrat.casecohort[-which((unstrat.subcohort == 0)), ]
    unstrat.casecohort <- as.data.frame(unstrat.casecohort)
    unstrat.casecohort$weights <- 1
    unstrat.casecohort$weights[unstrat.casecohort$status == 0] <- n / m # 1 / sampling fraction
    
    cohort <- cbind(cohort, subcohort, unstrat.subcohort)
    
    # --------------------------------------------------------------------------
    # Estimation using the whole cohort ----------------------------------------
    
    mod.cohort  <- coxph(Surv(times, status) ~ X1 + X2 + X3, data = cohort, 
                                       robust = TRUE)
    
    # parameters and influences estimation -------------------------------------
    estimation.cohort             <- influences(mod.cohort, Tau1 = Tau1, 
                                                Tau2 = Tau2, x = x1)
    
    beta.hat.cohort               <- estimation.cohort$beta.hat
    Lambda0.Tau1Tau2.hat.cohort   <- estimation.cohort$Lambda0.Tau1Tau2.hat
    Pi.x1.Tau1Tau2.hat.cohort     <- estimation.cohort$Pi.x.Tau1Tau2.hat
    
    infl.beta.cohort              <- estimation.cohort$infl.beta
    infl.Lambda0.Tau1Tau2.cohort  <- estimation.cohort$infl.Lambda0.Tau1Tau2
    infl.Pi.x1.Tau1Tau2.cohort    <- estimation.cohort$infl.Pi.x.Tau1Tau2
    
    estimation.cohort.x2 <- influences.PR(beta = beta.hat.cohort, 
                                          Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.cohort, 
                                          x = x2, infl.beta = infl.beta.cohort, 
                                          infl.Lambda0.Tau1Tau2 = infl.Lambda0.Tau1Tau2.cohort)
    estimation.cohort.x3 <- influences.PR(beta = beta.hat.cohort, 
                                          Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.cohort, 
                                          x = x3, infl.beta = infl.beta.cohort, 
                                          infl.Lambda0.Tau1Tau2 = infl.Lambda0.Tau1Tau2.cohort)
    
    Pi.x2.Tau1Tau2.hat.cohort     <- estimation.cohort.x2$Pi.x.Tau1Tau2.hat
    Pi.x3.Tau1Tau2.hat.cohort     <- estimation.cohort.x3$Pi.x.Tau1Tau2.hat
    
    infl.Pi.x2.Tau1Tau2.cohort    <- estimation.cohort.x2$infl.Pi.x.Tau1Tau2
    infl.Pi.x3.Tau1Tau2.cohort    <- estimation.cohort.x3$infl.Pi.x.Tau1Tau2
    
    # Robust variance estimation and corresponding CIs -------------------------
    coxph.var.beta.cohort               <- diag(mod.cohort$var) 
    robust.var.beta.cohort              <- robustvariance(infl.beta.cohort)
    robust.var.Lambda0.Tau1Tau2.cohort  <- robustvariance(infl.Lambda0.Tau1Tau2.cohort)
    robust.var.Pi.x1.Tau1Tau2.cohort    <- robustvariance(infl.Pi.x1.Tau1Tau2.cohort)
    robust.var.Pi.x2.Tau1Tau2.cohort    <- robustvariance(infl.Pi.x2.Tau1Tau2.cohort)
    robust.var.Pi.x3.Tau1Tau2.cohort    <- robustvariance(infl.Pi.x3.Tau1Tau2.cohort)
    
    CI.beta.cohort  <- conf.interval(beta.hat.cohort, robust.var.beta.cohort)
    CI.Lambda0.Tau1Tau2.cohort <- conf.interval(Lambda0.Tau1Tau2.hat.cohort, 
                                                robust.var.Lambda0.Tau1Tau2.cohort)
    CI.Pi.x1.Tau1Tau2.cohort <- conf.interval(Pi.x1.Tau1Tau2.hat.cohort,
                                            robust.var.Pi.x1.Tau1Tau2.cohort)
    CI.Pi.x2.Tau1Tau2.cohort <- conf.interval(Pi.x2.Tau1Tau2.hat.cohort,
                                             robust.var.Pi.x2.Tau1Tau2.cohort)
    CI.Pi.x3.Tau1Tau2.cohort <- conf.interval(Pi.x3.Tau1Tau2.hat.cohort,
                                             robust.var.Pi.x3.Tau1Tau2.cohort)
    
    res.cohort <- cbind("WholeCohort", t(beta.hat.cohort), 
                        Lambda0.Tau1Tau2.hat.cohort, Pi.x1.Tau1Tau2.hat.cohort, 
                        Pi.x2.Tau1Tau2.hat.cohort, Pi.x3.Tau1Tau2.hat.cohort, 
                        t(coxph.var.beta.cohort), t(robust.var.beta.cohort), 
                        robust.var.Lambda0.Tau1Tau2.cohort, 
                        robust.var.Pi.x1.Tau1Tau2.cohort,
                        robust.var.Pi.x2.Tau1Tau2.cohort,
                        robust.var.Pi.x3.Tau1Tau2.cohort,t(CI.beta.cohort), 
                        t(CI.Lambda0.Tau1Tau2.cohort), 
                        t(CI.Pi.x1.Tau1Tau2.cohort), t(CI.Pi.x2.Tau1Tau2.cohort),
                        t(CI.Pi.x3.Tau1Tau2.cohort))
    
    # --------------------------------------------------------------------------
    # Estimation using the stratified case cohort with design weights ----------
    
    mod                   <- coxph(Surv(times, status) ~ X1 + X2 + X3, 
                                   data = casecohort, weight = weights, 
                                   id = id, robust = TRUE)
    
    # parameters and influences estimation -------------------------------------
    estimation            <- influences(mod, Tau1 = Tau1, Tau2 = Tau2, x = x1)
    
    beta.hat              <- estimation$beta.hat
    Lambda0.Tau1Tau2.hat  <- estimation$Lambda0.Tau1Tau2.hat
    Pi.x1.Tau1Tau2.hat    <- estimation$Pi.x.Tau1Tau2.hat
    
    infl.beta             <- estimation$infl.beta
    infl.Lambda0.Tau1Tau2 <- estimation$infl.Lambda0.Tau1Tau2
    infl.Pi.x1.Tau1Tau2   <- estimation$infl.Pi.x.Tau1Tau2
    
    estimation.x2 <- influences.PR(beta = beta.hat, 
                                   Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat, 
                                   x = x2, infl.beta = infl.beta, 
                                   infl.Lambda0.Tau1Tau2 = infl.Lambda0.Tau1Tau2)
    estimation.x3 <- influences.PR(beta = beta.hat, 
                                   Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat, 
                                   x = x3, infl.beta = infl.beta, 
                                   infl.Lambda0.Tau1Tau2 = infl.Lambda0.Tau1Tau2)
    
    Pi.x2.Tau1Tau2.hat     <- estimation.x2$Pi.x.Tau1Tau2.hat
    Pi.x3.Tau1Tau2.hat     <- estimation.x3$Pi.x.Tau1Tau2.hat
    
    infl.Pi.x2.Tau1Tau2    <- estimation.x2$infl.Pi.x.Tau1Tau2
    infl.Pi.x3.Tau1Tau2    <- estimation.x3$infl.Pi.x.Tau1Tau2
    
    # Robust variance estimation and corresponding CIs -------------------------
    coxph.var.beta              <- diag(mod$var) 
    robust.var.beta             <- robustvariance(infl.beta) 
    robust.var.Lambda0.Tau1Tau2 <- robustvariance(infl.Lambda0.Tau1Tau2)
    robust.var.Pi.x1.Tau1Tau2   <- robustvariance(infl.Pi.x1.Tau1Tau2)
    robust.var.Pi.x2.Tau1Tau2   <- robustvariance(infl.Pi.x2.Tau1Tau2)
    robust.var.Pi.x3.Tau1Tau2   <- robustvariance(infl.Pi.x3.Tau1Tau2)

    CI.beta.robust  <- conf.interval(beta.hat, robust.var.beta)
    CI.Lambda0.Tau1Tau2.robust <- conf.interval(Lambda0.Tau1Tau2.hat, 
                                                robust.var.Lambda0.Tau1Tau2)
    CI.Pi.x1.Tau1Tau2.robust <- conf.interval(Pi.x1.Tau1Tau2.hat,
                                              robust.var.Pi.x1.Tau1Tau2)
    CI.Pi.x2.Tau1Tau2.robust <- conf.interval(Pi.x2.Tau1Tau2.hat,
                                             robust.var.Pi.x2.Tau1Tau2)
    CI.Pi.x3.Tau1Tau2.robust <- conf.interval(Pi.x3.Tau1Tau2.hat,
                                             robust.var.Pi.x3.Tau1Tau2)
    
    res.surveyrobust <- cbind("Case Cohort Robust", t(beta.hat), 
                              Lambda0.Tau1Tau2.hat, Pi.x1.Tau1Tau2.hat, 
                              Pi.x2.Tau1Tau2.hat, Pi.x3.Tau1Tau2.hat,
                              t(coxph.var.beta), t(robust.var.beta), 
                              robust.var.Lambda0.Tau1Tau2, 
                              robust.var.Pi.x1.Tau1Tau2, 
                              robust.var.Pi.x2.Tau1Tau2,
                              robust.var.Pi.x3.Tau1Tau2,
                              t(CI.beta.robust), 
                              t(CI.Lambda0.Tau1Tau2.robust), 
                              t(CI.Pi.x1.Tau1Tau2.robust),
                              t(CI.Pi.x2.Tau1Tau2.robust),
                              t(CI.Pi.x3.Tau1Tau2.robust))
    
    # Superpopulation and phase-two sampling variance estimation and 
    # corresponding CIs --------------------------------------------------------
    var.beta          <- variance(n = n, casecohort = casecohort, 
                                  infl = infl.beta, stratified = TRUE)
    var.Lambda0.Tau1Tau2 <- variance(n = n, casecohort = casecohort, 
                                     infl = infl.Lambda0.Tau1Tau2, 
                                     stratified = TRUE)
    var.Pi.x1.Tau1Tau2   <- variance(n = n, casecohort = casecohort, 
                                     infl = infl.Pi.x1.Tau1Tau2, 
                                     stratified = TRUE)
    var.Pi.x2.Tau1Tau2   <- variance(n = n, casecohort = casecohort, 
                                     infl = infl.Pi.x2.Tau1Tau2, 
                                     stratified = TRUE)
    var.Pi.x3.Tau1Tau2   <- variance(n = n, casecohort = casecohort, 
                                     infl = infl.Pi.x3.Tau1Tau2, 
                                     stratified = TRUE)
    
    CI.beta   <- conf.interval(beta.hat, var.beta)
 
    CI.Lambda0.Tau1Tau2   <- conf.interval(Lambda0.Tau1Tau2.hat, 
                                                 var.Lambda0.Tau1Tau2)
    CI.Pi.x1.Tau1Tau2  <- conf.interval(Pi.x1.Tau1Tau2.hat,
                                             var.Pi.x1.Tau1Tau2)
    CI.Pi.x2.Tau1Tau2  <- conf.interval(Pi.x2.Tau1Tau2.hat,
                                        var.Pi.x2.Tau1Tau2)
    CI.Pi.x3.Tau1Tau2  <- conf.interval(Pi.x3.Tau1Tau2.hat,
                                        var.Pi.x3.Tau1Tau2)
 
    res.survey <- cbind("Case Cohort", t(beta.hat), Lambda0.Tau1Tau2.hat, 
                        Pi.x1.Tau1Tau2.hat, Pi.x2.Tau1Tau2.hat, 
                        Pi.x3.Tau1Tau2.hat, t(coxph.var.beta), t(var.beta), 
                        var.Lambda0.Tau1Tau2, var.Pi.x1.Tau1Tau2,
                        var.Pi.x2.Tau1Tau2, var.Pi.x3.Tau1Tau2, t(CI.beta), 
                        t(CI.Lambda0.Tau1Tau2), t(CI.Pi.x1.Tau1Tau2),
                        t(CI.Pi.x2.Tau1Tau2), t(CI.Pi.x3.Tau1Tau2))

    # --------------------------------------------------------------------------
    # Estimation using the stratified case cohort with postratified weights ----
    
    # Poststratification of the weights ----------------------------------------
    casecohort$weights.poststrat <- 1
    casecohort$weights.poststrat[(casecohort$W == 0) & (casecohort$status == 0)] <- sum((cohort$W == 0) &
                                                                                           (cohort$status == 0)) / 
      sum((cohort$W == 0) & (cohort$status == 0) & (cohort$subcohort == 1))
    casecohort$weights.poststrat[(casecohort$W == 1) & (casecohort$status == 0)] <- sum((cohort$W == 1) &
                                                                                           (cohort$status == 0)) / 
      sum((cohort$W == 1) & (cohort$status == 0) & (cohort$subcohort == 1))
    casecohort$weights.poststrat[(casecohort$W == 2) & (casecohort$status == 0)] <- sum((cohort$W == 2) &
                                                                                           (cohort$status == 0)) / 
      sum((cohort$W == 2) & (cohort$status == 0) & (cohort$subcohort == 1))
    casecohort$weights.poststrat[(casecohort$W == 3) & (casecohort$status == 0)] <- sum((cohort$W == 3) &
                                                                                           (cohort$status == 0)) / 
      sum((cohort$W == 3) & (cohort$status == 0) & (cohort$subcohort == 1))
    
    mod.poststrat             <- coxph(Surv(times, status) ~ X1 + X2 + X3, 
                                       data = casecohort, 
                                       weight = weights.poststrat, id = id, 
                                       robust = TRUE)
    
    # Parameters and influences estimation -------------------------------------
    estimation.poststrat            <- influences(mod.poststrat, Tau1 = Tau1, 
                                                  Tau2 = Tau2, x = x1)
    
    beta.hat.poststrat              <- estimation.poststrat$beta.hat
    Lambda0.Tau1Tau2.hat.poststrat  <- estimation.poststrat$Lambda0.Tau1Tau2.hat
    Pi.x1.Tau1Tau2.hat.poststrat    <- estimation.poststrat$Pi.x.Tau1Tau2.hat
    
    infl.beta.poststrat             <- estimation.poststrat$infl.beta
    infl.Lambda0.Tau1Tau2.poststrat <- estimation.poststrat$infl.Lambda0.Tau1Tau2
    infl.Pi.x1.Tau1Tau2.poststrat   <- estimation.poststrat$infl.Pi.x.Tau1Tau2
    
    estimation.poststrat.x2 <- influences.PR(beta = beta.hat.poststrat,
                                             Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.poststrat, 
                                             x = x2, infl.beta = infl.beta.poststrat,
                                             infl.Lambda0.Tau1Tau2 = infl.Lambda0.Tau1Tau2.poststrat)
    
    estimation.poststrat.x3 <- influences.PR(beta = beta.hat.poststrat,
                                             Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.poststrat, 
                                             x = x3, infl.beta = infl.beta.poststrat,
                                             infl.Lambda0.Tau1Tau2 = infl.Lambda0.Tau1Tau2.poststrat)
    
    Pi.x2.Tau1Tau2.hat.poststrat    <- estimation.poststrat.x2$Pi.x.Tau1Tau2.hat
    Pi.x3.Tau1Tau2.hat.poststrat    <- estimation.poststrat.x3$Pi.x.Tau1Tau2.hat
    
    infl.Pi.x2.Tau1Tau2.poststrat   <- estimation.poststrat.x2$infl.Pi.x.Tau1Tau2
    infl.Pi.x3.Tau1Tau2.poststrat   <- estimation.poststrat.x3$infl.Pi.x.Tau1Tau2
    
    # Robust variance estimation and corresponding CIs -------------------------
    coxph.var.beta.poststrat            <- diag(mod.poststrat$var)
    robust.var.beta.poststrat           <- robustvariance(infl.beta.poststrat)
    robust.var.Lambda0.Tau1Tau2.poststrat <- robustvariance(infl.Lambda0.Tau1Tau2.poststrat)
    robust.var.Pi.x1.Tau1Tau2.poststrat  <- robustvariance(infl.Pi.x1.Tau1Tau2.poststrat)
    robust.var.Pi.x2.Tau1Tau2.poststrat  <- robustvariance(infl.Pi.x2.Tau1Tau2.poststrat)
    robust.var.Pi.x3.Tau1Tau2.poststrat  <- robustvariance(infl.Pi.x3.Tau1Tau2.poststrat)
    
    CI.beta.robust.poststrat <- conf.interval(beta.hat.poststrat,
                                              robust.var.beta.poststrat)
    CI.Lambda0.Tau1Tau2.robust.poststrat <- conf.interval(Lambda0.Tau1Tau2.hat.poststrat,
                                                          robust.var.Lambda0.Tau1Tau2.poststrat)
    CI.Pi.x1.Tau1Tau2.robust.poststrat <- conf.interval(Pi.x1.Tau1Tau2.hat.poststrat, 
                                                        robust.var.Pi.x1.Tau1Tau2.poststrat)
    CI.Pi.x2.Tau1Tau2.robust.poststrat <- conf.interval(Pi.x2.Tau1Tau2.hat.poststrat, 
                                                        robust.var.Pi.x2.Tau1Tau2.poststrat)
    CI.Pi.x3.Tau1Tau2.robust.poststrat <- conf.interval(Pi.x3.Tau1Tau2.hat.poststrat, 
                                                        robust.var.Pi.x3.Tau1Tau2.poststrat)
    
    res.surveyrobust.poststrat <- cbind("Case Cohort Poststratified Weights Robust", 
                                        t(beta.hat.poststrat), 
                                        Lambda0.Tau1Tau2.hat.poststrat, 
                                        Pi.x1.Tau1Tau2.hat.poststrat,
                                        Pi.x2.Tau1Tau2.hat.poststrat,
                                        Pi.x3.Tau1Tau2.hat.poststrat,
                                        t(coxph.var.beta.poststrat), 
                                        t(robust.var.beta.poststrat), 
                                        robust.var.Lambda0.Tau1Tau2.poststrat, 
                                        robust.var.Pi.x1.Tau1Tau2.poststrat, 
                                        robust.var.Pi.x2.Tau1Tau2.poststrat, 
                                        robust.var.Pi.x3.Tau1Tau2.poststrat, 
                                        t(CI.beta.robust.poststrat), 
                                        t(CI.Lambda0.Tau1Tau2.robust.poststrat), 
                                        t(CI.Pi.x1.Tau1Tau2.robust.poststrat),
                                        t(CI.Pi.x2.Tau1Tau2.robust.poststrat),
                                        t(CI.Pi.x3.Tau1Tau2.robust.poststrat))
    
    # Superpopulation and phase-two sampling variance estimation and 
    # corresponding CIs --------------------------------------------------------
    
    var.beta.poststrat <- variance(n = n, casecohort = casecohort, 
                                   infl = infl.beta.poststrat, 
                                   stratified = TRUE)
    
    var.Lambda0.Tau1Tau2.poststrat <- variance(n = n, casecohort = casecohort, 
                                               infl = infl.Lambda0.Tau1Tau2.poststrat, 
                                               stratified = TRUE)
    
    var.Pi.x1.Tau1Tau2.poststrat <- variance(n = n, casecohort = casecohort, 
                                               infl = infl.Pi.x1.Tau1Tau2.poststrat, 
                                               stratified = TRUE)
    
    var.Pi.x2.Tau1Tau2.poststrat <- variance(n = n, casecohort = casecohort, 
                                             infl = infl.Pi.x2.Tau1Tau2.poststrat, 
                                             stratified = TRUE)
    
    var.Pi.x3.Tau1Tau2.poststrat <- variance(n = n, casecohort = casecohort, 
                                             infl = infl.Pi.x3.Tau1Tau2.poststrat, 
                                             stratified = TRUE)

    CI.beta.poststrat <- conf.interval(beta.hat.poststrat,var.beta.poststrat)
    CI.Lambda0.Tau1Tau2.poststrat <- conf.interval(Lambda0.Tau1Tau2.hat.poststrat,
                                                   var.Lambda0.Tau1Tau2.poststrat)
    CI.Pi.x1.Tau1Tau2.poststrat <- conf.interval(Pi.x1.Tau1Tau2.hat.poststrat,
                                                 var.Pi.x1.Tau1Tau2.poststrat)
    CI.Pi.x2.Tau1Tau2.poststrat <- conf.interval(Pi.x2.Tau1Tau2.hat.poststrat,
                                                 var.Pi.x2.Tau1Tau2.poststrat)
    CI.Pi.x3.Tau1Tau2.poststrat <- conf.interval(Pi.x3.Tau1Tau2.hat.poststrat,
                                                 var.Pi.x3.Tau1Tau2.poststrat)
    
    res.survey.poststrat  <- cbind("Case Cohort Poststratified Weights", t(beta.hat.poststrat), 
                                   Lambda0.Tau1Tau2.hat.poststrat, 
                                   Pi.x1.Tau1Tau2.hat.poststrat,
                                   Pi.x2.Tau1Tau2.hat.poststrat,
                                   Pi.x3.Tau1Tau2.hat.poststrat, 
                                   t(coxph.var.beta.poststrat), t(var.beta.poststrat), 
                                   var.Lambda0.Tau1Tau2.poststrat, 
                                   var.Pi.x1.Tau1Tau2.poststrat,
                                   var.Pi.x2.Tau1Tau2.poststrat,
                                   var.Pi.x3.Tau1Tau2.poststrat, t(CI.beta.poststrat), 
                                   t(CI.Lambda0.Tau1Tau2.poststrat), 
                                   t(CI.Pi.x1.Tau1Tau2.poststrat), 
                                   t(CI.Pi.x2.Tau1Tau2.poststrat),
                                   t(CI.Pi.x3.Tau1Tau2.poststrat))

    # --------------------------------------------------------------------------
    # Estimation using the un-stratified case cohort with design weights -------
    
    mod.unstrat                 <- coxph(Surv(times, status) ~ X1 + X2 + X3, 
                                   data = unstrat.casecohort, weight = weights, 
                                   id = id, robust = TRUE)
    
    # Parameters and influences estimation -------------------------------------
    estimation.unstrat            <- influences(mod.unstrat, Tau1 = Tau1, 
                                                Tau2 = Tau2, x = x1)
    
    beta.hat.unstrat              <- estimation.unstrat$beta.hat
    Lambda0.Tau1Tau2.hat.unstrat  <- estimation.unstrat$Lambda0.Tau1Tau2.hat
    Pi.x1.Tau1Tau2.hat.unstrat    <- estimation.unstrat$Pi.x.Tau1Tau2.hat
    
    infl.beta.unstrat             <- estimation.unstrat$infl.beta
    infl.Lambda0.Tau1Tau2.unstrat <- estimation.unstrat$infl.Lambda0.Tau1Tau2
    infl.Pi.x1.Tau1Tau2.unstrat   <- estimation.unstrat$infl.Pi.x.Tau1Tau2
    
    estimation.x2.unstrat <- influences.PR(beta = beta.hat.unstrat, 
                                   Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.unstrat, 
                                   x = x2, infl.beta = infl.beta.unstrat, 
                                   infl.Lambda0.Tau1Tau2 = infl.Lambda0.Tau1Tau2.unstrat)
    estimation.x3.unstrat <- influences.PR(beta = beta.hat.unstrat, 
                                   Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.unstrat, 
                                   x = x3, infl.beta = infl.beta.unstrat, 
                                   infl.Lambda0.Tau1Tau2 = infl.Lambda0.Tau1Tau2.unstrat)
    
    Pi.x2.Tau1Tau2.hat.unstrat     <- estimation.x2.unstrat$Pi.x.Tau1Tau2.hat
    Pi.x3.Tau1Tau2.hat.unstrat     <- estimation.x3.unstrat$Pi.x.Tau1Tau2.hat
    
    infl.Pi.x2.Tau1Tau2.unstrat    <- estimation.x2.unstrat$infl.Pi.x.Tau1Tau2
    infl.Pi.x3.Tau1Tau2.unstrat    <- estimation.x3.unstrat$infl.Pi.x.Tau1Tau2
    
    # Robust variance estimation and corresponding CIs -------------------------
    coxph.var.beta.unstrat            <- diag(mod.unstrat$var)
    robust.var.beta.unstrat           <- robustvariance(infl.beta.unstrat)
    robust.var.Lambda0.Tau1Tau2.unstrat <- robustvariance(infl.Lambda0.Tau1Tau2.unstrat)
    robust.var.Pi.x1.Tau1Tau2.unstrat  <- robustvariance(infl.Pi.x1.Tau1Tau2.unstrat)
    robust.var.Pi.x2.Tau1Tau2.unstrat  <- robustvariance(infl.Pi.x2.Tau1Tau2.unstrat)
    robust.var.Pi.x3.Tau1Tau2.unstrat  <- robustvariance(infl.Pi.x3.Tau1Tau2.unstrat)
    
    CI.beta.robust.unstrat <- conf.interval(beta.hat.unstrat,
                                          robust.var.beta.unstrat)
    CI.Lambda0.Tau1Tau2.robust.unstrat <- conf.interval(Lambda0.Tau1Tau2.hat.unstrat,
                                                      robust.var.Lambda0.Tau1Tau2.unstrat)
    CI.Pi.x1.Tau1Tau2.robust.unstrat <- conf.interval(Pi.x1.Tau1Tau2.hat.unstrat, 
                                                    robust.var.Pi.x1.Tau1Tau2.unstrat)
    CI.Pi.x2.Tau1Tau2.robust.unstrat <- conf.interval(Pi.x2.Tau1Tau2.hat.unstrat, 
                                                    robust.var.Pi.x2.Tau1Tau2.unstrat)
    CI.Pi.x3.Tau1Tau2.robust.unstrat <- conf.interval(Pi.x3.Tau1Tau2.hat.unstrat, 
                                                    robust.var.Pi.x3.Tau1Tau2.unstrat)
    
    res.unstratrobust <- cbind("Unstrat Case Cohort Robust", 
                                    t(beta.hat.unstrat), 
                                    Lambda0.Tau1Tau2.hat.unstrat, 
                                    Pi.x1.Tau1Tau2.hat.unstrat,
                                    Pi.x2.Tau1Tau2.hat.unstrat,
                                    Pi.x3.Tau1Tau2.hat.unstrat,
                                    t(coxph.var.beta.unstrat), 
                                    t(robust.var.beta.unstrat), 
                                    robust.var.Lambda0.Tau1Tau2.unstrat, 
                                    robust.var.Pi.x1.Tau1Tau2.unstrat, 
                                    robust.var.Pi.x2.Tau1Tau2.unstrat, 
                                    robust.var.Pi.x3.Tau1Tau2.unstrat, 
                                    t(CI.beta.robust.unstrat), 
                                    t(CI.Lambda0.Tau1Tau2.robust.unstrat), 
                                    t(CI.Pi.x1.Tau1Tau2.robust.unstrat),
                                    t(CI.Pi.x2.Tau1Tau2.robust.unstrat),
                                    t(CI.Pi.x3.Tau1Tau2.robust.unstrat))
    
    # Superpopulation and phase-two sampling variance estimation and 
    # corresponding CIs --------------------------------------------------------
    var.beta.unstrat <- variance(n = n, casecohort = unstrat.casecohort, 
                                 infl = infl.beta.unstrat, stratified = FALSE)
    var.Lambda0.Tau1Tau2.unstrat <- variance(n = n, 
                                             casecohort = unstrat.casecohort, 
                                             infl = infl.Lambda0.Tau1Tau2.unstrat, 
                                             stratified = FALSE)
    var.Pi.x1.Tau1Tau2.unstrat <- variance(n = n, 
                                           casecohort = unstrat.casecohort, 
                                           infl = infl.Pi.x1.Tau1Tau2.unstrat,
                                           stratified = FALSE)
    var.Pi.x2.Tau1Tau2.unstrat <- variance(n = n, 
                                           casecohort = unstrat.casecohort, 
                                           infl = infl.Pi.x2.Tau1Tau2.unstrat, 
                                           stratified = FALSE)
    var.Pi.x3.Tau1Tau2.unstrat <- variance(n = n, 
                                           casecohort = unstrat.casecohort, 
                                           infl = infl.Pi.x3.Tau1Tau2.unstrat, 
                                           stratified = FALSE)
     
    CI.beta.unstrat             <- conf.interval(beta.hat.unstrat, 
                                                 var.beta.unstrat)
    CI.Lambda0.Tau1Tau2.unstrat <- conf.interval(Lambda0.Tau1Tau2.hat.unstrat, 
                                           var.Lambda0.Tau1Tau2.unstrat)
    CI.Pi.x1.Tau1Tau2.unstrat   <- conf.interval(Pi.x1.Tau1Tau2.hat.unstrat,
                                        var.Pi.x1.Tau1Tau2.unstrat)
    CI.Pi.x2.Tau1Tau2.unstrat   <- conf.interval(Pi.x2.Tau1Tau2.hat.unstrat,
                                                var.Pi.x2.Tau1Tau2.unstrat)
    CI.Pi.x3.Tau1Tau2.unstrat   <- conf.interval(Pi.x3.Tau1Tau2.hat.unstrat,
                                                var.Pi.x3.Tau1Tau2.unstrat)
    
    res.unstrat <- cbind("Unstrat Case Cohort", t(beta.hat.unstrat),
                         Lambda0.Tau1Tau2.hat.unstrat, 
                        Pi.x1.Tau1Tau2.hat.unstrat, Pi.x2.Tau1Tau2.hat.unstrat, 
                        Pi.x3.Tau1Tau2.hat.unstrat, t(coxph.var.beta.unstrat), 
                        t(var.beta.unstrat), var.Lambda0.Tau1Tau2.unstrat, 
                        var.Pi.x1.Tau1Tau2.unstrat, var.Pi.x2.Tau1Tau2.unstrat, 
                        var.Pi.x3.Tau1Tau2.unstrat, t(CI.beta.unstrat), 
                        t(CI.Lambda0.Tau1Tau2.unstrat), 
                        t(CI.Pi.x1.Tau1Tau2.unstrat),
                        t(CI.Pi.x2.Tau1Tau2.unstrat), 
                        t(CI.Pi.x3.Tau1Tau2.unstrat))

    # --------------------------------------------------------------------------
    # Estimation using the un-stratified case cohort with poststratified weights
    
    # Poststratification of the weights ----------------------------------------
    unstrat.casecohort$weights.poststrat <- 1
    unstrat.casecohort$weights.poststrat[unstrat.casecohort$status == 0] <- sum(cohort$status == 0) / 
      sum((cohort$status == 0) & (cohort$unstrat.subcohort == 1))
    
    mod.poststrat.unstrat           <- coxph(Surv(times, status) ~ X1 + X2 + X3, 
                                             data = unstrat.casecohort, 
                                             weight = weights.poststrat, id = id, 
                                             robust = TRUE)
    
    # parameters and influences estimation -------------------------------------
    estimation.poststrat.unstrat    <- influences(mod.poststrat.unstrat, 
                                                  Tau1 = Tau1, Tau2 = Tau2, 
                                                  x = x1)
    
    beta.hat.poststrat.unstrat              <- estimation.poststrat.unstrat$beta.hat
    Lambda0.Tau1Tau2.hat.poststrat.unstrat  <- estimation.poststrat.unstrat$Lambda0.Tau1Tau2.hat
    Pi.x1.Tau1Tau2.hat.poststrat.unstrat    <- estimation.poststrat.unstrat$Pi.x.Tau1Tau2.hat
    
    infl.beta.poststrat.unstrat             <- estimation.poststrat.unstrat$infl.beta
    infl.Lambda0.Tau1Tau2.poststrat.unstrat <- estimation.poststrat.unstrat$infl.Lambda0.Tau1Tau2
    infl.Pi.x1.Tau1Tau2.poststrat.unstrat    <- estimation.poststrat.unstrat$infl.Pi.x.Tau1Tau2
    
    estimation.poststrat.unstrat.x2 <- influences.PR(beta = beta.hat.poststrat.unstrat,
                                                     Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.poststrat.unstrat, 
                                                     x = x2, infl.beta = infl.beta.poststrat.unstrat,
                                                     infl.Lambda0.Tau1Tau2 = infl.Lambda0.Tau1Tau2.poststrat.unstrat)
    
    estimation.poststrat.unstrat.x3 <- influences.PR(beta = beta.hat.poststrat.unstrat,
                                                     Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.poststrat.unstrat, 
                                                     x = x3, infl.beta = infl.beta.poststrat.unstrat,
                                                     infl.Lambda0.Tau1Tau2 = infl.Lambda0.Tau1Tau2.poststrat.unstrat)
    
    Pi.x2.Tau1Tau2.hat.poststrat.unstrat    <- estimation.poststrat.unstrat.x2$Pi.x.Tau1Tau2.hat
    Pi.x3.Tau1Tau2.hat.poststrat.unstrat    <- estimation.poststrat.unstrat.x3$Pi.x.Tau1Tau2.hat
    
    infl.Pi.x2.Tau1Tau2.poststrat.unstrat   <- estimation.poststrat.unstrat.x2$infl.Pi.x.Tau1Tau2
    infl.Pi.x3.Tau1Tau2.poststrat.unstrat   <- estimation.poststrat.unstrat.x3$infl.Pi.x.Tau1Tau2
    
    # Robust variance estimation and corresponding CIs -------------------------
    coxph.var.beta.poststrat.unstrat            <- diag(mod.poststrat.unstrat$var)
    robust.var.beta.poststrat.unstrat           <- robustvariance(infl.beta.poststrat.unstrat)
    robust.var.Lambda0.Tau1Tau2.poststrat.unstrat <- robustvariance(infl.Lambda0.Tau1Tau2.poststrat.unstrat)
    robust.var.Pi.x1.Tau1Tau2.poststrat.unstrat  <- robustvariance(infl.Pi.x1.Tau1Tau2.poststrat.unstrat)
    robust.var.Pi.x2.Tau1Tau2.poststrat.unstrat  <- robustvariance(infl.Pi.x2.Tau1Tau2.poststrat.unstrat)
    robust.var.Pi.x3.Tau1Tau2.poststrat.unstrat  <- robustvariance(infl.Pi.x3.Tau1Tau2.poststrat.unstrat)
    
    CI.beta.robust.poststrat.unstrat <- conf.interval(beta.hat.poststrat.unstrat,
                                                      robust.var.beta.poststrat.unstrat)
    CI.Lambda0.Tau1Tau2.robust.poststrat.unstrat <- conf.interval(Lambda0.Tau1Tau2.hat.poststrat.unstrat,
                                                                  robust.var.Lambda0.Tau1Tau2.poststrat.unstrat)
    CI.Pi.x1.Tau1Tau2.robust.poststrat.unstrat <- conf.interval(Pi.x1.Tau1Tau2.hat.poststrat.unstrat, 
                                                                robust.var.Pi.x1.Tau1Tau2.poststrat.unstrat)
    CI.Pi.x2.Tau1Tau2.robust.poststrat.unstrat <- conf.interval(Pi.x2.Tau1Tau2.hat.poststrat.unstrat, 
                                                                robust.var.Pi.x2.Tau1Tau2.poststrat.unstrat)
    CI.Pi.x3.Tau1Tau2.robust.poststrat.unstrat <- conf.interval(Pi.x3.Tau1Tau2.hat.poststrat.unstrat, 
                                                                robust.var.Pi.x3.Tau1Tau2.poststrat.unstrat)
    
    res.unstratrobust.poststrat <- cbind("Unstrat Case Cohort Poststratified Weights Robust", 
                                         t(beta.hat.poststrat.unstrat), 
                                         Lambda0.Tau1Tau2.hat.poststrat.unstrat, 
                                         Pi.x1.Tau1Tau2.hat.poststrat.unstrat,
                                         Pi.x2.Tau1Tau2.hat.poststrat.unstrat,
                                         Pi.x3.Tau1Tau2.hat.poststrat.unstrat,
                                         t(coxph.var.beta.poststrat.unstrat), 
                                         t(robust.var.beta.poststrat.unstrat), 
                                         robust.var.Lambda0.Tau1Tau2.poststrat.unstrat, 
                                         robust.var.Pi.x1.Tau1Tau2.poststrat.unstrat, 
                                         robust.var.Pi.x2.Tau1Tau2.poststrat.unstrat, 
                                         robust.var.Pi.x3.Tau1Tau2.poststrat.unstrat, 
                                         t(CI.beta.robust.poststrat.unstrat), 
                                         t(CI.Lambda0.Tau1Tau2.robust.poststrat.unstrat), 
                                         t(CI.Pi.x1.Tau1Tau2.robust.poststrat.unstrat),
                                         t(CI.Pi.x2.Tau1Tau2.robust.poststrat.unstrat),
                                         t(CI.Pi.x3.Tau1Tau2.robust.poststrat.unstrat))
    
    # Superpopulation and phase-two sampling variance estimation and 
    # corresponding CIs --------------------------------------------------------
    var.beta.poststrat.unstrat <- variance(n = n, 
                                           casecohort = unstrat.casecohort, 
                                           infl = infl.beta.poststrat.unstrat, 
                                           stratified = FALSE)
    var.Lambda0.Tau1Tau2.poststrat.unstrat <- variance(n = n, 
                                                       casecohort = unstrat.casecohort, 
                                                       infl = infl.Lambda0.Tau1Tau2.poststrat.unstrat, 
                                                       stratified = FALSE)
    var.Pi.x1.Tau1Tau2.poststrat.unstrat <-  variance(n = n, 
                                                      casecohort = unstrat.casecohort, 
                                                      infl = infl.Pi.x1.Tau1Tau2.poststrat.unstrat, 
                                                      stratified = FALSE)
    var.Pi.x2.Tau1Tau2.poststrat.unstrat <-  variance(n = n, 
                                                      casecohort = unstrat.casecohort, 
                                                      infl = infl.Pi.x2.Tau1Tau2.poststrat.unstrat, 
                                                      stratified = FALSE)
    var.Pi.x3.Tau1Tau2.poststrat.unstrat <-  variance(n = n, 
                                                      casecohort = unstrat.casecohort, 
                                                      infl = infl.Pi.x3.Tau1Tau2.poststrat.unstrat, 
                                                      stratified = FALSE)

    CI.beta.poststrat.unstrat <- conf.interval(beta.hat.poststrat.unstrat,
                                               var.beta.poststrat.unstrat)
    CI.Lambda0.Tau1Tau2.poststrat.unstrat <- conf.interval(Lambda0.Tau1Tau2.hat.poststrat.unstrat,
                                                           var.Lambda0.Tau1Tau2.poststrat.unstrat)
    CI.Pi.x1.Tau1Tau2.poststrat.unstrat <- conf.interval(Pi.x1.Tau1Tau2.hat.poststrat.unstrat,
                                                         var.Pi.x1.Tau1Tau2.poststrat.unstrat)
    CI.Pi.x2.Tau1Tau2.poststrat.unstrat <- conf.interval(Pi.x2.Tau1Tau2.hat.poststrat.unstrat,
                                                         var.Pi.x2.Tau1Tau2.poststrat.unstrat)
    CI.Pi.x3.Tau1Tau2.poststrat.unstrat <- conf.interval(Pi.x3.Tau1Tau2.hat.poststrat.unstrat,
                                                         var.Pi.x3.Tau1Tau2.poststrat.unstrat)
    
    res.unstrat.poststrat  <- cbind("Unstrat Case Cohort Poststratified Weights", 
                                    t(beta.hat.poststrat.unstrat), 
                                    Lambda0.Tau1Tau2.hat.poststrat.unstrat, 
                                    Pi.x1.Tau1Tau2.hat.poststrat.unstrat,
                                    Pi.x2.Tau1Tau2.hat.poststrat.unstrat,
                                    Pi.x3.Tau1Tau2.hat.poststrat.unstrat, 
                                    t(coxph.var.beta.poststrat.unstrat), 
                                    t(var.beta.poststrat.unstrat), 
                                    var.Lambda0.Tau1Tau2.poststrat.unstrat, 
                                    var.Pi.x1.Tau1Tau2.poststrat.unstrat,
                                    var.Pi.x2.Tau1Tau2.poststrat.unstrat,
                                    var.Pi.x3.Tau1Tau2.poststrat.unstrat, 
                                    t(CI.beta.poststrat.unstrat), 
                                    t(CI.Lambda0.Tau1Tau2.poststrat.unstrat), 
                                    t(CI.Pi.x1.Tau1Tau2.poststrat.unstrat), 
                                    t(CI.Pi.x2.Tau1Tau2.poststrat.unstrat),
                                    t(CI.Pi.x3.Tau1Tau2.poststrat.unstrat))
    
    recap <- rbind(res.survey, res.surveyrobust, res.survey.poststrat, 
                   res.surveyrobust.poststrat, res.unstrat, res.unstratrobust, 
                   res.unstrat.poststrat, res.unstratrobust.poststrat, 
                   res.cohort)
    colnames(recap) <- c("Method", paste0(names.beta, ".hat"), 
                         "Lambda0.Tau1Tau2.hat", paste0(names.Pi, ".hat"),
                         paste0("coxph.var.", names.beta, ".hat"), 
                         paste0("var.", names.beta, ".hat"), 
                         "var.Lambda0.Tau1Tau2.hat", 
                         paste0("var.", names.Pi, ".hat"),
                         paste0("CI.left.", names.beta, ".hat"), 
                         paste0("CI.right.", names.beta, ".hat"), 
                         "CI.left.Lambda0.Tau1Tau2.hat", 
                         "CI.right.Lambda0.Tau1Tau2.hat", 
                         paste0("CI.left.", names.Pi, ".hat"), 
                         paste0("CI.right.", names.Pi, ".hat"))
    
   res <- rbind(res, cbind(recap, n, n.stratum.W0 = strata.n[1], 
                           n.stratum.W1 = strata.n[2], 
                           n.stratum.W2 = strata.n[3], 
                           n.stratum.W3 = strata.n[4], 
                           m.stratum.W0 = strata.m[1], 
                           m.stratum.W1 = strata.m[2], 
                           m.stratum.W2 = strata.m[3], 
                           m.stratum.W3 = strata.m[4], 
                           K, 
                           n.cases.stratum.W0 = strata.n.cases[1], 
                           n.cases.stratum.W1 = strata.n.cases[2], 
                           n.cases.stratum.W2 = strata.n.cases[3], 
                           n.cases.stratum.W3 =strata.n.cases[4],
                           m = m,
                           Prob.y, beta1, beta2, beta3, Lambda0.Tau1Tau2, 
                           Pi.x1.Tau1Tau2, Pi.x2.Tau1Tau2, Pi.x3.Tau1Tau2 ))
    row.names(res) <- NULL
     }
  myfile  <- paste0("SimulationResults_PoststratifiedWeights-n", n, "-Prob.y", 
                    Prob.y, "-K", K, "Part", part, ".Rdata")
  save(res, file = myfile)
}

# ------------------------------------------------------------------------------
# Running the simulations in parallel for all the scenarios --------------------

resultat <- mclapply(1:nrow(param), Onerun, mc.cores = 32)

# ------------------------------------------------------------------------------
# Combining the simulation results ---------------------------------------------

RES <- NULL
for (p in 1: nrow(param)) {
  
  Prob.y            <- param[p, 1] 
  n                 <- param[p, 2] 
  K                 <- param[p, 3] 
  part              <- param[p, 4]
  
  load(paste0("SimulationResults_PoststratifiedWeights-n", n, "-Prob.y", Prob.y,
              "-K", K, "Part", part, ".Rdata"))
  RES <- rbind(RES, res)
  
}
RECAP           <- as.data.frame(RES)
ColNames        <- colnames(RECAP[,c(2:55)])
RECAP[ColNames] <- sapply(RECAP[ColNames], as.numeric)
RECAP$Method    <- as.factor(RECAP$Method)
RECAP$n         <- as.factor(RECAP$n)
RECAP$K         <- as.factor(RECAP$K)
RECAP$Prob.y    <- as.factor(RECAP$Prob.y)
myfile  <- paste0("SimulationResults_PoststratifiedWeights.Rdata")
save(RECAP, file = myfile)

# ------------------------------------------------------------------------------
# Details of the results  ------------------------------------------------------

param <- scenarios(n = c(5000, 10000), prob.y = c(0.02, 0.05, 0.1), 
                   noncases.per.case = c(2, 4))

list.methods <- c("WholeCohort", "Case Cohort Robust", 
                  "Case Cohort", 
                  "Case Cohort Poststratified Weights Robust",  
                  "Case Cohort Poststratified Weights",
                  "Unstrat Case Cohort Robust", 
                  "Unstrat Case Cohort", 
                  "Unstrat Case Cohort Poststratified Weights Robust", 
                  "Unstrat Case Cohort Poststratified Weights")

list.methods.col <- c("Cohort", "SCC.Robust", "SCC", "SCC.Robust.Poststr",
                      "SCC.Poststr", "USCC.Robust", "USCC", 
                      "USCC.Robust.Poststr", "USCC.Poststr")

source("detailed.results.R")

detailed.results(RECAP = RECAP, param = param, Nreplic = 5000, 
                 list.methods = list.methods, 
                 list.methods.col = list.methods.col, 
                 simul.name = "_PoststratifiedWeights")
