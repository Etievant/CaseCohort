## -----------------------------------------------------------------------------
## Filename: simulation.R
## -----------------------------------------------------------------------------
## Description: This script compares the estimation of relative hazard, 
##              cumulative baseline and pure risks when using the whole cohort, 
##              a stratified case-cohort with design weights, and an 
##              un-stratified case cohort with design weights, when covariate 
##              information may be missing for individuals in phase-two sample.
##              When using the case-cohorts, we compare estimation with the 
##              estimated design weights, the true design weights, and the 
##              estimated design weights as if they were the true known ones.
##              We also compare the robust variance estimate and the estimate 
##              with two or three variance components.
##
##              This script gives the simulation results displayed in Web
##              Appendix I in Etievant and Gail (2022)
##
## -----------------------------------------------------------------------------
## Required Package: dplyr, ggplot2, grid, gridExtra, gtable, nnet, parallel, 
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
library(nnet)
library(parallel)
library(survival)
library(xtable)

source("help.functions.R")

# ------------------------------------------------------------------------------
# Different scenarios to be investigated over simulations ----------------------

param <- scenarios.missingdata(n = c(5000, 10000), prob.y = c(0.02, 0.05, 0.1), 
                               noncases.per.case = c(2, 4), 
                               strata.prob.missingness = rbind(c(0.9, 0.8), 
                                                               c(0.98, 0.9)), 
                               part = 4)

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
x1.cat    <- 1:3
w         <- 0:3
w.phase3  <- 0:1
Tau1      <- 0
Tau2      <- 8
x1        <- c(-1, 1, -0.6)  # given covariate profiles for the pure risk
x2        <- c(1, -1, 0.6)
x3        <- c(1, 1, 0.6)
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

Nreplic <- 5000 / 4

set.seed(12345)
seed <- round(abs(rnorm(4) * 10^5))

Onerun <- function (p) {
  
  Prob.y                  <- param[p, 1] 
  n                       <- param[p, 2] 
  K                       <- param[p, 3] # number of non-cases we wish to sample for each case
  strata.prob.missingness <- param[p, 4:5]
  part                    <- param[p, 6]
  
  set.seed(seed[part])
  
  # Time-fixed baseline hazard -------------------------------------------------
  bh                <- Prob.y / E.exp / time # (time-fixed) baseline hazard 
  names(bh)         <- NULL
  
  # Fixed numbers of individuals that will be sampled to obtain an stratified 
  # case-cohort ----------------------------------------------------------------
  E.strata.n        <- n * Prob.w # expected strata sizes
  strata.m          <- round(E.strata.n * K * bh * time * 
                               c(E.exp.w0, E.exp.w1, E.exp.w2, E.exp.w3) / 
                               (1 - bh * time * c(E.exp.w0, E.exp.w1, 
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
    cohort$n  <- n
    cohort$m  <- m
    cohort$W3 <- cohort$status # the J3 strata for the missingness are based on case status

    # --------------------------------------------------------------------------
    # Sampling of the stratified case cohort -----------------------------------
  
    strata.n  <- gener$strata.n
    strata.m  <- gener$strata.m # fixed numbers of individuals to be sampled in the strata
    strata.n.cases <- gener$strata.n.cases
    
    indiv.sampled.phase2 <- NULL
    for (ww in w) {
      indiv.sampled.phase2 <- c(indiv.sampled.phase2, 
                                sample(cohort$id[cohort$W == ww], 
                                       size = strata.m[ww + 1], 
                                       replace = FALSE))
    }
    indiv.sampled.phase2 <- unique(c(cohort$id[cohort$status == 1], 
                                     indiv.sampled.phase2))
    
    indiv.sampled.phase3 <- NULL
    for (ww in w.phase3) {
      indiv <- cohort$id[cohort$W3 == ww]
      v     <- rbinom(n = sum(cohort$W3 == ww), size = 1, 
                      prob = strata.prob.missingness[ww+1])
      indiv.sampled.phase3 <- c(indiv.sampled.phase3, indiv[which(v == 1)])
    }
    indiv.sampled <- c(indiv.sampled.phase2, 
                       indiv.sampled.phase3)[duplicated(c(indiv.sampled.phase2, 
                                                          indiv.sampled.phase3))]
    
    phase2            <- rep(0, n)
    phase2[indiv.sampled.phase2] <- 1
    
    phase3 <- rep(0, n)
    phase3[indiv.sampled] <- 1
    
    cohort$phase2     <- phase2
    cohort$phase3     <- phase3
    
    cohort$weights.true <- 1 / strata.prob.missingness[2]
    cohort$weights.true[cohort$status == 0] <- 1 / strata.prob.missingness[1] * 
      cohort$strata.n[cohort$status == 0] / 
      cohort$strata.m[cohort$status == 0] # 1 over the product of the phase-two and phase-three sampling probabilities
    cohort$weights.p2.true <- 1
    cohort$weights.p2.true[cohort$status == 0] <- cohort$strata.n[cohort$status == 0] / 
      cohort$strata.m[cohort$status == 0]
    cohort$weights.p3.true  <- 1 / strata.prob.missingness[cohort$W3 + 1]
    
    casecohort.phase2   <- cohort[-which(phase2 == 0), ]
    B.phase2            <- cbind(1 * (casecohort.phase2$W3 == 0), 
                                 1 * (casecohort.phase2$W3 == 1))
    rownames(B.phase2)  <- cohort[cohort$phase2 == 1, "id"]
    total.B.phase2      <- colSums(B.phase2)
    #casecohort.phase2$B <- B.phase2
    
    casecohort          <- cohort[-which(phase3 == 0), ]
    B.phase3            <- cbind(1 * (casecohort$W3 == 0), 
                                 1 * (casecohort$W3 == 1))
    rownames(B.phase3)  <- cohort[cohort$phase3 == 1, "id"]
    #casecohort$B        <- B.phase3
    
    estimation.weights.p3 <- weights.phase3(B.phase3, total.B.phase2, 
                                            gamma0 = rep(0, 2), 
                                            niter.max = 10^(4),
                                            epsilon.stop = 10^(-10))
    
    casecohort$weights.p3.est   <- estimation.weights.p3$estimated.weights
    casecohort$weights.est      <- casecohort$weights.p3.est * 
      casecohort$weights.p2.true
    
    casecohort.phase2$weights.p3.est  <- as.numeric(exp(B.phase2 %*% 
                                              estimation.weights.p3$gamma.hat))
    casecohort.phase2$weights.est     <- casecohort.phase2$weights.p3.est * 
      casecohort.phase2$weights.p2.true

    # --------------------------------------------------------------------------
    # Sampling of the unstratified case cohort ---------------------------------
    
    indiv.sampled.phase2.unstrat <- sample(cohort$id, size = m, replace = FALSE)
    indiv.sampled.phase2.unstrat <- unique(c(cohort$id[cohort$status == 1], 
                                     indiv.sampled.phase2.unstrat))
    
    indiv.sampled.phase3.unstrat <- NULL
    for (ww in w.phase3) {
      indiv <- cohort$id[cohort$W3 == ww]
      v     <- rbinom(n = sum(cohort$W3 == ww), size = 1, 
                      prob = strata.prob.missingness[ww+1])
      indiv.sampled.phase3.unstrat <- c(indiv.sampled.phase3.unstrat,
                                        indiv[which(v == 1)])
    }
    indiv.sampled.unstrat <- c(indiv.sampled.phase2.unstrat, 
                       indiv.sampled.phase3.unstrat)[duplicated(c(indiv.sampled.phase2.unstrat, 
                                                          indiv.sampled.phase3.unstrat))]
    
    unstrat.phase2            <- rep(0, n)
    unstrat.phase2[indiv.sampled.phase2.unstrat] <- 1
    
    unstrat.phase3 <- rep(0, n)
    unstrat.phase3[indiv.sampled.unstrat] <- 1
    
    cohort$unstrat.phase2   <- unstrat.phase2
    cohort$unstrat.phase3    <- unstrat.phase3
    
    cohort$unstrat.weights.true <- 1 / strata.prob.missingness[2]
    cohort$unstrat.weights.true[cohort$status == 0] <- 1 / strata.prob.missingness[1] * 
      cohort$n[cohort$status == 0] / 
      cohort$m[cohort$status == 0] # 1 over the product of the phase-two and phase-three sampling probabilities
    cohort$unstrat.weights.p2.true <- 1
    cohort$unstrat.weights.p2.true[cohort$status == 0] <- cohort$n[cohort$status == 0] / 
      cohort$m[cohort$status == 0]
    cohort$unstrat.weights.p3.true <- 1 / strata.prob.missingness[cohort$W3 + 1]
    
    unstrat.casecohort.phase2   <- cohort[-which(unstrat.phase2 == 0), ]
    B.phase2.unstrat            <- cbind(1 * (unstrat.casecohort.phase2$W3 == 0), 
                                 1 * (unstrat.casecohort.phase2$W3 == 1))
    rownames(B.phase2.unstrat)  <- cohort[cohort$unstrat.phase2 == 1, "id"]
    total.B.phase2.unstrat      <- colSums(B.phase2.unstrat)
    
    unstrat.casecohort          <- cohort[-which(unstrat.phase3 == 0), ]
    B.phase3.unstrat            <- cbind(1 * (unstrat.casecohort$W3 == 0), 
                                 1 * (unstrat.casecohort$W3 == 1))
    rownames(B.phase3.unstrat)  <- cohort[cohort$unstrat.phase3 == 1, "id"]
    
    estimation.unstrat.weights.p3 <- weights.phase3(B.phase3.unstrat, 
                                                    total.B.phase2.unstrat, 
                                                    gamma0 = rep(0, 2), 
                                                    niter.max = 10^(4),
                                                    epsilon.stop = 10^(-10))
    
    unstrat.casecohort$unstrat.weights.p3.est <- estimation.unstrat.weights.p3$estimated.weights
    unstrat.casecohort$unstrat.weights.est <- unstrat.casecohort$unstrat.weights.p3.est * 
      unstrat.casecohort$unstrat.weights.p2.true
    
    unstrat.casecohort.phase2$unstrat.weights.p3.est  <- as.numeric(exp(B.phase2.unstrat %*% 
                                               estimation.unstrat.weights.p3$gamma.hat))
    unstrat.casecohort.phase2$unstrat.weights.est     <- unstrat.casecohort.phase2$unstrat.weights.p3.est * 
      unstrat.casecohort.phase2$unstrat.weights.p2.true
    
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
    # Estimation using the stratified case cohort with true design weights -----
    
    mod                   <- coxph(Surv(times, status) ~ X1 + X2 + X3, 
                                   data = casecohort, weight = weights.true, 
                                   id = id, robust = TRUE)
    
    # At risk indicator matrix and counting process matrix for the phase-two 
    # data at all event times, even that with missing covariate data -----------
    mod.cohort.detail   <- coxph.detail(mod.cohort, riskmat = TRUE)   
    riskmat.phase2      <- with(cohort, mod.cohort.detail$riskmat[phase2 == 1,])
    rownames(riskmat.phase2)  <- cohort[cohort$phase2 == 1, "id"]
    observed.times.phase2 <- apply(riskmat.phase2, 1,
                                   function(v) {which.max(cumsum(v))})
    dNt.phase2            <- matrix(0, nrow(riskmat.phase2), 
                                    ncol(riskmat.phase2))
    dNt.phase2[cbind(1:nrow(riskmat.phase2), observed.times.phase2)] <- 1
    dNt.phase2            <- sweep(dNt.phase2, 1, casecohort.phase2$status, "*")
    colnames(dNt.phase2)  <- colnames(riskmat.phase2)
    rownames(dNt.phase2)  <- rownames(riskmat.phase2)
    
    # Parameters and influences estimation -------------------------------------
    estimation <- influences.missingdata(mod = mod, 
                                         riskmat.phase2 = riskmat.phase2, 
                                         dNt.phase2 = dNt.phase2, Tau1 = Tau1, 
                                         Tau2 = Tau2, x = x1)
    
    beta.hat                <- estimation$beta.hat
    Lambda0.Tau1Tau2.hat    <- estimation$Lambda0.Tau1Tau2.hat
    Pi.x1.Tau1Tau2.hat      <- estimation$Pi.x.Tau1Tau2.hat
    
    infl.beta              <- estimation$infl.beta
    infl.Lambda0.Tau1Tau2  <- estimation$infl.Lambda0.Tau1Tau2
    infl.Pi.x1.Tau1Tau2    <- estimation$infl.Pi.x.Tau1Tau2
    
    infl2.beta              <- estimation$infl2.beta
    infl2.Lambda0.Tau1Tau2  <- estimation$infl2.Lambda0.Tau1Tau2
    infl2.Pi.x1.Tau1Tau2    <- estimation$infl2.Pi.x.Tau1Tau2
    
    infl3.beta              <- estimation$infl3.beta
    infl3.Lambda0.Tau1Tau2  <- estimation$infl3.Lambda0.Tau1Tau2
    infl3.Pi.x1.Tau1Tau2    <- estimation$infl3.Pi.x.Tau1Tau2
    
    estimation.x2 <- influences.PR.missingdata(beta = beta.hat, 
                                   Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat, 
                                   x = x2, infl2.beta = infl2.beta, 
                                   infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.Tau1Tau2,
                                   infl3.beta = infl3.beta, 
                                   infl3.Lambda0.Tau1Tau2 = infl3.Lambda0.Tau1Tau2)
    estimation.x3 <- influences.PR.missingdata(beta = beta.hat, 
                                   Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat, 
                                   x = x3, infl2.beta = infl2.beta, 
                                   infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.Tau1Tau2,
                                   infl3.beta = infl3.beta, 
                                   infl3.Lambda0.Tau1Tau2 = infl3.Lambda0.Tau1Tau2)
    
    Pi.x2.Tau1Tau2.hat     <- estimation.x2$Pi.x.Tau1Tau2.hat
    Pi.x3.Tau1Tau2.hat     <- estimation.x3$Pi.x.Tau1Tau2.hat
    
    infl.Pi.x2.Tau1Tau2    <- estimation.x2$infl.Pi.x.Tau1Tau2
    infl.Pi.x3.Tau1Tau2    <- estimation.x3$infl.Pi.x.Tau1Tau2
    
    infl2.Pi.x2.Tau1Tau2    <- estimation.x2$infl2.Pi.x.Tau1Tau2
    infl2.Pi.x3.Tau1Tau2    <- estimation.x3$infl2.Pi.x.Tau1Tau2
    
    infl3.Pi.x2.Tau1Tau2    <- estimation.x2$infl3.Pi.x.Tau1Tau2
    infl3.Pi.x3.Tau1Tau2    <- estimation.x3$infl3.Pi.x.Tau1Tau2
    
    # Robust variance estimation and corresponding CIs -------------------------
    coxph.var.beta              <- diag(mod$var) 
    robust.var.beta             <- robustvariance(infl.beta) 
    robust.var.Lambda0.Tau1Tau2 <- robustvariance(infl.Lambda0.Tau1Tau2)
    robust.var.Pi.x1.Tau1Tau2    <- robustvariance(infl.Pi.x1.Tau1Tau2)
    robust.var.Pi.x2.Tau1Tau2    <- robustvariance(infl.Pi.x2.Tau1Tau2)
    robust.var.Pi.x3.Tau1Tau2    <- robustvariance(infl.Pi.x3.Tau1Tau2)

    CI.beta.robust  <- conf.interval(beta.hat, robust.var.beta)
    CI.Lambda0.Tau1Tau2.robust <- conf.interval(Lambda0.Tau1Tau2.hat, 
                                                robust.var.Lambda0.Tau1Tau2)
    CI.Pi.x1.Tau1Tau2.robust <- conf.interval(Pi.x1.Tau1Tau2.hat,
                                              robust.var.Pi.x1.Tau1Tau2)
    CI.Pi.x2.Tau1Tau2.robust <- conf.interval(Pi.x2.Tau1Tau2.hat,
                                             robust.var.Pi.x2.Tau1Tau2)
    CI.Pi.x3.Tau1Tau2.robust <- conf.interval(Pi.x3.Tau1Tau2.hat,
                                             robust.var.Pi.x3.Tau1Tau2)
    
    res.casecohort.robust <- cbind("Case Cohort Robust", t(beta.hat), 
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
    var.beta <- variance.missingdata(n = n, casecohort = casecohort, 
                                     casecohort.phase2 = casecohort.phase2, 
                                     weights = casecohort$weights.true, 
                                     weights.phase2 = casecohort.phase2$weights.true, 
                                     weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                                     infl2 = infl2.beta, infl3 = infl3.beta, 
                                     stratified.p2 = TRUE)
    var.Lambda0.Tau1Tau2 <- variance.missingdata(n = n, casecohort = casecohort, 
                                     casecohort.phase2 = casecohort.phase2, 
                                     weights = casecohort$weights.true, 
                                     weights.phase2 = casecohort.phase2$weights.true, 
                                     weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                                     infl2 = infl2.Lambda0.Tau1Tau2, 
                                     infl3 = infl3.Lambda0.Tau1Tau2, 
                                     stratified.p2 = TRUE)
    var.Pi.x1.Tau1Tau2 <- variance.missingdata(n = n, casecohort = casecohort, 
                                                 casecohort.phase2 = casecohort.phase2, 
                                                 weights = casecohort$weights.true, 
                                                 weights.phase2 = casecohort.phase2$weights.true, 
                                                 weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                                                 infl2 = infl2.Pi.x1.Tau1Tau2, 
                                                 infl3 = infl3.Pi.x1.Tau1Tau2, 
                                                 stratified.p2 = TRUE)
    var.Pi.x2.Tau1Tau2 <- variance.missingdata(n = n, casecohort = casecohort, 
                                               casecohort.phase2 = casecohort.phase2, 
                                               weights = casecohort$weights.true, 
                                               weights.phase2 = casecohort.phase2$weights.true, 
                                               weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                                               infl2 = infl2.Pi.x2.Tau1Tau2, 
                                               infl3 = infl3.Pi.x2.Tau1Tau2, 
                                               stratified.p2 = TRUE)
    var.Pi.x3.Tau1Tau2 <- variance.missingdata(n = n, casecohort = casecohort, 
                                               casecohort.phase2 = casecohort.phase2, 
                                               weights = casecohort$weights.true, 
                                               weights.phase2 = casecohort.phase2$weights.true, 
                                               weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                                               infl2 = infl2.Pi.x3.Tau1Tau2, 
                                               infl3 = infl3.Pi.x3.Tau1Tau2, 
                                               stratified.p2 = TRUE)
    
    CI.beta   <- conf.interval(beta.hat, var.beta)
    CI.Lambda0.Tau1Tau2 <- conf.interval(Lambda0.Tau1Tau2.hat, 
                                           var.Lambda0.Tau1Tau2)
    CI.Pi.x1.Tau1Tau2   <- conf.interval(Pi.x1.Tau1Tau2.hat,
                                        var.Pi.x1.Tau1Tau2)
    CI.Pi.x2.Tau1Tau2   <- conf.interval(Pi.x2.Tau1Tau2.hat,
                                        var.Pi.x2.Tau1Tau2)
    CI.Pi.x3.Tau1Tau2   <- conf.interval(Pi.x3.Tau1Tau2.hat,
                                        var.Pi.x3.Tau1Tau2)
    
    res.casecohort <- cbind("Case Cohort", t(beta.hat), Lambda0.Tau1Tau2.hat, 
                        Pi.x1.Tau1Tau2.hat, Pi.x2.Tau1Tau2.hat, 
                        Pi.x3.Tau1Tau2.hat, t(coxph.var.beta), t(var.beta), 
                        var.Lambda0.Tau1Tau2, var.Pi.x1.Tau1Tau2,
                        var.Pi.x2.Tau1Tau2, var.Pi.x3.Tau1Tau2, t(CI.beta), 
                        t(CI.Lambda0.Tau1Tau2), t(CI.Pi.x1.Tau1Tau2),
                        t(CI.Pi.x2.Tau1Tau2), t(CI.Pi.x3.Tau1Tau2))

    # --------------------------------------------------------------------------
    # Estimation using the stratified case cohort with estimated weights as if 
    # they were the true weights -----------------------------------------------
    
    mod.est <- coxph(Surv(times, status) ~ X1 + X2 + X3, 
                           data = casecohort, weight = weights.est, id = id, 
                           robust = TRUE)
    
    # parameters and influences estimation -------------------------------------
    estimation.est.naive  <- influences.missingdata(mod.est, 
                                                    riskmat.phase2 = riskmat.phase2, 
                                                    dNt.phase2 = dNt.phase2, 
                                                    Tau1 = Tau1, 
                                                    Tau2 = Tau2, x = x1)
    
    beta.hat.est.naive              <- estimation.est.naive$beta.hat
    Lambda0.Tau1Tau2.hat.est.naive  <- estimation.est.naive$Lambda0.Tau1Tau2.hat
    Pi.x1.Tau1Tau2.hat.est.naive    <- estimation.est.naive$Pi.x.Tau1Tau2.hat
    
    infl.beta.est.naive             <- estimation.est.naive$infl.beta
    infl.Lambda0.Tau1Tau2.est.naive <- estimation.est.naive$infl.Lambda0.Tau1Tau2
    infl.Pi.x1.Tau1Tau2.est.naive   <- estimation.est.naive$infl.Pi.x.Tau1Tau2
    
    infl2.beta.est.naive            <- estimation.est.naive$infl2.beta
    infl2.Lambda0.Tau1Tau2.est.naive <- estimation.est.naive$infl2.Lambda0.Tau1Tau2
    infl2.Pi.x1.Tau1Tau2.est.naive  <- estimation.est.naive$infl2.Pi.x.Tau1Tau2
    
    infl3.beta.est.naive            <- estimation.est.naive$infl3.beta
    infl3.Lambda0.Tau1Tau2.est.naive <- estimation.est.naive$infl3.Lambda0.Tau1Tau2
    infl3.Pi.x1.Tau1Tau2.est.naive  <- estimation.est.naive$infl3.Pi.x.Tau1Tau2
    
    estimation.est.naive.x2 <- influences.PR.missingdata(beta = beta.hat.est.naive,
                                         Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.est.naive, 
                                         x = x2, infl2.beta = infl2.beta.est.naive,
                                         infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.Tau1Tau2.est.naive,
                                         infl3.beta = infl3.beta.est.naive,
                                         infl3.Lambda0.Tau1Tau2 = infl3.Lambda0.Tau1Tau2.est.naive)

    estimation.est.naive.x3 <- influences.PR.missingdata(beta = beta.hat.est.naive,
                                                         Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.est.naive, 
                                                         x = x3, infl2.beta = infl2.beta.est.naive,
                                                         infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.Tau1Tau2.est.naive,
                                                         infl3.beta = infl3.beta.est.naive,
                                                         infl3.Lambda0.Tau1Tau2 = infl3.Lambda0.Tau1Tau2.est.naive)
    
    Pi.x2.Tau1Tau2.hat.est.naive    <- estimation.est.naive.x2$Pi.x.Tau1Tau2.hat
    Pi.x3.Tau1Tau2.hat.est.naive    <- estimation.est.naive.x3$Pi.x.Tau1Tau2.hat
    
    infl.Pi.x2.Tau1Tau2.est.naive   <- estimation.est.naive.x2$infl.Pi.x.Tau1Tau2
    infl.Pi.x3.Tau1Tau2.est.naive   <- estimation.est.naive.x3$infl.Pi.x.Tau1Tau2
    
    infl2.Pi.x2.Tau1Tau2.est.naive  <- estimation.est.naive.x2$infl2.Pi.x.Tau1Tau2
    infl2.Pi.x3.Tau1Tau2.est.naive  <- estimation.est.naive.x3$infl2.Pi.x.Tau1Tau2
    
    infl3.Pi.x2.Tau1Tau2.est.naive  <- estimation.est.naive.x2$infl3.Pi.x.Tau1Tau2
    infl3.Pi.x3.Tau1Tau2.est.naive  <- estimation.est.naive.x3$infl3.Pi.x.Tau1Tau2
    
    # Robust variance estimation and corresponding CIs -------------------------
    coxph.var.beta.est.naive            <- diag(mod.est$var)
    robust.var.beta.est.naive           <- robustvariance(infl.beta.est.naive)
    robust.var.Lambda0.Tau1Tau2.est.naive <- robustvariance(infl.Lambda0.Tau1Tau2.est.naive)
    robust.var.Pi.x1.Tau1Tau2.est.naive  <- robustvariance(infl.Pi.x1.Tau1Tau2.est.naive)
    robust.var.Pi.x2.Tau1Tau2.est.naive  <- robustvariance(infl.Pi.x2.Tau1Tau2.est.naive)
    robust.var.Pi.x3.Tau1Tau2.est.naive  <- robustvariance(infl.Pi.x3.Tau1Tau2.est.naive)
    
    CI.beta.robust.est.naive <- conf.interval(beta.hat.est.naive,
                                          robust.var.beta.est.naive)
    CI.Lambda0.Tau1Tau2.robust.est.naive <- conf.interval(Lambda0.Tau1Tau2.hat.est.naive,
                                                      robust.var.Lambda0.Tau1Tau2.est.naive)
    CI.Pi.x1.Tau1Tau2.robust.est.naive <- conf.interval(Pi.x1.Tau1Tau2.hat.est.naive, 
                                                  robust.var.Pi.x1.Tau1Tau2.est.naive)
    CI.Pi.x2.Tau1Tau2.robust.est.naive <- conf.interval(Pi.x2.Tau1Tau2.hat.est.naive, 
                                                    robust.var.Pi.x2.Tau1Tau2.est.naive)
    CI.Pi.x3.Tau1Tau2.robust.est.naive <- conf.interval(Pi.x3.Tau1Tau2.hat.est.naive, 
                                                    robust.var.Pi.x3.Tau1Tau2.est.naive)
    
    res.casecohort.est.robust.naive <- cbind("Case Cohort Robust Estimation Naive", 
                                    t(beta.hat.est.naive), 
                                    Lambda0.Tau1Tau2.hat.est.naive, 
                                    Pi.x1.Tau1Tau2.hat.est.naive,
                                    Pi.x2.Tau1Tau2.hat.est.naive,
                                    Pi.x3.Tau1Tau2.hat.est.naive,
                                    t(coxph.var.beta.est.naive), 
                                    t(robust.var.beta.est.naive), 
                                    robust.var.Lambda0.Tau1Tau2.est.naive, 
                                    robust.var.Pi.x1.Tau1Tau2.est.naive, 
                                    robust.var.Pi.x2.Tau1Tau2.est.naive, 
                                    robust.var.Pi.x3.Tau1Tau2.est.naive, 
                                    t(CI.beta.robust.est.naive), 
                                    t(CI.Lambda0.Tau1Tau2.robust.est.naive), 
                                    t(CI.Pi.x1.Tau1Tau2.robust.est.naive),
                                    t(CI.Pi.x2.Tau1Tau2.robust.est.naive),
                                    t(CI.Pi.x3.Tau1Tau2.robust.est.naive))
    
    # Superpopulation and phase-two sampling variance estimation and 
    # corresponding CIs --------------------------------------------------------
    var.beta.est.naive <- variance.missingdata(n = n, casecohort = casecohort, 
                                               casecohort.phase2 = casecohort.phase2, 
                                               weights = casecohort$weights.est, 
                                               weights.phase2 = casecohort.phase2$weights.est, 
                                               weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                                               infl2 = infl2.beta.est.naive, 
                                               infl3 = infl3.beta.est.naive, 
                                               stratified.p2 = TRUE)
    var.Lambda0.Tau1Tau2.est.naive <- variance.missingdata(n = n, casecohort = casecohort, 
                                               casecohort.phase2 = casecohort.phase2, 
                                               weights = casecohort$weights.est, 
                                               weights.phase2 = casecohort.phase2$weights.est, 
                                               weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                                               infl2 = infl2.Lambda0.Tau1Tau2.est.naive, 
                                               infl3 = infl3.Lambda0.Tau1Tau2.est.naive, 
                                               stratified.p2 = TRUE)
    var.Pi.x1.Tau1Tau2.est.naive <- variance.missingdata(n = n, casecohort = casecohort, 
                                                           casecohort.phase2 = casecohort.phase2, 
                                                           weights = casecohort$weights.est, 
                                                           weights.phase2 = casecohort.phase2$weights.est, 
                                                           weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                                                           infl2 = infl2.Pi.x1.Tau1Tau2.est.naive, 
                                                           infl3 = infl3.Pi.x1.Tau1Tau2.est.naive, 
                                                           stratified.p2 = TRUE)
    var.Pi.x2.Tau1Tau2.est.naive <- variance.missingdata(n = n, casecohort = casecohort, 
                                                           casecohort.phase2 = casecohort.phase2, 
                                                           weights = casecohort$weights.est, 
                                                           weights.phase2 = casecohort.phase2$weights.est, 
                                                           weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                                                           infl2 = infl2.Pi.x2.Tau1Tau2.est.naive, 
                                                           infl3 = infl3.Pi.x2.Tau1Tau2.est.naive, 
                                                           stratified.p2 = TRUE)
    var.Pi.x3.Tau1Tau2.est.naive <- variance.missingdata(n = n, casecohort = casecohort, 
                                                           casecohort.phase2 = casecohort.phase2, 
                                                           weights = casecohort$weights.est, 
                                                           weights.phase2 = casecohort.phase2$weights.est, 
                                                           weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                                                           infl2 = infl2.Pi.x3.Tau1Tau2.est.naive, 
                                                           infl3 = infl3.Pi.x3.Tau1Tau2.est.naive, 
                                                           stratified.p2 = TRUE)
 
    CI.beta.est.naive <- conf.interval(beta.hat.est.naive,var.beta.est.naive)
    CI.Lambda0.Tau1Tau2.est.naive <- conf.interval(Lambda0.Tau1Tau2.hat.est.naive,
                                                var.Lambda0.Tau1Tau2.est.naive)
    CI.Pi.x1.Tau1Tau2.est.naive <- conf.interval(Pi.x1.Tau1Tau2.hat.est.naive,
                                            var.Pi.x1.Tau1Tau2.est.naive)
    CI.Pi.x2.Tau1Tau2.est.naive <- conf.interval(Pi.x2.Tau1Tau2.hat.est.naive,
                                            var.Pi.x2.Tau1Tau2.est.naive)
    CI.Pi.x3.Tau1Tau2.est.naive <- conf.interval(Pi.x3.Tau1Tau2.hat.est.naive,
                                            var.Pi.x3.Tau1Tau2.est.naive)

    res.casecohort.est.naive  <- cbind("Case Cohort Estimation Naive", t(beta.hat.est.naive), 
                               Lambda0.Tau1Tau2.hat.est.naive, 
                               Pi.x1.Tau1Tau2.hat.est.naive,
                               Pi.x2.Tau1Tau2.hat.est.naive,
                               Pi.x3.Tau1Tau2.hat.est.naive, 
                               t(coxph.var.beta.est.naive), t(var.beta.est.naive), 
                               var.Lambda0.Tau1Tau2.est.naive, 
                               var.Pi.x1.Tau1Tau2.est.naive,
                               var.Pi.x2.Tau1Tau2.est.naive,
                               var.Pi.x3.Tau1Tau2.est.naive, t(CI.beta.est.naive), 
                               t(CI.Lambda0.Tau1Tau2.est.naive), 
                               t(CI.Pi.x1.Tau1Tau2.est.naive), 
                               t(CI.Pi.x2.Tau1Tau2.est.naive),
                               t(CI.Pi.x3.Tau1Tau2.est.naive))

    # --------------------------------------------------------------------------
    # Estimation using the stratified case cohort with estimated weights, and
    # accounting for the estimation through the influences ---------------------
    
    # Using the same model as when the estimated weights are used as if they 
    # were true; only the variance estimation changes --------------------------
    
    # Parameters and influences estimation -------------------------------------
    
    estimation.est  <- influences.missingdata(mod.est, riskmat.phase2 = riskmat.phase2, 
                                              dNt.phase2 = dNt.phase2, 
                                              estimated.weights = TRUE,
                                              B.phase2 = B.phase2, Tau1 = Tau1, 
                                              Tau2 = Tau2, x = x1)
    
    beta.hat.est              <- estimation.est$beta.hat
    Lambda0.Tau1Tau2.hat.est  <- estimation.est$Lambda0.Tau1Tau2.hat
    Pi.x1.Tau1Tau2.hat.est    <- estimation.est$Pi.x.Tau1Tau2.hat
    
    infl.beta.est             <- estimation.est$infl.beta
    infl.Lambda0.Tau1Tau2.est <- estimation.est$infl.Lambda0.Tau1Tau2
    infl.Pi.x1.Tau1Tau2.est   <- estimation.est$infl.Pi.x.Tau1Tau2
    
    infl2.beta.est            <- estimation.est$infl2.beta
    infl2.Lambda0.Tau1Tau2.est <- estimation.est$infl2.Lambda0.Tau1Tau2
    infl2.Pi.x1.Tau1Tau2.est  <- estimation.est$infl2.Pi.x.Tau1Tau2
    
    infl3.beta.est            <- estimation.est$infl3.beta
    infl3.Lambda0.Tau1Tau2.est <- estimation.est$infl3.Lambda0.Tau1Tau2
    infl3.Pi.x1.Tau1Tau2.est  <- estimation.est$infl3.Pi.x.Tau1Tau2
    
    estimation.est.x2 <- influences.PR.missingdata(beta = beta.hat.est,
                                                         Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.est, 
                                                         x = x2, infl2.beta = infl2.beta.est,
                                                         infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.Tau1Tau2.est,
                                                         infl3.beta = infl3.beta.est,
                                                         infl3.Lambda0.Tau1Tau2 = infl3.Lambda0.Tau1Tau2.est)
    
    estimation.est.x3 <- influences.PR.missingdata(beta = beta.hat.est,
                                                         Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.est, 
                                                         x = x3, infl2.beta = infl2.beta.est,
                                                         infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.Tau1Tau2.est,
                                                         infl3.beta = infl3.beta.est,
                                                         infl3.Lambda0.Tau1Tau2 = infl3.Lambda0.Tau1Tau2.est)
    
    Pi.x2.Tau1Tau2.hat.est    <- estimation.est.x2$Pi.x.Tau1Tau2.hat
    Pi.x3.Tau1Tau2.hat.est    <- estimation.est.x3$Pi.x.Tau1Tau2.hat
    
    infl.Pi.x2.Tau1Tau2.est   <- estimation.est.x2$infl.Pi.x.Tau1Tau2
    infl.Pi.x3.Tau1Tau2.est   <- estimation.est.x3$infl.Pi.x.Tau1Tau2
    
    infl2.Pi.x2.Tau1Tau2.est  <- estimation.est.x2$infl2.Pi.x.Tau1Tau2
    infl2.Pi.x3.Tau1Tau2.est  <- estimation.est.x3$infl2.Pi.x.Tau1Tau2
    
    infl3.Pi.x2.Tau1Tau2.est  <- estimation.est.x2$infl3.Pi.x.Tau1Tau2
    infl3.Pi.x3.Tau1Tau2.est  <- estimation.est.x3$infl3.Pi.x.Tau1Tau2
    
    # Robust variance estimation and corresponding CIs -------------------------
    coxph.var.beta.est            <- diag(mod.est$var)
    robust.var.beta.est           <- robustvariance(infl.beta.est)
    robust.var.Lambda0.Tau1Tau2.est <- robustvariance(infl.Lambda0.Tau1Tau2.est)
    robust.var.Pi.x1.Tau1Tau2.est  <- robustvariance(infl.Pi.x1.Tau1Tau2.est)
    robust.var.Pi.x2.Tau1Tau2.est  <- robustvariance(infl.Pi.x2.Tau1Tau2.est)
    robust.var.Pi.x3.Tau1Tau2.est  <- robustvariance(infl.Pi.x3.Tau1Tau2.est)
    
    CI.beta.robust.est <- conf.interval(beta.hat.est,
                                              robust.var.beta.est)
    CI.Lambda0.Tau1Tau2.robust.est <- conf.interval(Lambda0.Tau1Tau2.hat.est,
                                                          robust.var.Lambda0.Tau1Tau2.est)
    CI.Pi.x1.Tau1Tau2.robust.est <- conf.interval(Pi.x1.Tau1Tau2.hat.est, 
                                                        robust.var.Pi.x1.Tau1Tau2.est)
    CI.Pi.x2.Tau1Tau2.robust.est <- conf.interval(Pi.x2.Tau1Tau2.hat.est, 
                                                        robust.var.Pi.x2.Tau1Tau2.est)
    CI.Pi.x3.Tau1Tau2.robust.est <- conf.interval(Pi.x3.Tau1Tau2.hat.est, 
                                                        robust.var.Pi.x3.Tau1Tau2.est)
    
    res.casecohort.robust.est <- cbind("Case Cohort Robust Estimation", 
                                      t(beta.hat.est), 
                                      Lambda0.Tau1Tau2.hat.est, 
                                      Pi.x1.Tau1Tau2.hat.est,
                                      Pi.x2.Tau1Tau2.hat.est,
                                      Pi.x3.Tau1Tau2.hat.est,
                                      t(coxph.var.beta.est), 
                                      t(robust.var.beta.est), 
                                      robust.var.Lambda0.Tau1Tau2.est, 
                                      robust.var.Pi.x1.Tau1Tau2.est, 
                                      robust.var.Pi.x2.Tau1Tau2.est, 
                                      robust.var.Pi.x3.Tau1Tau2.est, 
                                      t(CI.beta.robust.est), 
                                      t(CI.Lambda0.Tau1Tau2.robust.est), 
                                      t(CI.Pi.x1.Tau1Tau2.robust.est),
                                      t(CI.Pi.x2.Tau1Tau2.robust.est),
                                      t(CI.Pi.x3.Tau1Tau2.robust.est))
    
    # Superpopulation and phase-two sampling variance estimation and 
    # corresponding CIs --------------------------------------------------------
    var.beta.est <- variance.missingdata(n = n, casecohort = casecohort, 
                                         casecohort.phase2 = casecohort.phase2, 
                                         weights = casecohort$weights.est, 
                                         weights.phase2 = casecohort.phase2$weights.est, 
                                         weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                                         infl2 = infl2.beta.est, 
                                         infl3 = infl3.beta.est, 
                                         stratified.p2 = TRUE,
                                         estimated.weights = TRUE)
    var.Lambda0.Tau1Tau2.est <- variance.missingdata(n = n, casecohort = casecohort, 
                                         casecohort.phase2 = casecohort.phase2, 
                                         weights = casecohort$weights.est, 
                                         weights.phase2 = casecohort.phase2$weights.est, 
                                         weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                                         infl2 = infl2.Lambda0.Tau1Tau2.est, 
                                         infl3 = infl3.Lambda0.Tau1Tau2.est, 
                                         stratified.p2 = TRUE,
                                         estimated.weights = TRUE)
    var.Pi.x1.Tau1Tau2.est <- variance.missingdata(n = n, casecohort = casecohort, 
                                         casecohort.phase2 = casecohort.phase2, 
                                         weights = casecohort$weights.est, 
                                         weights.phase2 = casecohort.phase2$weights.est, 
                                         weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                                         infl2 = infl2.Pi.x1.Tau1Tau2.est, 
                                         infl3 = infl3.Pi.x1.Tau1Tau2.est, 
                                         stratified.p2 = TRUE,
                                         estimated.weights = TRUE)
    var.Pi.x2.Tau1Tau2.est <- variance.missingdata(n = n, casecohort = casecohort, 
                                         casecohort.phase2 = casecohort.phase2, 
                                         weights = casecohort$weights.est, 
                                         weights.phase2 = casecohort.phase2$weights.est, 
                                         weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                                         infl2 = infl2.Pi.x2.Tau1Tau2.est, 
                                         infl3 = infl3.Pi.x2.Tau1Tau2.est, 
                                         stratified.p2 = TRUE,
                                         estimated.weights = TRUE)
    var.Pi.x3.Tau1Tau2.est <- variance.missingdata(n = n, casecohort = casecohort, 
                                         casecohort.phase2 = casecohort.phase2, 
                                         weights = casecohort$weights.est, 
                                         weights.phase2 = casecohort.phase2$weights.est, 
                                         weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                                         infl2 = infl2.Pi.x3.Tau1Tau2.est, 
                                         infl3 = infl3.Pi.x3.Tau1Tau2.est, 
                                         stratified.p2 = TRUE,
                                         estimated.weights = TRUE)
                         
    CI.beta.est <- conf.interval(beta.hat.est,var.beta.est)
    CI.Lambda0.Tau1Tau2.est <- conf.interval(Lambda0.Tau1Tau2.hat.est,
                                                   var.Lambda0.Tau1Tau2.est)
    CI.Pi.x1.Tau1Tau2.est <- conf.interval(Pi.x1.Tau1Tau2.hat.est,
                                                 var.Pi.x1.Tau1Tau2.est)
    CI.Pi.x2.Tau1Tau2.est <- conf.interval(Pi.x2.Tau1Tau2.hat.est,
                                                 var.Pi.x2.Tau1Tau2.est)
    CI.Pi.x3.Tau1Tau2.est <- conf.interval(Pi.x3.Tau1Tau2.hat.est,
                                                 var.Pi.x3.Tau1Tau2.est)

    res.casecohort.est  <- cbind("Case Cohort Estimation", t(beta.hat.est), 
                                       Lambda0.Tau1Tau2.hat.est, 
                                       Pi.x1.Tau1Tau2.hat.est,
                                       Pi.x2.Tau1Tau2.hat.est,
                                       Pi.x3.Tau1Tau2.hat.est, 
                                       t(coxph.var.beta.est), t(var.beta.est), 
                                       var.Lambda0.Tau1Tau2.est, 
                                       var.Pi.x1.Tau1Tau2.est,
                                       var.Pi.x2.Tau1Tau2.est,
                                       var.Pi.x3.Tau1Tau2.est, t(CI.beta.est), 
                                       t(CI.Lambda0.Tau1Tau2.est), 
                                       t(CI.Pi.x1.Tau1Tau2.est), 
                                       t(CI.Pi.x2.Tau1Tau2.est),
                                       t(CI.Pi.x3.Tau1Tau2.est))
    
    # --------------------------------------------------------------------------
    # Estimation using the un-stratified case cohort with design weights -------
    
    mod.unstrat                 <- coxph(Surv(times, status) ~ X1 + X2 + X3, 
                                   data = unstrat.casecohort,
                                   weight = unstrat.weights.true, 
                                   id = id, robust = TRUE)
    
    # At risk indicator matrix and counting process matrix for the phase-two 
    # data at all event times, even that with missing covariate data -----------
    mod.cohort.detail   <- coxph.detail(mod.cohort, riskmat = TRUE)   
    riskmat.unstrat.phase2  <- with(cohort, 
                                    mod.cohort.detail$riskmat[unstrat.phase2 == 1,])
    rownames(riskmat.unstrat.phase2)  <- cohort[cohort$unstrat.phase2 == 1, "id"]
    observed.times.unstrat.phase2 <- apply(riskmat.unstrat.phase2, 1,
                                   function(v) {which.max(cumsum(v))})
    dNt.unstrat.phase2            <- matrix(0, nrow(riskmat.unstrat.phase2), 
                                    ncol(riskmat.unstrat.phase2))
    dNt.unstrat.phase2[cbind(1:nrow(riskmat.unstrat.phase2), 
                             observed.times.unstrat.phase2)] <- 1
    dNt.unstrat.phase2            <- sweep(dNt.unstrat.phase2, 1, 
                                           unstrat.casecohort.phase2$status, 
                                           "*")
    colnames(dNt.unstrat.phase2)  <- colnames(riskmat.unstrat.phase2)
    rownames(dNt.unstrat.phase2)  <- rownames(riskmat.unstrat.phase2)
    
    # Parameters and influences estimation -------------------------------------
    estimation.unstrat <- influences.missingdata(mod = mod.unstrat, 
                                         riskmat.phase2 = riskmat.unstrat.phase2, 
                                         dNt.phase2 = dNt.unstrat.phase2, 
                                         Tau1 = Tau1, 
                                         Tau2 = Tau2, x = x1)
    
    beta.hat.unstrat                <- estimation.unstrat$beta.hat
    Lambda0.Tau1Tau2.hat.unstrat    <- estimation.unstrat$Lambda0.Tau1Tau2.hat
    Pi.x1.Tau1Tau2.hat.unstrat      <- estimation.unstrat$Pi.x.Tau1Tau2.hat
    
    infl.beta.unstrat              <- estimation.unstrat$infl.beta
    infl.Lambda0.Tau1Tau2.unstrat  <- estimation.unstrat$infl.Lambda0.Tau1Tau2
    infl.Pi.x1.Tau1Tau2.unstrat    <- estimation.unstrat$infl.Pi.x.Tau1Tau2
    
    infl2.beta.unstrat              <- estimation.unstrat$infl2.beta
    infl2.Lambda0.Tau1Tau2.unstrat  <- estimation.unstrat$infl2.Lambda0.Tau1Tau2
    infl2.Pi.x1.Tau1Tau2.unstrat    <- estimation.unstrat$infl2.Pi.x.Tau1Tau2
    
    infl3.beta.unstrat              <- estimation.unstrat$infl3.beta
    infl3.Lambda0.Tau1Tau2.unstrat  <- estimation.unstrat$infl3.Lambda0.Tau1Tau2
    infl3.Pi.x1.Tau1Tau2.unstrat    <- estimation.unstrat$infl3.Pi.x.Tau1Tau2
    
    estimation.x2.unstrat <- influences.PR.missingdata(beta = beta.hat.unstrat, 
                                               Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.unstrat, 
                                               x = x2, infl2.beta = infl2.beta.unstrat, 
                                               infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.Tau1Tau2.unstrat,
                                               infl3.beta = infl3.beta.unstrat, 
                                               infl3.Lambda0.Tau1Tau2 = infl3.Lambda0.Tau1Tau2.unstrat)
    estimation.x3.unstrat <- influences.PR.missingdata(beta = beta.hat.unstrat, 
                                               Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.unstrat, 
                                               x = x3, infl2.beta = infl2.beta.unstrat, 
                                               infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.Tau1Tau2.unstrat,
                                               infl3.beta = infl3.beta.unstrat, 
                                               infl3.Lambda0.Tau1Tau2 = infl3.Lambda0.Tau1Tau2.unstrat)
    
    Pi.x2.Tau1Tau2.hat.unstrat     <- estimation.x2.unstrat$Pi.x.Tau1Tau2.hat
    Pi.x3.Tau1Tau2.hat.unstrat     <- estimation.x3.unstrat$Pi.x.Tau1Tau2.hat
    
    infl.Pi.x2.Tau1Tau2.unstrat    <- estimation.x2.unstrat$infl.Pi.x.Tau1Tau2
    infl.Pi.x3.Tau1Tau2.unstrat    <- estimation.x3.unstrat$infl.Pi.x.Tau1Tau2
    
    infl2.Pi.x2.Tau1Tau2.unstrat    <- estimation.x2.unstrat$infl2.Pi.x.Tau1Tau2
    infl2.Pi.x3.Tau1Tau2.unstrat    <- estimation.x3.unstrat$infl2.Pi.x.Tau1Tau2
    
    infl3.Pi.x2.Tau1Tau2.unstrat    <- estimation.x2.unstrat$infl3.Pi.x.Tau1Tau2
    infl3.Pi.x3.Tau1Tau2.unstrat    <- estimation.x3.unstrat$infl3.Pi.x.Tau1Tau2
    
    # Robust variance estimation and corresponding CIs -------------------------
    coxph.var.beta.unstrat              <- diag(mod.unstrat$var) 
    robust.var.beta.unstrat             <- robustvariance(infl.beta.unstrat) 
    robust.var.Lambda0.Tau1Tau2.unstrat <- robustvariance(infl.Lambda0.Tau1Tau2.unstrat)
    robust.var.Pi.x1.Tau1Tau2.unstrat    <- robustvariance(infl.Pi.x1.Tau1Tau2.unstrat)
    robust.var.Pi.x2.Tau1Tau2.unstrat    <- robustvariance(infl.Pi.x2.Tau1Tau2.unstrat)
    robust.var.Pi.x3.Tau1Tau2.unstrat    <- robustvariance(infl.Pi.x3.Tau1Tau2.unstrat)
    
    CI.beta.robust.unstrat  <- conf.interval(beta.hat.unstrat, 
                                             robust.var.beta.unstrat)
    CI.Lambda0.Tau1Tau2.robust.unstrat <- conf.interval(Lambda0.Tau1Tau2.hat.unstrat, 
                                                robust.var.Lambda0.Tau1Tau2.unstrat)
    CI.Pi.x1.Tau1Tau2.robust.unstrat <- conf.interval(Pi.x1.Tau1Tau2.hat.unstrat,
                                              robust.var.Pi.x1.Tau1Tau2.unstrat)
    CI.Pi.x2.Tau1Tau2.robust.unstrat <- conf.interval(Pi.x2.Tau1Tau2.hat.unstrat,
                                              robust.var.Pi.x2.Tau1Tau2.unstrat)
    CI.Pi.x3.Tau1Tau2.robust.unstrat <- conf.interval(Pi.x3.Tau1Tau2.hat.unstrat,
                                              robust.var.Pi.x3.Tau1Tau2.unstrat)
    
    res.unstrat.robust <- cbind("Unstrat Case Cohort Robust", 
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
    
    var.beta.unstrat <- variance.missingdata(n = n, 
                                             casecohort = unstrat.casecohort, 
                                             casecohort.phase2 = unstrat.casecohort.phase2, 
                                             weights = unstrat.casecohort$unstrat.weights.true, 
                                             weights.phase2 = unstrat.casecohort.phase2$unstrat.weights.true, 
                                             weights.p2.phase2 = unstrat.casecohort.phase2$unstrat.weights.p2.true,
                                             infl2 = infl2.beta.unstrat,
                                             infl3 = infl3.beta.unstrat, 
                                             stratified.p2 = FALSE)
    var.Lambda0.Tau1Tau2.unstrat <- variance.missingdata(n = n, 
                                             casecohort = unstrat.casecohort, 
                                             casecohort.phase2 = unstrat.casecohort.phase2, 
                                             weights = unstrat.casecohort$unstrat.weights.true, 
                                             weights.phase2 = unstrat.casecohort.phase2$unstrat.weights.true, 
                                             weights.p2.phase2 = unstrat.casecohort.phase2$unstrat.weights.p2.true,
                                             infl2 = infl2.Lambda0.Tau1Tau2.unstrat,
                                             infl3 = infl3.Lambda0.Tau1Tau2.unstrat, 
                                             stratified.p2 = FALSE)
    var.Pi.x1.Tau1Tau2.unstrat <- variance.missingdata(n = n, 
                                             casecohort = unstrat.casecohort, 
                                             casecohort.phase2 = unstrat.casecohort.phase2, 
                                             weights = unstrat.casecohort$unstrat.weights.true, 
                                             weights.phase2 = unstrat.casecohort.phase2$unstrat.weights.true, 
                                             weights.p2.phase2 = unstrat.casecohort.phase2$unstrat.weights.p2.true,
                                             infl2 = infl2.Pi.x1.Tau1Tau2.unstrat,
                                             infl3 = infl3.Pi.x1.Tau1Tau2.unstrat, 
                                             stratified.p2 = FALSE)
    var.Pi.x2.Tau1Tau2.unstrat <- variance.missingdata(n = n, 
                                             casecohort = unstrat.casecohort, 
                                             casecohort.phase2 = unstrat.casecohort.phase2, 
                                             weights = unstrat.casecohort$unstrat.weights.true, 
                                             weights.phase2 = unstrat.casecohort.phase2$unstrat.weights.true, 
                                             weights.p2.phase2 = unstrat.casecohort.phase2$unstrat.weights.p2.true,
                                             infl2 = infl2.Pi.x2.Tau1Tau2.unstrat,
                                             infl3 = infl3.Pi.x2.Tau1Tau2.unstrat, 
                                             stratified.p2 = FALSE)
    var.Pi.x3.Tau1Tau2.unstrat <- variance.missingdata(n = n, 
                                             casecohort = unstrat.casecohort, 
                                             casecohort.phase2 = unstrat.casecohort.phase2, 
                                             weights = unstrat.casecohort$unstrat.weights.true, 
                                             weights.phase2 = unstrat.casecohort.phase2$unstrat.weights.true, 
                                             weights.p2.phase2 = unstrat.casecohort.phase2$unstrat.weights.p2.true,
                                             infl2 = infl2.Pi.x3.Tau1Tau2.unstrat,
                                             infl3 = infl3.Pi.x3.Tau1Tau2.unstrat, 
                                             stratified.p2 = FALSE)

    CI.beta.unstrat   <- conf.interval(beta.hat.unstrat, var.beta.unstrat)
    CI.Lambda0.Tau1Tau2.unstrat   <- conf.interval(Lambda0.Tau1Tau2.hat.unstrat, 
                                           var.Lambda0.Tau1Tau2.unstrat)
    CI.Pi.x1.Tau1Tau2.unstrat  <- conf.interval(Pi.x1.Tau1Tau2.hat.unstrat,
                                        var.Pi.x1.Tau1Tau2.unstrat)
    CI.Pi.x2.Tau1Tau2.unstrat  <- conf.interval(Pi.x2.Tau1Tau2.hat.unstrat,
                                                var.Pi.x2.Tau1Tau2.unstrat)
    CI.Pi.x3.Tau1Tau2.unstrat  <- conf.interval(Pi.x3.Tau1Tau2.hat.unstrat,
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
    # Estimation using the un-stratified case cohort with estimated weights as 
    # if they were the true weights --------------------------------------------
    
    mod.unstrat.est <- coxph(Surv(times, status) ~ X1 + X2 + X3, 
                     data = unstrat.casecohort, weight = unstrat.weights.est, 
                     id = id, robust = TRUE)
    
    # parameters and influences estimation -------------------------------------
    
    estimation.unstrat.est.naive  <- influences.missingdata(mod.unstrat.est, 
                                                    riskmat.phase2 = riskmat.unstrat.phase2, 
                                                    dNt.phase2 = dNt.unstrat.phase2, 
                                                    Tau1 = Tau1, 
                                                    Tau2 = Tau2, x = x1)
    
    beta.hat.unstrat.est.naive              <- estimation.unstrat.est.naive$beta.hat
    Lambda0.Tau1Tau2.hat.unstrat.est.naive  <- estimation.unstrat.est.naive$Lambda0.Tau1Tau2.hat
    Pi.x1.Tau1Tau2.hat.unstrat.est.naive    <- estimation.unstrat.est.naive$Pi.x.Tau1Tau2.hat
    
    infl.beta.unstrat.est.naive             <- estimation.unstrat.est.naive$infl.beta
    infl.Lambda0.Tau1Tau2.unstrat.est.naive <- estimation.unstrat.est.naive$infl.Lambda0.Tau1Tau2
    infl.Pi.x1.Tau1Tau2.unstrat.est.naive   <- estimation.unstrat.est.naive$infl.Pi.x.Tau1Tau2
    
    infl2.beta.unstrat.est.naive            <- estimation.unstrat.est.naive$infl2.beta
    infl2.Lambda0.Tau1Tau2.unstrat.est.naive <- estimation.unstrat.est.naive$infl2.Lambda0.Tau1Tau2
    infl2.Pi.x1.Tau1Tau2.unstrat.est.naive  <- estimation.unstrat.est.naive$infl2.Pi.x.Tau1Tau2
    
    infl3.beta.unstrat.est.naive            <- estimation.unstrat.est.naive$infl3.beta
    infl3.Lambda0.Tau1Tau2.unstrat.est.naive <- estimation.unstrat.est.naive$infl3.Lambda0.Tau1Tau2
    infl3.Pi.x1.Tau1Tau2.unstrat.est.naive  <- estimation.unstrat.est.naive$infl3.Pi.x.Tau1Tau2
    
    estimation.unstrat.est.naive.x2 <- influences.PR.missingdata(beta = beta.hat.unstrat.est.naive,
                                                         Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.unstrat.est.naive, 
                                                         x = x2, infl2.beta = infl2.beta.unstrat.est.naive,
                                                         infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.Tau1Tau2.unstrat.est.naive,
                                                         infl3.beta = infl3.beta.unstrat.est.naive,
                                                         infl3.Lambda0.Tau1Tau2 = infl3.Lambda0.Tau1Tau2.unstrat.est.naive)
    
    estimation.unstrat.est.naive.x3 <- influences.PR.missingdata(beta = beta.hat.unstrat.est.naive,
                                                         Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.unstrat.est.naive, 
                                                         x = x3, infl2.beta = infl2.beta.unstrat.est.naive,
                                                         infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.Tau1Tau2.unstrat.est.naive,
                                                         infl3.beta = infl3.beta.unstrat.est.naive,
                                                         infl3.Lambda0.Tau1Tau2 = infl3.Lambda0.Tau1Tau2.unstrat.est.naive)
    
    Pi.x2.Tau1Tau2.hat.unstrat.est.naive    <- estimation.unstrat.est.naive.x2$Pi.x.Tau1Tau2.hat
    Pi.x3.Tau1Tau2.hat.unstrat.est.naive    <- estimation.unstrat.est.naive.x3$Pi.x.Tau1Tau2.hat
    
    infl.Pi.x2.Tau1Tau2.unstrat.est.naive   <- estimation.unstrat.est.naive.x2$infl.Pi.x.Tau1Tau2
    infl.Pi.x3.Tau1Tau2.unstrat.est.naive   <- estimation.unstrat.est.naive.x3$infl.Pi.x.Tau1Tau2
    
    infl2.Pi.x2.Tau1Tau2.unstrat.est.naive  <- estimation.unstrat.est.naive.x2$infl2.Pi.x.Tau1Tau2
    infl2.Pi.x3.Tau1Tau2.unstrat.est.naive  <- estimation.unstrat.est.naive.x3$infl2.Pi.x.Tau1Tau2
    
    infl3.Pi.x2.Tau1Tau2.unstrat.est.naive  <- estimation.unstrat.est.naive.x2$infl3.Pi.x.Tau1Tau2
    infl3.Pi.x3.Tau1Tau2.unstrat.est.naive  <- estimation.unstrat.est.naive.x3$infl3.Pi.x.Tau1Tau2
    
    # Robust variance estimation and corresponding CIs -------------------------
    coxph.var.beta.unstrat.est.naive            <- diag(mod.unstrat.est$var)
    robust.var.beta.unstrat.est.naive           <- robustvariance(infl.beta.unstrat.est.naive)
    robust.var.Lambda0.Tau1Tau2.unstrat.est.naive <- robustvariance(infl.Lambda0.Tau1Tau2.unstrat.est.naive)
    robust.var.Pi.x1.Tau1Tau2.unstrat.est.naive  <- robustvariance(infl.Pi.x1.Tau1Tau2.unstrat.est.naive)
    robust.var.Pi.x2.Tau1Tau2.unstrat.est.naive  <- robustvariance(infl.Pi.x2.Tau1Tau2.unstrat.est.naive)
    robust.var.Pi.x3.Tau1Tau2.unstrat.est.naive  <- robustvariance(infl.Pi.x3.Tau1Tau2.unstrat.est.naive)
    
    CI.beta.robust.unstrat.est.naive <- conf.interval(beta.hat.unstrat.est.naive,
                                              robust.var.beta.unstrat.est.naive)
    CI.Lambda0.Tau1Tau2.robust.unstrat.est.naive <- conf.interval(Lambda0.Tau1Tau2.hat.unstrat.est.naive,
                                                          robust.var.Lambda0.Tau1Tau2.unstrat.est.naive)
    CI.Pi.x1.Tau1Tau2.robust.unstrat.est.naive <- conf.interval(Pi.x1.Tau1Tau2.hat.unstrat.est.naive, 
                                                        robust.var.Pi.x1.Tau1Tau2.unstrat.est.naive)
    CI.Pi.x2.Tau1Tau2.robust.unstrat.est.naive <- conf.interval(Pi.x2.Tau1Tau2.hat.unstrat.est.naive, 
                                                        robust.var.Pi.x2.Tau1Tau2.unstrat.est.naive)
    CI.Pi.x3.Tau1Tau2.robust.unstrat.est.naive <- conf.interval(Pi.x3.Tau1Tau2.hat.unstrat.est.naive, 
                                                        robust.var.Pi.x3.Tau1Tau2.unstrat.est.naive)
    
    res.unstrat.est.robust.naive <- cbind("Unstrat Case Cohort Robust Estimation Naive", 
                                      t(beta.hat.unstrat.est.naive), 
                                      Lambda0.Tau1Tau2.hat.unstrat.est.naive, 
                                      Pi.x1.Tau1Tau2.hat.unstrat.est.naive,
                                      Pi.x2.Tau1Tau2.hat.unstrat.est.naive,
                                      Pi.x3.Tau1Tau2.hat.unstrat.est.naive,
                                      t(coxph.var.beta.unstrat.est.naive), 
                                      t(robust.var.beta.unstrat.est.naive), 
                                      robust.var.Lambda0.Tau1Tau2.unstrat.est.naive, 
                                      robust.var.Pi.x1.Tau1Tau2.unstrat.est.naive, 
                                      robust.var.Pi.x2.Tau1Tau2.unstrat.est.naive, 
                                      robust.var.Pi.x3.Tau1Tau2.unstrat.est.naive, 
                                      t(CI.beta.robust.unstrat.est.naive), 
                                      t(CI.Lambda0.Tau1Tau2.robust.unstrat.est.naive), 
                                      t(CI.Pi.x1.Tau1Tau2.robust.unstrat.est.naive),
                                      t(CI.Pi.x2.Tau1Tau2.robust.unstrat.est.naive),
                                      t(CI.Pi.x3.Tau1Tau2.robust.unstrat.est.naive))
    
    # Superpopulation and phase-two sampling variance estimation and 
    # corresponding CIs --------------------------------------------------------
    var.beta.unstrat.est.naive <- variance.missingdata(n = n, casecohort = unstrat.casecohort, 
                                                       casecohort.phase2 = unstrat.casecohort.phase2, 
                                                       weights = unstrat.casecohort$unstrat.weights.est, 
                                                       weights.phase2 = unstrat.casecohort.phase2$unstrat.weights.est, 
                                                       weights.p2.phase2 = unstrat.casecohort.phase2$unstrat.weights.p2.true,
                                                       infl2 = infl2.beta.unstrat.est.naive, 
                                                       infl3 = infl3.beta.unstrat.est.naive, 
                                                       stratified.p2 = FALSE)
    var.Lambda0.Tau1Tau2.unstrat.est.naive <- variance.missingdata(n = n, casecohort = unstrat.casecohort, 
                                                       casecohort.phase2 = unstrat.casecohort.phase2, 
                                                       weights = unstrat.casecohort$unstrat.weights.est, 
                                                       weights.phase2 = unstrat.casecohort.phase2$unstrat.weights.est, 
                                                       weights.p2.phase2 = unstrat.casecohort.phase2$unstrat.weights.p2.true,
                                                       infl2 = infl2.Lambda0.Tau1Tau2.unstrat.est.naive, 
                                                       infl3 = infl3.Lambda0.Tau1Tau2.unstrat.est.naive, 
                                                       stratified.p2 = FALSE)
    var.Pi.x1.Tau1Tau2.unstrat.est.naive <- variance.missingdata(n = n, casecohort = unstrat.casecohort, 
                                                       casecohort.phase2 = unstrat.casecohort.phase2, 
                                                       weights = unstrat.casecohort$unstrat.weights.est, 
                                                       weights.phase2 = unstrat.casecohort.phase2$unstrat.weights.est, 
                                                       weights.p2.phase2 = unstrat.casecohort.phase2$unstrat.weights.p2.true,
                                                       infl2 = infl2.Pi.x1.Tau1Tau2.unstrat.est.naive, 
                                                       infl3 = infl3.Pi.x1.Tau1Tau2.unstrat.est.naive, 
                                                       stratified.p2 = FALSE)
    var.Pi.x2.Tau1Tau2.unstrat.est.naive <- variance.missingdata(n = n, casecohort = unstrat.casecohort, 
                                                       casecohort.phase2 = unstrat.casecohort.phase2, 
                                                       weights = unstrat.casecohort$unstrat.weights.est, 
                                                       weights.phase2 = unstrat.casecohort.phase2$unstrat.weights.est, 
                                                       weights.p2.phase2 = unstrat.casecohort.phase2$unstrat.weights.p2.true,
                                                       infl2 = infl2.Pi.x2.Tau1Tau2.unstrat.est.naive, 
                                                       infl3 = infl3.Pi.x2.Tau1Tau2.unstrat.est.naive, 
                                                       stratified.p2 = FALSE)
    var.Pi.x3.Tau1Tau2.unstrat.est.naive <- variance.missingdata(n = n, casecohort = unstrat.casecohort, 
                                                       casecohort.phase2 = unstrat.casecohort.phase2, 
                                                       weights = unstrat.casecohort$unstrat.weights.est, 
                                                       weights.phase2 = unstrat.casecohort.phase2$unstrat.weights.est, 
                                                       weights.p2.phase2 = unstrat.casecohort.phase2$unstrat.weights.p2.true,
                                                       infl2 = infl2.Pi.x3.Tau1Tau2.unstrat.est.naive, 
                                                       infl3 = infl3.Pi.x3.Tau1Tau2.unstrat.est.naive, 
                                                       stratified.p2 = FALSE)

    CI.beta.unstrat.est.naive <- conf.interval(beta.hat.unstrat.est.naive,var.beta.unstrat.est.naive)
    CI.Lambda0.Tau1Tau2.unstrat.est.naive <- conf.interval(Lambda0.Tau1Tau2.hat.unstrat.est.naive,
                                                   var.Lambda0.Tau1Tau2.unstrat.est.naive)
    CI.Pi.x1.Tau1Tau2.unstrat.est.naive <- conf.interval(Pi.x1.Tau1Tau2.hat.unstrat.est.naive,
                                                 var.Pi.x1.Tau1Tau2.unstrat.est.naive)
    CI.Pi.x2.Tau1Tau2.unstrat.est.naive <- conf.interval(Pi.x2.Tau1Tau2.hat.unstrat.est.naive,
                                                 var.Pi.x2.Tau1Tau2.unstrat.est.naive)
    CI.Pi.x3.Tau1Tau2.unstrat.est.naive <- conf.interval(Pi.x3.Tau1Tau2.hat.unstrat.est.naive,
                                                 var.Pi.x3.Tau1Tau2.unstrat.est.naive)
    
    res.unstrat.est.naive  <- cbind("Unstrat Case Cohort Estimation Naive", t(beta.hat.unstrat.est.naive), 
                                       Lambda0.Tau1Tau2.hat.unstrat.est.naive, 
                                       Pi.x1.Tau1Tau2.hat.unstrat.est.naive,
                                       Pi.x2.Tau1Tau2.hat.unstrat.est.naive,
                                       Pi.x3.Tau1Tau2.hat.unstrat.est.naive, 
                                       t(coxph.var.beta.unstrat.est.naive), t(var.beta.unstrat.est.naive), 
                                       var.Lambda0.Tau1Tau2.unstrat.est.naive, 
                                       var.Pi.x1.Tau1Tau2.unstrat.est.naive,
                                       var.Pi.x2.Tau1Tau2.unstrat.est.naive,
                                       var.Pi.x3.Tau1Tau2.unstrat.est.naive, t(CI.beta.unstrat.est.naive), 
                                       t(CI.Lambda0.Tau1Tau2.unstrat.est.naive), 
                                       t(CI.Pi.x1.Tau1Tau2.unstrat.est.naive), 
                                       t(CI.Pi.x2.Tau1Tau2.unstrat.est.naive),
                                       t(CI.Pi.x3.Tau1Tau2.unstrat.est.naive))
    
    # --------------------------------------------------------------------------
    # Estimation using the un-stratified case cohort with estimated weights, and
    # accounting for the estimation through the influences ---------------------
    
    # Using the same model as when the estimated weights are used as if they 
    # were true; only the variance estimation changes --------------------------
    
    # Parameters and influences estimation -------------------------------------
    
    estimation.unstrat.est  <- influences.missingdata(mod.unstrat.est, 
                                                      riskmat.phase2 = riskmat.unstrat.phase2, 
                                              dNt.phase2 = dNt.unstrat.phase2, 
                                              estimated.weights = TRUE,
                                              B.phase2 = B.phase2.unstrat, Tau1 = Tau1, 
                                              Tau2 = Tau2, x = x1)
    
    beta.hat.unstrat.est              <- estimation.unstrat.est$beta.hat
    Lambda0.Tau1Tau2.hat.unstrat.est  <- estimation.unstrat.est$Lambda0.Tau1Tau2.hat
    Pi.x1.Tau1Tau2.hat.unstrat.est    <- estimation.unstrat.est$Pi.x.Tau1Tau2.hat
    
    infl.beta.unstrat.est             <- estimation.unstrat.est$infl.beta
    infl.Lambda0.Tau1Tau2.unstrat.est <- estimation.unstrat.est$infl.Lambda0.Tau1Tau2
    infl.Pi.x1.Tau1Tau2.unstrat.est   <- estimation.unstrat.est$infl.Pi.x.Tau1Tau2
    
    infl2.beta.unstrat.est            <- estimation.unstrat.est$infl2.beta
    infl2.Lambda0.Tau1Tau2.unstrat.est <- estimation.unstrat.est$infl2.Lambda0.Tau1Tau2
    infl2.Pi.x1.Tau1Tau2.unstrat.est  <- estimation.unstrat.est$infl2.Pi.x.Tau1Tau2
    
    infl3.beta.unstrat.est            <- estimation.unstrat.est$infl3.beta
    infl3.Lambda0.Tau1Tau2.unstrat.est <- estimation.unstrat.est$infl3.Lambda0.Tau1Tau2
    infl3.Pi.x1.Tau1Tau2.unstrat.est  <- estimation.unstrat.est$infl3.Pi.x.Tau1Tau2
    
    estimation.unstrat.est.x2 <- influences.PR.missingdata(beta = beta.hat.unstrat.est,
                                                   Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.unstrat.est, 
                                                   x = x2, infl2.beta = infl2.beta.unstrat.est,
                                                   infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.Tau1Tau2.unstrat.est,
                                                   infl3.beta = infl3.beta.unstrat.est,
                                                   infl3.Lambda0.Tau1Tau2 = infl3.Lambda0.Tau1Tau2.unstrat.est)
    
    estimation.unstrat.est.x3 <- influences.PR.missingdata(beta = beta.hat.unstrat.est,
                                                   Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.unstrat.est, 
                                                   x = x3, infl2.beta = infl2.beta.unstrat.est,
                                                   infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.Tau1Tau2.unstrat.est,
                                                   infl3.beta = infl3.beta.unstrat.est,
                                                   infl3.Lambda0.Tau1Tau2 = infl3.Lambda0.Tau1Tau2.unstrat.est)
    
    Pi.x2.Tau1Tau2.hat.unstrat.est    <- estimation.unstrat.est.x2$Pi.x.Tau1Tau2.hat
    Pi.x3.Tau1Tau2.hat.unstrat.est    <- estimation.unstrat.est.x3$Pi.x.Tau1Tau2.hat
    
    infl.Pi.x2.Tau1Tau2.unstrat.est   <- estimation.unstrat.est.x2$infl.Pi.x.Tau1Tau2
    infl.Pi.x3.Tau1Tau2.unstrat.est   <- estimation.unstrat.est.x3$infl.Pi.x.Tau1Tau2
    
    infl2.Pi.x2.Tau1Tau2.unstrat.est  <- estimation.unstrat.est.x2$infl2.Pi.x.Tau1Tau2
    infl2.Pi.x3.Tau1Tau2.unstrat.est  <- estimation.unstrat.est.x3$infl2.Pi.x.Tau1Tau2
    
    infl3.Pi.x2.Tau1Tau2.unstrat.est  <- estimation.unstrat.est.x2$infl3.Pi.x.Tau1Tau2
    infl3.Pi.x3.Tau1Tau2.unstrat.est  <- estimation.unstrat.est.x3$infl3.Pi.x.Tau1Tau2
    
    # Robust variance estimation and corresponding CIs -------------------------
    coxph.var.beta.unstrat.est            <- diag(mod.unstrat.est$var)
    robust.var.beta.unstrat.est           <- robustvariance(infl.beta.unstrat.est)
    robust.var.Lambda0.Tau1Tau2.unstrat.est <- robustvariance(infl.Lambda0.Tau1Tau2.unstrat.est)
    robust.var.Pi.x1.Tau1Tau2.unstrat.est  <- robustvariance(infl.Pi.x1.Tau1Tau2.unstrat.est)
    robust.var.Pi.x2.Tau1Tau2.unstrat.est  <- robustvariance(infl.Pi.x2.Tau1Tau2.unstrat.est)
    robust.var.Pi.x3.Tau1Tau2.unstrat.est  <- robustvariance(infl.Pi.x3.Tau1Tau2.unstrat.est)
    
    CI.beta.robust.unstrat.est <- conf.interval(beta.hat.unstrat.est,
                                        robust.var.beta.unstrat.est)
    CI.Lambda0.Tau1Tau2.robust.unstrat.est <- conf.interval(Lambda0.Tau1Tau2.hat.unstrat.est,
                                                    robust.var.Lambda0.Tau1Tau2.unstrat.est)
    CI.Pi.x1.Tau1Tau2.robust.unstrat.est <- conf.interval(Pi.x1.Tau1Tau2.hat.unstrat.est, 
                                                  robust.var.Pi.x1.Tau1Tau2.unstrat.est)
    CI.Pi.x2.Tau1Tau2.robust.unstrat.est <- conf.interval(Pi.x2.Tau1Tau2.hat.unstrat.est, 
                                                  robust.var.Pi.x2.Tau1Tau2.unstrat.est)
    CI.Pi.x3.Tau1Tau2.robust.unstrat.est <- conf.interval(Pi.x3.Tau1Tau2.hat.unstrat.est, 
                                                  robust.var.Pi.x3.Tau1Tau2.unstrat.est)
    
    res.unstrat.est.robust <- cbind("Unstrat Case Cohort Robust Estimation", 
                                t(beta.hat.unstrat.est), 
                                Lambda0.Tau1Tau2.hat.unstrat.est, 
                                Pi.x1.Tau1Tau2.hat.unstrat.est,
                                Pi.x2.Tau1Tau2.hat.unstrat.est,
                                Pi.x3.Tau1Tau2.hat.unstrat.est,
                                t(coxph.var.beta.unstrat.est), 
                                t(robust.var.beta.unstrat.est), 
                                robust.var.Lambda0.Tau1Tau2.unstrat.est, 
                                robust.var.Pi.x1.Tau1Tau2.unstrat.est, 
                                robust.var.Pi.x2.Tau1Tau2.unstrat.est, 
                                robust.var.Pi.x3.Tau1Tau2.unstrat.est, 
                                t(CI.beta.robust.unstrat.est), 
                                t(CI.Lambda0.Tau1Tau2.robust.unstrat.est), 
                                t(CI.Pi.x1.Tau1Tau2.robust.unstrat.est),
                                t(CI.Pi.x2.Tau1Tau2.robust.unstrat.est),
                                t(CI.Pi.x3.Tau1Tau2.robust.unstrat.est))
    

    # Superpopulation and phase-two sampling variance estimation and 
    # corresponding CIs --------------------------------------------------------
    var.beta.unstrat.est <- variance.missingdata(n = n, 
                                                 casecohort = unstrat.casecohort, 
                                                 casecohort.phase2 = unstrat.casecohort.phase2, 
                                                 weights = unstrat.casecohort$unstrat.weights.est, 
                                                 weights.phase2 = unstrat.casecohort.phase2$unstrat.weights.est, 
                                                 weights.p2.phase2 = unstrat.casecohort.phase2$unstrat.weights.p2.true,
                                                 infl2 = infl2.beta.unstrat.est, 
                                                 infl3 = infl3.beta.unstrat.est, 
                                                 stratified.p2 = FALSE,
                                                 estimated.weights = TRUE)
    var.Lambda0.Tau1Tau2.unstrat.est <- variance.missingdata(n = n, 
                                                 casecohort = unstrat.casecohort, 
                                                 casecohort.phase2 = unstrat.casecohort.phase2, 
                                                 weights = unstrat.casecohort$unstrat.weights.est, 
                                                 weights.phase2 = unstrat.casecohort.phase2$unstrat.weights.est, 
                                                 weights.p2.phase2 = unstrat.casecohort.phase2$unstrat.weights.p2.true,
                                                 infl2 = infl2.Lambda0.Tau1Tau2.unstrat.est, 
                                                 infl3 = infl3.Lambda0.Tau1Tau2.unstrat.est, 
                                                 stratified.p2 = FALSE,
                                                 estimated.weights = TRUE)
    var.Pi.x1.Tau1Tau2.unstrat.est <- variance.missingdata(n = n, 
                                                 casecohort = unstrat.casecohort, 
                                                 casecohort.phase2 = unstrat.casecohort.phase2, 
                                                 weights = unstrat.casecohort$unstrat.weights.est, 
                                                 weights.phase2 = unstrat.casecohort.phase2$unstrat.weights.est, 
                                                 weights.p2.phase2 = unstrat.casecohort.phase2$unstrat.weights.p2.true,
                                                 infl2 = infl2.Pi.x1.Tau1Tau2.unstrat.est, 
                                                 infl3 = infl3.Pi.x1.Tau1Tau2.unstrat.est, 
                                                 stratified.p2 = FALSE,
                                                 estimated.weights = TRUE)
    var.Pi.x2.Tau1Tau2.unstrat.est <- variance.missingdata(n = n, 
                                                 casecohort = unstrat.casecohort, 
                                                 casecohort.phase2 = unstrat.casecohort.phase2, 
                                                 weights = unstrat.casecohort$unstrat.weights.est, 
                                                 weights.phase2 = unstrat.casecohort.phase2$unstrat.weights.est, 
                                                 weights.p2.phase2 = unstrat.casecohort.phase2$unstrat.weights.p2.true,
                                                 infl2 = infl2.Pi.x2.Tau1Tau2.unstrat.est, 
                                                 infl3 = infl3.Pi.x2.Tau1Tau2.unstrat.est, 
                                                 stratified.p2 = FALSE,
                                                 estimated.weights = TRUE)
    var.Pi.x3.Tau1Tau2.unstrat.est <- variance.missingdata(n = n, 
                                                 casecohort = unstrat.casecohort, 
                                                 casecohort.phase2 = unstrat.casecohort.phase2, 
                                                 weights = unstrat.casecohort$unstrat.weights.est, 
                                                 weights.phase2 = unstrat.casecohort.phase2$unstrat.weights.est, 
                                                 weights.p2.phase2 = unstrat.casecohort.phase2$unstrat.weights.p2.true,
                                                 infl2 = infl2.Pi.x3.Tau1Tau2.unstrat.est, 
                                                 infl3 = infl3.Pi.x3.Tau1Tau2.unstrat.est, 
                                                 stratified.p2 = FALSE,
                                                 estimated.weights = TRUE)
                                                          
    CI.beta.unstrat.est <- conf.interval(beta.hat.unstrat.est, var.beta.unstrat.est)
    CI.Lambda0.Tau1Tau2.unstrat.est <- conf.interval(Lambda0.Tau1Tau2.hat.unstrat.est,
                                             var.Lambda0.Tau1Tau2.unstrat.est)
    CI.Pi.x1.Tau1Tau2.unstrat.est <- conf.interval(Pi.x1.Tau1Tau2.hat.unstrat.est,
                                           var.Pi.x1.Tau1Tau2.unstrat.est)
    CI.Pi.x2.Tau1Tau2.unstrat.est <- conf.interval(Pi.x2.Tau1Tau2.hat.unstrat.est,
                                           var.Pi.x2.Tau1Tau2.unstrat.est)
    CI.Pi.x3.Tau1Tau2.unstrat.est <- conf.interval(Pi.x3.Tau1Tau2.hat.unstrat.est,
                                           var.Pi.x3.Tau1Tau2.unstrat.est)
    
    res.unstrat.est  <- cbind("Unstrat Case Cohort Estimation", 
                                         t(beta.hat.unstrat.est), 
                                 Lambda0.Tau1Tau2.hat.unstrat.est, 
                                 Pi.x1.Tau1Tau2.hat.unstrat.est,
                                 Pi.x2.Tau1Tau2.hat.unstrat.est,
                                 Pi.x3.Tau1Tau2.hat.unstrat.est, 
                                 t(coxph.var.beta.unstrat.est), 
                              t(var.beta.unstrat.est), 
                                 var.Lambda0.Tau1Tau2.unstrat.est, 
                                 var.Pi.x1.Tau1Tau2.unstrat.est,
                                 var.Pi.x2.Tau1Tau2.unstrat.est,
                                 var.Pi.x3.Tau1Tau2.unstrat.est, 
                              t(CI.beta.unstrat.est), 
                                 t(CI.Lambda0.Tau1Tau2.unstrat.est), 
                                 t(CI.Pi.x1.Tau1Tau2.unstrat.est), 
                                 t(CI.Pi.x2.Tau1Tau2.unstrat.est),
                                 t(CI.Pi.x3.Tau1Tau2.unstrat.est))

    recap <- rbind(res.cohort, res.casecohort.robust, res.casecohort, 
                   res.casecohort.robust.est,
                   res.casecohort.est, 
                   res.casecohort.est.robust.naive, res.casecohort.est.naive,
                   res.unstrat.robust, res.unstrat,
                   res.unstrat.est.robust, res.unstrat.est,
                   res.unstrat.est.robust.naive, 
                   res.unstrat.est.naive)
    
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
                           n.cases.stratum.W3 = strata.n.cases[4], 
                           m = m,
                           Prob.y, beta1, beta2, beta3, Lambda0.Tau1Tau2, 
                           Pi.x1.Tau1Tau2, Pi.x2.Tau1Tau2, Pi.x3.Tau1Tau2,
                           strata.prob.missing1 = strata.prob.missingness[1],
                           strata.prob.missing2 = strata.prob.missingness[2]))
    row.names(res) <- NULL
     }
  myfile  <- paste0("SimulationResults_MissingData-n", n, "-Prob.y", Prob.y, 
                    "-K", K, 
                    "strata.prob.missingness", 
                    paste(strata.prob.missingness, collapse = "_"),
                    "Part", part, ".Rdata")
  save(res, file = myfile)
}

# ------------------------------------------------------------------------------
# Running the simulations in parallel for all the scenarios --------------------

resultat <- mclapply(1:80, Onerun, mc.cores = 32)

# ------------------------------------------------------------------------------
# Combining the simulation results ---------------------------------------------

param <- param[1:80, ]
RES <- NULL
for (p in 1: nrow(param)) {
  
  Prob.y                  <- param[p, 1] 
  n                       <- param[p, 2] 
  K       <- param[p, 3] 
  strata.prob.missingness <- param[p, 4:5]
  part                    <- param[p, 6]
  
  load(paste0("SimulationResults_MissingData-n", n, "-Prob.y", Prob.y, 
              "-K", K, 
              "strata.prob.missingness", 
              paste(strata.prob.missingness, collapse = "_"),
              "Part", part, 
              ".Rdata"))
  RES <- rbind(RES, res)
}

RECAP <- as.data.frame(RES)
ColNames <- colnames(RECAP[,c(2:57)])
RECAP[ColNames] <- sapply(RECAP[ColNames], as.numeric)
RECAP$Method <- as.factor(RECAP$Method)
RECAP$n <- as.factor(RECAP$n)
RECAP$K <- as.factor(RECAP$K)
RECAP$Prob.y <- as.factor(RECAP$Prob.y)

myfile  <- paste0("SimulationResults_MissingData.Rdata")
save(RECAP, file = myfile)

param <- scenarios.missingdata(n = c(5000, 10000), prob.y = c(0.02, 0.05, 0.1), 
                               noncases.per.case = c(2, 4), 
                               strata.prob.missingness = rbind(c(0.9, 0.8),
                                                               c(0.98, 0.9)))

param <- param[-c(21:24), ]

# ------------------------------------------------------------------------------
# Details of the results -------------------------------------------------------

details.beta1   <- NULL
details.beta2   <- NULL
details.beta3   <- NULL
details.Lambda0 <- NULL
details.Pi.x1   <- NULL
details.Pi.x2   <- NULL
details.Pi.x3   <- NULL
relative.efficiency <- NULL
Nreplic <- 5000
for (i in 1:nrow(param)) {
  
  RECAP1 <- RECAP[((i-1) * (13 * Nreplic) + 1):(i * (13 * Nreplic)), ]
  
  Lambda0 <- RECAP1[1,]$Lambda0.Tau1Tau2
  Pi.x1   <- RECAP1[1,]$Pi.x1.Tau1Tau2
  Pi.x2   <- RECAP1[1,]$Pi.x2.Tau1Tau2
  Pi.x3   <- RECAP1[1,]$Pi.x3.Tau1Tau2
  
    
  list.methods <- c("WholeCohort", "Case Cohort Robust", 
                    "Case Cohort", 
                    "Case Cohort Robust Estimation", 
                    "Case Cohort Estimation", 
                    "Case Cohort Robust Estimation Naive",
                    "Case Cohort Estimation Naive",             
                    "Unstrat Case Cohort Robust",
                    "Unstrat Case Cohort",
                    "Unstrat Case Cohort Robust Estimation",      
                   "Unstrat Case Cohort Estimation",            
                    "Unstrat Case Cohort Robust Estimation Naive", 
                   "Unstrat Case Cohort Estimation Naive")
    
  list.methods.col <- c("Cohort", "SCC.Robust.True", 
                        "SCC.True", 
                        "SCC.Robust.Est", 
                        "SCC.Est", 
                        "SCC.Robust.Naive",
                        "SCC.Naive",             
                        "USCC.Robust.True",
                        "USCC.True",
                        "USCC.Robust.Est",      
                        "USCC.Est",            
                        "USCC.Robust.Naive", 
                        "USCC.Naive")
  
  nameColResults <- function (name) {
    a <- NULL
    for (j in 1:length(list.methods)) {
      a <- c(a, paste0(name, ".", list.methods.col[j]))
    }
    return(a)
  }
  
  # Coverage of the confidence intervals 
  coverage.beta1 <- NULL
  coverage.beta2 <- NULL
  coverage.beta3 <- NULL
  coverage.Lambda0 <- NULL
  coverage.Pi.x1 <- NULL
  coverage.Pi.x2 <- NULL
  coverage.Pi.x3 <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    coverage.beta1 <- cbind(coverage.beta1,
                            mean((RECAP1$CI.left.beta1.hat[which(RECAP1$Method == method)] < beta1) & 
                                   (RECAP1$CI.right.beta1.hat[which(RECAP1$Method == method)] > beta1)) )
    coverage.beta2 <- cbind(coverage.beta2,
                            mean((RECAP1$CI.left.beta2.hat[which(RECAP1$Method == method)] < beta2) &
                                   (RECAP1$CI.right.beta2.hat[which(RECAP1$Method == method)] > beta2)) )
    coverage.beta3 <- cbind(coverage.beta3,
                            mean((RECAP1$CI.left.beta3.hat[which(RECAP1$Method == method)] < beta3) &
                                   (RECAP1$CI.right.beta3.hat[which(RECAP1$Method == method)] > beta3)) )
    coverage.Lambda0 <- cbind(coverage.Lambda0, mean((RECAP1$CI.left.Lambda0.Tau1Tau2.hat[which(RECAP1$Method == method)] < Lambda0) & (RECAP1$CI.right.Lambda0.Tau1Tau2.hat[which(RECAP1$Method == method)] > Lambda0)) )
    coverage.Pi.x1 <- cbind(coverage.Pi.x1,
                            mean((RECAP1$CI.left.Pi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)] < Pi.x1) &
                                   (RECAP1$CI.right.Pi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)] > Pi.x1)) )
    coverage.Pi.x2 <- cbind(coverage.Pi.x2,
                            mean((RECAP1$CI.left.Pi.x2.Tau1Tau2.hat[which(RECAP1$Method == method)] < Pi.x2) &
                                   (RECAP1$CI.right.Pi.x2.Tau1Tau2.hat[which(RECAP1$Method == method)] > Pi.x2)) )
    coverage.Pi.x3 <- cbind(coverage.Pi.x3,
                            mean((RECAP1$CI.left.Pi.x3.Tau1Tau2.hat[which(RECAP1$Method == method)] < Pi.x3) &
                                   (RECAP1$CI.right.Pi.x3.Tau1Tau2.hat[which(RECAP1$Method == method)] > Pi.x3)) )
      }
  coverage.beta1 <- as.data.frame(coverage.beta1)
  colnames(coverage.beta1) <- nameColResults("coverage.beta1")
  coverage.beta2 <- as.data.frame(coverage.beta2)
  colnames(coverage.beta2) <- nameColResults("coverage.beta2")
  coverage.beta3 <- as.data.frame(coverage.beta3)
  colnames(coverage.beta3) <- nameColResults("coverage.beta3")
  coverage.Lambda0 <- as.data.frame(coverage.Lambda0)
  colnames(coverage.Lambda0) <- nameColResults("coverage.Lambda0")
  coverage.Pi.x1 <- as.data.frame(coverage.Pi.x1)
  colnames(coverage.Pi.x1) <- nameColResults("coverage.Pi.x1")
  coverage.Pi.x2 <- as.data.frame(coverage.Pi.x2)
  colnames(coverage.Pi.x2) <- nameColResults("coverage.Pi.x2")
  coverage.Pi.x3 <- as.data.frame(coverage.Pi.x3)
  colnames(coverage.Pi.x3) <- nameColResults("coverage.Pi.x3")
  
  # Empirical variance
  empir.var.beta1 <- NULL
  empir.var.beta2 <- NULL
  empir.var.beta3 <- NULL
  empir.var.Lambda0 <- NULL
  empir.var.Pi.x1 <- NULL
  empir.var.Pi.x2 <- NULL
  empir.var.Pi.x3 <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    empir.var.beta1 <- cbind(empir.var.beta1, 
                             var(RECAP1$beta1.hat[which(RECAP1$Method == method)]))
    empir.var.beta2 <- cbind(empir.var.beta2, 
                             var(RECAP1$beta2.hat[which(RECAP1$Method == method)]))
    empir.var.beta3 <- cbind(empir.var.beta3, 
                             var(RECAP1$beta3.hat[which(RECAP1$Method == method)]))
    empir.var.Lambda0 <- cbind(empir.var.Lambda0, 
                               var(RECAP1$Lambda0.Tau1Tau2.hat[which(RECAP1$Method == method)]))
    empir.var.Pi.x1 <- cbind(empir.var.Pi.x1, 
                             var(RECAP1$Pi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)]))
    empir.var.Pi.x2 <- cbind(empir.var.Pi.x2, 
                             var(RECAP1$Pi.x2.Tau1Tau2.hat[which(RECAP1$Method == method)]))
    empir.var.Pi.x3 <- cbind(empir.var.Pi.x3, 
                             var(RECAP1$Pi.x3.Tau1Tau2.hat[which(RECAP1$Method == method)]))
  }
  empir.var.beta1 <- as.data.frame(empir.var.beta1)
  colnames(empir.var.beta1) <- nameColResults("empir.var.beta1")
  empir.var.beta2 <- as.data.frame(empir.var.beta2)
  colnames(empir.var.beta2) <- nameColResults("empir.var.beta2")
  empir.var.beta3 <- as.data.frame(empir.var.beta3)
  colnames(empir.var.beta3) <- nameColResults("empir.var.beta3")
  empir.var.Lambda0 <- as.data.frame(empir.var.Lambda0)
  colnames(empir.var.Lambda0) <- nameColResults("empir.var.Lambda0")
  empir.var.Pi.x1 <- as.data.frame(empir.var.Pi.x1)
  colnames(empir.var.Pi.x1) <- nameColResults("empir.var.Pi.x1")
  empir.var.Pi.x2 <- as.data.frame(empir.var.Pi.x2)
  colnames(empir.var.Pi.x2) <- nameColResults("empir.var.Pi.x2")
  empir.var.Pi.x3 <- as.data.frame(empir.var.Pi.x3)
  colnames(empir.var.Pi.x3) <- nameColResults("empir.var.Pi.x3")
  
  # Mean variance
  var.beta1 <- NULL
  var.beta2 <- NULL
  var.beta3 <- NULL
  var.Lambda0 <- NULL
  var.Pi.x1 <- NULL
  var.Pi.x2 <- NULL
  var.Pi.x3 <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    var.beta1 <- cbind(var.beta1, mean((RECAP1$var.beta1.hat[which(RECAP1$Method == method)])))
    var.beta2 <- cbind(var.beta2, mean((RECAP1$var.beta2.hat[which(RECAP1$Method == method)])))
    var.beta3 <- cbind(var.beta3, mean((RECAP1$var.beta3.hat[which(RECAP1$Method == method)])))
    var.Lambda0 <- cbind(var.Lambda0, mean((RECAP1$var.Lambda0.Tau1Tau2.hat[which(RECAP1$Method == method)])))
    var.Pi.x1 <- cbind(var.Pi.x1, mean((RECAP1$var.Pi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)])))
    var.Pi.x2 <- cbind(var.Pi.x2, mean((RECAP1$var.Pi.x2.Tau1Tau2.hat[which(RECAP1$Method == method)])))
    var.Pi.x3 <- cbind(var.Pi.x3, mean((RECAP1$var.Pi.x2.Tau1Tau2.hat[which(RECAP1$Method == method)])))
  }
  var.beta1 <- as.data.frame(var.beta1)
  colnames(var.beta1) <- nameColResults("var.beta1")
  var.beta2 <- as.data.frame(var.beta2)
  colnames(var.beta2) <- nameColResults("var.beta2")
  var.beta3 <- as.data.frame(var.beta3)
  colnames(var.beta3) <- nameColResults("var.beta3")
  var.Lambda0 <- as.data.frame(var.Lambda0)
  colnames(var.Lambda0) <- nameColResults("var.Lambda0")
  var.Pi.x1 <- as.data.frame(var.Pi.x1)
  colnames(var.Pi.x1) <- nameColResults("var.Pi.x1")
  var.Pi.x2 <- as.data.frame(var.Pi.x2)
  colnames(var.Pi.x2) <- nameColResults("var.Pi.x2")
  var.Pi.x3 <- as.data.frame(var.Pi.x3)
  colnames(var.Pi.x3) <- nameColResults("var.Pi.x3")
  
  coxph.var.beta1 <- NULL
  coxph.var.beta2 <- NULL
  coxph.var.beta3 <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    coxph.var.beta1 <- cbind(coxph.var.beta1, mean((RECAP1$coxph.var.beta1.hat[which(RECAP1$Method == method)])))
    coxph.var.beta2 <- cbind(coxph.var.beta2, mean((RECAP1$coxph.var.beta2.hat[which(RECAP1$Method == method)])))
    coxph.var.beta3 <- cbind(coxph.var.beta3, mean((RECAP1$coxph.var.beta3.hat[which(RECAP1$Method == method)])))
  }
  coxph.var.beta1 <- as.data.frame(coxph.var.beta1)
  colnames(coxph.var.beta1) <- nameColResults("coxph.var.beta1")
  coxph.var.beta2 <- as.data.frame(coxph.var.beta2)
  colnames(coxph.var.beta2) <- nameColResults("coxph.var.beta2")
  coxph.var.beta3 <- as.data.frame(coxph.var.beta3)
  colnames(coxph.var.beta3) <- nameColResults("coxph.var.beta3")
  
  
  # Bias/mean estimate
  mean.beta1 <- NULL
  mean.beta2 <- NULL
  mean.beta3 <- NULL
  mean.Lambda0 <- NULL
  mean.Pi.x1 <- NULL
  mean.Pi.x2 <- NULL
  mean.Pi.x3 <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    mean.beta1 <- cbind(mean.beta1, mean((RECAP1$beta1.hat[which(RECAP1$Method == method)])))
    mean.beta2 <- cbind(mean.beta2, mean((RECAP1$beta2.hat[which(RECAP1$Method == method)])))
    mean.beta3 <- cbind(mean.beta3, mean((RECAP1$beta3.hat[which(RECAP1$Method == method)])))
    mean.Lambda0 <- cbind(mean.Lambda0, mean((RECAP1$Lambda0.Tau1Tau2.hat[which(RECAP1$Method == method)])))
    mean.Pi.x1 <- cbind(mean.Pi.x1, mean((RECAP1$Pi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)])))
    mean.Pi.x2 <- cbind(mean.Pi.x2, mean((RECAP1$Pi.x2.Tau1Tau2.hat[which(RECAP1$Method == method)])))
    mean.Pi.x3 <- cbind(mean.Pi.x3, mean((RECAP1$Pi.x3.Tau1Tau2.hat[which(RECAP1$Method == method)])))
  }
  mean.beta1 <- as.data.frame(mean.beta1)
  colnames(mean.beta1) <- nameColResults("mean.beta1")
  mean.beta2 <- as.data.frame(mean.beta2)
  colnames(mean.beta2) <- nameColResults("mean.beta2")
  mean.beta3 <- as.data.frame(mean.beta3)
  colnames(mean.beta3) <- nameColResults("mean.beta3")
  mean.Lambda0 <- as.data.frame(mean.Lambda0)
  colnames(mean.Lambda0) <- nameColResults("mean.Lambda0")
  mean.Pi.x1 <- as.data.frame(mean.Pi.x1)
  colnames(mean.Pi.x1) <- nameColResults("mean.Pi.x1")
  mean.Pi.x2 <- as.data.frame(mean.Pi.x2)
  colnames(mean.Pi.x2) <- nameColResults("mean.Pi.x2")
  mean.Pi.x3 <- as.data.frame(mean.Pi.x3)
  colnames(mean.Pi.x3) <- nameColResults("mean.Pi.x3")
  
  bias.beta1 <- mean.beta1 - beta1
  colnames(bias.beta1) <- nameColResults("bias.beta1")
  bias.beta2 <- mean.beta2 - beta2
  colnames(bias.beta2) <- nameColResults("bias.beta2")
  bias.beta3 <- mean.beta3 - beta3
  colnames(bias.beta3) <- nameColResults("bias.beta3")
  bias.Lambda0 <- mean.Lambda0 - Lambda0
  colnames(bias.Lambda0) <- nameColResults("bias.Lambda0")
  bias.Pi.x1 <- mean.Pi.x1 - Pi.x1
  colnames(bias.Pi.x1) <- nameColResults("bias.Pi.x1")
  bias.Pi.x2 <- mean.Pi.x2 - Pi.x2
  colnames(bias.Pi.x2) <- nameColResults("bias.Pi.x2")
  bias.Pi.x3 <- mean.Pi.x3 - Pi.x3
  colnames(bias.Pi.x3) <- nameColResults("bias.Pi.x3")
  
  # relative efficiency
  eff.beta1 <- NULL
  eff.beta2 <- NULL
  eff.beta3 <- NULL
  eff.Lambda0 <- NULL
  eff.Pi.x1 <- NULL
  eff.Pi.x2 <- NULL
  eff.Pi.x3 <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    eff.beta1 <- cbind(eff.beta1, as.numeric(empir.var.beta1[1] / 
                                               empir.var.beta1[j]))
    eff.beta2 <- cbind(eff.beta2, as.numeric(empir.var.beta2[1] / 
                                               empir.var.beta2[j]))
    eff.beta3 <- cbind(eff.beta3, as.numeric(empir.var.beta3[1] / 
                                               empir.var.beta3[j]))
    eff.Lambda0 <- cbind(eff.Lambda0, as.numeric(empir.var.Lambda0[1] / 
                                                   empir.var.Lambda0[j]))
    eff.Pi.x1 <- cbind(eff.Pi.x1, as.numeric(empir.var.Pi.x1[1] / 
                                               empir.var.Pi.x1[j]))
    eff.Pi.x2 <- cbind(eff.Pi.x2, as.numeric(empir.var.Pi.x2[1] / 
                                               empir.var.Pi.x2[j]))
    eff.Pi.x3 <- cbind(eff.Pi.x3, as.numeric(empir.var.Pi.x3[1] / 
                                               empir.var.Pi.x3[j]))
  }
  eff.beta1 <- as.data.frame(eff.beta1)
  colnames(eff.beta1) <- nameColResults("eff.beta1")
  eff.beta2 <- as.data.frame(eff.beta2)
  colnames(eff.beta2) <- nameColResults("eff.beta2")
  eff.beta3 <- as.data.frame(eff.beta3)
  colnames(eff.beta3) <- nameColResults("eff.beta3")
  eff.Lambda0 <- as.data.frame(eff.Lambda0)
  colnames(eff.Lambda0) <- nameColResults("eff.Lambda0")
  eff.Pi.x1 <- as.data.frame(eff.Pi.x1)
  colnames(eff.Pi.x1) <- nameColResults("eff.Pi.x1")
  eff.Pi.x2 <- as.data.frame(eff.Pi.x2)
  colnames(eff.Pi.x2) <- nameColResults("eff.Pi.x2")
  eff.Pi.x3 <- as.data.frame(eff.Pi.x3)
  colnames(eff.Pi.x3) <- nameColResults("eff.Pi.x3")

  details.beta1 <- rbind(details.beta1, c(bias.beta1, empir.var.beta1, 
                                          var.beta1, coxph.var.beta1, 
                                          coverage.beta1, eff.beta1, 
                                          beta1 = beta1, beta2 = beta2,  
                                          beta3 = beta3, Lambda0 = Lambda0, 
                                          Pi.x1 = Pi.x1, Pi.x2 = Pi.x2,
                                          Pi.x3 = Pi.x3,
                                          n = as.character(RECAP1[1,]$n), 
                                          K = as.character(RECAP1[1,]$K), 
                                          sampl.fracW0 = as.character(RECAP1[1,]$sampl.fracW0),
                                          sampl.fracW1 = as.character(RECAP1[1,]$sampl.fracW1),
                                          sampl.fracW2 = as.character(RECAP1[1,]$sampl.fracW2),
                                          sampl.fracW3 = as.character(RECAP1[1,]$sampl.fracW3), 
                                          Prob.y = as.character(RECAP1[1,]$Prob.y), 
                                          m = as.character(RECAP1[1,]$m),
                                          strata.prob.missing1 = RECAP1$strata.prob.missing1[1],
                                          strata.prob.missing2 = RECAP1$strata.prob.missing2[1]))
  
  details.beta2 <- rbind(details.beta2, c(bias.beta2, empir.var.beta2, 
                                          var.beta2, coxph.var.beta2, 
                                          coverage.beta2, eff.beta2,
                                         beta1 = beta1, beta2 = beta2,  
                                         beta3 = beta3, Lambda0 = Lambda0, 
                                         Pi.x1 = Pi.x1, Pi.x2 = Pi.x2,
                                         Pi.x3 = Pi.x3,
                                         n = as.character(RECAP1[1,]$n), 
                                         K = as.character(RECAP1[1,]$K), 
                                         sampl.fracW0 = as.character(RECAP1[1,]$sampl.fracW0), 
                                         sampl.fracW1 = as.character(RECAP1[1,]$sampl.fracW1), 
                                         sampl.fracW2 = as.character(RECAP1[1,]$sampl.fracW2), 
                                         sampl.fracW3 = as.character(RECAP1[1,]$sampl.fracW3), 
                                         Prob.y = as.character(RECAP1[1,]$Prob.y), 
                                         m = as.character(RECAP1[1,]$m),
                                         strata.prob.missing1 = RECAP1$strata.prob.missing1[1],
                                         strata.prob.missing2 = RECAP1$strata.prob.missing2[1]))
  
  details.beta3 <- rbind(details.beta3, c(bias.beta3, empir.var.beta3, 
                                          var.beta3, coxph.var.beta3, 
                                          coverage.beta3, eff.beta3,
                                         beta1 = beta1, beta2 = beta2,  
                                         beta3 = beta3, Lambda0 = Lambda0, 
                                         Pi.x1 = Pi.x1, Pi.x2 = Pi.x2,
                                         Pi.x3 = Pi.x3,
                                         n = as.character(RECAP1[1,]$n), 
                                         K = as.character(RECAP1[1,]$K), 
                                         sampl.fracW0 = as.character(RECAP1[1,]$sampl.fracW0), 
                                         sampl.fracW1 = as.character(RECAP1[1,]$sampl.fracW1), 
                                         sampl.fracW2 = as.character(RECAP1[1,]$sampl.fracW2), 
                                         sampl.fracW3 = as.character(RECAP1[1,]$sampl.fracW3), 
                                         Prob.y = as.character(RECAP1[1,]$Prob.y), 
                                         m = as.character(RECAP1[1,]$m),
                                         strata.prob.missing1 = RECAP1$strata.prob.missing1[1],
                                         strata.prob.missing2 = RECAP1$strata.prob.missing2[1]))
  
  details.Lambda0 <- rbind(details.Lambda0, c(bias.Lambda0, empir.var.Lambda0, 
                                              var.Lambda0, 
                                              coverage.Lambda0, eff.Lambda0,
                                             beta1 = beta1, beta2 = beta2,  
                                             beta3 = beta3, Lambda0 = Lambda0, 
                                             Pi.x1 = Pi.x1, Pi.x2 = Pi.x2,
                                             Pi.x3 = Pi.x3,
                                             n = as.character(RECAP1[1,]$n), 
                                             K = as.character(RECAP1[1,]$K), 
                                             sampl.fracW0 = as.character(RECAP1[1,]$sampl.fracW0), 
                                             sampl.fracW1 = as.character(RECAP1[1,]$sampl.fracW1), 
                                             sampl.fracW2 = as.character(RECAP1[1,]$sampl.fracW2), 
                                             sampl.fracW3 = as.character(RECAP1[1,]$sampl.fracW3), 
                                             Prob.y = as.character(RECAP1[1,]$Prob.y), 
                                             m = as.character(RECAP1[1,]$m),
                                             strata.prob.missing1 = RECAP1$strata.prob.missing1[1],
                                             strata.prob.missing2 = RECAP1$strata.prob.missing2[1]))
  
  details.Pi.x1 <- rbind(details.Pi.x1, c(bias.Pi.x1, empir.var.Pi.x1, 
                                          var.Pi.x1, coverage.Pi.x1, eff.Pi.x1,
                                          beta1 = beta1, beta2 = beta2,  
                                          beta3 = beta3,
                                          Lambda0 = Lambda0, Pi.x1 = Pi.x1, 
                                          Pi.x2 = Pi.x2, Pi.x3 = Pi.x3,  
                                          n = as.character(RECAP1[1,]$n), 
                                          K = as.character(RECAP1[1,]$K), 
                                          sampl.fracW0 = as.character(RECAP1[1,]$sampl.fracW0), 
                                          sampl.fracW1 = as.character(RECAP1[1,]$sampl.fracW1), 
                                          sampl.fracW2 = as.character(RECAP1[1,]$sampl.fracW2), 
                                          sampl.fracW3 = as.character(RECAP1[1,]$sampl.fracW3), 
                                          Prob.y = as.character(RECAP1[1,]$Prob.y), 
                                          m = as.character(RECAP1[1,]$m),
                                          strata.prob.missing1 = RECAP1$strata.prob.missing1[1],
                                          strata.prob.missing2 = RECAP1$strata.prob.missing2[1]))

  details.Pi.x2 <- rbind(details.Pi.x2, c(bias.Pi.x2, empir.var.Pi.x2, 
                                          var.Pi.x2,  coverage.Pi.x2, eff.Pi.x2,
                                          beta1 = beta1, beta2 = beta2,  
                                          beta3 = beta3, Lambda0 = Lambda0, 
                                          Pi.x1 = Pi.x1, Pi.x2 = Pi.x2,
                                          Pi.x3 = Pi.x3,  
                                          n = as.character(RECAP1[1,]$n), 
                                          K = as.character(RECAP1[1,]$K), 
                                          sampl.fracW0 = as.numeric(as.character(RECAP1[1,]$sampl.fracW0)), 
                                          sampl.fracW1 = as.character(RECAP1[1,]$sampl.fracW1), 
                                          sampl.fracW2 = as.character(RECAP1[1,]$sampl.fracW2), 
                                          sampl.fracW3 = as.character(RECAP1[1,]$sampl.fracW3), 
                                          Prob.y = as.character(RECAP1[1,]$Prob.y), 
                                          m = as.character(RECAP1[1,]$m),
                                          strata.prob.missing1 = RECAP1$strata.prob.missing1[1],
                                          strata.prob.missing2 = RECAP1$strata.prob.missing2[1]))

  details.Pi.x3 <- rbind(details.Pi.x3, c(bias.Pi.x3, empir.var.Pi.x3, 
                                          var.Pi.x3, coverage.Pi.x3, eff.Pi.x3,
                                          beta1 = beta1, beta2 = beta2, 
                                          beta3 = beta3, Lambda0 = Lambda0, 
                                          Pi.x1 = Pi.x1, Pi.x2 = Pi.x2,
                                          Pi.x3 = Pi.x3,  
                                          n = as.character(RECAP1[1,]$n), 
                                          K = as.character(RECAP1[1,]$K), 
                                          sampl.fracW0 = as.numeric(as.character(RECAP1[1,]$sampl.fracW0)), 
                                          sampl.fracW1 = as.character(RECAP1[1,]$sampl.fracW1), 
                                          sampl.fracW2 = as.character(RECAP1[1,]$sampl.fracW2), 
                                          sampl.fracW3 = as.character(RECAP1[1,]$sampl.fracW3), 
                                          Prob.y = as.character(RECAP1[1,]$Prob.y), 
                                          m = as.character(RECAP1[1,]$m),
                                          strata.prob.missing1 = RECAP1$strata.prob.missing1[1],
                                          strata.prob.missing2 = RECAP1$strata.prob.missing2[1]))
  
  relative.efficiency <- rbind(relative.efficiency, c(eff.beta1, eff.beta2, 
                                                      eff.beta3, eff.Lambda0, 
                                                      eff.Pi.x1, eff.Pi.x2,
                                                      eff.Pi.x3, beta1 = beta1, 
                                                      beta2 = beta2,
                                                      beta3 = beta3, 
                                                      Lambda0 = Lambda0, 
                                                      Pi.x1 = Pi.x1, 
                                                      Pi.x2 = Pi.x2,
                                                      Pi.x3 = Pi.x3,
                                                      n = as.character(RECAP1[1,]$n), 
                                                      K = as.character(RECAP1[1,]$K), 
                                                      sampl.fracW0 = as.numeric(as.character(RECAP1[1,]$sampl.fracW0)), 
                                                      sampl.fracW1 = as.character(RECAP1[1,]$sampl.fracW1), 
                                                      sampl.fracW2 = as.character(RECAP1[1,]$sampl.fracW2), 
                                                      sampl.fracW3 = as.character(RECAP1[1,]$sampl.fracW3), 
                                                   Prob.y = as.character(RECAP1[1,]$Prob.y), 
                                                   m = as.character(RECAP1[1,]$m),
                                                   strata.prob.missing1 = RECAP1$strata.prob.missing1[1],
                                                   strata.prob.missing2 = RECAP1$strata.prob.missing2[1]))
}

details.beta1 <- as.data.frame(details.beta1)
details.beta2 <- as.data.frame(details.beta2)
details.beta3 <- as.data.frame(details.beta3)
details.Lambda0 <- as.data.frame(details.Lambda0)
details.Pi.x1 <- as.data.frame(details.Pi.x1)
details.Pi.x2 <- as.data.frame(details.Pi.x2)
details.Pi.x3 <- as.data.frame(details.Pi.x3)
relative.efficiency <- as.data.frame(relative.efficiency)

col.names.beta1 <- colnames(details.beta1[,c(1:91)])
details.beta1[col.names.beta1] <- sapply(details.beta1[col.names.beta1], 
                                         as.numeric)
myfile  <- paste0("details.beta1_MissingData.Rdata")
save(details.beta1, file = myfile)
details.beta1[col.names.beta1] <- round(details.beta1[col.names.beta1], 
                                        digits = 4)
write.csv(details.beta1, file = "details.beta1_MissingData.csv")

col.names.beta2 <- colnames(details.beta2[,c(1:91)])
details.beta2[col.names.beta2] <- sapply(details.beta2[col.names.beta2], 
                                         as.numeric)
myfile  <- paste0("details.beta2_MissingData.Rdata")
save(details.beta2, file = myfile)
details.beta2[col.names.beta2] <- round(details.beta2[col.names.beta2], 
                                        digits = 4)
write.csv(details.beta2, file = "details.beta2_MissingData.csv")

col.names.beta3 <- colnames(details.beta3[,c(1:91)])
details.beta3[col.names.beta3] <- sapply(details.beta3[col.names.beta3], 
                                         as.numeric)
myfile  <- paste0("details.beta3_MissingData.Rdata")
save(details.beta3, file = myfile)
details.beta3[col.names.beta3] <- round(details.beta3[col.names.beta3], 
                                        digits = 4)
write.csv(details.beta3, file = "details.beta3_MissingData.csv")

col.names.Lambda0 <- colnames(details.Lambda0[,c(1:78)])
details.Lambda0[col.names.Lambda0] <- sapply(details.Lambda0[col.names.Lambda0], 
                                             as.numeric)
myfile  <- paste0("details.Lambda0_MissingData.Rdata")
save(details.Lambda0, file = myfile)
write.csv(details.Lambda0, file = "details.Lambda0_MissingData.csv")

col.names.Pi.x1 <- colnames(details.Pi.x1[,c(1:78)])
details.Pi.x1[col.names.Pi.x1] <- sapply(details.Pi.x1[col.names.Pi.x1], 
                                         as.numeric)
myfile  <- paste0("details.Pi.x1_MissingData.Rdata")
save(details.Pi.x1, file = myfile)
write.csv(details.Pi.x1, file = "details.Pi.x1_MissingData.csv")

col.names.Pi.x2 <- colnames(details.Pi.x2[,c(1:78)])
details.Pi.x2[col.names.Pi.x2] <- sapply(details.Pi.x2[col.names.Pi.x2], 
                                         as.numeric)
myfile  <- paste0("details.Pi.x2_MissingData.Rdata")
save(details.Pi.x2, file = myfile)
write.csv(details.Pi.x2, file = "details.Pi.x2_MissingData.csv")

col.names.Pi.x3 <- colnames(details.Pi.x3[,c(1:78)])
details.Pi.x3[col.names.Pi.x3] <- sapply(details.Pi.x3[col.names.Pi.x3], 
                                         as.numeric)
myfile  <- paste0("details.Pi.x3_MissingData.Rdata")
save(details.Pi.x3, file = myfile)
write.csv(details.Pi.x3, file = "details.Pi.x3_MissingData.csv")

col.names.Eff <- colnames(relative.efficiency[,c(1:104)])
relative.efficiency[col.names.Eff] <- sapply(relative.efficiency[col.names.Eff], 
                                             as.numeric)
myfile  <- paste0("relative.efficiency_MissingData.Rdata")
save(relative.efficiency, file = myfile)
relative.efficiency[col.names.Eff] <- round(relative.efficiency[col.names.Eff], 
                                            digits = 3)
write.csv(relative.efficiency, file <- "relative.efficiency_MissingData.csv")

# ------------------------------------------------------------------------------
# Details of the results after log-transformation of the cumulative baseline 
# hazard and covariate specific pure risks estimates ---------------------------

details.transformedLambda0 <- NULL
details.transformedPi.x1 <- NULL
details.transformedPi.x2 <- NULL
details.transformedPi.x3 <- NULL
Nreplic <- 5000
for (i in 1:nrow(param)) {
  
  RECAP1 <- RECAP[((i-1) * (13 * Nreplic) + 1):(i * (13 * Nreplic)),]
  
  Lambda0 <- RECAP1[1,]$Lambda0.Tau1Tau2
  Pi.x1   <- RECAP1[1,]$Pi.x1.Tau1Tau2
  Pi.x2   <- RECAP1[1,]$Pi.x2.Tau1Tau2
  Pi.x3   <- RECAP1[1,]$Pi.x3.Tau1Tau2
  
  logLambda0 <- log(Lambda0)
  RECAP1$logLambda0.Tau1Tau2.hat <- log(RECAP1$Lambda0.Tau1Tau2.hat)
  logPi.x1 <- log(Pi.x1)
  RECAP1$logPi.x1.Tau1Tau2.hat <- log(RECAP1$Pi.x1.Tau1Tau2.hat)
  logPi.x2 <- log(Pi.x2)
  RECAP1$logPi.x2.Tau1Tau2.hat <- log(RECAP1$Pi.x2.Tau1Tau2.hat)
  logPi.x3 <- log(Pi.x3)
  RECAP1$logPi.x3.Tau1Tau2.hat <- log(RECAP1$Pi.x3.Tau1Tau2.hat)
  
  list.methods <- c("WholeCohort", "Case Cohort Robust", 
                    "Case Cohort", 
                    "Case Cohort Robust Estimation", 
                    "Case Cohort Estimation", 
                    "Case Cohort Robust Estimation Naive",
                    "Case Cohort Estimation Naive",             
                    "Unstrat Case Cohort Robust",
                    "Unstrat Case Cohort",
                    "Unstrat Case Cohort Robust Estimation",      
                    "Unstrat Case Cohort Estimation",            
                    "Unstrat Case Cohort Robust Estimation Naive", 
                    "Unstrat Case Cohort Estimation Naive")
  
  list.methods.col <- c("Cohort", "SCC.Robust.True", 
                        "SCC.True", 
                        "SCC.Robust.Est", 
                        "SCC.Est", 
                        "SCC.Robust.Naive",
                        "SCC.Naive",             
                        "USCC.Robust.True",
                        "USCC.True",
                        "USCC.Robust.Est",      
                        "USCC.Est",            
                        "USCC.Robust.Naive", 
                        "USCC.Naive")
  
  nameColResults <- function (name) {
    a <- NULL
    for (j in 1:length(list.methods)) {
      a <- c(a, paste0(name, ".", list.methods.col[j]))
    }
    return(a)
  }
  
  # Empirical variance
  empir.var.Lambda0 <- NULL
  empir.var.Pi.x1 <- NULL
  empir.var.Pi.x2 <- NULL
  empir.var.Pi.x3 <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    empir.var.Lambda0 <- cbind(empir.var.Lambda0, var(RECAP1$Lambda0.Tau1Tau2.hat[which(RECAP1$Method == method)]))
    empir.var.Pi.x1 <- cbind(empir.var.Pi.x1, var(RECAP1$Pi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)]))
    empir.var.Pi.x2 <- cbind(empir.var.Pi.x2, var(RECAP1$Pi.x2.Tau1Tau2.hat[which(RECAP1$Method == method)]))
    empir.var.Pi.x3 <- cbind(empir.var.Pi.x3, var(RECAP1$Pi.x3.Tau1Tau2.hat[which(RECAP1$Method == method)]))
  }
  empir.var.Lambda0 <- as.data.frame(empir.var.Lambda0)
  colnames(empir.var.Lambda0) <- nameColResults("empir.var.Lambda0")
  empir.var.Pi.x1 <- as.data.frame(empir.var.Pi.x1)
  colnames(empir.var.Pi.x1) <- nameColResults("empir.var.Pi.x1")
  empir.var.Pi.x2 <- as.data.frame(empir.var.Pi.x2)
  colnames(empir.var.Pi.x2) <- nameColResults("empir.var.Pi.x2")
  empir.var.Pi.x3 <- as.data.frame(empir.var.Pi.x3)
  colnames(empir.var.Pi.x3) <- nameColResults("empir.var.Pi.x3")
  
  ## Empirical variance for the transformed estimates
  empir.var.logLambda0 <- NULL
  empir.var.logPi.x1 <- NULL
  empir.var.logPi.x2 <- NULL
  empir.var.logPi.x3 <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    empir.var.logLambda0 <- cbind(empir.var.logLambda0, var(RECAP1$logLambda0.Tau1Tau2.hat[which(RECAP1$Method == method)]))
    empir.var.logPi.x1 <- cbind(empir.var.logPi.x1, var(RECAP1$logPi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)]))
    empir.var.logPi.x2 <- cbind(empir.var.logPi.x2, var(RECAP1$logPi.x2.Tau1Tau2.hat[which(RECAP1$Method == method)]))
    empir.var.logPi.x3 <- cbind(empir.var.logPi.x3, var(RECAP1$logPi.x3.Tau1Tau2.hat[which(RECAP1$Method == method)]))
  }
  empir.var.logLambda0 <- as.data.frame(empir.var.logLambda0)
  colnames(empir.var.logLambda0) <- nameColResults("empir.var.logLambda0")
  empir.var.logPi.x1 <- as.data.frame(empir.var.logPi.x1)
  colnames(empir.var.logPi.x1) <- nameColResults("empir.var.logPi.x1")
  empir.var.logPi.x2 <- as.data.frame(empir.var.logPi.x2)
  colnames(empir.var.logPi.x2) <- nameColResults("empir.var.logPi.x2")
  empir.var.logPi.x3 <- as.data.frame(empir.var.logPi.x3)
  colnames(empir.var.logPi.x3) <- nameColResults("empir.var.logPi.x3")
  

  ## Coverage of the confidence intervals without transformation
  coverage.Lambda0 <- NULL
  coverage.Pi.x1 <- NULL
  coverage.Pi.x2 <- NULL
  coverage.Pi.x3 <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    coverage.Lambda0 <- cbind(coverage.Lambda0, mean((RECAP1$CI.left.Lambda0.Tau1Tau2.hat[which(RECAP1$Method == method)] < Lambda0)&(RECAP1$CI.right.Lambda0.Tau1Tau2.hat[which(RECAP1$Method == method)] > Lambda0)) )
    coverage.Pi.x1 <- cbind(coverage.Pi.x1, mean((RECAP1$CI.left.Pi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)] < Pi.x1)&(RECAP1$CI.right.Pi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)] > Pi.x1)) )
    coverage.Pi.x2 <- cbind(coverage.Pi.x2, mean((RECAP1$CI.left.Pi.x2.Tau1Tau2.hat[which(RECAP1$Method == method)] < Pi.x2)&(RECAP1$CI.right.Pi.x2.Tau1Tau2.hat[which(RECAP1$Method == method)] > Pi.x2)) )
    coverage.Pi.x3 <- cbind(coverage.Pi.x3, mean((RECAP1$CI.left.Pi.x3.Tau1Tau2.hat[which(RECAP1$Method == method)] < Pi.x3)&(RECAP1$CI.right.Pi.x3.Tau1Tau2.hat[which(RECAP1$Method == method)] > Pi.x3)) )
  }
  coverage.Lambda0 <- as.data.frame(coverage.Lambda0)
  colnames(coverage.Lambda0) <- nameColResults("coverage.Lambda0")
  coverage.Pi.x1 <- as.data.frame(coverage.Pi.x1)
  colnames(coverage.Pi.x1) <- nameColResults("coverage.Pi.x1")
  coverage.Pi.x2 <- as.data.frame(coverage.Pi.x2)
  colnames(coverage.Pi.x2) <- nameColResults("coverage.Pi.x2")
  coverage.Pi.x3 <- as.data.frame(coverage.Pi.x3)
  colnames(coverage.Pi.x3) <- nameColResults("coverage.Pi.x3")
  
  ## Confidences intervals for the transformed estimates
  
  ## se estimation for the transformed estimates
  se.est.logLambda0 <- NULL
  se.est.logPi.x1 <- NULL
  se.est.logPi.x2 <- NULL
  se.est.logPi.x3 <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    se.est.logLambda0 <- cbind(se.est.logLambda0, sqrt((RECAP1$var.Lambda0.Tau1Tau2.hat[which(RECAP1$Method == method)]) / (RECAP1$Lambda0.Tau1Tau2.hat[which(RECAP1$Method == method)])^2))
    
    se.est.logPi.x1 <- cbind(se.est.logPi.x1, sqrt((RECAP1$var.Pi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)]) / (RECAP1$Pi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)])^2))
    se.est.logPi.x2 <- cbind(se.est.logPi.x2, sqrt((RECAP1$var.Pi.x2.Tau1Tau2.hat[which(RECAP1$Method == method)]) / (RECAP1$Pi.x2.Tau1Tau2.hat[which(RECAP1$Method == method)])^2))
    se.est.logPi.x3 <- cbind(se.est.logPi.x3, sqrt((RECAP1$var.Pi.x3.Tau1Tau2.hat[which(RECAP1$Method == method)]) / (RECAP1$Pi.x3.Tau1Tau2.hat[which(RECAP1$Method == method)])^2))
  }
  se.est.logLambda0 <- as.data.frame(se.est.logLambda0)
  colnames(se.est.logLambda0) <- nameColResults("se.est.logLambda0")
  se.est.logPi.x1 <- as.data.frame(se.est.logPi.x1)
  colnames(se.est.logPi.x1) <- nameColResults("se.est.logPi.x1")
  se.est.logPi.x2 <- as.data.frame(se.est.logPi.x2)
  colnames(se.est.logPi.x2) <- nameColResults("se.est.logPi.x2")
  se.est.logPi.x3 <- as.data.frame(se.est.logPi.x3)
  colnames(se.est.logPi.x3) <- nameColResults("se.est.logPi.x3")
  
  ## variance estimation for the transformed estimates
  var.est.logLambda0 <- NULL
  var.est.logPi.x1 <- NULL
  var.est.logPi.x2 <- NULL
  var.est.logPi.x3 <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    var.est.logLambda0 <- cbind(var.est.logLambda0, ((RECAP1$var.Lambda0.Tau1Tau2.hat[which(RECAP1$Method == method)]) / (RECAP1$Lambda0.Tau1Tau2.hat[which(RECAP1$Method == method)])^2))
    
    var.est.logPi.x1 <- cbind(var.est.logPi.x1, ((RECAP1$var.Pi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)]) / (RECAP1$Pi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)])^2))
    var.est.logPi.x2 <- cbind(var.est.logPi.x2, ((RECAP1$var.Pi.x2.Tau1Tau2.hat[which(RECAP1$Method == method)]) / (RECAP1$Pi.x2.Tau1Tau2.hat[which(RECAP1$Method == method)])^2))
    var.est.logPi.x3 <- cbind(var.est.logPi.x3, ((RECAP1$var.Pi.x3.Tau1Tau2.hat[which(RECAP1$Method == method)]) / (RECAP1$Pi.x3.Tau1Tau2.hat[which(RECAP1$Method == method)])^2))
  }
  var.est.logLambda0 <- as.data.frame(var.est.logLambda0)
  colnames(var.est.logLambda0) <- nameColResults("var.est.logLambda0")
  var.est.logPi.x1 <- as.data.frame(var.est.logPi.x1)
  colnames(var.est.logPi.x1) <- nameColResults("var.est.logPi.x1")
  var.est.logPi.x2 <- as.data.frame(var.est.logPi.x2)
  colnames(var.est.logPi.x2) <- nameColResults("var.est.logPi.x2")
  var.est.logPi.x3 <- as.data.frame(var.est.logPi.x3)
  colnames(var.est.logPi.x3) <- nameColResults("var.est.logPi.x3")
  
  # Construction of the CIs
  RECAP1$CI.left.logLambda0.Tau1Tau2.hat <-  RECAP1$CI.left.Lambda0.Tau1Tau2.hat
  RECAP1$CI.right.logLambda0.Tau1Tau2.hat <- RECAP1$CI.right.Lambda0.Tau1Tau2.hat
  RECAP1$CI.left.logPi.x1.Tau1Tau2.hat <-  RECAP1$CI.left.Pi.x1.Tau1Tau2.hat
  RECAP1$CI.right.logPi.x1.Tau1Tau2.hat <- RECAP1$CI.right.Pi.x1.Tau1Tau2.hat
  RECAP1$CI.left.logPi.x2.Tau1Tau2.hat <-  RECAP1$CI.left.Pi.x2.Tau1Tau2.hat
  RECAP1$CI.right.logPi.x2.Tau1Tau2.hat <- RECAP1$CI.right.Pi.x2.Tau1Tau2.hat
  RECAP1$CI.left.logPi.x3.Tau1Tau2.hat <-  RECAP1$CI.left.Pi.x3.Tau1Tau2.hat
  RECAP1$CI.right.logPi.x3.Tau1Tau2.hat <- RECAP1$CI.right.Pi.x3.Tau1Tau2.hat    
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    RECAP1$CI.left.logLambda0.Tau1Tau2.hat[which(RECAP1$Method == method)] <-    RECAP1$logLambda0[which(RECAP1$Method == method)] - qnorm(0.975) * (se.est.logLambda0[,j])
    RECAP1$CI.right.logLambda0.Tau1Tau2.hat[which(RECAP1$Method == method)] <-    RECAP1$logLambda0[which(RECAP1$Method == method)] + qnorm(0.975) * (se.est.logLambda0[,j])
    RECAP1$CI.left.logPi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)] <-    RECAP1$logPi.x1[which(RECAP1$Method == method)] - qnorm(0.975) * (se.est.logPi.x1[,j])
    RECAP1$CI.right.logPi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)] <-    RECAP1$logPi.x1[which(RECAP1$Method == method)] + qnorm(0.975) * (se.est.logPi.x1[,j])
    RECAP1$CI.left.logPi.x2.Tau1Tau2.hat[which(RECAP1$Method == method)] <-    RECAP1$logPi.x2[which(RECAP1$Method == method)] - qnorm(0.975) * (se.est.logPi.x2[,j])
    RECAP1$CI.right.logPi.x2.Tau1Tau2.hat[which(RECAP1$Method == method)] <-    RECAP1$logPi.x2[which(RECAP1$Method == method)] + qnorm(0.975) * (se.est.logPi.x2[,j])
    RECAP1$CI.left.logPi.x3.Tau1Tau2.hat[which(RECAP1$Method == method)] <-    RECAP1$logPi.x3[which(RECAP1$Method == method)] - qnorm(0.975) * (se.est.logPi.x3[,j])
    RECAP1$CI.right.logPi.x3.Tau1Tau2.hat[which(RECAP1$Method == method)] <-    RECAP1$logPi.x3[which(RECAP1$Method == method)] + qnorm(0.975) * (se.est.logPi.x3[,j])
  }
  
  # coverage of these CIs
  coverage.logLambda0 <- NULL
  coverage.logPi.x1 <- NULL
  coverage.logPi.x2 <- NULL
  coverage.logPi.x3 <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    coverage.logLambda0 <- cbind(coverage.logLambda0, mean((RECAP1$CI.left.logLambda0.Tau1Tau2.hat[which(RECAP1$Method == method)] < logLambda0)&(RECAP1$CI.right.logLambda0.Tau1Tau2.hat[which(RECAP1$Method == method)] > logLambda0)))
    coverage.logPi.x1 <- cbind(coverage.logPi.x1, mean((RECAP1$CI.left.logPi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)] < logPi.x1)&(RECAP1$CI.right.logPi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)] > logPi.x1)))
    coverage.logPi.x2 <- cbind(coverage.logPi.x2, mean((RECAP1$CI.left.logPi.x2.Tau1Tau2.hat[which(RECAP1$Method == method)] < logPi.x2)&(RECAP1$CI.right.logPi.x2.Tau1Tau2.hat[which(RECAP1$Method == method)] > logPi.x2)))
    coverage.logPi.x3 <- cbind(coverage.logPi.x3, mean((RECAP1$CI.left.logPi.x3.Tau1Tau2.hat[which(RECAP1$Method == method)] < logPi.x3)&(RECAP1$CI.right.logPi.x3.Tau1Tau2.hat[which(RECAP1$Method == method)] > logPi.x3)))
     }
  coverage.logLambda0 <- as.data.frame(coverage.logLambda0)
  colnames(coverage.logLambda0) <- nameColResults("coverage.logLambda0")
  coverage.logPi.x1 <- as.data.frame(coverage.logPi.x1)
  colnames(coverage.logPi.x1) <- nameColResults("coverage.logPi.x1")
  coverage.logPi.x2 <- as.data.frame(coverage.logPi.x2)
  colnames(coverage.logPi.x2) <- nameColResults("coverage.logPi.x2")
  coverage.logPi.x3 <- as.data.frame(coverage.logPi.x3)
  colnames(coverage.logPi.x3) <- nameColResults("coverage.logPi.x3")
  
  # Mean variance of untransformed estimates
  var.Lambda0 <- NULL
  var.Pi.x1 <- NULL
  var.Pi.x2 <- NULL
  var.Pi.x3 <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    var.Lambda0 <- cbind(var.Lambda0, mean((RECAP1$var.Lambda0.Tau1Tau2.hat[which(RECAP1$Method == method)])))
    var.Pi.x1 <- cbind(var.Pi.x1, mean((RECAP1$var.Pi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)])))
    var.Pi.x2 <- cbind(var.Pi.x2, mean((RECAP1$var.Pi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)])))
    var.Pi.x3 <- cbind(var.Pi.x3, mean((RECAP1$var.Pi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)])))
  }
  var.Lambda0 <- as.data.frame(var.Lambda0)
  colnames(var.Lambda0) <- nameColResults("var.Lambda0")
  var.Pi.x1 <- as.data.frame(var.Pi.x1)
  colnames(var.Pi.x1) <- nameColResults("var.Pi.x1")
  var.Pi.x2 <- as.data.frame(var.Pi.x2)
  colnames(var.Pi.x2) <- nameColResults("var.Pi.x2")
  var.Pi.x3 <- as.data.frame(var.Pi.x3)
  colnames(var.Pi.x3) <- nameColResults("var.Pi.x3")
  
  # Mean variance of transformed estimates
  var.logLambda0 <- apply(var.est.logLambda0, 2, mean)
  var.logPi.x1 <- apply(var.est.logPi.x1, 2, mean)
  var.logPi.x2 <- apply(var.est.logPi.x2, 2, mean)
  var.logPi.x3 <- apply(var.est.logPi.x3, 2, mean)
  names(var.logLambda0) <- nameColResults("var.logLambda0")
  names(var.logPi.x1) <- nameColResults("var.logPi.x1")
  names(var.logPi.x2) <- nameColResults("var.logPi.x2")
  names(var.logPi.x3) <- nameColResults("var.logPi.x3")
  
  # Bias/mean untransformed estimate
  mean.Lambda0 <- NULL
  mean.Pi.x1 <- NULL
  mean.Pi.x2 <- NULL
  mean.Pi.x3 <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    mean.Lambda0 <- cbind(mean.Lambda0, mean((RECAP1$Lambda0.Tau1Tau2.hat[which(RECAP1$Method == method)])))
    mean.Pi.x1 <- cbind(mean.Pi.x1, mean((RECAP1$Pi.x1.Tau1Tau2.hat[which(RECAP1$Method == method)])))
    mean.Pi.x2 <- cbind(mean.Pi.x2, mean((RECAP1$Pi.x2.Tau1Tau2.hat[which(RECAP1$Method == method)])))
    mean.Pi.x3 <- cbind(mean.Pi.x3, mean((RECAP1$Pi.x3.Tau1Tau2.hat[which(RECAP1$Method == method)])))
  }
  mean.Lambda0 <- as.data.frame(mean.Lambda0)
  colnames(mean.Lambda0) <- nameColResults("mean.Lambda0")
  mean.Pi.x1 <- as.data.frame(mean.Pi.x1)
  colnames(mean.Pi.x1) <- nameColResults("mean.Pi.x1")
  mean.Pi.x2 <- as.data.frame(mean.Pi.x2)
  colnames(mean.Pi.x2) <- nameColResults("mean.Pi.x2")
  mean.Pi.x3 <- as.data.frame(mean.Pi.x3)
  colnames(mean.Pi.x3) <- nameColResults("mean.Pi.x3")
  
  bias.Lambda0 <- mean.Lambda0 - Lambda0
  colnames(bias.Lambda0) <- nameColResults("bias.Lambda0")
  bias.Pi.x1 <- mean.Pi.x1 - Pi.x1
  colnames(bias.Pi.x1) <- nameColResults("bias.Pi.x1")
  bias.Pi.x2 <- mean.Pi.x2 - Pi.x2
  colnames(bias.Pi.x2) <- nameColResults("bias.Pi.x2")
  bias.Pi.x3 <- mean.Pi.x3 - Pi.x3
  colnames(bias.Pi.x3) <- nameColResults("bias.Pi.x3")
  
  # Bias/mean transformed estimate
  mean.logLambda0 <- NULL
  mean.logPi.x1 <- NULL
  mean.logPi.x2 <- NULL
  mean.logPi.x3 <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    mean.logLambda0 <- cbind(mean.logLambda0, mean((RECAP1$logLambda0[which(RECAP1$Method == method)])))
    mean.logPi.x1 <- cbind(mean.logPi.x1, mean((RECAP1$logPi.x1[which(RECAP1$Method == method)])))
    mean.logPi.x2 <- cbind(mean.logPi.x2, mean((RECAP1$logPi.x2[which(RECAP1$Method == method)])))
    mean.logPi.x3 <- cbind(mean.logPi.x3, mean((RECAP1$logPi.x3[which(RECAP1$Method == method)])))
  }
  mean.logLambda0 <- as.data.frame(mean.logLambda0)
  colnames(mean.logLambda0) <- nameColResults("mean.logLambda0")
  mean.logPi.x1 <- as.data.frame(mean.logPi.x1)
  colnames(mean.logPi.x1) <- nameColResults("mean.logPi.x1")
  mean.logPi.x2 <- as.data.frame(mean.logPi.x2)
  colnames(mean.logPi.x2) <- nameColResults("mean.logPi.x2")
  mean.logPi.x3 <- as.data.frame(mean.logPi.x3)
  colnames(mean.logPi.x3) <- nameColResults("mean.logPi.x3")
  
  bias.logLambda0 <- mean.logLambda0 - logLambda0
  colnames(bias.logLambda0) <- nameColResults("bias.logLambda0")
  bias.logPi.x1 <- mean.logPi.x1 - logPi.x1
  colnames(bias.logPi.x1) <- nameColResults("bias.logPi.x1")
  bias.logPi.x2 <- mean.logPi.x2 - logPi.x2
  colnames(bias.logPi.x2) <- nameColResults("bias.logPi.x2")
  bias.logPi.x3 <- mean.logPi.x3 - logPi.x3
  colnames(bias.logPi.x3) <- nameColResults("bias.logPi.x3")
  
  # saving the results 
  
  details.transformedLambda0 <- rbind(details.transformedLambda0, cbind(empir.var.Lambda0, var.Lambda0, empir.var.logLambda0, 
                                                                       t(var.logLambda0),
                                                                       coverage.Lambda0, coverage.logLambda0, bias.Lambda0, bias.logLambda0, 
                                                                       beta1 = beta1, beta2 = beta2, beta3 = beta3, Lambda0 = Lambda0, Pi.x1 = Pi.x1, Pi.x2 = Pi.x2, Pi.x3 = Pi.x3, logLambda0 = logLambda0, logPi.x1 = logPi.x1,  logPi.x2 = logPi.x2, logPi.x3 = logPi.x3,
                                                                       n = as.character(RECAP1[1,]$n), K = as.character(RECAP1[1,]$K), Prob.y = as.character(RECAP1[1,]$Prob.y) ,
                                                                       strata.prob.missing1 = RECAP1$strata.prob.missing1[1],
                                                                       strata.prob.missing2 = RECAP1$strata.prob.missing1[2]))
  
  details.transformedPi.x1 <- rbind(details.transformedPi.x1, cbind(empir.var.Pi.x1, var.Pi.x1, empir.var.logPi.x1, 
                                                             t(var.logPi.x1),
                                                             coverage.Pi.x1, coverage.logPi.x1, bias.Pi.x1, bias.logPi.x1, 
                                                             beta1 = beta1, beta2 = beta2, beta3 = beta3, Lambda0 = Lambda0, Pi.x1 = Pi.x1, Pi.x2 = Pi.x2, Pi.x3 = Pi.x3, logLambda0 = logLambda0, logPi.x1 = logPi.x1,  logPi.x2 = logPi.x2, logPi.x3 = logPi.x3, 
                                                             n = as.character(RECAP1[1,]$n), K = as.character(RECAP1[1,]$K), Prob.y = as.character(RECAP1[1,]$Prob.y) ,
                                                             strata.prob.missing1 = RECAP1$strata.prob.missing1[1],
                                                             strata.prob.missing2 = RECAP1$strata.prob.missing1[2]))

  details.transformedPi.x2 <- rbind(details.transformedPi.x2, cbind(empir.var.Pi.x2, var.Pi.x2, empir.var.logPi.x2, 
                                                                    t(var.logPi.x2),
                                                                    coverage.Pi.x2, coverage.logPi.x2, bias.Pi.x2, bias.logPi.x2, 
                                                                    beta1 = beta1, beta2 = beta2, beta3 = beta3, Lambda0 = Lambda0, Pi.x1 = Pi.x1, Pi.x2 = Pi.x2, Pi.x3 = Pi.x3, logLambda0 = logLambda0, logPi.x1 = logPi.x1,  logPi.x2 = logPi.x2, logPi.x3 = logPi.x3,
                                                                    n = as.character(RECAP1[1,]$n), K = as.character(RECAP1[1,]$K), Prob.y = as.character(RECAP1[1,]$Prob.y) ,
                                                                    strata.prob.missing1 = RECAP1$strata.prob.missing1[1],
                                                                    strata.prob.missing2 = RECAP1$strata.prob.missing1[2]))

  details.transformedPi.x3 <- rbind(details.transformedPi.x3, cbind(empir.var.Pi.x3, var.Pi.x3, empir.var.logPi.x3, 
                                                                    t(var.logPi.x3),
                                                                    coverage.Pi.x3, coverage.logPi.x3, bias.Pi.x3, bias.logPi.x3, 
                                                                    beta1 = beta1, beta2 = beta2, beta3 = beta3, Lambda0 = Lambda0, Pi.x1 = Pi.x1, Pi.x2 = Pi.x2, Pi.x3 = Pi.x3, logLambda0 = logLambda0,logPi.x1 = logPi.x1,  logPi.x2 = logPi.x2, logPi.x3 = logPi.x3, 
                                                                    n = as.character(RECAP1[1,]$n), K = as.character(RECAP1[1,]$K), Prob.y = as.character(RECAP1[1,]$Prob.y) ,
                                                                    strata.prob.missing1 = RECAP1$strata.prob.missing1[1],
                                                                    strata.prob.missing2 = RECAP1$strata.prob.missing1[2]))
}

details.transformedLambda0 <- as.data.frame(details.transformedLambda0)
details.transformedPi.x1 <- as.data.frame(details.transformedPi.x1)
details.transformedPi.x2 <- as.data.frame(details.transformedPi.x2)
details.transformedPi.x3 <- as.data.frame(details.transformedPi.x3)

col.names.transformedLambda0 <- colnames(details.transformedLambda0[,c(1:120)])
details.transformedLambda0[col.names.transformedLambda0] <- sapply(details.transformedLambda0[col.names.transformedLambda0], as.numeric)
myfile  <- paste0("details.transformedLambda0_sd_MissingData.Rdata")
save(details.transformedLambda0, file = myfile)
details.transformedLambda0[col.names.transformedLambda0[-c(52:78)]] <- round(details.transformedLambda0[col.names.transformedLambda0[-c(52:78)]], digits = 4)
details.transformedLambda0[col.names.transformedLambda0[c(52:78)]] <- round(details.transformedLambda0[col.names.transformedLambda0[c(52:78)]], digits = 8)
write.csv(details.transformedLambda0, file = "details.transformedLambda0_sd_MissingData.csv")

col.names.transformedPi.x1 <- colnames(details.transformedPi.x1[,c(1:120)])
details.transformedPi.x1[col.names.transformedPi.x1] <- sapply(details.transformedPi.x1[col.names.transformedPi.x1], as.numeric)
myfile  <- paste0("details.transformedPi.x1_sd_MissingData.Rdata")
save(details.transformedPi.x1, file = myfile)
details.transformedPi.x1[col.names.transformedPi.x1[-c(52:78)]] <- round(details.transformedPi.x1[col.names.transformedPi.x1[-c(52:78)]], digits = 4)
details.transformedPi.x1[col.names.transformedPi.x1[c(52:78)]] <- round(details.transformedPi.x1[col.names.transformedPi.x1[c(52:78)]], digits = 8)
write.csv(details.transformedPi.x1, file = "details.transformedPi.x1_sd_MissingData.csv")

col.names.transformedPi.x2 <- colnames(details.transformedPi.x2[,c(1:120)])
details.transformedPi.x2[col.names.transformedPi.x2] <- sapply(details.transformedPi.x2[col.names.transformedPi.x2], as.numeric)
myfile  <- paste0("details.transformedPi.x2_sd_MissingData.Rdata")
save(details.transformedPi.x2, file = myfile)
details.transformedPi.x2[col.names.transformedPi.x2[-c(52:78)]] <- round(details.transformedPi.x2[col.names.transformedPi.x2[-c(52:78)]], digits = 4)
details.transformedPi.x2[col.names.transformedPi.x2[c(52:78)]] <- round(details.transformedPi.x2[col.names.transformedPi.x2[c(52:78)]], digits = 8)
write.csv(details.transformedPi.x2, file = "details.transformedPi.x2_sd_MissingData.csv")

col.names.transformedPi.x3 <- colnames(details.transformedPi.x3[,c(1:120)])
details.transformedPi.x3[col.names.transformedPi.x3] <- sapply(details.transformedPi.x3[col.names.transformedPi.x3], as.numeric)
myfile  <- paste0("details.transformedPi.x3_sd_MissingData.Rdata")
save(details.transformedPi.x3, file = myfile)
details.transformedPi.x3[col.names.transformedPi.x3[-c(52:78)]] <- round(details.transformedPi.x3[col.names.transformedPi.x3[-c(52:78)]], digits = 4)
details.transformedPi.x3[col.names.transformedPi.x3[c(52:78)]] <- round(details.transformedPi.x3[col.names.transformedPi.x3[c(52:78)]], digits = 8)
write.csv(details.transformedPi.x3, file = "details.transformedPi.x3_sd_MissingData.csv")



