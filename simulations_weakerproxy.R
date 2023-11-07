## -----------------------------------------------------------------------------
## Filename: simulation_weakerproxy.R
## -----------------------------------------------------------------------------
## Description: This script compares the estimation of relative hazard, 
##              cumulative baseline and pure risks when using the whole cohort, 
##              a stratified case-cohort with design weights, a stratified
##              case-cohort with weights calibrated as proposed in Breslow et
##              al. (2009) and in Shin and al. (2020), an un-stratified case 
##              cohort with design weights, and an un-stratified case-cohort
##              with weights calibrated as proposed in Breslow et al. (2009) and 
##              in Shin and al. (2020). When using the case-cohorts, we also 
##              compare the robust variance estimate and the estimate with two 
##              components, one for the variation due to sampling from the 
##              superpopulation, and one due to sampling from the cohort, as
##              proposed in Etievant and Gail (2022).
##
##              This script gives the simulation results displayed in Web
##              Appendix D.6 in Etievant and Gail (2022), i.e., the same 
##              simulations as in Section 7 and Web Appendix D.2, but with a
##              weaker proxy for X1 (correlation between X1 and its proxy of 
##              0.64 instead of 0.8)
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
    
    # Generation of the covariates ---------------------------------------------
    X1        <- rnorm(n)
    X1.cat    <- 1 * (X1 < -2) + 2 *((X1 >= -2)&(X1 < 1)) + 3 * (X1 > 1)
    X2        <- rep(NA, n)
    for (i in 1:length(x1.cat)) {
      X2[which(X1.cat == x1.cat[i])] <- sample(0:2, sum(X1.cat == x1.cat[i]), 
                                               replace = TRUE, prob = Prob.x2[i,])
    }
    mu3       <- alpha1 * X1 + alpha2 * X2
    X3        <- mapply(mu3, FUN = X.generation)
    
    # Generation of the strata -------------------------------------------------
    W <- 1 * ((X1 < 0)&(X2 < 2)) + 2 * ((X1 >= 0)&(X2 > 0)) +
      3 * ((X1 < 0)&(X2 == 2)) + 0 *((X1 >= 0)&(X2 == 0))
    
    # Generation of the proxy variables for the covariates ---------------------
    X1.proxy  <- X1 + rnorm(n, sd = 1.2) # instead of sd = 0.75 in the previous simulations
    X3.proxy  <- X3 + rnorm(n, sd = 0.75)
    eps       <- runif(n)
    X2.proxy  <- X2
    X2.proxy[eps < 0.1] <- 2*(X2[eps < 0.1] == 0) + 1*(X2[eps < 0.1] == 2) + 0*(X2[eps < 0.1] == 1)
    X2.proxy[eps > 0.9] <- 2*(X2[eps > 0.9] == 1) + 1*(X2[eps > 0.9] == 0) + 0*(X2[eps > 0.9] == 2)
    
    # Generation of the entry times, event times, and times to independent 
    # censoring ----------------------------------------------------------------
    entry.times     <- runif(n, min = 0, max = 5) 
    scale.outcome   <- 1 / (bh * exp(beta1 * X1 + beta2 * X2 + beta3 * X3))
    outcome.times   <- rweibull(n, shape = 1, scale = scale.outcome)
    U               <- runif(n)
    death.times     <- (-log(U) / (0.02 / 10))
    
    # Observed times and case indicators ---------------------------------------
    times           <- pmin(death.times, outcome.times, 10 - entry.times) 
    status          <- 1 * (times == outcome.times) 
    
    # Numbers of individuals, of cases, and of individuals to be sampled, in 
    # each strata --------------------------------------------------------------
    strata.n.cases  <- c(sum(W == 0 & status == 1), sum(W == 1 & status == 1), 
                         sum(W == 2 & status == 1), sum(W == 3 & status == 1))
    strata.n        <- c(sum(W == 0), sum(W == 1), sum(W == 2), sum(W == 3))
    
    cohort <- cbind(id = 1:n, X1, X2, X3, W,X1.proxy, X2.proxy, X3.proxy, times, 
                    status, n = n, m = m, strata.n = strata.n[W+1], 
                    strata.m = strata.m[W+1], 
                    strata.n.cases = strata.n.cases[W+1])
    cohort <- as.data.frame(cohort)

    # --------------------------------------------------------------------------
    # Sampling of the stratified case cohort -----------------------------------
  
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
    # Estimation using the stratified case cohort with calibrated weights ------
    
    # Imputation of the covariate values for the individuals who are not in the 
    # stratified case cohort, from the proxy variables -------------------------
    mod.pred.X2     <- multinom(X2 ~ W + X2.proxy + X1.proxy, 
                                data = casecohort, weights = weights)
    cohort$X2.pred  <-  predict(mod.pred.X2, 
                                newdata = as.data.frame(cohort %>% 
                                                          select(W, X2.proxy, 
                                                                 X1.proxy)))
    mod.pred.X3     <- lm(X3 ~ X3.proxy + X1.proxy, data = casecohort, 
                          weights = weights)
    cohort$X3.pred  <- predict(mod.pred.X3, 
                               newdata = as.data.frame(cohort %>% 
                                                         select(X3.proxy, 
                                                                X1.proxy)), 
                               type = "response")
    mod.pred.X1     <- lm(X1 ~ X1.proxy, data = casecohort, weights = weights)
    cohort$X1.pred  <- predict(mod.pred.X1, 
                               newdata = as.data.frame(cohort %>% 
                                                         select(X1.proxy)), 
                               type = "response")
    
    # Running the coxph model on the imputed cohort data -----------------------
    mod.imputedcohort <- coxph(Surv(times, status) ~ X1.pred + X2 + 
                                 X3.pred, data = cohort, robust = TRUE)
    
    # Building the auxiliary variables proposed by Breslow et al. (2019) -------
    A.Breslow <- auxiliary.construction(mod.imputedcohort, Tau1, Tau2) 
    
    # Building the auxiliary variables proposed by Shin et al. (2020) ----------
    time.on.study <- pmax(pmin(Tau2, cohort$Time) - Tau1, 0)
    A             <- as.matrix(cbind(1, A.Breslow$A.RH))
    total         <- colSums(A)
    indiv.phase2  <- rownames(casecohort)
    calibration1  <- calibration(A.phase2 = A[indiv.phase2,], 
                                 design.weights = casecohort$weights, 
                                 total = total, eta0 = rep(0, ncol(A)), 
                                 niter.max = 10^3, epsilon.stop = 10^(-10))
    casecohort$weights.calib1 <- c(calibration1$calibrated.weights)
    mod.calib1    <- coxph(Surv(times, status) ~ X1 + X2 + X3, 
                           data = casecohort, weight = weights.calib1, 
                           id = id, robust = TRUE)
    A.Shin <- c(time.on.study * exp(mod.calib1$coefficients %*% 
                                      t(cbind(cohort$X1.pred, cohort$X2, 
                                              cohort$X3.pred))))
    
    # Calibrating the design weights against the auxiliary variables proposed by
    # Breslow et al. (2009) and Shin et al. (2020) -----------------------------
    A             <- as.matrix(cbind(1, A.Breslow$A.RH, A.Shin))
    total         <- colSums(A)
    calibration2  <- calibration(A.phase2 = A[indiv.phase2,], 
                                 design.weights = casecohort$weights, 
                                 total = total, eta0 = rep(0, ncol(A)), 
                                 niter.max = 10^3, epsilon.stop = 10^(-10))
    casecohort$weights.calib <- c(calibration2$calibrated.weights)
    
    mod.calib             <- coxph(Surv(times, status) ~ X1 + X2 + X3, 
                                   data = casecohort, 
                                   weight = weights.calib, id = id, 
                                   robust = TRUE)
    
    # Parameters and influences estimation -------------------------------------
    estimation.calib            <- influences(mod.calib, A = A, 
                                              calibrated = TRUE, Tau1 = Tau1, 
                                              Tau2 = Tau2, x = x1)
    
    beta.hat.calib              <- estimation.calib$beta.hat
    Lambda0.Tau1Tau2.hat.calib  <- estimation.calib$Lambda0.Tau1Tau2.hat
    Pi.x1.Tau1Tau2.hat.calib    <- estimation.calib$Pi.x.Tau1Tau2.hat
    
    infl.beta.calib             <- estimation.calib$infl.beta
    infl.Lambda0.Tau1Tau2.calib <- estimation.calib$infl.Lambda0.Tau1Tau2
    infl.Pi.x1.Tau1Tau2.calib   <- estimation.calib$infl.Pi.x.Tau1Tau2
    
    infl2.beta.calib            <- estimation.calib$infl2.beta
    infl2.Lambda0.Tau1Tau2.calib<- estimation.calib$infl2.Lambda0.Tau1Tau2
    infl2.Pi.x1.Tau1Tau2.calib  <- estimation.calib$infl2.Pi.x.Tau1Tau2
    
    estimation.calib.x2 <- influences.PR(beta = beta.hat.calib,
                                         Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.calib, 
                                         x = x2, infl.beta = infl.beta.calib,
                                         infl.Lambda0.Tau1Tau2 = infl.Lambda0.Tau1Tau2.calib,
                                         calibrated = TRUE, 
                                         infl2.beta = infl2.beta.calib,
                                         infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.Tau1Tau2.calib)
    
    estimation.calib.x3 <- influences.PR(beta = beta.hat.calib,
                                         Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.calib, 
                                         x = x3, infl.beta = infl.beta.calib,
                                         infl.Lambda0.Tau1Tau2 = infl.Lambda0.Tau1Tau2.calib,
                                         calibrated = TRUE, 
                                         infl2.beta = infl2.beta.calib,
                                         infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.Tau1Tau2.calib)
    
    Pi.x2.Tau1Tau2.hat.calib    <- estimation.calib.x2$Pi.x.Tau1Tau2.hat
    Pi.x3.Tau1Tau2.hat.calib    <- estimation.calib.x3$Pi.x.Tau1Tau2.hat
    
    infl.Pi.x2.Tau1Tau2.calib   <- estimation.calib.x2$infl.Pi.x.Tau1Tau2
    infl.Pi.x3.Tau1Tau2.calib   <- estimation.calib.x3$infl.Pi.x.Tau1Tau2
    
    infl2.Pi.x2.Tau1Tau2.calib  <- estimation.calib.x2$infl2.Pi.x.Tau1Tau2
    infl2.Pi.x3.Tau1Tau2.calib  <- estimation.calib.x3$infl2.Pi.x.Tau1Tau2
    
    # Robust variance estimation and corresponding CIs -------------------------
    coxph.var.beta.calib            <- diag(mod.calib$var)
    robust.var.beta.calib           <- robustvariance(infl.beta.calib)
    robust.var.Lambda0.Tau1Tau2.calib <- robustvariance(infl.Lambda0.Tau1Tau2.calib)
    robust.var.Pi.x1.Tau1Tau2.calib  <- robustvariance(infl.Pi.x1.Tau1Tau2.calib)
    robust.var.Pi.x2.Tau1Tau2.calib  <- robustvariance(infl.Pi.x2.Tau1Tau2.calib)
    robust.var.Pi.x3.Tau1Tau2.calib  <- robustvariance(infl.Pi.x3.Tau1Tau2.calib)
    
    CI.beta.robust.calib <- conf.interval(beta.hat.calib,
                                          robust.var.beta.calib)
    CI.Lambda0.Tau1Tau2.robust.calib <- conf.interval(Lambda0.Tau1Tau2.hat.calib,
                                                      robust.var.Lambda0.Tau1Tau2.calib)
    CI.Pi.x1.Tau1Tau2.robust.calib <- conf.interval(Pi.x1.Tau1Tau2.hat.calib, 
                                                    robust.var.Pi.x1.Tau1Tau2.calib)
    CI.Pi.x2.Tau1Tau2.robust.calib <- conf.interval(Pi.x2.Tau1Tau2.hat.calib, 
                                                    robust.var.Pi.x2.Tau1Tau2.calib)
    CI.Pi.x3.Tau1Tau2.robust.calib <- conf.interval(Pi.x3.Tau1Tau2.hat.calib, 
                                                    robust.var.Pi.x3.Tau1Tau2.calib)
    
    res.surveyrobust.calib <- cbind("Case Cohort Calibrated Weights Robust", 
                                    t(beta.hat.calib), 
                                    Lambda0.Tau1Tau2.hat.calib, 
                                    Pi.x1.Tau1Tau2.hat.calib,
                                    Pi.x2.Tau1Tau2.hat.calib,
                                    Pi.x3.Tau1Tau2.hat.calib,
                                    t(coxph.var.beta.calib), 
                                    t(robust.var.beta.calib), 
                                    robust.var.Lambda0.Tau1Tau2.calib, 
                                    robust.var.Pi.x1.Tau1Tau2.calib, 
                                    robust.var.Pi.x2.Tau1Tau2.calib, 
                                    robust.var.Pi.x3.Tau1Tau2.calib, 
                                    t(CI.beta.robust.calib), 
                                    t(CI.Lambda0.Tau1Tau2.robust.calib), 
                                    t(CI.Pi.x1.Tau1Tau2.robust.calib),
                                    t(CI.Pi.x2.Tau1Tau2.robust.calib),
                                    t(CI.Pi.x3.Tau1Tau2.robust.calib))
    
    # Superpopulation and phase-two sampling variance estimation and 
    # corresponding CIs --------------------------------------------------------
    var.beta.calib          <- variance(n = n, casecohort = casecohort, 
                                        cohort = cohort, calibrated = TRUE, 
                                        stratified = TRUE, 
                                        infl = infl.beta.calib,
                                        infl2 = infl2.beta.calib)
    var.Lambda0.Tau1Tau2.calib <- variance(n = n, casecohort = casecohort, 
                                           cohort = cohort, calibrated = TRUE, 
                                           stratified = TRUE, 
                                           infl = infl.Lambda0.Tau1Tau2.calib,
                                           infl2 = infl2.Lambda0.Tau1Tau2.calib)
    var.Pi.x1.Tau1Tau2.calib <- variance(n = n, casecohort = casecohort, 
                                         cohort = cohort, calibrated = TRUE, 
                                         stratified = TRUE, 
                                         infl = infl.Pi.x1.Tau1Tau2.calib,
                                         infl2 = infl2.Pi.x1.Tau1Tau2.calib)
    var.Pi.x2.Tau1Tau2.calib <- variance(n = n, casecohort = casecohort, 
                                         cohort = cohort, calibrated = TRUE, 
                                         stratified = TRUE, 
                                         infl = infl.Pi.x2.Tau1Tau2.calib,
                                         infl2 = infl2.Pi.x2.Tau1Tau2.calib)
    var.Pi.x3.Tau1Tau2.calib <- variance(n = n, casecohort = casecohort, 
                                         cohort = cohort, calibrated = TRUE, 
                                         stratified = TRUE, 
                                         infl = infl.Pi.x3.Tau1Tau2.calib,
                                         infl2 = infl2.Pi.x3.Tau1Tau2.calib)
    
    CI.beta.calib <- conf.interval(beta.hat.calib,var.beta.calib)
    CI.Lambda0.Tau1Tau2.calib <- conf.interval(Lambda0.Tau1Tau2.hat.calib,
                                               var.Lambda0.Tau1Tau2.calib)
    CI.Pi.x1.Tau1Tau2.calib <- conf.interval(Pi.x1.Tau1Tau2.hat.calib,
                                             var.Pi.x1.Tau1Tau2.calib)
    CI.Pi.x2.Tau1Tau2.calib <- conf.interval(Pi.x2.Tau1Tau2.hat.calib,
                                             var.Pi.x2.Tau1Tau2.calib)
    CI.Pi.x3.Tau1Tau2.calib <- conf.interval(Pi.x3.Tau1Tau2.hat.calib,
                                             var.Pi.x3.Tau1Tau2.calib)
    
    res.survey.calib  <- cbind("Case Cohort Calibrated Weights", 
                               t(beta.hat.calib), Lambda0.Tau1Tau2.hat.calib, 
                               Pi.x1.Tau1Tau2.hat.calib,
                               Pi.x2.Tau1Tau2.hat.calib,
                               Pi.x3.Tau1Tau2.hat.calib, 
                               t(coxph.var.beta.calib), t(var.beta.calib), 
                               var.Lambda0.Tau1Tau2.calib, 
                               var.Pi.x1.Tau1Tau2.calib,
                               var.Pi.x2.Tau1Tau2.calib,
                               var.Pi.x3.Tau1Tau2.calib, t(CI.beta.calib), 
                               t(CI.Lambda0.Tau1Tau2.calib), 
                               t(CI.Pi.x1.Tau1Tau2.calib), 
                               t(CI.Pi.x2.Tau1Tau2.calib),
                               t(CI.Pi.x3.Tau1Tau2.calib))
    
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
    # Estimation using the un-stratified case cohort with calibrated weights ---
    
    # Imputation of the covariate values for the individuals who are not in the 
    # stratified case cohort, from the proxy variables -------------------------
    mod.pred.X2.unstrat     <- multinom(X2 ~ W + X2.proxy + X1.proxy, 
                                        data = unstrat.casecohort, weights = weights)
    cohort$X2.pred.unstrat  <-  predict(mod.pred.X2.unstrat, 
                                        newdata = as.data.frame(cohort %>% 
                                                                  select(W, X2.proxy, 
                                                                         X1.proxy)))
    mod.pred.X3.unstrat     <- lm(X3 ~ X3.proxy + X1.proxy, 
                                  data = unstrat.casecohort, weights = weights)
    cohort$X3.pred.unstrat  <- predict(mod.pred.X3.unstrat, 
                                       newdata = as.data.frame(cohort %>% 
                                                                 select(X3.proxy, 
                                                                        X1.proxy)), 
                                       type = "response")
    mod.pred.X1.unstrat     <- lm(X1 ~ X1.proxy, data = unstrat.casecohort, 
                                  weights = weights)
    cohort$X1.pred.unstrat  <- predict(mod.pred.X1.unstrat, 
                                       newdata = as.data.frame(cohort %>% 
                                                                 select(X1.proxy)), 
                                       type = "response")
    
    # Running the coxph model on the imputed cohort data -----------------------
    mod.imputedcohort.unstrat <- coxph(Surv(times, status) ~ X1.pred.unstrat + 
                                         X2 + X3.pred.unstrat, 
                                       data = cohort, robust = TRUE)
    
    # Building the auxiliary variables proposed by Breslow et al. (2019) -------
    A.Breslow <- auxiliary.construction(mod.imputedcohort.unstrat, Tau1, Tau2) 
    
    # Building the auxiliary variables proposed by Shin et al. (2020) ----------
    time.on.study <- pmax(pmin(Tau2, cohort$Time) - Tau1, 0)
    A             <- as.matrix(cbind(1, A.Breslow$A.RH))
    total         <- colSums(A)
    indiv.phase2  <- rownames(unstrat.casecohort)
    calibration1  <- calibration(A.phase2 = A[indiv.phase2,], 
                                 design.weights = unstrat.casecohort$weights, 
                                 total = total, eta0 = rep(0, ncol(A)), 
                                 niter.max = 10^3, epsilon.stop = 10^(-10))
    unstrat.casecohort$weights.calib1 <- c(calibration1$calibrated.weights)
    mod.calib1.unstrat    <- coxph(Surv(times, status) ~ X1 + X2 + X3, 
                                   data = unstrat.casecohort, weight = weights.calib1, 
                                   id = id, robust = TRUE)
    A.Shin <- c(time.on.study * exp(mod.calib1.unstrat$coefficients %*% 
                                      t(cbind(cohort$X1.pred.unstrat, 
                                              cohort$X2, 
                                              cohort$X3.pred.unstrat))))
    
    # Calibrating the design weights against the auxiliary variables proposed by
    # Breslow et al. (2009) and Shin et al. (2020) -----------------------------
    A             <- as.matrix(cbind(1, A.Breslow$A.RH, A.Shin))
    total         <- colSums(A)
    calibration2  <- calibration(A.phase2 = A[indiv.phase2,], 
                                 design.weights = unstrat.casecohort$weights, 
                                 total = total, eta0 = rep(0, ncol(A)), 
                                 niter.max = 10^3, epsilon.stop = 10^(-10))
    unstrat.casecohort$weights.calib <- c(calibration2$calibrated.weights)
    
    mod.calib.unstrat           <- coxph(Surv(times, status) ~ X1 + X2 + X3, 
                                         data = unstrat.casecohort, 
                                         weight = weights.calib, id = id, 
                                         robust = TRUE)
    
    # Parameters and influences estimation -------------------------------------
    estimation.calib.unstrat    <- influences(mod.calib.unstrat, A = A, 
                                              calibrated = TRUE, Tau1 = Tau1, 
                                              Tau2 = Tau2, x = x1)
    
    beta.hat.calib.unstrat              <- estimation.calib.unstrat$beta.hat
    Lambda0.Tau1Tau2.hat.calib.unstrat  <- estimation.calib.unstrat$Lambda0.Tau1Tau2.hat
    Pi.x1.Tau1Tau2.hat.calib.unstrat    <- estimation.calib.unstrat$Pi.x.Tau1Tau2.hat
    
    infl.beta.calib.unstrat             <- estimation.calib.unstrat$infl.beta
    infl.Lambda0.Tau1Tau2.calib.unstrat <- estimation.calib.unstrat$infl.Lambda0.Tau1Tau2
    infl.Pi.x1.Tau1Tau2.calib.unstrat    <- estimation.calib.unstrat$infl.Pi.x.Tau1Tau2
    
    infl2.beta.calib.unstrat            <- estimation.calib.unstrat$infl2.beta
    infl2.Lambda0.Tau1Tau2.calib.unstrat <- estimation.calib.unstrat$infl2.Lambda0.Tau1Tau2
    infl2.Pi.x1.Tau1Tau2.calib.unstrat  <- estimation.calib.unstrat$infl2.Pi.x.Tau1Tau2
    
    estimation.calib.unstrat.x2 <- influences.PR(beta = beta.hat.calib.unstrat,
                                                 Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.calib.unstrat, 
                                                 x = x2, infl.beta = infl.beta.calib.unstrat,
                                                 infl.Lambda0.Tau1Tau2 = infl.Lambda0.Tau1Tau2.calib.unstrat,
                                                 calibrated = TRUE, 
                                                 infl2.beta = infl2.beta.calib.unstrat,
                                                 infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.Tau1Tau2.calib.unstrat)
    
    estimation.calib.unstrat.x3 <- influences.PR(beta = beta.hat.calib.unstrat,
                                                 Lambda0.Tau1Tau2 = Lambda0.Tau1Tau2.hat.calib.unstrat, 
                                                 x = x3, infl.beta = infl.beta.calib.unstrat,
                                                 infl.Lambda0.Tau1Tau2 = infl.Lambda0.Tau1Tau2.calib.unstrat,
                                                 calibrated = TRUE, 
                                                 infl2.beta = infl2.beta.calib.unstrat,
                                                 infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.Tau1Tau2.calib.unstrat)
    
    Pi.x2.Tau1Tau2.hat.calib.unstrat    <- estimation.calib.unstrat.x2$Pi.x.Tau1Tau2.hat
    Pi.x3.Tau1Tau2.hat.calib.unstrat    <- estimation.calib.unstrat.x3$Pi.x.Tau1Tau2.hat
    
    infl.Pi.x2.Tau1Tau2.calib.unstrat   <- estimation.calib.unstrat.x2$infl.Pi.x.Tau1Tau2
    infl.Pi.x3.Tau1Tau2.calib.unstrat   <- estimation.calib.unstrat.x3$infl.Pi.x.Tau1Tau2
    
    infl2.Pi.x2.Tau1Tau2.calib.unstrat  <- estimation.calib.unstrat.x2$infl2.Pi.x.Tau1Tau2
    infl2.Pi.x3.Tau1Tau2.calib.unstrat  <- estimation.calib.unstrat.x3$infl2.Pi.x.Tau1Tau2
    
    # Robust variance estimation and corresponding CIs -------------------------
    coxph.var.beta.calib.unstrat            <- diag(mod.calib.unstrat$var)
    robust.var.beta.calib.unstrat           <- robustvariance(infl.beta.calib.unstrat)
    robust.var.Lambda0.Tau1Tau2.calib.unstrat <- robustvariance(infl.Lambda0.Tau1Tau2.calib.unstrat)
    robust.var.Pi.x1.Tau1Tau2.calib.unstrat  <- robustvariance(infl.Pi.x1.Tau1Tau2.calib.unstrat)
    robust.var.Pi.x2.Tau1Tau2.calib.unstrat  <- robustvariance(infl.Pi.x2.Tau1Tau2.calib.unstrat)
    robust.var.Pi.x3.Tau1Tau2.calib.unstrat  <- robustvariance(infl.Pi.x3.Tau1Tau2.calib.unstrat)
    
    CI.beta.robust.calib.unstrat <- conf.interval(beta.hat.calib.unstrat,
                                                  robust.var.beta.calib.unstrat)
    CI.Lambda0.Tau1Tau2.robust.calib.unstrat <- conf.interval(Lambda0.Tau1Tau2.hat.calib.unstrat,
                                                              robust.var.Lambda0.Tau1Tau2.calib.unstrat)
    CI.Pi.x1.Tau1Tau2.robust.calib.unstrat <- conf.interval(Pi.x1.Tau1Tau2.hat.calib.unstrat, 
                                                            robust.var.Pi.x1.Tau1Tau2.calib.unstrat)
    CI.Pi.x2.Tau1Tau2.robust.calib.unstrat <- conf.interval(Pi.x2.Tau1Tau2.hat.calib.unstrat, 
                                                            robust.var.Pi.x2.Tau1Tau2.calib.unstrat)
    CI.Pi.x3.Tau1Tau2.robust.calib.unstrat <- conf.interval(Pi.x3.Tau1Tau2.hat.calib.unstrat, 
                                                            robust.var.Pi.x3.Tau1Tau2.calib.unstrat)
    
    res.unstratrobust.calib <- cbind("Unstrat Case Cohort Calibrated Weights Robust", 
                                     t(beta.hat.calib.unstrat), 
                                     Lambda0.Tau1Tau2.hat.calib.unstrat, 
                                     Pi.x1.Tau1Tau2.hat.calib.unstrat,
                                     Pi.x2.Tau1Tau2.hat.calib.unstrat,
                                     Pi.x3.Tau1Tau2.hat.calib.unstrat,
                                     t(coxph.var.beta.calib.unstrat), 
                                     t(robust.var.beta.calib.unstrat), 
                                     robust.var.Lambda0.Tau1Tau2.calib.unstrat, 
                                     robust.var.Pi.x1.Tau1Tau2.calib.unstrat, 
                                     robust.var.Pi.x2.Tau1Tau2.calib.unstrat, 
                                     robust.var.Pi.x3.Tau1Tau2.calib.unstrat, 
                                     t(CI.beta.robust.calib.unstrat), 
                                     t(CI.Lambda0.Tau1Tau2.robust.calib.unstrat), 
                                     t(CI.Pi.x1.Tau1Tau2.robust.calib.unstrat),
                                     t(CI.Pi.x2.Tau1Tau2.robust.calib.unstrat),
                                     t(CI.Pi.x3.Tau1Tau2.robust.calib.unstrat))
    
    # Superpopulation and phase-two sampling variance estimation and 
    # corresponding CIs --------------------------------------------------------
    var.beta.calib.unstrat <- variance(n = n, casecohort = unstrat.casecohort, 
                                       cohort = cohort, calibrated = TRUE, 
                                       stratified = FALSE, 
                                       infl = infl.beta.calib.unstrat,
                                       infl2 = infl2.beta.calib.unstrat)
    var.Lambda0.Tau1Tau2.calib.unstrat <- variance(n = n, 
                                                   casecohort = unstrat.casecohort, 
                                                   cohort = cohort, 
                                                   calibrated = TRUE, 
                                                   stratified = FALSE, 
                                                   infl = infl.Lambda0.Tau1Tau2.calib.unstrat,
                                                   infl2 = infl2.Lambda0.Tau1Tau2.calib.unstrat)
    var.Pi.x1.Tau1Tau2.calib.unstrat <- variance(n = n, 
                                                 casecohort = unstrat.casecohort, 
                                                 cohort = cohort, 
                                                 calibrated = TRUE, 
                                                 stratified = FALSE, 
                                                 infl = infl.Pi.x1.Tau1Tau2.calib.unstrat,
                                                 infl2 = infl2.Pi.x1.Tau1Tau2.calib.unstrat)
    var.Pi.x2.Tau1Tau2.calib.unstrat <- variance(n = n, 
                                                 casecohort = unstrat.casecohort, 
                                                 cohort = cohort, 
                                                 calibrated = TRUE, 
                                                 stratified = FALSE, 
                                                 infl = infl.Pi.x2.Tau1Tau2.calib.unstrat,
                                                 infl2 = infl2.Pi.x2.Tau1Tau2.calib.unstrat)
    var.Pi.x3.Tau1Tau2.calib.unstrat <- variance(n = n, 
                                                 casecohort = unstrat.casecohort, 
                                                 cohort = cohort, 
                                                 calibrated = TRUE, 
                                                 stratified = FALSE, 
                                                 infl = infl.Pi.x3.Tau1Tau2.calib.unstrat,
                                                 infl2 = infl2.Pi.x3.Tau1Tau2.calib.unstrat)
    
    CI.beta.calib.unstrat <- conf.interval(beta.hat.calib.unstrat,
                                           var.beta.calib.unstrat)
    CI.Lambda0.Tau1Tau2.calib.unstrat <- conf.interval(Lambda0.Tau1Tau2.hat.calib.unstrat,
                                                       var.Lambda0.Tau1Tau2.calib.unstrat)
    CI.Pi.x1.Tau1Tau2.calib.unstrat <- conf.interval(Pi.x1.Tau1Tau2.hat.calib.unstrat,
                                                     var.Pi.x1.Tau1Tau2.calib.unstrat)
    CI.Pi.x2.Tau1Tau2.calib.unstrat <- conf.interval(Pi.x2.Tau1Tau2.hat.calib.unstrat,
                                                     var.Pi.x2.Tau1Tau2.calib.unstrat)
    CI.Pi.x3.Tau1Tau2.calib.unstrat <- conf.interval(Pi.x3.Tau1Tau2.hat.calib.unstrat,
                                                     var.Pi.x3.Tau1Tau2.calib.unstrat)
    
    res.unstrat.calib  <- cbind("Unstrat Case Cohort Calibrated Weights", 
                                t(beta.hat.calib.unstrat), 
                                Lambda0.Tau1Tau2.hat.calib.unstrat, 
                                Pi.x1.Tau1Tau2.hat.calib.unstrat,
                                Pi.x2.Tau1Tau2.hat.calib.unstrat,
                                Pi.x3.Tau1Tau2.hat.calib.unstrat, 
                                t(coxph.var.beta.calib.unstrat), 
                                t(var.beta.calib.unstrat), 
                                var.Lambda0.Tau1Tau2.calib.unstrat, 
                                var.Pi.x1.Tau1Tau2.calib.unstrat,
                                var.Pi.x2.Tau1Tau2.calib.unstrat,
                                var.Pi.x3.Tau1Tau2.calib.unstrat, 
                                t(CI.beta.calib.unstrat), 
                                t(CI.Lambda0.Tau1Tau2.calib.unstrat), 
                                t(CI.Pi.x1.Tau1Tau2.calib.unstrat), 
                                t(CI.Pi.x2.Tau1Tau2.calib.unstrat),
                                t(CI.Pi.x3.Tau1Tau2.calib.unstrat))
    
    recap <- rbind(res.survey, res.surveyrobust, res.survey.calib, 
                   res.surveyrobust.calib, res.unstrat, res.unstratrobust, 
                   res.unstrat.calib, res.unstratrobust.calib, res.cohort)
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
  myfile  <- paste0("SimulationResults_WeakerProxy-n", n, "-Prob.y", Prob.y, 
                    "-K", K, "Part", part, 
                    ".Rdata")
  save(res, file = myfile)
}

# ------------------------------------------------------------------------------
# Running the simulations in parallel for all the scenarios --------------------

resultat <- mclapply(41:48, Onerun, mc.cores = 8)

# ------------------------------------------------------------------------------
# Combining the simulation results ---------------------------------------------

RES <- NULL
for (p in 1: nrow(param)) {
  
  Prob.y            <- param[p, 1] 
  n                 <- param[p, 2] 
  K                 <- param[p, 3] 
  part              <- param[p, 4]
  
  load(paste0("SimulationResults_WeakerProxy-n", n, "-Prob.y", Prob.y, "-K", 
              K, "Part", part, ".Rdata"))
  RES <- rbind(RES, res)
  
}
RECAP           <- as.data.frame(RES)
ColNames        <- colnames(RECAP[,c(2:55)])
RECAP[ColNames] <- sapply(RECAP[ColNames], as.numeric)
RECAP$Method    <- as.factor(RECAP$Method)
RECAP$n         <- as.factor(RECAP$n)
RECAP$K         <- as.factor(RECAP$K)
RECAP$Prob.y    <- as.factor(RECAP$Prob.y)
myfile  <- paste0("SimulationResults_WeakerProxy.Rdata")
save(RECAP, file = myfile)

# ------------------------------------------------------------------------------
# Details of the results  ------------------------------------------------------

param <- scenarios(n = c(5000, 10000), prob.y = c(0.02, 0.05, 0.1), 
                   noncases.per.case = c(2, 4))

list.methods <- c("WholeCohort", "Case Cohort Robust", 
                  "Case Cohort", 
                  "Case Cohort Calibrated Weights Robust",  
                  "Case Cohort Calibrated Weights",
                  "Unstrat Case Cohort Robust", 
                  "Unstrat Case Cohort", 
                  "Unstrat Case Cohort Calibrated Weights Robust", 
                  "Unstrat Case Cohort Calibrated Weights")

list.methods.col <- c("Cohort", "SCC.Robust", "SCC", "SCC.Robust.Cal",
                      "SCC.Cal", "USCC.Robust", "USCC", "USCC.Robust.Cal", 
                      "USCC.Cal")

source("detailed.results.R")

detailed.results(RECAP = RECAP, param = param, Nreplic = 5000, 
                 list.methods = list.methods, 
                 list.methods.col = list.methods.col, 
                 simul.name = "_WeakerProxy")
