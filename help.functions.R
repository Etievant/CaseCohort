## -----------------------------------------------------------------------------
## Function: estimation()
## -----------------------------------------------------------------------------
## Description: This function estimates the log-relative hazard, baseline 
##              hazards at each unique event time, cumulative baseline hazard in 
##              the time interval [Tau1, Tau2] and pure risk in [Tau1, Tau2] and
##              for covariate profile x
## -----------------------------------------------------------------------------
## Arguments:
##
##  mod             a cox model object, result of function coxph 
##
##  Tau1            left bound of the time interval considered for the 
##                  cumulative baseline hazard and pure risk. Default is the 
##                  first event time
##
##  Tau2            right bound of the time interval considered for the 
##                  cumulative baseline hazard and pure risk. Default is the 
##                  last event time
##
##  x               vector of length p, specifying the covariate profile 
##                  considered for the pure risk. Default is (0,...,0).
##
##  missing.data    is data on the p covariates missing for certain individuals
##                  in the stratified case cohort? If missing.data = TRUE, the
##                  arguments below need to be provided
##
##  riskmat.phase2  at risk matrix for the stratified case cohort (phase-two 
##                  data) at all of the cases event times, even that with 
##                  missing covariate data. If missing.data = TRUE, this 
##                  argument needs to be provided
##
##  dNt.phase2      counting process matrix for failures in the stratified
##                  case cohort (phase-two data). 
##
##  status.phase2   vector indicating the case status in the stratified case 
##                  cohort (phase-two data). If missing.data = TRUE, this 
##                  argument needs to be provided
## -----------------------------------------------------------------------------

estimation <- function (mod, Tau1 = NULL, Tau2 = NULL, x = NULL, 
                        missing.data = NULL, riskmat.phase2 = NULL, 
                        dNt.phase2 = NULL, status.phase2 = NULL) {
  
  if (is.null(missing.data)) {
    missing.data <- FALSE
    }
  
  # ----------------------------------------------------------------------------
  # If missing covariate data --------------------------------------------------
  
  if (missing.data == TRUE) {
    if (is.null(riskmat.phase2)) {
      stop("The at risk matrix for the phase-two data at all of the cases event 
           times need to be provided, even that with missing covariate data")
    } else {
      if (is.null(dNt.phase2)&is.null(status.phase2)) {
        stop("The status for the phase-two data, or the counting process matrix 
             for failures in the phase-two data, need to be provided")
      } else {
        if (is.null(dNt.phase2)) { 
          
          # --------------------------------------------------------------------
          # If not provided, compute the counting process matrix for failures in 
          # the phase-two data -------------------------------------------------
          
          observed.times.phase2 <- apply(riskmat.phase2, 1,
                                         function(v) {which.max(cumsum(v))})
          dNt.phase2            <- matrix(0, nrow(riskmat.phase2), 
                                          ncol(riskmat.phase2))
          dNt.phase2[cbind(1:nrow(riskmat.phase2), observed.times.phase2)] <- 1
          dNt.phase2            <- sweep(dNt.phase2, 1, status.phase2, "*") 
        }
        
        # ----------------------------------------------------------------------
        # Quantities needed for the estimation ---------------------------------
        
        mod.detail      <- coxph.detail(mod, riskmat = TRUE)
        riskmat         <- mod.detail$riskmat
        X               <- model.matrix(mod)
        indiv.phase3    <- row.names(X) 
        riskmat.phase3  <- riskmat.phase2[indiv.phase3, ] # use the cases' 
        # actual failure time, known even for cases with missing covariate data
        weights         <- mod$weights # weights used for the fit, omega_i
        if (is.null(weights)) {
          weights       <- 1
          } # when using the whole cohort, omega_i = 1
        
        # ----------------------------------------------------------------------
        # Estimation of the log-relative hazard --------------------------------
        
        beta.hat        <- mod$coefficients 
        
        exp.X.weighted  <- weights * exp(X %*% beta.hat)
        S0t             <- t(riskmat.phase3) %*% (exp.X.weighted) # for any
        # event time t (all of the cases failure times, even the one of cases 
        # who have missing covariate data)
        cases.times     <- colnames(riskmat.phase2) # cases failure times, even 
        # the one of cases who have missing covariate data
        if (is.null(Tau1)) {
          Tau1          <- floor(min(mod.detail$time)) 
        }
        if (is.null(Tau2)) {
          Tau2          <- floor(max(mod.detail$time))
        }
        Tau1Tau2.times  <- which((cases.times > Tau1 & cases.times <= Tau2)) 

        # ----------------------------------------------------------------------
        # Estimation of the baseline hazards at each unique event time and 
        # cumulative baseline hazard in time interval [Tau1, Tau2] -------------
        
        lambda0.t.hat   <- t(colSums(dNt.phase2) / S0t) # we do not use weights
        # in the numerator of the Breslow estimator

        Lambda0.Tau1Tau2.hat  <- sum(lambda0.t.hat[Tau1Tau2.times]) 
        
        # ----------------------------------------------------------------------
        # Estimation of the pure risk in [Tau1, Tau2] and for covariate profile 
        # x --------------------------------------------------------------------
        
        if ((length(x) != length(beta.hat)) | (is.null(x))) {
          x             <- rep(0, length(beta.hat))
          } # if  no covariate profile provided, use 0 as reference level
        exp.x           <- c(exp(x %*% beta.hat))
        Pi.x.hat        <- 1 - exp(- exp.x * Lambda0.Tau1Tau2.hat)
        
        return(list(beta.hat = beta.hat, lambda0.t.hat = lambda0.t.hat, 
                    Lambda0.Tau1Tau2.hat = Lambda0.Tau1Tau2.hat, 
                    Pi.x.Tau1Tau2.hat = Pi.x.hat))
      }
    }
    
  } else {
    # --------------------------------------------------------------------------
    # If no missing covariate data ---------------------------------------------
    
    # --------------------------------------------------------------------------
    # Quantities needed for the estimation -------------------------------------
    
    mod.detail            <- coxph.detail(mod, riskmat = TRUE)
    riskmat               <- mod.detail$riskmat 
    X                     <- model.matrix(mod) 
    weights               <- mod$weights # weights used for the fit, omega_i
    if (is.null(weights)) {
      weights             <- 1
      } # when using the whole cohort, omega_i = 1
    
    # --------------------------------------------------------------------------
    # Estimation of the log-relative hazard ------------------------------------
    
    beta.hat              <- mod$coefficients 
    
    exp.X.weighted        <- weights * exp(X %*% beta.hat)
    S0t                   <- t(riskmat) %*% (exp.X.weighted)
    observed.times        <- apply(riskmat, 1,
                                       function(v) {which.max(cumsum(v))})
    dNt                   <- matrix(0, nrow(riskmat), ncol(riskmat)) 
    dNt[cbind(1:nrow(riskmat), observed.times)] <- 1 
    dNt                   <- sweep(dNt, 1, mod$y[, ncol(mod$y)], "*") 
    if (is.null(Tau1)) {
      Tau1                <- floor(min(mod.detail$time))
    }
    if (is.null(Tau2)) {
      Tau2                <- floor(max(mod.detail$time))
    }
    Tau1Tau2.times        <- which((mod.detail$time > Tau1 &
                                      mod.detail$time <= Tau2)) 
    
    # --------------------------------------------------------------------------
    # Estimation of the baseline hazards at each unique event time and 
    # cumulative baseline hazard in time interval [Tau1, Tau2] -----------------
    
    lambda0.t.hat         <- t(colSums(dNt) / S0t) # we do not use weights in 
    # the numerator of the Breslow estimator

    Lambda0.Tau1Tau2.hat  <- sum(lambda0.t.hat[Tau1Tau2.times]) 
    
    # --------------------------------------------------------------------------
    # Estimation of the pure risk in [Tau1, Tau2] and for covariate profile x
    
    if ((length(x) != length(beta.hat)) | (is.null(x))) {
      x                   <- rep(0, length(beta.hat))
      } # if no covariate profile provided, use 0 as reference level
    exp.x                 <- c(exp(x %*% beta.hat))
    Pi.x.hat              <- 1 - exp(- exp.x * Lambda0.Tau1Tau2.hat) 
    
    return(list(beta.hat = beta.hat, lambda0.t.hat = lambda0.t.hat, 
                Lambda0.Tau1Tau2.hat = Lambda0.Tau1Tau2.hat, 
                Pi.x.Tau1Tau2.hat = Pi.x.hat))
  }
}


## -----------------------------------------------------------------------------
## Function: estimation.CumBH()
## -----------------------------------------------------------------------------
## Description: This function estimates the log-relative hazard, baseline 
##              hazards at each unique event time and cumulative baseline hazard 
##              in the time interval [Tau1, Tau2]
## -----------------------------------------------------------------------------
## Arguments:
##
##  mod             a cox model object, result of function coxph 
##
##  Tau1            left bound of the time interval considered for the 
##                  cumulative baseline hazard. Default is the first event time
##
##  Tau2            right bound of the time interval considered for the 
##                  cumulative baseline hazard. Default is the last event time
##
##
##  missing.data    is data on the p covariates missing for certain individuals
##                  in the stratified case cohort? If missing.data = TRUE, the
##                  arguments below need to be provided
##
##  riskmat.phase2  at risk matrix for the stratified case cohort (phase-two 
##                  data) at all of the cases event times, even that with 
##                  missing covariate data. If missing.data = TRUE, this 
##                  argument needs to be provided
##
##  dNt.phase2      counting process matrix for failures in the stratified
##                  case cohort (phase-two data). 
##
##  status.phase2   vector indicating the case status in the stratified case 
##                  cohort (phase-two data). If missing.data = TRUE, this 
##                  argument needs to be provided
## -----------------------------------------------------------------------------

estimation.CumBH  <- function (mod, Tau1 = NULL, Tau2 = NULL, 
                               missing.data = FALSE, riskmat.phase2 = NULL, 
                               dNt.phase2 = NULL, status.phase2 = NULL) {
  
  if (is.null(missing.data)) {
    missing.data = FALSE
  }
  
  # ----------------------------------------------------------------------------
  # If missing covariate data --------------------------------------------------
  
  if (missing.data == TRUE) {
    if (is.null(riskmat.phase2)) {
      stop("The at risk matrix for the phase-two data at all of the cases event 
           times need to be provided, even that with missing covariate data")
    } else {
      if (is.null(dNt.phase2)&is.null(status.phase2)) {
        stop("The status for the phase-two data, or the counting process matrix 
             for failures in the phase-two data, need to be provided")
      } else {
        if (is.null(dNt.phase2)) {
          
          # --------------------------------------------------------------------
          # If not provided, compute the counting process matrix for failures in 
          # the phase-two data -------------------------------------------------
          
          observed.times.phase2 <- apply(riskmat.phase2, 1, 
                                         function(v) {which.max(cumsum(v))}) 
          dNt.phase2            <- matrix(0, nrow(riskmat.phase2), 
                                          ncol(riskmat.phase2)) 
          dNt.phase2[cbind(1:nrow(riskmat.phase2), observed.times.phase2)] <-1 
          dNt.phase2            <- sweep(dNt.phase2, 1, status.phase2, "*") 
        }
        
        # ----------------------------------------------------------------------
        # Quantities needed for the estimation ---------------------------------
        
        mod.detail            <- coxph.detail(mod, riskmat = TRUE)
        riskmat               <- mod.detail$riskmat 
        X                     <- model.matrix(mod) 
        indiv.phase3          <- row.names(X)
        riskmat.phase3        <- riskmat.phase2[indiv.phase3, ] # use the cases' 
        # actual failure time, known even for cases with missing covariate data
        weights               <- mod$weights # weights used for the fit, omega_i
        if (is.null(weights)) {
          weights             <- 1 # when using the whole cohort, omega_i = 1
        } 
        
        # ----------------------------------------------------------------------
        # Estimation of the log relative hazard --------------------------------
        
        beta.hat              <- mod$coefficients 
        
        exp.X.weighted  <- weights * exp(X %*% beta.hat)
        S0t             <- t(riskmat.phase3) %*% (exp.X.weighted) # for any 
        # event time t (all of the cases failure times, even the one of cases 
        # who have missing covariate data)
        cases.times     <- colnames(riskmat.phase2) # cases failure times, even 
        # the one of cases who have missing covariate data
        if (is.null(Tau1)) {
          Tau1          <- floor(min(mod.detail$time))
        }
        if (is.null(Tau2)) {
          Tau2          <- floor(max(mod.detail$time))
        }
        Tau1Tau2.times  <- which((cases.times > Tau1 & cases.times <= Tau2))
        
        # ----------------------------------------------------------------------
        # Estimation of the baseline hazards at each unique event time and 
        # cumulative baseline hazard in time interval [Tau1, Tau2] -------------
        
        lambda0.t.hat   <- t(colSums(dNt.phase2) / S0t) # we do not use weights 
        # in the numerator of the Breslow estimator

        Lambda0.Tau1Tau2.hat  <- sum(lambda0.t.hat[Tau1Tau2.times])
        
        return(list(beta.hat = beta.hat, lambda0.t.hat = lambda0.t.hat, 
                    Lambda0.Tau1Tau2.hat = Lambda0.Tau1Tau2.hat))
      }
    }
  } else {
    
    # --------------------------------------------------------------------------
    # Quantities needed for the estimation -------------------------------------
    
    mod.detail            <- coxph.detail(mod, riskmat = TRUE)
    riskmat               <- mod.detail$riskmat 
    X                     <- model.matrix(mod)
    weights               <- mod$weights # weights used for the fit, omega_i
    if (is.null(weights)) {
      weights             <- 1
      } # when using the whole cohort, omega_i = 1
    
    # --------------------------------------------------------------------------
    # Estimation of the log-relative hazard ------------------------------------
    
    beta.hat              <- mod$coefficients 
    
    # --------------------------------------------------------------------------
    # Estimation of the baseline hazards at each unique event time and 
    # cumulative baseline hazard in time interval [Tau1, Tau2] -----------------
    
    exp.X.weighted        <- weights * exp(X %*% beta.hat)
    S0t                   <- t(riskmat) %*% (exp.X.weighted) 
    observed.times        <- apply(riskmat, 1, 
                                   function(v) {which.max(cumsum(v))}) 
    dNt                   <- matrix(0, nrow(riskmat), ncol(riskmat)) 
    dNt[cbind(1:nrow(riskmat), observed.times)] <- 1 
    dNt                   <- sweep(dNt, 1, mod$y[, ncol(mod$y)], "*") 
    
    lambda0.t.hat         <- t(colSums(dNt) / S0t) # we do not use weights in
    # the numerator of the Breslow estimator
    if (is.null(Tau1)) {
      Tau1                <- floor(min(mod.detail$time))
      }
    if (is.null(Tau2)) {
      Tau2                <- floor(max(mod.detail$time))
      }
    Tau1Tau2.times        <- which((mod.detail$time > Tau1 &
                                      mod.detail$time <= Tau2)) 
    Lambda0.Tau1Tau2.hat  <- sum(lambda0.t.hat[Tau1Tau2.times]) 
    
    return(list(beta.hat = beta.hat, lambda0.t.hat = lambda0.t.hat, 
                Lambda0.Tau1Tau2.hat = Lambda0.Tau1Tau2.hat))
  }
}


## -----------------------------------------------------------------------------
## Function: estimation.PR()
## -----------------------------------------------------------------------------
## Description: This function estimates the pure risk in the time interval 
##              [Tau1, Tau2] and for covariate profile x, from the values of the
##              log-relative hazard and cumulative baseline hazard in 
##              [Tau1, Tau2]
## -----------------------------------------------------------------------------
##  beta             log-relative hazard 
##
##  Lambda0.Tau1Tau2  cumulative baseline hazard in the time interval that is   
##                    considered for the pure risk
##
##  x                 vector of length p, specifying the covariate profile 
##                    considered for the pure risk. Default is (0,...,0)
## -----------------------------------------------------------------------------

estimation.PR     <- function (beta, Lambda0.Tau1Tau2, x = NULL) {
  
  if ((length(x) != length(beta)) | (is.null(x))) {
    x     <- rep(0, length(beta))
    } # if no covariate profile provided, use 0 as reference level

  exp.x   <- c(exp(x %*% beta))
  Pi.x    <- 1 - exp(- exp.x * Lambda0.Tau1Tau2) 
  
  return(Pi.x.Tau1Tau2.hat = Pi.x)
}


## -----------------------------------------------------------------------------
## Function: influences()
## -----------------------------------------------------------------------------
## Description: This function estimates the influences of the log-relative 
##              hazard, baseline hazards at each unique event time, cumulative 
##              baseline hazard in the time interval [Tau1, Tau2] and pure risk 
##              in [Tau1, Tau2] and for covariate profile x. It also returns the
##              parameters estimates. If calibrated weights were used and the
##              auxiliary variables observations are provided, the function
##              returns the overall influences, as well as the phase-two
##              influences
## -----------------------------------------------------------------------------
## Arguments:
##
##  mod             a cox model object, result of function coxph 
##
##  Tau1            left bound of the time interval considered for the 
##                  cumulative baseline hazard and pure risk. Default is the 
##                  first event time
##
##  Tau2            right bound of the time interval considered for the 
##                  cumulative baseline hazard and pure risk. Default is the 
##                  last event time
##
##  x               vector of length p, specifying the covariate profile 
##                  considered for the pure risk. Default is (0,...,0)
##
##  calibrated      are calibrated weights used for the estimation of the log-  
##                  relative hazard, cumulative baseline hazard and pure risk?
##                  If calibrated = TRUE, the argument below needs to be 
##                  provided
##
##  A               matrix with the values of the auxiliary variables used for 
##                  the calibration of the weights in the whole cohort. If 
##                  calibrated = TRUE, this argument needs to be provided
## -----------------------------------------------------------------------------

influences  <- function (mod, Tau1 = NULL, Tau2 = NULL, x = NULL,
                         calibrated = NULL, A = NULL) {
  
  # ----------------------------------------------------------------------------
  # Quantities needed for the influences ---------------------------------------
  
  mod.detail        <- coxph.detail(mod, riskmat = T)
  riskmat           <- mod.detail$riskmat 
  number.times      <- ncol(riskmat) 
  n                 <- nrow(riskmat) 
  X                 <- model.matrix(mod) 
  p                 <- ncol(X) 
  weights           <- mod$weights # weights used for the fit, omega_i
  if (is.null(weights)) {
    weights         <- 1
    } # when using the whole cohort, omega_i = 1

  beta.hat          <- mod$coefficients
  
  exp.X.weighted    <- weights * exp(X %*% beta.hat)
  Y.exp.X.weighted  <- riskmat * matrix(exp.X.weighted, nrow = n, 
                                        ncol = number.times, byrow = FALSE) 
  S0t               <- t(riskmat) %*% (exp.X.weighted) 
  S1t               <- t(riskmat) %*% (X * matrix(exp.X.weighted, nrow = n,
                                                  ncol = p, byrow = FALSE)) 
  observed.times    <- apply(riskmat, 1, function(v) {which.max(cumsum(v))}) 
  dNt               <- matrix(0, n, number.times) 
  dNt[cbind(1:nrow(riskmat), observed.times)] <- 1 
  dNt               <- dNt * matrix(mod$y[, ncol(mod$y)], nrow(riskmat), 
                                    number.times, byrow = FALSE) 
  dNt.weighted      <- dNt * matrix(weights, nrow = n, ncol = number.times, 
                                    byrow = FALSE) 
  infomat.indiv     <- mod.detail$imat  
  infomat           <- apply(X = infomat.indiv, MARGIN = c(1,2), FUN = sum)
  if (is.null(Tau1)) {
    Tau1            <- floor(min(mod.detail$time))
  }
  if (is.null(Tau2)) {
    Tau2            <- floor(max(mod.detail$time))
  }
  Tau1Tau2.times    <- which((Tau1 < mod.detail$time) & 
                               (mod.detail$time <= Tau2)) 
  Tau1Times         <- which((mod.detail$time <= Tau1)) 
  Tau2Times         <- which((mod.detail$time <= Tau2))
  
  lambda0.t.hat     <- t(colSums(dNt) / S0t) 
  
  Lambda0.Tau1Tau2.hat <- sum(lambda0.t.hat[Tau1Tau2.times])
  
  if ((length(x) != p) | (is.null(x))) {
    x               <- rep(0, p)
  } # if no covariate profile provided, use 0 as reference level
  exp.x             <- c(exp(x %*% beta.hat))
  
  Pi.x.hat          <- 1 - exp(- exp.x * Lambda0.Tau1Tau2.hat)

  if (is.null(calibrated)) {
    calibrated <- FALSE
    }
  
  # ----------------------------------------------------------------------------
  # If calibrated weights ------------------------------------------------------
  
  if (calibrated == TRUE) {
    if (is.null(A)) {
      print("A matrix with auxiliary variables values need to be provided")
      break}
  
    # --------------------------------------------------------------------------
    # Computation of the influences for the Lagrangian multipliers, eta --------
    
    indiv.phase2      <- row.names(X)
    A.phase2          <- A[indiv.phase2, ]
    q                 <- ncol(A.phase2)
    
    A.weighted         <- sweep(A.phase2, 1, weights, "*")
    AA.weighted        <- array(NA, dim = c(q, q, n))
    for(i in 1:nrow(A.phase2)) {
      AA.weighted[,, i] <- tcrossprod(A.phase2[i,], A.weighted[i,])
      }
    sum.AA.weighted    <- apply(X = AA.weighted, MARGIN = c(1,2), FUN = sum)
    sum.AA.weighted.inv <- solve(sum.AA.weighted)
    
    infl1.eta         <- A %*% sum.AA.weighted.inv
    infl2.eta         <- 0 * A
    infl2.eta[indiv.phase2,] <- - A.weighted %*% sum.AA.weighted.inv
    infl.eta          <- infl1.eta + infl2.eta 
    
    # --------------------------------------------------------------------------
    # Computation of the influences for log-relative hazard, beta --------------

    XA              <- array(NA, dim = c(p, q, n))
    for(i in 1:n) {
      XA[,, i]      <- tcrossprod(X[i,], A.phase2[i,])
      }
    drond.G1t.eta   <- array(NA, dim = c(p, q, number.times))
    drond.G0t.eta   <- array(NA, dim = c(q, number.times))
    drond.S1t.eta   <- array(NA, dim = c(p, q, number.times))
    drond.S0t.eta   <- array(NA, dim = c(q, number.times))
    for(t in 1:number.times) { 
      drond.G1t.eta[,, t] <- apply(X = sweep(XA, 3, dNt.weighted[, t], "*"), 
                                   MARGIN = c(1,2), FUN = sum)
      drond.G0t.eta[, t]  <- colSums(A.phase2 * matrix(dNt.weighted[, t], 
                                                       nrow = n, ncol = q, 
                                                       byrow = FALSE))
      drond.S1t.eta[,, t] <- apply(X = sweep(XA, 3, Y.exp.X.weighted[, t],
                                               "*"), MARGIN = c(1,2), FUN = sum)
      drond.S0t.eta[, t]  <- colSums(A.phase2 * matrix(Y.exp.X.weighted[, t], 
                                                       nrow = n, ncol = q, 
                                                       byrow = FALSE))
    }
    drond.Ut.eta    <- array(NA, dim = c(p, q, number.times))
    X.expect        <- S1t / matrix(S0t, nrow = number.times, ncol = p, 
                                    byrow = FALSE)
    for(t in 1:number.times) { 
      drond.Ut.eta[,, t] <- drond.G1t.eta[,, t] - tcrossprod(X.expect[t, ],
                                                             drond.G0t.eta[, t]) -
                            (colSums(dNt.weighted) / S0t)[t] * drond.S1t.eta[,, t] +
                            (colSums(dNt.weighted) / S0t ^ 2)[t] * 
                            tcrossprod(S1t[t, ], drond.S0t.eta[, t])
    }
    drond.U.eta     <- apply(X = drond.Ut.eta, MARGIN = c(1,2), FUN = sum)
    infl1.beta      <- infl1.eta %*% t(drond.U.eta) %*% solve(infomat)
    score           <- array(NA, dim = c(n, p, number.times))
    for(i in 1:n) {
      score[i,,]    <- - sweep(t(X.expect), 1, X[i,], "-")
      }
    for(t in 1:number.times) {
      dMt           <- dNt.weighted[, t] - riskmat[, t] *
        exp.X.weighted * (colSums(dNt.weighted)/S0t)[t]
      score[,, t]   <- score[,, t] * matrix(dMt, nrow = n, ncol = p, 
                                            byrow = FALSE)
    }
    score.beta      <- apply(X = score, MARGIN = c(1,2), FUN = sum) 
    infl2.beta      <- infl2.eta %*% t(drond.U.eta) %*% solve(infomat)
    infl2.beta[indiv.phase2,] <- infl2.beta[indiv.phase2,] +
      score.beta %*% solve(infomat) 
    infl.beta       <- infl1.beta + infl2.beta
    
    # --------------------------------------------------------------------------
    # Computation of the influences for the baseline hazards at each unique 
    # event time and for the cumulative baseline hazard in time interval 
    # [Tau1, Tau2] -------------------------------------------------------------
    
    infl1.lambda0.t <- (- (infl1.beta %*% t(S1t) + infl1.eta %*% drond.S0t.eta) *
                        matrix(lambda0.t.hat, nrow = nrow(A), 
                                 ncol = number.times, byrow = TRUE)) /
                        matrix(S0t, nrow = nrow(A), ncol = number.times, 
                               byrow = TRUE) 
    infl2.lambda0.t <- 0 * infl1.lambda0.t
    infl2.lambda0.t[indiv.phase2,] <- (dNt - (Y.exp.X.weighted + 
                                              infl2.beta[indiv.phase2,] %*% 
                                              t(S1t) + infl2.eta[indiv.phase2,] %*%
                                              drond.S0t.eta) * 
                                      matrix(lambda0.t.hat, nrow = n, 
                                             ncol = number.times, byrow = TRUE)) /
                                      matrix(S0t, nrow = n, ncol = number.times,
                                             byrow = TRUE) 
    infl.lambda0.t  <- infl1.lambda0.t + infl2.lambda0.t

    infl1.Lambda0.Tau1Tau2  <- rowSums(infl1.lambda0.t[, Tau1Tau2.times])
    infl2.Lambda0.Tau1Tau2  <- rowSums(infl2.lambda0.t[, Tau1Tau2.times])
    infl.Lambda0.Tau1Tau2   <- infl1.Lambda0.Tau1Tau2 + infl2.Lambda0.Tau1Tau2 
  
    # --------------------------------------------------------------------------
    # Computation of the influences for the pure risk in [Tau1, Tau2] and for 
    # covariate profile x ------------------------------------------------------
    
    infl1.Pi.x  <- c((1 - Pi.x.hat) * exp.x * Lambda0.Tau1Tau2.hat) * 
                  infl1.beta %*% x + c((1 - Pi.x.hat) * exp.x) * 
                  infl1.Lambda0.Tau1Tau2
    infl2.Pi.x  <- c((1 - Pi.x.hat) * exp.x * Lambda0.Tau1Tau2.hat) * 
                  infl2.beta %*% x + c((1 - Pi.x.hat) * exp.x) * 
                  infl2.Lambda0.Tau1Tau2
    infl.Pi.x   <- infl1.Pi.x + infl2.Pi.x
    
    return(list(infl.beta = infl.beta, infl.lambda0.t = infl.lambda0.t, 
                infl.Lambda0.Tau1Tau2 = as.matrix(infl.Lambda0.Tau1Tau2), 
                infl.Pi.x.Tau1Tau2 = infl.Pi.x, infl2.beta = infl2.beta, 
                infl2.lambda0.t = infl2.lambda0.t, 
                infl2.Lambda0.Tau1Tau2 = as.matrix(infl2.Lambda0.Tau1Tau2), 
                infl2.Pi.x.Tau1Tau2 = infl2.Pi.x, beta.hat = beta.hat, 
                lambda0.t.hat = lambda0.t.hat, 
                Lambda0.Tau1Tau2.hat = Lambda0.Tau1Tau2.hat, 
                Pi.x.Tau1Tau2.hat = Pi.x.hat))
  }else{
    # --------------------------------------------------------------------------
    # If design weights --------------------------------------------------------
    
    # --------------------------------------------------------------------------
    # Computation of the influences for log-relative hazard, beta --------------
    
    infl.beta <- residuals(mod, type = "dfbeta", weighted = T)

    # --------------------------------------------------------------------------
    # Computation of the influences for the baseline hazards at each unique 
    # event time and for the cumulative baseline hazard in time interval 
    # [Tau1, Tau2] -------------------------------------------------------------
    
    infl.lambda0.t <- (dNt - (Y.exp.X.weighted + infl.beta %*% t(S1t)) * 
                      matrix(lambda0.t.hat, nrow = n, ncol = number.times,
                                   byrow = TRUE)) /
                      matrix(S0t, nrow = n, ncol = number.times, 
                                 byrow = TRUE) 
    infl.Lambda0.Tau1Tau2 <- rowSums(infl.lambda0.t[, Tau1Tau2.times])
    
    # --------------------------------------------------------------------------
    # Computation of the influences for the pure risk in [Tau1, Tau2] and for 
    # covariate profile x ------------------------------------------------------
    
    infl.Pi.x <- c((1 - Pi.x.hat) * exp.x * Lambda0.Tau1Tau2.hat) *
                infl.beta %*% x + c((1 - Pi.x.hat) * exp.x) * 
                infl.Lambda0.Tau1Tau2
    
    return(list(infl.beta = infl.beta, infl.lambda0.t = infl.lambda0.t, 
                infl.Lambda0.Tau1Tau2 = as.matrix(infl.Lambda0.Tau1Tau2), 
                infl.Pi.x.Tau1Tau2 = infl.Pi.x, beta.hat = beta.hat, 
                lambda0.t.hat = lambda0.t.hat, 
                Lambda0.Tau1Tau2.hat = Lambda0.Tau1Tau2.hat, 
                Pi.x.Tau1Tau2.hat = Pi.x.hat))
  }
}


## -----------------------------------------------------------------------------
## Function: influences.RH()
## -----------------------------------------------------------------------------
## Description: This function estimates the influences of the log-relative 
##              hazard and also returns the parameter estimate. If calibrated 
##              weights were used and the auxiliary variables observations are 
##              provided, the function returns the overall influences, as well 
##              as the phase-two influences
## -----------------------------------------------------------------------------
## Arguments:
##
##  mod             a cox model object, result of function coxph 
##
##  calibrated      are calibrated weights used for the estimation of the log-  
##                  relative hazard, cumulative baseline hazard and pure risk?
##                  If calibrated = TRUE, the argument below needs to be 
##                  provided
##
##  A               matrix with the values of the auxiliary variables used for 
##                  the calibration of the weights in the whole cohort. If 
##                  calibrated = TRUE, this argument needs to be provided
## -----------------------------------------------------------------------------

influences.RH   <- function (mod, calibrated = NULL, A = NULL) {
  
  if (is.null(calibrated)) {
    calibrated <- FALSE
  }
  
  # ----------------------------------------------------------------------------
  # If calibrated weights ------------------------------------------------------
  
  if (calibrated == TRUE) {
    if (is.null(A)) {
      print("A matrix with auxiliary variables values need to be provided")
      break}
    
    # --------------------------------------------------------------------------
    # Quantities needed for the influences -------------------------------------
    
    mod.detail        <- coxph.detail(mod, riskmat = T)
    riskmat           <- mod.detail$riskmat 
    number.times      <- ncol(riskmat) 
    n                 <- nrow(riskmat) 
    beta.hat          <- mod$coefficients 
    X                 <- model.matrix(mod) 
    p                 <- ncol(X)
    weights           <- mod$weights # weights used for the fit, omega_i
    if (is.null(weights)) {
      weights           <- 1
    } # when using the whole cohort, omega_i = 1
    exp.X.weighted    <- weights * exp(X %*% beta.hat)
    Y.exp.X.weighted  <- riskmat * matrix(exp.X.weighted, nrow = n, 
                                          ncol = number.times, byrow = FALSE)  
    S0t               <- t(riskmat) %*% (exp.X.weighted) 
    S1t               <- t(riskmat) %*% (X * matrix(exp.X.weighted, nrow = n,
                                                    ncol = p, byrow = FALSE))  
    observed.times    <- apply(riskmat, 1, function(v) {which.max(cumsum(v))}) 
    dNt               <- matrix(0, n, number.times) 
    dNt[cbind(1:nrow(riskmat), observed.times)] <- 1 
    dNt               <- dNt * matrix(mod$y[, ncol(mod$y)], nrow(riskmat), 
                                      number.times, byrow = FALSE) 
    dNt.weighted      <- dNt * matrix(weights, nrow = n, ncol = number.times, 
                                      byrow = FALSE) 
    infomat.indiv     <- mod.detail$imat
    infomat           <- apply(X = infomat.indiv, MARGIN = c(1,2), FUN = sum) 
    
    # --------------------------------------------------------------------------
    # Computation of the influences for the Lagrangian multipliers, eta --------
    
    indiv.phase2      <- row.names(X)
    A.phase2          <- A[indiv.phase2, ]
    q                 <- ncol(A.phase2)
    
    A.weighted        <- sweep(A.phase2, 1, weights, "*")
    AA.weighted       <- array(NA, dim = c(q, q, n))
    for(i in 1:nrow(A.phase2)) {
      AA.weighted[,, i] <- tcrossprod(A.phase2[i,], A.weighted[i,])
    }
    sum.AA.weighted     <- apply(X = AA.weighted, MARGIN = c(1,2), FUN = sum)
    sum.AA.weighted.inv <- solve(sum.AA.weighted)
    
    infl1.eta         <- A %*% sum.AA.weighted.inv
    infl2.eta         <- 0 * A
    infl2.eta[indiv.phase2,] <- - A.weighted %*% sum.AA.weighted.inv
    infl.eta          <- infl1.eta + infl2.eta 
    
    # --------------------------------------------------------------------------
    # Computation of the influences for log-relative hazard, beta --------------
    
    XA <- array(NA, dim = c(p, q, n))
    for(i in 1:n) {
      XA[,, i] <- tcrossprod(X[i,], A.phase2[i,])
    }
    drond.G1t.eta   <- array(NA, dim = c(p, q, number.times))
    drond.G0t.eta   <- array(NA, dim = c(q, number.times))
    drond.S1t.eta   <- array(NA, dim = c(p, q, number.times))
    drond.S0t.eta   <- array(NA, dim = c(q, number.times))
    for(t in 1:number.times) { 
      drond.G1t.eta[,, t] <- apply(X = sweep(XA, 3, dNt.weighted[, t], "*"), 
                                   MARGIN = c(1,2), FUN = sum)
      drond.G0t.eta[, t]  <- colSums(A.phase2 * matrix(dNt.weighted[, t], 
                                                       nrow = n, ncol = q, 
                                                       byrow = FALSE))
      drond.S1t.eta[,, t] <- apply(X = sweep(XA, 3, Y.exp.X.weighted[, t], "*"), 
                                   MARGIN = c(1,2), FUN = sum)
      drond.S0t.eta[, t]  <- colSums(A.phase2 * matrix(Y.exp.X.weighted[, t], 
                                                       nrow = n, ncol = q, 
                                                       byrow = FALSE))
    }
    drond.Ut.eta  <- array(NA, dim = c(p, q, number.times))
    X.expect      <- S1t / matrix(S0t, nrow = number.times, ncol = p, 
                             byrow = FALSE)
    for(t in 1:number.times) { 
      drond.Ut.eta[,, t] <- drond.G1t.eta[,, t] - tcrossprod(X.expect[t, ], 
                                                             drond.G0t.eta[, t]) -
                            (colSums(dNt.weighted)/S0t)[t] * drond.S1t.eta[,, t] +
                            (colSums(dNt.weighted)/S0t^2)[t] * 
                            tcrossprod(S1t[t, ], drond.S0t.eta[, t])
    }
    drond.U.eta   <- apply(X = drond.Ut.eta, MARGIN = c(1,2), FUN = sum)
    infl1.beta    <- infl1.eta %*% t(drond.U.eta) %*% solve(infomat)
    score         <- array(NA, dim = c(n, p, number.times))
    for(i in 1:n) {
      score[i,,]  <- - sweep(t(X.expect), 1, X[i,], "-")
    }
    for(t in 1:number.times) {
      dMt           <- dNt.weighted[, t] - riskmat[, t] * exp.X.weighted * 
                      (colSums(dNt.weighted)/S0t)[t]
      score[,, t]   <- score[,, t] * matrix(dMt, nrow = n, ncol = p, 
                                            byrow = FALSE)
    }
    score.beta        <- apply(X = score, MARGIN = c(1,2), FUN = sum)
    infl2.beta        <- infl2.eta %*% t(drond.U.eta) %*% solve(infomat)
    infl2.beta[indiv.phase2,] <- infl2.beta[indiv.phase2,] + score.beta %*%
                                  solve(infomat) 
    infl.beta         <- infl1.beta + infl2.beta
    
    return(list(infl.beta = infl.beta, infl2.beta = infl2.beta, 
                beta.hat = beta.hat))
  }else{
    # --------------------------------------------------------------------------
    # If design weights --------------------------------------------------------
    
    infl.beta <- residuals(mod, type = "dfbeta", weighted = T)
  
    return(list(infl.beta = infl.beta, beta.hat = beta.hat))
  }
}


## -----------------------------------------------------------------------------
## Function: influences.CumBH()
## -----------------------------------------------------------------------------
## Description: This function estimates the influences of the log-relative 
##              hazard, baseline hazards at each unique event time, and 
##              cumulative baseline hazard in the time interval [Tau1, Tau2]. It 
##              also returns the parameters estimates. If calibrated weights 
##              were used and the auxiliary variables observations are provided, 
##              the function returns the overall influences, as well as the 
##              phase-two influences
## -----------------------------------------------------------------------------
## Arguments:
##
##  mod             a cox model object, result of function coxph 
##
##  Tau1            left bound of the time interval considered for the 
##                  cumulative baseline hazard. Default is the first event time
##
##  Tau2            right bound of the time interval considered for the 
##                  cumulative baseline hazard. Default is the last event time
##
##  calibrated      are calibrated weights used for the estimation of the log-  
##                  relative hazard, cumulative baseline hazard and pure risk?
##                  If calibrated = TRUE, the argument below needs to be 
##                  provided
##
##  A               matrix with the values of the auxiliary variables used for 
##                  the calibration of the weights in the whole cohort. If 
##                  calibrated = TRUE, this argument needs to be provided
## -----------------------------------------------------------------------------

influences.CumBH  <- function (mod, Tau1 = NULL, Tau2 = NULL, A = NULL, 
                               calibrated = NULL) {
  
  # ----------------------------------------------------------------------------
  # Quantities needed for the influences ---------------------------------------
  
  mod.detail        <- coxph.detail(mod, riskmat = T)
  riskmat           <- mod.detail$riskmat 
  number.times      <- ncol(riskmat) 
  n                 <- nrow(riskmat) 
  X                 <- model.matrix(mod) 
  p                 <- ncol(X) 
  weights           <- mod$weights # weights used for the fit, omega_i
  if (is.null(weights)) {
    weights         <- 1
  } # when using the whole cohort, omega_i = 1
  
  beta.hat          <- mod$coefficients
  
  exp.X.weighted    <- weights * exp(X %*% beta.hat)
  Y.exp.X.weighted  <- riskmat * matrix(exp.X.weighted, nrow = n, 
                                        ncol = number.times, byrow = FALSE) 
  S0t               <- t(riskmat) %*% (exp.X.weighted) 
  S1t               <- t(riskmat) %*% (X * matrix(exp.X.weighted, nrow = n,
                                                  ncol = p, byrow = FALSE)) 
  observed.times    <- apply(riskmat, 1, function(v) {which.max(cumsum(v))}) 
  dNt               <- matrix(0, n, number.times) 
  dNt[cbind(1:nrow(riskmat), observed.times)] <- 1 
  dNt               <- dNt * matrix(mod$y[, ncol(mod$y)], nrow(riskmat), 
                                    number.times, byrow = FALSE) 
  dNt.weighted      <- dNt * matrix(weights, nrow = n, ncol = number.times, 
                                    byrow = FALSE) 
  infomat.indiv     <- mod.detail$imat  
  infomat           <- apply(X = infomat.indiv, MARGIN = c(1,2), FUN = sum)
  if (is.null(Tau1)) {
    Tau1            <- floor(min(mod.detail$time))
  }
  if (is.null(Tau2)) {
    Tau2            <- floor(max(mod.detail$time))
  }
  Tau1Tau2.times    <- which((Tau1 < mod.detail$time) & 
                               (mod.detail$time <= Tau2)) 
  Tau1Times         <- which((mod.detail$time <= Tau1)) 
  Tau2Times         <- which((mod.detail$time <= Tau2))
  
  lambda0.t.hat     <- t(colSums(dNt) / S0t) 
  
  Lambda0.Tau1Tau2.hat <- sum(lambda0.t.hat[Tau1Tau2.times])
  
  # ----------------------------------------------------------------------------
  # If calibrated weights ------------------------------------------------------
  
  if (calibrated == TRUE) {
    if (is.null(A)) {
      print("A matrix with auxiliary variables values need to be provided")
      break}
    
    # --------------------------------------------------------------------------
    # Computation of the influences for the Lagrangian multipliers, eta --------
    
    indiv.phase2      <- row.names(X)
    A.phase2          <- A[indiv.phase2, ]
    q                 <- ncol(A.phase2)
    
    A.weighted         <- sweep(A.phase2, 1, weights, "*")
    AA.weighted        <- array(NA, dim = c(q, q, n))
    for(i in 1:nrow(A.phase2)) {
      AA.weighted[,, i] <- tcrossprod(A.phase2[i,], A.weighted[i,])
    }
    sum.AA.weighted    <- apply(X = AA.weighted, MARGIN = c(1,2), FUN = sum)
    sum.AA.weighted.inv <- solve(sum.AA.weighted)
    
    infl1.eta         <- A %*% sum.AA.weighted.inv
    infl2.eta         <- 0 * A
    infl2.eta[indiv.phase2,] <- - A.weighted %*% sum.AA.weighted.inv
    infl.eta          <- infl1.eta + infl2.eta 
    
    # --------------------------------------------------------------------------
    # Computation of the influences for log-relative hazard, beta --------------
    
    XA              <- array(NA, dim = c(p, q, n))
    for(i in 1:n) {
      XA[,, i]      <- tcrossprod(X[i,], A.phase2[i,])
    }
    drond.G1t.eta   <- array(NA, dim = c(p, q, number.times))
    drond.G0t.eta   <- array(NA, dim = c(q, number.times))
    drond.S1t.eta   <- array(NA, dim = c(p, q, number.times))
    drond.S0t.eta   <- array(NA, dim = c(q, number.times))
    for(t in 1:number.times) { 
      drond.G1t.eta[,, t] <- apply(X = sweep(XA, 3, dNt.weighted[, t], "*"), 
                                   MARGIN = c(1,2), FUN = sum)
      drond.G0t.eta[, t]  <- colSums(A.phase2 * matrix(dNt.weighted[, t], 
                                                       nrow = n, ncol = q, 
                                                       byrow = FALSE))
      drond.S1t.eta[,, t] <- apply(X = sweep(XA, 3, Y.exp.X.weighted[, t],
                                             "*"), MARGIN = c(1,2), FUN = sum)
      drond.S0t.eta[, t]  <- colSums(A.phase2 * matrix(Y.exp.X.weighted[, t], 
                                                       nrow = n, ncol = q, 
                                                       byrow = FALSE))
    }
    drond.Ut.eta    <- array(NA, dim = c(p, q, number.times))
    X.expect        <- S1t / matrix(S0t, nrow = number.times, ncol = p, 
                                    byrow = FALSE)
    for(t in 1:number.times) { 
      drond.Ut.eta[,, t] <- drond.G1t.eta[,, t] - tcrossprod(X.expect[t, ],
                                                             drond.G0t.eta[, t]) -
        (colSums(dNt.weighted) / S0t)[t] * drond.S1t.eta[,, t] +
        (colSums(dNt.weighted) / S0t ^ 2)[t] * 
        tcrossprod(S1t[t, ], drond.S0t.eta[, t])
    }
    drond.U.eta     <- apply(X = drond.Ut.eta, MARGIN = c(1,2), FUN = sum)
    infl1.beta      <- infl1.eta %*% t(drond.U.eta) %*% solve(infomat)
    score           <- array(NA, dim = c(n, p, number.times))
    for(i in 1:n) {
      score[i,,]    <- - sweep(t(X.expect), 1, X[i,], "-")
    }
    for(t in 1:number.times) {
      dMt           <- dNt.weighted[, t] - riskmat[, t] *
        exp.X.weighted * (colSums(dNt.weighted)/S0t)[t]
      score[,, t]   <- score[,, t] * matrix(dMt, nrow = n, ncol = p, 
                                            byrow = FALSE)
    }
    score.beta      <- apply(X = score, MARGIN = c(1,2), FUN = sum) 
    infl2.beta      <- infl2.eta %*% t(drond.U.eta) %*% solve(infomat)
    infl2.beta[indiv.phase2,] <- infl2.beta[indiv.phase2,] +
      score.beta %*% solve(infomat) 
    infl.beta       <- infl1.beta + infl2.beta
    
    # --------------------------------------------------------------------------
    # Computation of the influences for the baseline hazards at each unique 
    # event time and for the cumulative baseline hazard in time interval 
    # [Tau1, Tau2] -------------------------------------------------------------
    
    infl1.lambda0.t <- (- (infl1.beta %*% t(S1t) + infl1.eta %*% drond.S0t.eta) *
                          matrix(lambda0.t.hat, nrow = nrow(A), 
                                 ncol = number.times, byrow = TRUE)) /
      matrix(S0t, nrow = nrow(A), ncol = number.times, 
             byrow = TRUE) 
    infl2.lambda0.t <- 0 * infl1.lambda0.t
    infl2.lambda0.t[indiv.phase2,] <- (dNt - (Y.exp.X.weighted + 
                                              infl2.beta[indiv.phase2,] %*% 
                                              t(S1t) +
                                              infl2.eta[indiv.phase2,] %*%
                                              drond.S0t.eta) * 
                                        matrix(lambda0.t.hat, nrow = n, 
                                                ncol = number.times, 
                                               byrow = TRUE)) /
                                        matrix(S0t, nrow = n, 
                                               ncol = number.times, 
                                               byrow = TRUE) 
    infl.lambda0.t  <- infl1.lambda0.t + infl2.lambda0.t
    infl1.Lambda0.Tau1Tau2  <- rowSums(infl1.lambda0.t[, Tau1Tau2.times])
    infl2.Lambda0.Tau1Tau2  <- rowSums(infl2.lambda0.t[, Tau1Tau2.times])
    infl.Lambda0.Tau1Tau2   <- infl1.Lambda0.Tau1Tau2 + infl2.Lambda0.Tau1Tau2 
    
    return(list(infl.beta = infl.beta, infl.lambda0.t = infl.lambda0.t, 
                infl.Lambda0.Tau1Tau2 = as.matrix(infl.Lambda0.Tau1Tau2), 
                infl2.beta = infl2.beta, infl2.lambda0.t = infl2.lambda0.t, 
                infl2.Lambda0.Tau1Tau2 = as.matrix(infl2.Lambda0.Tau1Tau2), 
                beta.hat = beta.hat, lambda0.t.hat = lambda0.t.hat, 
                Lambda0.Tau1Tau2.hat = Lambda0.Tau1Tau2.hat))
  }else{
    # --------------------------------------------------------------------------
    # If design weights --------------------------------------------------------
    
    # --------------------------------------------------------------------------
    # Computation of the influences for log-relative hazard, beta --------------
    
    infl.beta <- residuals(mod, type = "dfbeta", weighted = T)
    
    # --------------------------------------------------------------------------
    # Computation of the influences for the baseline hazards at each unique 
    # event time and for the cumulative baseline hazard in time interval 
    # [Tau1, Tau2] -------------------------------------------------------------
    
    infl.lambda0.t <- (dNt - (Y.exp.X.weighted + infl.beta %*% t(S1t)) * 
                         matrix(lambda0.t.hat, nrow = n, ncol = number.times,
                                byrow = TRUE)) /
      matrix(S0t, nrow = n, ncol = number.times, 
             byrow = TRUE) 
    infl.Lambda0.Tau1Tau2 <- rowSums(infl.lambda0.t[, Tau1Tau2.times])
    
    return(list(infl.beta = infl.beta, infl.lambda0.t = infl.lambda0.t, 
                infl.Lambda0.Tau1Tau2 = as.matrix(infl.Lambda0.Tau1Tau2), 
                beta.hat = beta.hat, lambda0.t.hat = lambda0.t.hat, 
                Lambda0.Tau1Tau2.hat = Lambda0.Tau1Tau2.hat))
  }
}


## -----------------------------------------------------------------------------
## Function: influences.PR()
## -----------------------------------------------------------------------------
## Description: This function estimates the influences of the pure risk in the 
##              time interval [Tau1, Tau2] and for covariate profile x, from
##              that of the log-relative hazard and cumulative baseline hazard 
##              [Tau1, Tau2]. It also returns the pure risk estimate. If
##              calibrated weights were used and the phase-two influences are
##              provided, the function returns the overall influences, as well 
##              as the phase-two influences
## -----------------------------------------------------------------------------
## Arguments:
##
##  beta                    log-relative hazard 
##
##  Lambda0.Tau1Tau2        cumulative baseline hazard in the time interval that 
##                          is considered for the pure risk
##
##  x                       vector of length p, specifying the covariate profile 
##                          considered for the pure risk. Default is (0,...,0)
##
##  infl.beta               log-relative hazard influences
##
##  infl.Lambda0.Tau1Tau2   influences for the cumulative baseline hazard in the  
##                          same time interval as that of the pure risk 
##
##  calibrated              are calibrated weights used for the estimation of 
##                          the log-relative hazard, cumulative baseline hazard 
##                          and pure risk? If calibrated = TRUE, the arguments 
##                          below need to be provided
##
##  infl2.beta              phase-two influences for the log-relative hazard
##                  
##  infl2.Lambda0.Tau1Tau2  phase-two influences for the cumulative baseline
##                          hazard in the same time interval as that of the pure 
##                          risk 
## -----------------------------------------------------------------------------

influences.PR  <- function (beta, Lambda0.Tau1Tau2, x = NULL, infl.beta, 
                            infl.Lambda0.Tau1Tau2, calibrated = NULL,
                            infl2.beta = NULL, infl2.Lambda0.Tau1Tau2 = NULL) {
  
  # ----------------------------------------------------------------------------
  # Quantities needed for the influences ---------------------------------------
  p                 <- length(beta)
  if ((length(x) != p) | (is.null(x))) {
    x               <- rep(0, p)
  } # if no covariate profile provided, use 0 as reference level
  exp.x             <- c(exp(x %*% beta))
  Pi.x              <- 1 - exp(- exp.x * Lambda0.Tau1Tau2)
  
  if (is.null(calibrated)) {
    calibrated <- FALSE
  }
  
  # ----------------------------------------------------------------------------
  # If calibrated weights ------------------------------------------------------
  
  if (calibrated == TRUE) {
    if (is.null(infl2.beta) | is.null(infl2.Lambda0.Tau1Tau2)) {
      print("Phase-two influences must be provided for the log-RH and CumBH")
      break}
    
    infl1.beta  <- infl.beta - infl2.beta
    infl1.Lambda0.Tau1Tau2 <- infl.Lambda0.Tau1Tau2 - infl2.Lambda0.Tau1Tau2
    
    infl1.Pi.x  <- c((1 - Pi.x) * exp.x * Lambda0.Tau1Tau2) * 
                    infl1.beta %*% x + c((1 - Pi.x) * exp.x) * 
                    infl1.Lambda0.Tau1Tau2
    infl2.Pi.x  <- c((1 - Pi.x) * exp.x * Lambda0.Tau1Tau2) * 
                    infl2.beta %*% x + c((1 - Pi.x) * exp.x) * 
                    infl2.Lambda0.Tau1Tau2
    infl.Pi.x   <- infl1.Pi.x + infl2.Pi.x
    
    return(list(infl.Pi.x.Tau1Tau2 = infl.Pi.x, infl2.Pi.x.Tau1Tau2 = infl2.Pi.x, 
                Pi.x.Tau1Tau2.hat = Pi.x))
  }else{
    # --------------------------------------------------------------------------
    # If design weights --------------------------------------------------------
    
    infl.Pi.x <- c((1 - Pi.x) * exp.x * Lambda0.Tau1Tau2) * infl.beta %*% x +
      c((1 - Pi.x) * exp.x) * infl.Lambda0.Tau1Tau2
    
    return(list(infl.Pi.x.Tau1Tau2 = infl.Pi.x, Pi.x.Tau1Tau2.hat = Pi.x))
  }
}


## -----------------------------------------------------------------------------
## Function: calibration()
## -----------------------------------------------------------------------------
## Description: This function calibrates the design weights using the raking 
##              procedure. It matches the weighted totals of the auxiliary 
##              variables in the stratified case-cohort (with calibrated 
##              weights), to the un-weighted auxiliary variables totals in the 
##              whole cohort. The Newton Raphson method is used to solve the 
##              optimization problem
## -----------------------------------------------------------------------------
## Arguments:
##
##  A               matrix with the values of the auxiliary variables to be used 
##                  for the calibration of the weights in the stratified 
##                  case-cohort (phase-two data)
##
##  design.weights  design weights to be calibrated
##
##  total           un-weighted auxiliary variables totals in the whole cohort
##
##  eta0            vector of initial/possible values for eta, to be used as 
##                  seed in the iterative procedure. Default is (0,...,0)
##
##  niter.max       maximum number of iterations for the Newton Raphson method. 
##                  Default is 10^4 iterations
##
##  epsilon.stop    threshold for the difference between the estimated weighted 
##                  total and the total in the whole cohort. If this difference 
##                  is less than the value of epsilon.stop, no more iterations 
##                  will be performed. Default is 10^(-10)
## -----------------------------------------------------------------------------

calibration <- function (A.phase2, design.weights, total, eta0 = NULL, 
                         niter.max = NULL, epsilon.stop = NULL) {
  
  A.phase2            <- as.matrix(A.phase2)
  
  n                   <- nrow(A.phase2)
  q                   <- ncol(A.phase2)
  
  if ((is.null(eta0)) || (length(eta0) != q)) {
    eta0          <- rep(0, q)
    print("eta0 has not been filled out. Default is 0. ")
  }
  if ((is.null(epsilon.stop)) || (epsilon.stop <= 0)) {
    epsilon.stop  <- 10^(-10)
    print("epsilon.stop has not been filled out. Default is 10^(-10). ")
  }
  if ((is.null(niter.max)) || (niter.max <= 0)) {
    niter.max     <- 10^4
    print("niter.max has not been filled out. Default is 10^4 iterations.")
  }
  
  for (niter in 1:niter.max) {
    calibrated.weights  <- design.weights * exp(A.phase2 %*% eta0)
    A.calweighted       <- sweep(A.phase2, 1, calibrated.weights, "*")
    f                   <- colSums(A.calweighted) - total
    
    epsilon             <- max(abs(f))
    if(epsilon < epsilon.stop){
      break
      }
    
    AA.calweighted      <- array(NA, dim = c(q, q, n))
    for (i in 1:n) {
      AA.calweighted[,, i] <- tcrossprod(A.phase2[i,], A.calweighted[i,])
      }
    drond.f.eta0        <- apply(X = AA.calweighted, MARGIN = c(1,2), FUN = sum)
    eta0                <- eta0 - c(f %*% solve(drond.f.eta0))
  }
  
  calibrated.weights    <- design.weights * exp(A.phase2 %*% eta0)
  A.calweighted         <- sweep(A.phase2, 1, calibrated.weights, "*")
  
  return(list(eta.hat = eta0, calibrated.weights = calibrated.weights, 
              estimated.total = colSums(A.calweighted)))
}


## -----------------------------------------------------------------------------
## Function: auxiliary.construction()
## -----------------------------------------------------------------------------
## Description: This function returns the auxiliary variables proposed by 
##              Breslow et al. (Stat. Biosci., 2009) and Breslow et al. (IMS,
##              2013). They are based on the relative hazard and cumulative 
##              baseline hazard influences, when using the cohort data with 
##              imputed covariates values
## -----------------------------------------------------------------------------
## Arguments:
##
##  mod             a cox model object, result of function coxph run on the 
##                  cohort data with imputed covariates values
##
##  Tau1            left bound of the time interval considered for the 
##                  cumulative baseline hazard. Default is the first event time
##
##  Tau2            right bound of the time interval considered for the 
##                  cumulative baseline hazard. Default is the last event time
## -----------------------------------------------------------------------------

auxiliary.construction <- function (mod, Tau1 = NULL, Tau2 = NULL) {
  
  # ----------------------------------------------------------------------------
  # Quantities needed for the influences ---------------------------------------
  
  mod.detail        <- coxph.detail(mod, riskmat = T)
  riskmat           <- mod.detail$riskmat 
  number.times      <- ncol(riskmat) 
  n                 <- nrow(riskmat) 
  X                 <- model.matrix(mod) 
  p                 <- ncol(X) 
  weights           <- mod$weights # weights used for the fit, omega_i
  if (is.null(weights)) {
    weights         <- 1
  } # when using the whole cohort, omega_i = 1. Should be the case here
  
  beta.hat          <- mod$coefficients
  
  exp.X.weighted    <- weights * exp(X %*% beta.hat)
  Y.exp.X.weighted  <- riskmat * matrix(exp.X.weighted, nrow = n, 
                                        ncol = number.times, byrow = FALSE) 
  S0t               <- t(riskmat) %*% (exp.X.weighted) 
  S1t               <- t(riskmat) %*% (X * matrix(exp.X.weighted, nrow = n,
                                                  ncol = p, byrow = FALSE)) 
  observed.times    <- apply(riskmat, 1, function(v) {which.max(cumsum(v))}) 
  dNt               <- matrix(0, n, number.times) 
  dNt[cbind(1:nrow(riskmat), observed.times)] <- 1 
  dNt               <- dNt * matrix(mod$y[, ncol(mod$y)], nrow(riskmat), 
                                    number.times, byrow = FALSE) 
  lambda0.t.hat     <- t(colSums(dNt) / S0t)
  if (is.null(Tau1)) {
    Tau1            <- floor(min(mod.detail$time))
  }
  if (is.null(Tau2)) {
    Tau2            <- floor(max(mod.detail$time))
  }
  Tau1Tau2.times  <- which((Tau1 < mod.detail$time) & (mod.detail$time <= Tau2)) 
  
  # ----------------------------------------------------------------------------
  # Computation of the influences for log-relative hazard, beta ----------------
  
  infl.beta         <- residuals(mod, type = "dfbeta", weighted = T)
  
  # ----------------------------------------------------------------------------
  # Computation of the influences for the baseline hazards at each unique 
  # event time and for the cumulative baseline hazard in time interval 
  # [Tau1, Tau2] ---------------------------------------------------------------
  
  infl.lambda0.t    <- (dNt - (Y.exp.X.weighted + infl.beta %*% t(S1t)) * 
                        matrix(lambda0.t.hat, nrow = n, ncol = number.times,
                                byrow = TRUE)) /
                        matrix(S0t, nrow = n, ncol = number.times,
                                byrow = TRUE) 
  infl.Lambda0.Tau1Tau2 <- rowSums(infl.lambda0.t[, Tau1Tau2.times])
 
  return(list(A.RH = infl.beta, A.CumBH = infl.Lambda0.Tau1Tau2))
}


## -----------------------------------------------------------------------------
## Function: product.covar.weight()
## -----------------------------------------------------------------------------
## Description: This function computes the product of the joint design weights
##              and joint sampling indicators covariances, needed for the 
##              phase-two component of the variance
## -----------------------------------------------------------------------------
## Argument:
##
##  stratified      was the sampling of the case-cohort stratified on W?
##
##  casecohort      if stratified = TRUE, data frame with status (case status), 
##                  W (strata), strata.m (number of sampled individuals in the
##                  stratum) and strata.n (stratum size), for each individual
##                  in the stratified case cohort data. If stratified = FALSE,
##                  data frame with status (case status), m (number of sampled
##                  individuals) and n (cohort size), for each individual in
##                  the un-stratified case cohort data
## -----------------------------------------------------------------------------

product.covar.weight <- function (casecohort, stratified = NULL) {
  
  if (is.null(stratified)) {
    stratified <- FALSE
  }
  
  if(stratified == TRUE){
    
    noncases              <- which(casecohort$status == 0)
    n.noncases            <- length(noncases)
    casecohort.W          <- casecohort$W[noncases]
    casecohort.m          <- casecohort$strata.m[noncases]
    casecohort.n          <- casecohort$strata.n[noncases]
    casecohort.m.n        <- casecohort.m / casecohort.n
    casecohort.n.m        <- casecohort.n / casecohort.m
    
    #---------------------------------------------------------------------------
    # Computation of the joint design weights, given by 1 over the sampling 
    # indicators probabilities -------------------------------------------------
    
    omega <- function (indiv) { 
      return((casecohort.W == casecohort.W[indiv]) * casecohort.n.m * 
               (casecohort.n - 1) / (casecohort.m - 1))
    }
    Omega                      <- do.call(rbind, lapply(1:n.noncases, omega))
    diag(Omega)                <- casecohort.n.m # omega_ij = 1 / P(Xi_i = Xi_j = 1), i \neq j, and the omega_ii = omega_i = 1 / P(Xi_i = 1)
    
    #---------------------------------------------------------------------------
    # Computation of the joint sampling indicators covariances -----------------
    
    covar                       <- function (indiv) {
      return((casecohort.W == casecohort.W[indiv]) * 
               (casecohort.m.n * (casecohort.m - 1) / (casecohort.n - 1) - 
                  casecohort.m.n^2))
    }
    Covar                       <- do.call(rbind, lapply(1:n.noncases, covar))
    diag(Covar)                 <- casecohort.m.n * (1 - casecohort.m.n)

    
  } else {
    
    noncases              <- which(casecohort$status == 0)
    n.noncases            <- length(noncases)
    casecohort.m          <- casecohort$m[noncases]
    casecohort.n          <- casecohort$n[noncases]
    casecohort.m.n        <- casecohort.m / casecohort.n
    casecohort.n.m        <- casecohort.n / casecohort.m
    
    #---------------------------------------------------------------------------
    # Computation of the joint design weights, given by 1 over the sampling 
    # indicators probabilities -------------------------------------------------
    
    Omega <- matrix(casecohort.n.m * (casecohort.n - 1) / (casecohort.m - 1),
                    nrow = n.noncases, ncol = n.noncases)
    diag(Omega) <- casecohort.n.m
    
    #---------------------------------------------------------------------------
    # Computation of the joint sampling indicators covariances -----------------
    
    Covar <- matrix(casecohort.m.n * (casecohort.m - 1) / (casecohort.n - 1) - 
                      casecohort.m.n ^ 2, nrow = n.noncases, ncol = n.noncases)
    diag(Covar) <- casecohort.m.n * (1 - casecohort.m.n)

  }
  
  #-----------------------------------------------------------------------------
  # Computation of the product -------------------------------------------------
  
  product.covar.weight  <- Covar * Omega # sigma_ij * omega_ij
  product.covar.weight[which(is.na(product.covar.weight))] <- 0
  
  return(product.covar.weight)
  
}

## -----------------------------------------------------------------------------
## Function: conf.interval()
## -----------------------------------------------------------------------------
## Description: This function computes the 95% confidence interval of a 
##              parameter estimate and its variance estimate, assuming normality
## -----------------------------------------------------------------------------
## Arguments:
##
##  param.estimate  parameter estimate
##                    
##  var.estimate    variance estimate for the parameter
## -----------------------------------------------------------------------------

conf.interval <- function (param.estimate, var.estimate) {
  
  return(c(param.estimate - sqrt(var.estimate) * qnorm(0.975), 
           param.estimate + sqrt(var.estimate) * qnorm(0.975)))
  
}


## -----------------------------------------------------------------------------
## Function: cohort.generation()
## -----------------------------------------------------------------------------
## Description: This function generates a cohort following the design used in 
##              the simulation study in Section 5 and Web Appendix E
## -----------------------------------------------------------------------------
## Arguments:
##
##  n                     cohort size
##
##  time                  study duration
##
##  bh                    time-fixed baseline hazard
##
##  alpha1                effect of covariate X1 on X3
##
##  alpha2                effect of covariate X2 on X3
##
##  beta1                 effect of covariate X1 on the outcome
##
##  beta2                 effect of covariate X2 on the outcome
##
##  beta3                 effect of covariate X3 on the outcome
##
##  Prob.x2               probabilities for the three possible values of 
##                        covariate X2, conditional on the value of covariate X1
##
##  strata.samplingfrac   stratum-specific sampling fractions
##
##  strata.m              fixed numbers of individuals to be sampled in the 
##                        strata
## -----------------------------------------------------------------------------

cohort.generation <- function (n, time = 10, bh, alpha1 = 0.05, alpha2 = -0.35, 
                               beta1 = - 0.2, beta2 = 0.25, beta3 = - 0.3,
                               Prob.x2 = matrix(c(0.7, 0.05, 0.25, 0.45, 0.2, 
                                                0.35, 0.40, 0.3, 0.3), nrow = 3, 
                                              byrow = TRUE), 
                               strata.samplingfrac = NULL,
                               strata.m = NULL) {
  
  # ----------------------------------------------------------------------------
  # Generation of the covariates -----------------------------------------------
  
  X1        <- rnorm(n)
  X1.cat    <- 1 * (X1 < -2) + 2 *((X1 >= -2)&(X1 < 1)) + 3 * (X1 > 1)
  X2        <- rep(NA, n)
  for (i in 1:length(x1.cat)) {
    X2[which(X1.cat == x1.cat[i])] <- sample(0:2, sum(X1.cat == x1.cat[i]), 
                                             replace = TRUE, prob = Prob.x2[i,])
  }
  mu3       <- alpha1 * X1 + alpha2 * X2
  X3        <- mapply(mu3, FUN = X.generation)
  
  # ----------------------------------------------------------------------------
  # Generation of the strata ---------------------------------------------------
  
  W <- 1 * ((X1 < 0)&(X2 < 2)) + 2 * ((X1 >= 0)&(X2 > 0)) +
    3 * ((X1 < 0)&(X2 == 2)) + 0 *((X1 >= 0)&(X2 == 0))
  
  # ----------------------------------------------------------------------------
  # Generation of the proxy variables for the covariates -----------------------
  
  X1.proxy  <- X1 + rnorm(n, sd = 0.75)
  X3.proxy  <- X3 + rnorm(n, sd = 0.75)
  eps       <- runif(n)
  X2.proxy  <- X2
  X2.proxy[eps < 0.1] <- 2*(X2[eps < 0.1] == 0) + 1*(X2[eps < 0.1] == 2) + 0*(X2[eps < 0.1] == 1)
  X2.proxy[eps > 0.9] <- 2*(X2[eps > 0.9] == 1) + 1*(X2[eps > 0.9] == 0) + 0*(X2[eps > 0.9] == 2)
  
  # ----------------------------------------------------------------------------
  # Generation of the entry times, event times, and times to independent 
  # censoring ------------------------------------------------------------------
  
  entry.times     <- runif(n, min = 0, max = 5) 
  
  scale.outcome   <- 1 / (bh * exp(beta1 * X1 + beta2 * X2 + beta3 * X3))
  outcome.times   <- rweibull(n, shape = 1, scale = scale.outcome)
  
  U               <- runif(n)
  death.times     <- (-log(U) / (0.02 / 10)) # 2% censoring by death at 10 years
  
  # Observed times and case indicators -----------------------------------------
  times           <- pmin(death.times, outcome.times, 10 - entry.times) 
  status          <- 1 * (times == outcome.times) 
  
  # Number of cases
  n.cases         <- sum(status == 1)
  
  if (!is.null(strata.samplingfrac)) {
    # Numbers of individuals, of cases, and of individuals to be sampled, in each
    # strata ---------------------------------------------------------------------
    strata.n.cases  <- c(sum(W == 0 & status == 1), sum(W == 1 & status == 1), 
                         sum(W == 2 & status == 1), sum(W == 3 & status == 1))
    strata.n        <- c(sum(W == 0), sum(W == 1), sum(W == 2), sum(W == 3))
    strata.m        <- floor(strata.n * strata.samplingfrac) # fixed number of individuals that will be sampled in each stratum to obtain a stratified case-cohort
    
    cohort <- cbind(id = 1:n, X1, X2, X3, W,X1.proxy, X2.proxy, X3.proxy, times, 
                    status, strata.n = strata.n[W+1], 
                    strata.m = strata.m[W+1], 
                    strata.n.cases = strata.n.cases[W+1], n.cases = n.cases)
    
    cohort <- as.data.frame(cohort)
    
    return(list(cohort = cohort, strata.n = strata.n, strata.m = strata.m, 
                strata.n.cases = strata.n.cases, n.cases = n.cases))
  }
  
  if (!is.null(strata.m)) {
    # Numbers of individuals, of cases, and of individuals to be sampled, in each
    # strata ---------------------------------------------------------------------
    strata.n.cases  <- c(sum(W == 0 & status == 1), sum(W == 1 & status == 1), 
                         sum(W == 2 & status == 1), sum(W == 3 & status == 1))
    strata.n        <- c(sum(W == 0), sum(W == 1), sum(W == 2), sum(W == 3))

    cohort <- cbind(id = 1:n, X1, X2, X3, W,X1.proxy, X2.proxy, X3.proxy, times, 
                    status, strata.n = strata.n[W+1], 
                    strata.m = strata.m[W+1], 
                    strata.n.cases = strata.n.cases[W+1], n.cases = n.cases)
    cohort <- as.data.frame(cohort)
    
    return(list(cohort = cohort, strata.n = strata.n, strata.m = strata.m, 
                strata.n.cases = strata.n.cases, n.cases = n.cases))
  }
  
  if(is.null(strata.m) & is.null(strata.samplingfrac)){
    # Numbers of individuals, of cases, and of individuals to be sampled, in each
    # strata ---------------------------------------------------------------------
    strata.n.cases  <- c(sum(W == 0 & status == 1), sum(W == 1 & status == 1), 
                         sum(W == 2 & status == 1), sum(W == 3 & status == 1))
    strata.n        <- c(sum(W == 0), sum(W == 1), sum(W == 2), sum(W == 3))
    
    cohort <- cbind(id = 1:n, X1, X2, X3, W,X1.proxy, X2.proxy, X3.proxy, times, 
                    status, strata.n = strata.n[W+1], 
                    strata.n.cases = strata.n.cases[W+1], n.cases = n.cases)
    
    cohort <- as.data.frame(cohort)
    
    return(list(cohort = cohort, strata.n = strata.n, 
                strata.n.cases = strata.n.cases, n.cases = n.cases))
  }
}


## -----------------------------------------------------------------------------
## Function: influences.missingdata()
## -----------------------------------------------------------------------------
## Description: This function estimates the influences of the log-relative 
##              hazard, baseline hazards at each unique event time, cumulative 
##              baseline hazard in the time interval [Tau1, Tau2] and pure risk 
##              in [Tau1, Tau2] and for covariate profile x, when covariate data
##              is missing for individuals in the case cohort. It returns the
##              overall influences, as well as the phase-two and phase-three 
##              influences. It also returns the parameters estimates
## -----------------------------------------------------------------------------
## Arguments:
##
##  mod                 a cox model object, result of function coxph 
##
##  Tau1                left bound of the time interval considered for the 
##                      cumulative baseline hazard and pure risk. Default is the 
##                      first event time
##
##  Tau2                right bound of the time interval considered for the 
##                      cumulative baseline hazard and pure risk. Default is the 
##                      last event time
##
##  x                   vector of length p, specifying the covariate profile 
##                      considered for the pure risk. Default is (0,...,0)
##
##  estimated.weights   are the weights for the third phase of sampling (due to 
##                      missingness) estimated ? If estimated.weights = TRUE, 
##                      B needs to be provided
##
##  B.phase2            matrix for the phase.two data with as many columns as 
##                      phase-three strata, with one 1 per row, to indicate the 
##                      phase-three stratum position of the whole cohort. This 
##                      argument needs to be provided if 
##                      estimated.weights = TRUE
##
##  riskmat.phase2      at risk matrix for the stratified case cohort (phase-two 
##                      data) at all of the cases event times, even that with 
##                      missing covariate data
##
##  dNt.phase2          counting process matrix for failures in the stratified
##                      case cohort (phase-two data). If status.phase2 = NULL, 
##                      this argument needs to be provided
##
##  status.phase2       vector indicating the case status in the stratified case 
##                      cohort (phase-two data). If dNt.phase2 = NULL, this 
##                      argument needs to be provided
## -----------------------------------------------------------------------------

influences.missingdata <- function (mod, riskmat.phase2, dNt.phase2 = NULL, 
                                    status.phase2 = NULL, 
                                    Tau1 = NULL, Tau2 = NULL, x = NULL, 
                                    estimated.weights = FALSE, B.phase2 = NULL) {
  
  if (is.null(estimated.weights)) {
    estimated.weights <- FALSE
  }
  
  if (is.null(dNt.phase2)&is.null(status.phase2)) {
    stop("The status for the phase-two data, or the counting process matrix 
             for failures in the phase-two data, need to be provided")
  } else {
    if (is.null(dNt.phase2)) { 
      # ------------------------------------------------------------------------
      # If not provided, compute the counting process matrix for failures in the
      # phase-two data ---------------------------------------------------------
      
      observed.times.phase2 <- apply(riskmat.phase2, 1,
                                     function(v) {which.max(cumsum(v))})
      dNt.phase2            <- matrix(0, nrow(riskmat.phase2), 
                                      ncol(riskmat.phase2))
      dNt.phase2[cbind(1:nrow(riskmat.phase2), observed.times.phase2)] <- 1
      dNt.phase2            <- sweep(dNt.phase2, 1, status.phase2, "*") 
    }
  }
  
  # ----------------------------------------------------------------------------
  # Quantities needed for the estimation ---------------------------------------
  
  mod.detail        <- coxph.detail(mod, riskmat = TRUE)
  riskmat           <- mod.detail$riskmat 
  number.times      <- ncol(riskmat)
  number.times.phase2 <- ncol(riskmat.phase2)
  eventtimes.phase2  <- colnames(riskmat.phase2)
  X                 <- model.matrix(mod)
  p                 <- ncol(X)
  n.phase3          <- nrow(X)
  n.phase2          <- nrow(riskmat.phase2)
  indiv.phase3      <- row.names(X)
  riskmat.phase3    <- riskmat.phase2[indiv.phase3, ] 
  beta.hat          <- mod$coefficients 
  weights           <- mod$weights # use the cases' actual failure times, known
  # even for cases with missing covariate data
  if (is.null(weights)) {
    weights <- 1
  } # when using the whole cohort, omega_i = 1
  exp.X.weighted    <- weights * exp(X %*% beta.hat)
  
  observed.times    <- apply(riskmat, 1, function(v) {which.max(cumsum(v))}) 
  dNt               <- matrix(0, n.phase3, number.times) 
  dNt[cbind(1:nrow(riskmat), observed.times)] <- 1 
  dNt               <- dNt * matrix(mod$y[, ncol(mod$y)], nrow(riskmat), 
                                    number.times, byrow = FALSE) 
  dNt.weighted      <- dNt * matrix(weights, nrow = n.phase3, 
                                    ncol = number.times,  byrow = FALSE) 
  infomat.indiv     <- mod.detail$imat  
  infomat           <- apply(X = infomat.indiv, MARGIN = c(1,2), FUN = sum)
  # can also be obtained from the following code:
  #XX                <- array(NA, dim = c(p, p, n.phase3))
  #for(i in 1:n.phase3){
  #  XX[,, i] <- tcrossprod(X[i,])
  #}
  #S2t               <- array(NA, dim = c(p, p, number.times))
  #S1S1t             <- array(NA, dim = c(p, p, number.times))
  #for(t in 1:number.times){ 
  #  S1S1t[,, t]     <- tcrossprod(S1t[t, ])
  #  S2t[,, t]       <- apply(X = sweep(XX, 3, riskmat[, t] * exp.X.weighted, "*"), MARGIN = c(1,2), FUN = sum)
  #}
  #infomat           <- rowSums((sweep(S2t, 3, colSums(dNt.weighted)/S0t, "*") - sweep(S1S1t, 3, colSums(dNt.weighted)/(S0t ^ 2), "*")), dims = 2) 
  
  Y.exp.X.weighted  <- sweep(riskmat, 1, exp.X.weighted, "*") 
  Y.exp.X.weighted.casestimes   <- sweep(riskmat.phase3, 1, exp.X.weighted, "*")
  
  S0t   <- t(riskmat) %*% (exp.X.weighted) 
  S1t   <- t(riskmat) %*% (sweep(X, 1, exp.X.weighted, "*")) 
  X.expect <- S1t / matrix(S0t, nrow = number.times, ncol = p, byrow = FALSE)
  S0t.casestimes  <- t(riskmat.phase3) %*% (exp.X.weighted)
  S1t.casestimes  <- t(riskmat.phase3) %*% (sweep(X, 1, exp.X.weighted, "*"))
  
  score         <- array(NA, dim = c(n.phase3, p, number.times))
  for(i in 1:n.phase3){
    score[i,,]  <- - sweep(t(X.expect), 1, X[i,], "-")
  }
  for(t in 1:number.times) {
    dMt         <- dNt.weighted[, t] - riskmat[, t] * exp.X.weighted * 
      (colSums(dNt.weighted) / S0t)[t]
    score[,, t] <- score[,, t] * matrix(dMt, nrow = n.phase3, ncol = p, 
                                        byrow = FALSE)
  }
  score.beta        <- apply(X = score, MARGIN = c(1,2), FUN = sum)
  
  if (is.null(Tau1)) {
    Tau1            <- floor(min(mod.detail$time))
  }
  if (is.null(Tau2)) {
    Tau2            <- floor(max(mod.detail$time))
  }
  Tau1Tau2.times    <- which((Tau1 < eventtimes.phase2) & 
                               (eventtimes.phase2 <= Tau2)) 
  Tau1Times         <- which((eventtimes.phase2 <= Tau1)) 
  Tau2Times         <- which((eventtimes.phase2 <= Tau2))
  
  lambda0.t.hat     <- t(colSums(dNt.phase2) / S0t.casestimes)
  
  Lambda0.Tau1Tau2.hat <- sum(lambda0.t.hat[Tau1Tau2.times])
  
  if ((length(x) != p) | (is.null(x))) {
    x               <- rep(0, p)
  } # if no covariate profile provided, use 0 as reference level
  exp.x             <- c(exp(x %*% beta.hat))
  
  Pi.x.hat          <- 1 - exp(- exp.x * Lambda0.Tau1Tau2.hat)
  
  # ----------------------------------------------------------------------------
  # If missingness weights estimated -------------------------------------------
  
  if (estimated.weights == TRUE) {
    if (is.null(B.phase2)) {
      print("Strata for the third phase of sampling (missingness) 
            must be provided")
      break
    }
    
    # --------------------------------------------------------------------------
    # Computation of the influences for gamma ----------------------------------
    
    J3              <- ncol(B.phase2)
    
    B.phase3        <- B.phase2[indiv.phase3, ]
    weights.p3.est  <- weights.phase3(B.phase3 = B.phase3, 
                                      total.phase2 = colSums(B.phase2), 
                                      gamma0 = rep(0, J3), niter.max = 10^3, 
                                      epsilon.stop = 10^(-10))$estimated.weights
    
    B.p3weights   <- sweep(B.phase3, 1, weights.p3.est, "*")
    BB.p3weights  <- array(NA, dim = c(J3, J3, n.phase3))
    
    for (i in 1:n.phase3) {
      BB.p3weights[,, i] <- tcrossprod(B.phase3[i,], B.p3weights[i,])
    }
    sum.BB.p3weights <- apply(X = BB.p3weights, MARGIN = c(1,2), FUN = sum)
    sum.BB.p3weights.inv <- solve(sum.BB.p3weights)
    
    infl2.gamma   <- B.phase2 %*% sum.BB.p3weights.inv
    infl3.gamma   <- 0 * B.phase2
    infl3.gamma[indiv.phase3,] <- - B.p3weights %*% sum.BB.p3weights.inv
    infl.gamma    <- infl2.gamma + infl3.gamma 
    
    # --------------------------------------------------------------------------
    # Computation of the influences for log-relative hazard, beta --------------
    
    XB <- array(NA, dim = c(p, J3, n.phase3))
    for (i in 1:n.phase3) {
      XB[,, i] = tcrossprod(X[i,], B.phase3[i,])
    }
    drond.G1t.gamma <- array(NA, dim = c(p, J3, number.times))
    drond.G0t.gamma <- array(NA, dim = c(J3, number.times))
    drondS1.gamma   <- array(NA, dim = c(p, J3, number.times))
    drond.S0t.gamma <- array(NA, dim = c(J3, number.times))
    for (t in 1:number.times) { 
      drond.G1t.gamma[,, t] <- apply(X = sweep(XB, 3, dNt.weighted[, t], "*"),
                                     MARGIN = c(1,2), FUN = sum)
      drond.G0t.gamma[, t]  <- colSums(B.phase3 * matrix(dNt.weighted[, t], 
                                                         nrow = n.phase3, 
                                                         ncol = J3, 
                                                         byrow = FALSE))
      drondS1.gamma[,, t]   <- apply(X = sweep(XB, 3,
                                               Y.exp.X.weighted[, t], "*"),
                                     MARGIN = c(1,2), FUN = sum)
      drond.S0t.gamma[, t]  <- colSums(B.phase3 * matrix(Y.exp.X.weighted[, t], 
                                                         nrow = n.phase3, 
                                                         ncol = J3, 
                                                         byrow = FALSE))
    }
    drond.Ut.gamma <- array(NA, dim = c(p, J3, number.times))
    for (t in 1:number.times) { 
      drond.Ut.gamma[,, t] <- drond.G1t.gamma[,, t] - 
        tcrossprod(X.expect[t, ], drond.G0t.gamma[, t]) - 
        (colSums(dNt.weighted) / S0t)[t] * drondS1.gamma[,, t] +
        (colSums(dNt.weighted) / S0t ^ 2)[t] * tcrossprod(S1t[t, ], 
                                                          drond.S0t.gamma[, t])
    }
    drond.U.gamma <- apply(X = drond.Ut.gamma, MARGIN = c(1,2), FUN = sum)
    
    infl2.beta    <- infl2.gamma %*% t(drond.U.gamma) %*% solve(infomat)
    
    infl3.beta        <- infl3.gamma %*% t(drond.U.gamma) %*% solve(infomat)
    infl3.beta[indiv.phase3, ] <- infl3.beta[indiv.phase3,] + score.beta %*% 
      solve(infomat) 
    
    infl.beta         <- infl2.beta + infl3.beta
    
    # --------------------------------------------------------------------------
    # Computation of the influences for the baseline hazards at each unique 
    # event time and for the cumulative baseline hazard in time interval 
    # [Tau1, Tau2] -------------------------------------------------------------
    
    drondS1.gamma.casestimes  <- array(NA, dim = c(p, J3, number.times.phase2))
    drond.S0t.gamma.casestimes  <- array(NA, dim = c(J3, number.times.phase2))
    for (t in 1:number.times.phase2) { 
      drondS1.gamma.casestimes[,, t]  <- apply(X = sweep(XB, 3, 
                                                         Y.exp.X.weighted.casestimes[, t], 
                                                         "*"), 
                                               MARGIN = c(1,2), FUN = sum)
      drond.S0t.gamma.casestimes[, t] <- colSums(B.phase3 * matrix(Y.exp.X.weighted.casestimes[, t], 
                                                                   nrow = n.phase3, 
                                                                   ncol = J3, 
                                                                   byrow = FALSE))
    }
    
    infl2.lambda0.t <- (dNt.phase2 - (infl2.beta %*% t(S1t.casestimes) + 
                                        infl2.gamma %*% drond.S0t.gamma.casestimes) * 
                          matrix(lambda0.t.hat, nrow = n.phase2, 
                                 ncol = number.times.phase2, byrow = TRUE)) / 
      matrix(S0t.casestimes, nrow = n.phase2, 
             ncol = number.times.phase2, byrow = TRUE)
    infl3.lambda0.t <- 0 * infl2.lambda0.t
    infl3.lambda0.t[indiv.phase3,] <- - (Y.exp.X.weighted.casestimes + 
                                           infl3.beta[indiv.phase3,] %*% 
                                           t(S1t.casestimes) + 
                                           infl3.gamma[indiv.phase3,] %*% 
                                           drond.S0t.gamma.casestimes) * 
      matrix(lambda0.t.hat,
             nrow = n.phase3, 
             ncol = number.times.phase2, 
             byrow = TRUE) / 
      matrix(S0t.casestimes, 
             nrow = n.phase3, 
             ncol = number.times.phase2, 
             byrow = TRUE)
    infl.lambda0.t    <- infl2.lambda0.t + infl3.lambda0.t
    
    infl2.Lambda0.Tau1Tau2 <- rowSums(infl2.lambda0.t[, Tau1Tau2.times])
    infl3.Lambda0.Tau1Tau2 <- rowSums(infl3.lambda0.t[, Tau1Tau2.times])
    infl.Lambda0.Tau1Tau2 <- infl2.Lambda0.Tau1Tau2 + 
      infl3.Lambda0.Tau1Tau2
    
    # --------------------------------------------------------------------------
    # Computation of the influences for the pure risk in [Tau1, Tau2] and for 
    # covariate profile x ------------------------------------------------------
    
    infl2.Pi.x  <- c((1 - Pi.x.hat) * exp.x * Lambda0.Tau1Tau2.hat) * 
      infl2.beta %*% x + c((1 - Pi.x.hat) * exp.x) * 
      infl2.Lambda0.Tau1Tau2
    infl3.Pi.x  <- c((1 - Pi.x.hat) * exp.x * Lambda0.Tau1Tau2.hat) * 
      infl3.beta %*% x + c((1 - Pi.x.hat) * exp.x) * 
      infl3.Lambda0.Tau1Tau2
    infl.Pi.x   <- infl2.Pi.x + infl3.Pi.x
    
    return(list(infl.beta = infl.beta, infl.lambda0.t = infl.lambda0.t, 
                infl.Lambda0.Tau1Tau2 = as.matrix(infl.Lambda0.Tau1Tau2), 
                infl.Pi.x.Tau1Tau2 = infl.Pi.x, infl2.beta = infl2.beta, 
                infl2.lambda0.t = infl2.lambda0.t, 
                infl2.Lambda0.Tau1Tau2 = as.matrix(infl2.Lambda0.Tau1Tau2), 
                infl2.Pi.x.Tau1Tau2 = infl2.Pi.x, infl3.beta = infl3.beta, 
                infl3.lambda0.t = infl3.lambda0.t, 
                infl3.Lambda0.Tau1Tau2 = as.matrix(infl3.Lambda0.Tau1Tau2), 
                infl3.Pi.x.Tau1Tau2 = infl3.Pi.x, beta.hat = beta.hat, 
                lambda0.t.hat = lambda0.t.hat, 
                Lambda0.Tau1Tau2.hat = Lambda0.Tau1Tau2.hat, 
                Pi.x.Tau1Tau2.hat = Pi.x.hat))
  } else {
    # --------------------------------------------------------------------------
    # If missingness weights known ---------------------------------------------
    
    # --------------------------------------------------------------------------
    # Computation of the influences for log-relative hazard, beta --------------
    
    infl2.beta <- matrix(0, nrow = nrow(riskmat.phase2), ncol = p)
    infl3.beta <- infl2.beta
    rownames(infl3.beta) <- rownames(riskmat.phase2)
    infl3.beta[indiv.phase3,] <- score.beta %*% solve(infomat) 
    infl.beta <- infl2.beta + infl3.beta
    
    # --------------------------------------------------------------------------
    # Computation of the influences for the baseline hazards at each unique 
    # event time and for the cumulative baseline hazard in time interval 
    # [Tau1, Tau2] -------------------------------------------------------------
    
    infl2.lambda0.t <- (dNt.phase2) / matrix(S0t.casestimes, 
                                             nrow = nrow(dNt.phase2), 
                                             ncol = number.times.phase2, 
                                             byrow = TRUE)
    
    infl3.lambda0.t    <- 0 * infl2.lambda0.t
    infl3.lambda0.t[indiv.phase3,] <- - (Y.exp.X.weighted.casestimes + 
                                           infl3.beta[indiv.phase3,] %*% 
                                           t(S1t.casestimes)) * 
      matrix(lambda0.t.hat, 
             nrow = n.phase3, 
             ncol = number.times.phase2, 
             byrow = TRUE) / 
      matrix(S0t.casestimes, 
             nrow = n.phase3, 
             ncol = number.times.phase2, 
             byrow = TRUE)
    infl.lambda0.t            <- infl2.lambda0.t + infl3.lambda0.t
    
    infl2.Lambda0.Tau1Tau2     <- rowSums(infl2.lambda0.t[, Tau1Tau2.times])
    infl3.Lambda0.Tau1Tau2     <- rowSums(infl3.lambda0.t[, Tau1Tau2.times])
    infl.Lambda0.Tau1Tau2      <- infl2.Lambda0.Tau1Tau2 + infl3.Lambda0.Tau1Tau2
    
    # --------------------------------------------------------------------------
    # Computation of the influences for the pure risk in [Tau1, Tau2] and for 
    # covariate profile x ------------------------------------------------------
    
    infl2.Pi.x  <- c((1 - Pi.x.hat) * exp.x * Lambda0.Tau1Tau2.hat) * 
      infl2.beta %*% x + c((1 - Pi.x.hat) * exp.x) * 
      infl2.Lambda0.Tau1Tau2
    infl3.Pi.x  <- c((1 - Pi.x.hat) * exp.x * Lambda0.Tau1Tau2.hat) * 
      infl3.beta %*% x + c((1 - Pi.x.hat) * exp.x) * 
      infl3.Lambda0.Tau1Tau2
    infl.Pi.x   <- infl2.Pi.x + infl3.Pi.x
    
    return(list(infl.beta = infl.beta, infl.lambda0.t = infl.lambda0.t, 
                infl.Lambda0.Tau1Tau2 = as.matrix(infl.Lambda0.Tau1Tau2), 
                infl.Pi.x.Tau1Tau2 = infl.Pi.x, infl2.beta = infl2.beta, 
                infl2.lambda0.t = infl2.lambda0.t, 
                infl2.Lambda0.Tau1Tau2 = as.matrix(infl2.Lambda0.Tau1Tau2), 
                infl2.Pi.x.Tau1Tau2 = infl2.Pi.x, infl3.beta = infl3.beta, 
                infl3.lambda0.t = infl3.lambda0.t, 
                infl3.Lambda0.Tau1Tau2 = as.matrix(infl3.Lambda0.Tau1Tau2), 
                infl3.Pi.x.Tau1Tau2 = infl3.Pi.x, beta.hat = beta.hat, 
                lambda0.t.hat = lambda0.t.hat, 
                Lambda0.Tau1Tau2.hat = Lambda0.Tau1Tau2.hat, 
                Pi.x.Tau1Tau2.hat = Pi.x.hat))
  }
}


## -----------------------------------------------------------------------------
## Function: influences.RH.missingdata()
## -----------------------------------------------------------------------------
## Description: This function estimates the influences of the log-relative 
##              hazard, when covariate data is missing for individuals in the 
##              case cohort. It returns the overall influences, as well as the
##              phase-two and phase-three influences and parameter estimate
## -----------------------------------------------------------------------------
## Arguments:
##
##  mod                 a cox model object, result of function coxph 
##
##  estimated.weights   are the weights for the third phase of sampling (due to 
##                      missingness) estimated ? If estimated.weights = TRUE, 
##                      B needs to be provided
##
##  B.phase2            matrix for the phase.two data with as many columns as 
##                      phase-three strata, with one 1 per row, to indicate the 
##                      phase-three stratum position of the whole cohort. This 
##                      argument needs to be provided if 
##                      estimated.weights = TRUE
##
##  riskmat.phase2      at risk matrix for the stratified case cohort (phase-two 
##                      data) at all of the cases event times, even that with 
##                      missing covariate data
##
##  dNt.phase2          counting process matrix for failures in the stratified
##                      case cohort (phase-two data). If status.phase2 = NULL, 
##                      this argument needs to be provided
##
##  status.phase2       vector indicating the case status in the stratified case 
##                      cohort (phase-two data). If dNt.phase2 = NULL, this 
##                      argument needs to be provided
## -----------------------------------------------------------------------------

influences.RH.missingdata <- function (mod, riskmat.phase2, dNt.phase2 = NULL, 
                                       status.phase2 = NULL, 
                                       estimated.weights = FALSE, 
                                       B.phase2 = NULL) {
  
  if (is.null(estimated.weights)) {
    estimated.weights <- FALSE
  }
  
  if (is.null(dNt.phase2)&is.null(status.phase2)) {
    stop("The status for the phase-two data, or the counting process matrix 
             for failures in the phase-two data, need to be provided")
  } else {
    if (is.null(dNt.phase2)) { 
      # ------------------------------------------------------------------------
      # If not provided, compute the counting process matrix for failures in the
      # phase-two data ---------------------------------------------------------
      
      observed.times.phase2 <- apply(riskmat.phase2, 1,
                                     function(v) {which.max(cumsum(v))})
      dNt.phase2            <- matrix(0, nrow(riskmat.phase2), 
                                      ncol(riskmat.phase2))
      dNt.phase2[cbind(1:nrow(riskmat.phase2), observed.times.phase2)] <- 1
      dNt.phase2            <- sweep(dNt.phase2, 1, status.phase2, "*") 
    }
  }
  
  # ----------------------------------------------------------------------------
  # Quantities needed for the estimation ---------------------------------------
  
  mod.detail        <- coxph.detail(mod, riskmat = TRUE)
  riskmat           <- mod.detail$riskmat 
  number.times      <- ncol(riskmat)
  number.times.phase2 <- ncol(riskmat.phase2)
  X                 <- model.matrix(mod)
  p                 <- ncol(X)
  n.phase3          <- nrow(X)
  n.phase2          <- nrow(riskmat.phase2)
  indiv.phase3      <- row.names(X)
  riskmat.phase3    <- riskmat.phase2[indiv.phase3, ] 
  beta.hat          <- mod$coefficients 
  weights           <- mod$weights # use the cases' actual failure times, known
  # even for cases with missing covariate data
  if (is.null(weights)) {
    weights <- 1
  } # when using the whole cohort, omega_i = 1
  exp.X.weighted    <- weights * exp(X %*% beta.hat)
  
  observed.times    <- apply(riskmat, 1, function(v) {which.max(cumsum(v))}) 
  dNt               <- matrix(0, n.phase3, number.times) 
  dNt[cbind(1:nrow(riskmat), observed.times)] <- 1 
  dNt               <- dNt * matrix(mod$y[, ncol(mod$y)], nrow(riskmat), 
                                    number.times, byrow = FALSE) 
  dNt.weighted      <- dNt * matrix(weights, nrow = n.phase3, 
                                    ncol = number.times,  byrow = FALSE) 
  infomat.indiv     <- mod.detail$imat  
  infomat           <- apply(X = infomat.indiv, MARGIN = c(1,2), FUN = sum)
  
  indiv.phase3      <- row.names(X)
  riskmat.phase3    <- riskmat.phase2[indiv.phase3, ] 
  
  Y.exp.X.weighted  <- sweep(riskmat, 1, exp.X.weighted, "*") 
  Y.exp.X.weighted.casestimes   <- sweep(riskmat.phase3, 1, exp.X.weighted, "*")
  
  S0t   <- t(riskmat) %*% (exp.X.weighted) 
  S1t   <- t(riskmat) %*% (sweep(X, 1, exp.X.weighted, "*")) 
  X.expect <- S1t / matrix(S0t, nrow = number.times, ncol = p, byrow = FALSE)
  S0t.casestimes  <- t(riskmat.phase3) %*% (exp.X.weighted)
  S1t.casestimes  <- t(riskmat.phase3) %*% (sweep(X, 1, exp.X.weighted, "*"))
  
  score         <- array(NA, dim = c(n.phase3, p, number.times))
  for(i in 1:n.phase3){
    score[i,,]  <- - sweep(t(X.expect), 1, X[i,], "-")
  }
  for(t in 1:number.times) {
    dMt         <- dNt.weighted[, t] - riskmat[, t] * exp.X.weighted * 
      (colSums(dNt.weighted) / S0t)[t]
    score[,, t] <- score[,, t] * matrix(dMt, nrow = n.phase3, ncol = p, 
                                        byrow = FALSE)
  }
  score.beta        <- apply(X = score, MARGIN = c(1,2), FUN = sum)
  
  # ----------------------------------------------------------------------------
  # If missingness weights estimated -------------------------------------------
  
  if (estimated.weights == TRUE) {
    if (is.null(B.phase2)) {
      print("Strata for the third phase of sampling (missingness) 
            must be provided")
      break
    }

    # --------------------------------------------------------------------------
    # Computation of the influences for gamma ----------------------------------
    
    J3              <- ncol(B.phase2)
    
    B.phase3        <- B.phase2[indiv.phase3, ]
    weights.p3.est  <- weights.phase3(B.phase3 = B.phase3, 
                                      total.phase2 = colSums(B.phase2), 
                                      gamma0 = rep(0, J3), niter.max = 10^3, 
                                      epsilon.stop = 10^(-10))$estimated.weights
    
    B.p3weights   <- sweep(B.phase3, 1, weights.p3.est, "*")
    BB.p3weights  <- array(NA, dim = c(J3, J3, n.phase3))
    
    for (i in 1:n.phase3) {
      BB.p3weights[,, i] <- tcrossprod(B.phase3[i,], B.p3weights[i,])
    }
    sum.BB.p3weights <- apply(X = BB.p3weights, MARGIN = c(1,2), FUN = sum)
    sum.BB.p3weights.inv <- solve(sum.BB.p3weights)
    
    infl2.gamma   <- B.phase2 %*% sum.BB.p3weights.inv
    infl3.gamma   <- 0 * B.phase2
    infl3.gamma[indiv.phase3,] <- - B.p3weights %*% sum.BB.p3weights.inv
    infl.gamma    <- infl2.gamma + infl3.gamma 
    
    # --------------------------------------------------------------------------
    # Computation of the influences for log-relative hazard, beta --------------
    
    XB <- array(NA, dim = c(p, J3, n.phase3))
    for (i in 1:n.phase3) {
      XB[,, i] = tcrossprod(X[i,], B.phase3[i,])
    }
    drond.G1t.gamma <- array(NA, dim = c(p, J3, number.times))
    drond.G0t.gamma <- array(NA, dim = c(J3, number.times))
    drondS1.gamma   <- array(NA, dim = c(p, J3, number.times))
    drond.S0t.gamma <- array(NA, dim = c(J3, number.times))
    for (t in 1:number.times) { 
      drond.G1t.gamma[,, t] <- apply(X = sweep(XB, 3, dNt.weighted[, t], "*"),
                                     MARGIN = c(1,2), FUN = sum)
      drond.G0t.gamma[, t]  <- colSums(B.phase3 * matrix(dNt.weighted[, t], 
                                                         nrow = n.phase3, 
                                                         ncol = J3, 
                                                         byrow = FALSE))
      drondS1.gamma[,, t]   <- apply(X = sweep(XB, 3,
                                               Y.exp.X.weighted[, t], "*"),
                                     MARGIN = c(1,2), FUN = sum)
      drond.S0t.gamma[, t]  <- colSums(B.phase3 * matrix(Y.exp.X.weighted[, t], 
                                                         nrow = n.phase3, 
                                                         ncol = J3, 
                                                         byrow = FALSE))
    }
    drond.Ut.gamma <- array(NA, dim = c(p, J3, number.times))
    for (t in 1:number.times) { 
      drond.Ut.gamma[,, t] <- drond.G1t.gamma[,, t] - 
        tcrossprod(X.expect[t, ], drond.G0t.gamma[, t]) - 
        (colSums(dNt.weighted) / S0t)[t] * drondS1.gamma[,, t] +
        (colSums(dNt.weighted) / S0t ^ 2)[t] * tcrossprod(S1t[t, ], 
                                                          drond.S0t.gamma[, t])
    }
    drond.U.gamma <- apply(X = drond.Ut.gamma, MARGIN = c(1,2), FUN = sum)
    
    infl2.beta    <- infl2.gamma %*% t(drond.U.gamma) %*% solve(infomat)
    
    infl3.beta        <- infl3.gamma %*% t(drond.U.gamma) %*% solve(infomat)
    infl3.beta[indiv.phase3, ] <- infl3.beta[indiv.phase3,] + score.beta %*% 
      solve(infomat) 
    
    infl.beta         <- infl2.beta + infl3.beta
    
    return(list(infl.beta = infl.beta, infl2.beta = infl2.beta, 
                infl3.beta = infl3.beta, beta.hat = beta.hat))
  } else {
    # --------------------------------------------------------------------------
    # If missingness weights known ---------------------------------------------
    
    # --------------------------------------------------------------------------
    # Computation of the influences for log-relative hazard, beta --------------
    
    infl2.beta <- matrix(0, nrow = nrow(riskmat.phase2), ncol = p)
    infl3.beta <- infl2.beta
    rownames(infl3.beta) <- rownames(riskmat.phase2)
    infl3.beta[indiv.phase3,] <- score.beta %*% solve(infomat) 
    infl.beta <- infl2.beta + infl3.beta
    
    return(list(infl.beta = infl.beta, infl2.beta = infl2.beta, 
                infl3.beta = infl3.beta, beta.hat = beta.hat))
  }
}


## -----------------------------------------------------------------------------
## Function: influences.CumBH.missingdata()
## -----------------------------------------------------------------------------
## Description: This function estimates the influences of the log-relative 
##              hazard, baseline hazards at each unique event time and 
##              cumulative baseline hazard in the time interval [Tau1, Tau2], 
##              when covariate data is missing for individuals in the case 
##              cohort. It returns the overall influences, as well as the 
##              phase-two and phase-three influences. It also returns the 
##              parameters estimates
## -----------------------------------------------------------------------------
## Arguments:
##
##  mod                 a cox model object, result of function coxph 
##
##  Tau1                left bound of the time interval considered for the 
##                      cumulative baseline hazard. Default is the first event
##                      time
##
##  Tau2                right bound of the time interval considered for the 
##                      cumulative baseline hazard. Default is the last event 
##                      time
##
##  estimated.weights   are the weights for the third phase of sampling (due to 
##                      missingness) estimated ? If estimated.weights = TRUE, 
##                      B needs to be provided
##
##  B.phase2            matrix for the phase.two data with as many columns as 
##                      phase-three strata, with one 1 per row, to indicate the 
##                      phase-three stratum position of the whole cohort. This 
##                      argument needs to be provided if 
##                      estimated.weights = TRUE
##
##  riskmat.phase2      at risk matrix for the stratified case cohort (phase-two 
##                      data) at all of the cases event times, even that with 
##                      missing covariate data
##
##  dNt.phase2          counting process matrix for failures in the stratified
##                      case cohort (phase-two data). If status.phase2 = NULL, 
##                      this argument needs to be provided
##
##  status.phase2       vector indicating the case status in the stratified case 
##                      cohort (phase-two data). If dNt.phase2 = NULL, this 
##                      argument needs to be provided
## -----------------------------------------------------------------------------

influences.CumBH.missingdata <- function (mod, riskmat.phase2, 
                                          dNt.phase2 = NULL, 
                                          status.phase2 = NULL, Tau1 = NULL, 
                                          Tau2 = NULL, x = NULL, 
                                          estimated.weights = FALSE, 
                                          B.phase2 = NULL) {
  
  if (is.null(estimated.weights)) {
    estimated.weights <- FALSE
  }
  
  if (is.null(dNt.phase2)&is.null(status.phase2)) {
    stop("The status for the phase-two data, or the counting process matrix 
             for failures in the phase-two data, need to be provided")
  } else {
    if (is.null(dNt.phase2)) { 
      # ------------------------------------------------------------------------
      # If not provided, compute the counting process matrix for failures in the
      # phase-two data ---------------------------------------------------------
      
      observed.times.phase2 <- apply(riskmat.phase2, 1,
                                     function(v) {which.max(cumsum(v))})
      dNt.phase2            <- matrix(0, nrow(riskmat.phase2), 
                                      ncol(riskmat.phase2))
      dNt.phase2[cbind(1:nrow(riskmat.phase2), observed.times.phase2)] <- 1
      dNt.phase2            <- sweep(dNt.phase2, 1, status.phase2, "*") 
    }
  }
  
  # ----------------------------------------------------------------------------
  # Quantities needed for the estimation ---------------------------------------
  
  mod.detail        <- coxph.detail(mod, riskmat = TRUE)
  riskmat           <- mod.detail$riskmat 
  number.times      <- ncol(riskmat)
  number.times.phase2 <- ncol(riskmat.phase2)
  eventtimes.phase2  <- colnames(riskmat.phase2)
  X                 <- model.matrix(mod)
  p                 <- ncol(X)
  n.phase3          <- nrow(X)
  n.phase2          <- nrow(riskmat.phase2)
  indiv.phase3      <- row.names(X)
  riskmat.phase3    <- riskmat.phase2[indiv.phase3, ] 
  beta.hat          <- mod$coefficients 
  weights           <- mod$weights # use the cases' actual failure times, known
  # even for cases with missing covariate data
  if (is.null(weights)) {
    weights <- 1
  } # when using the whole cohort, omega_i = 1
  exp.X.weighted    <- weights * exp(X %*% beta.hat)
  
  observed.times    <- apply(riskmat, 1, function(v) {which.max(cumsum(v))}) 
  dNt               <- matrix(0, n.phase3, number.times) 
  dNt[cbind(1:nrow(riskmat), observed.times)] <- 1 
  dNt               <- dNt * matrix(mod$y[, ncol(mod$y)], nrow(riskmat), 
                                    number.times, byrow = FALSE) 
  dNt.weighted      <- dNt * matrix(weights, nrow = n.phase3, 
                                    ncol = number.times,  byrow = FALSE) 
  infomat.indiv     <- mod.detail$imat  
  infomat           <- apply(X = infomat.indiv, MARGIN = c(1,2), FUN = sum)
  
  indiv.phase3      <- row.names(X)
  riskmat.phase3    <- riskmat.phase2[indiv.phase3, ] 
  
  Y.exp.X.weighted  <- sweep(riskmat, 1, exp.X.weighted, "*") 
  Y.exp.X.weighted.casestimes   <- sweep(riskmat.phase3, 1, exp.X.weighted, "*")
  
  S0t   <- t(riskmat) %*% (exp.X.weighted) 
  S1t   <- t(riskmat) %*% (sweep(X, 1, exp.X.weighted, "*")) 
  X.expect <- S1t / matrix(S0t, nrow = number.times, ncol = p, byrow = FALSE)
  S0t.casestimes  <- t(riskmat.phase3) %*% (exp.X.weighted)
  S1t.casestimes  <- t(riskmat.phase3) %*% (sweep(X, 1, exp.X.weighted, "*"))
  
  score         <- array(NA, dim = c(n.phase3, p, number.times))
  for(i in 1:n.phase3){
    score[i,,]  <- - sweep(t(X.expect), 1, X[i,], "-")
  }
  for(t in 1:number.times) {
    dMt         <- dNt.weighted[, t] - riskmat[, t] * exp.X.weighted * 
      (colSums(dNt.weighted) / S0t)[t]
    score[,, t] <- score[,, t] * matrix(dMt, nrow = n.phase3, ncol = p, 
                                        byrow = FALSE)
  }
  score.beta        <- apply(X = score, MARGIN = c(1,2), FUN = sum)
  
  if (is.null(Tau1)) {
    Tau1            <- floor(min(mod.detail$time))
  }
  if (is.null(Tau2)) {
    Tau2            <- floor(max(mod.detail$time))
  }
  Tau1Tau2.times    <- which((Tau1 < eventtimes.phase2) & 
                               (eventtimes.phase2 <= Tau2)) 
  Tau1Times         <- which((eventtimes.phase2 <= Tau1)) 
  Tau2Times         <- which((eventtimes.phase2 <= Tau2))
  
  lambda0.t.hat     <- t(colSums(dNt.phase2) / S0t.casestimes)
  
  Lambda0.Tau1Tau2.hat <- sum(lambda0.t.hat[Tau1Tau2.times])
  
  # ----------------------------------------------------------------------------
  # If missingness weights estimated -------------------------------------------
  
  if (estimated.weights == TRUE) {
    if (is.null(B.phase2)) {
      print("Strata for the third phase of sampling (missingness) 
            must be provided")
      break
    }
    
    # --------------------------------------------------------------------------
    # Computation of the influences for gamma ----------------------------------
    
    J3              <- ncol(B.phase2)
    
    B.phase3        <- B.phase2[indiv.phase3, ]
    weights.p3.est  <- weights.phase3(B.phase3 = B.phase3, 
                                      total.phase2 = colSums(B.phase2), 
                                      gamma0 = rep(0, J3), niter.max = 10^3, 
                                      epsilon.stop = 10^(-10))$estimated.weights
    
    B.p3weights   <- sweep(B.phase3, 1, weights.p3.est, "*")
    BB.p3weights  <- array(NA, dim = c(J3, J3, n.phase3))
    
    for (i in 1:n.phase3) {
      BB.p3weights[,, i] <- tcrossprod(B.phase3[i,], B.p3weights[i,])
    }
    sum.BB.p3weights <- apply(X = BB.p3weights, MARGIN = c(1,2), FUN = sum)
    sum.BB.p3weights.inv <- solve(sum.BB.p3weights)
    
    infl2.gamma   <- B.phase2 %*% sum.BB.p3weights.inv
    infl3.gamma   <- 0 * B.phase2
    infl3.gamma[indiv.phase3,] <- - B.p3weights %*% sum.BB.p3weights.inv
    infl.gamma    <- infl2.gamma + infl3.gamma 
    
    # --------------------------------------------------------------------------
    # Computation of the influences for log-relative hazard, beta --------------
    
    XB <- array(NA, dim = c(p, J3, n.phase3))
    for (i in 1:n.phase3) {
      XB[,, i] = tcrossprod(X[i,], B.phase3[i,])
    }
    drond.G1t.gamma <- array(NA, dim = c(p, J3, number.times))
    drond.G0t.gamma <- array(NA, dim = c(J3, number.times))
    drondS1.gamma   <- array(NA, dim = c(p, J3, number.times))
    drond.S0t.gamma <- array(NA, dim = c(J3, number.times))
    for (t in 1:number.times) { 
      drond.G1t.gamma[,, t] <- apply(X = sweep(XB, 3, dNt.weighted[, t], "*"),
                                     MARGIN = c(1,2), FUN = sum)
      drond.G0t.gamma[, t]  <- colSums(B.phase3 * matrix(dNt.weighted[, t], 
                                                         nrow = n.phase3, 
                                                         ncol = J3, 
                                                         byrow = FALSE))
      drondS1.gamma[,, t]   <- apply(X = sweep(XB, 3,
                                               Y.exp.X.weighted[, t], "*"),
                                     MARGIN = c(1,2), FUN = sum)
      drond.S0t.gamma[, t]  <- colSums(B.phase3 * matrix(Y.exp.X.weighted[, t], 
                                                         nrow = n.phase3, 
                                                         ncol = J3, 
                                                         byrow = FALSE))
    }
    drond.Ut.gamma <- array(NA, dim = c(p, J3, number.times))
    for (t in 1:number.times) { 
      drond.Ut.gamma[,, t] <- drond.G1t.gamma[,, t] - 
        tcrossprod(X.expect[t, ], drond.G0t.gamma[, t]) - 
        (colSums(dNt.weighted) / S0t)[t] * drondS1.gamma[,, t] +
        (colSums(dNt.weighted) / S0t ^ 2)[t] * tcrossprod(S1t[t, ], 
                                                          drond.S0t.gamma[, t])
    }
    drond.U.gamma <- apply(X = drond.Ut.gamma, MARGIN = c(1,2), FUN = sum)
    
    infl2.beta    <- infl2.gamma %*% t(drond.U.gamma) %*% solve(infomat)
    
    infl3.beta        <- infl3.gamma %*% t(drond.U.gamma) %*% solve(infomat)
    infl3.beta[indiv.phase3, ] <- infl3.beta[indiv.phase3,] + score.beta %*% 
      solve(infomat) 
    
    infl.beta         <- infl2.beta + infl3.beta
    
    # --------------------------------------------------------------------------
    # Computation of the influences for the baseline hazards at each unique 
    # event time and for the cumulative baseline hazard in time interval 
    # [Tau1, Tau2] -------------------------------------------------------------
    
    drondS1.gamma.casestimes  <- array(NA, dim = c(p, J3, number.times.phase2))
    drond.S0t.gamma.casestimes  <- array(NA, dim = c(J3, number.times.phase2))
    for (t in 1:number.times.phase2) { 
      drondS1.gamma.casestimes[,, t]  <- apply(X = sweep(XB, 3, 
                                                         Y.exp.X.weighted.casestimes[, t], 
                                                         "*"), 
                                               MARGIN = c(1,2), FUN = sum)
      drond.S0t.gamma.casestimes[, t] <- colSums(B.phase3 * matrix(Y.exp.X.weighted.casestimes[, t], 
                                                                   nrow = n.phase3, 
                                                                   ncol = J3, 
                                                                   byrow = FALSE))
    }
    
    infl2.lambda0.t <- (dNt.phase2 - (infl2.beta %*% t(S1t.casestimes) + 
                                        infl2.gamma %*% drond.S0t.gamma.casestimes) * 
                          matrix(lambda0.t.hat, nrow = n.phase2, 
                                 ncol = number.times.phase2, byrow = TRUE)) / 
      matrix(S0t.casestimes, nrow = n.phase2, 
             ncol = number.times.phase2, byrow = TRUE)
    infl3.lambda0.t <- 0 * infl2.lambda0.t
    infl3.lambda0.t[indiv.phase3,] <- - (Y.exp.X.weighted.casestimes + 
                                           infl3.beta[indiv.phase3,] %*% 
                                           t(S1t.casestimes) + 
                                           infl3.gamma[indiv.phase3,] %*% 
                                           drond.S0t.gamma.casestimes) * 
      matrix(lambda0.t.hat,
             nrow = n.phase3, 
             ncol = number.times.phase2, 
             byrow = TRUE) / 
      matrix(S0t.casestimes, 
             nrow = n.phase3, 
             ncol = number.times.phase2, 
             byrow = TRUE)
    infl.lambda0.t    <- infl2.lambda0.t + infl3.lambda0.t
    
    infl2.Lambda0.Tau1Tau2 <- rowSums(infl2.lambda0.t[, Tau1Tau2.times])
    infl3.Lambda0.Tau1Tau2 <- rowSums(infl3.lambda0.t[, Tau1Tau2.times])
    infl.Lambda0.Tau1Tau2 <- infl2.Lambda0.Tau1Tau2 + 
      infl3.Lambda0.Tau1Tau2
    
    return(list(infl.beta = infl.beta, infl.lambda0.t = infl.lambda0.t, 
                infl.Lambda0.Tau1Tau2 = as.matrix(infl.Lambda0.Tau1Tau2), 
                infl2.beta = infl2.beta, infl2.lambda0.t = infl2.lambda0.t, 
                infl2.Lambda0.Tau1Tau2 = as.matrix(infl2.Lambda0.Tau1Tau2), 
                infl3.beta = infl3.beta, infl3.lambda0.t = infl3.lambda0.t, 
                infl3.Lambda0.Tau1Tau2 = as.matrix(infl3.Lambda0.Tau1Tau2), 
                beta.hat = beta.hat, lambda0.t.hat = lambda0.t.hat, 
                Lambda0.Tau1Tau2.hat = Lambda0.Tau1Tau2.hat))
  } else {
    # --------------------------------------------------------------------------
    # If missingness weights known ---------------------------------------------
    
    # --------------------------------------------------------------------------
    # Computation of the influences for log-relative hazard, beta --------------
    
    infl2.beta <- matrix(0, nrow = nrow(riskmat.phase2), ncol = p)
    infl3.beta <- infl2.beta
    rownames(infl3.beta) <- rownames(riskmat.phase2)
    infl3.beta[indiv.phase3,] <- score.beta %*% solve(infomat) 
    infl.beta <- infl2.beta + infl3.beta
    
    # --------------------------------------------------------------------------
    # Computation of the influences for the baseline hazards at each unique 
    # event time and for the cumulative baseline hazard in time interval 
    # [Tau1, Tau2] -------------------------------------------------------------
    
    infl2.lambda0.t <- (dNt.phase2) / matrix(S0t.casestimes, 
                                             nrow = nrow(dNt.phase2), 
                                             ncol = number.times.phase2, 
                                             byrow = TRUE)
    
    infl3.lambda0.t    <- 0 * infl2.lambda0.t
    infl3.lambda0.t[indiv.phase3,] <- - (Y.exp.X.weighted.casestimes + 
                                           infl3.beta[indiv.phase3,] %*% 
                                           t(S1t.casestimes)) * 
      matrix(lambda0.t.hat, 
             nrow = n.phase3, 
             ncol = number.times.phase2, 
             byrow = TRUE) / 
      matrix(S0t.casestimes, 
             nrow = n.phase3, 
             ncol = number.times.phase2, 
             byrow = TRUE)
    infl.lambda0.t            <- infl2.lambda0.t + infl3.lambda0.t
    
    infl2.Lambda0.Tau1Tau2    <- rowSums(infl2.lambda0.t[, Tau1Tau2.times])
    infl3.Lambda0.Tau1Tau2    <- rowSums(infl3.lambda0.t[, Tau1Tau2.times])
    infl.Lambda0.Tau1Tau2     <- infl2.Lambda0.Tau1Tau2 + infl3.Lambda0.Tau1Tau2
    
    return(list(infl.beta = infl.beta, infl.lambda0.t = infl.lambda0.t, 
                infl.Lambda0.Tau1Tau2 = as.matrix(infl.Lambda0.Tau1Tau2), 
                infl2.beta = infl2.beta, infl2.lambda0.t = infl2.lambda0.t, 
                infl2.Lambda0.Tau1Tau2 = as.matrix(infl2.Lambda0.Tau1Tau2), 
                infl3.beta = infl3.beta, infl3.lambda0.t = infl3.lambda0.t, 
                infl3.Lambda0.Tau1Tau2 = as.matrix(infl3.Lambda0.Tau1Tau2), 
                beta.hat = beta.hat, lambda0.t.hat = lambda0.t.hat, 
                Lambda0.Tau1Tau2.hat = Lambda0.Tau1Tau2.hat))
  }
}


## -----------------------------------------------------------------------------
## Function: influences.PR.missingdata()
## -----------------------------------------------------------------------------
## Description: This function estimates the influences of the pure risk in the 
##              time interval [Tau1, Tau2] and for covariate profile x, from
##              that of the log-relative hazard and cumulative baseline hazard 
##              [Tau1, Tau2], when covariate data is missing for individuals in 
##              the case cohort. It also returns the pure risk estimate. It
##              returns the overall influences, as well as the phase-two and
##              phase-three influences
## -----------------------------------------------------------------------------
## Arguments:
##
##  beta                    log-relative hazard 
##
##  Lambda0.Tau1Tau2        cumulative baseline hazard in the time interval that 
##                          is considered for the pure risk
##
##  x                       vector of length p, specifying the covariate profile 
##                          considered for the pure risk. Default is (0,...,0)
##
##  infl2.beta              phase-two influences for the log-relative hazard
##                  
##  infl2.Lambda0.Tau1Tau2  phase-two influences for the cumulative baseline
##                          hazard in the same time interval as that of the pure 
##                          risk
##
##  infl3.beta              phase-three influences for the log-relative hazard
##                  
##  infl3.Lambda0.Tau1Tau2  phase-three influences for the cumulative baseline
##                          hazard in the same time interval as that of the pure 
##                          risk 
## -----------------------------------------------------------------------------

influences.PR.missingdata  <- function (beta, Lambda0.Tau1Tau2, x = NULL, 
                                        infl2.beta, infl2.Lambda0.Tau1Tau2,
                                        infl3.beta, infl3.Lambda0.Tau1Tau2) {
  
  p <- length(beta)
  
  # ----------------------------------------------------------------------------
  # Quantities needed for the influences ---------------------------------------
  
  if ((length(x) != p) | (is.null(x))) {
    x               <- rep(0, p)
  } # if no covariate profile provided, use 0 as reference level
  exp.x             <- c(exp(x %*% beta))
  Pi.x              <- 1 - exp(- exp.x * Lambda0.Tau1Tau2)
  
  infl2.Pi.x  <- c((1 - Pi.x) * exp.x * Lambda0.Tau1Tau2) * 
    infl2.beta %*% x + c((1 - Pi.x) * exp.x) * 
    infl2.Lambda0.Tau1Tau2
  infl3.Pi.x  <- c((1 - Pi.x) * exp.x * Lambda0.Tau1Tau2) * 
    infl3.beta %*% x + c((1 - Pi.x) * exp.x) * 
    infl3.Lambda0.Tau1Tau2
  infl.Pi.x   <- infl2.Pi.x + infl3.Pi.x
  
  return(list(infl.Pi.x.Tau1Tau2 = infl.Pi.x, infl2.Pi.x.Tau1Tau2 = infl2.Pi.x, 
              infl3.Pi.x.Tau1Tau2 = infl3.Pi.x, Pi.x.Tau1Tau2.hat = Pi.x))
}



## -----------------------------------------------------------------------------
## Function: weights.phase3()
## -----------------------------------------------------------------------------
## Description: This function estimates the weights for the third phase of 
##              sampling (due to missingness)
## -----------------------------------------------------------------------------
## Arguments:
##
##  B.phase3        matrix with as many columns as phase-three strata, with one 
##                  1 per row, to indicate the phase-three stratum position 
##
##  total.phase2    un-weighted totals in the stratified case-cohort, using 
##                  all the individuals, even the ones with missing covariate
##                  data (phase-three data)
##
##  gamma0          vector of initial/possible values for gamma, to be used as 
##                  seed in the iterative procedure. Default is (0,...,0)
##
##  niter.max       maximum number of iterations for the Newton Raphson method. 
##                  Default is 10^4 iterations
##
##  epsilon.stop    threshold for the difference between the estimated weighted 
##                  total and the total in the whole cohort. If this difference 
##                  is less than the value of epsilon.stop, no more iterations 
##                  will be performed. Default is 10^(-10)
## -----------------------------------------------------------------------------

weights.phase3 <- function (B.phase3, total.phase2, gamma0 = NULL,
                            niter.max = NULL, epsilon.stop = NULL) {
  
  if ((is.null(gamma0)) || (length(gamma0) != ncol(B.phase3))) {
    gamma0        <- rep(0, ncol(B.phase3))
    print("gamma0 has not been filled out. Default is 0. ")
  }
  if ((is.null(epsilon.stop)) || (epsilon.stop <= 0)) {
    epsilon.stop  <- 10^(-10)
    print("epsilon.stop has not been filled out. Default is 10^(-10). ")
  }
  if ((is.null(niter.max)) || (niter.max <= 0)) {
    niter.max     <- 10^4
    print("niter.max has not been filled out. Default is 10^4 iterations.")
  }
  
  for (niter in 1:niter.max) {
    estimated.weights <- exp(B.phase3 %*% gamma0)
    B.estweighted     <- sweep(B.phase3, 1, estimated.weights, "*")
    f                 <- colSums(B.estweighted) - total.phase2
    
    epsilon           <- max(abs(f))
    if(epsilon < epsilon.stop){break}
    
    BB.estweighted    <- array(NA, dim = c(ncol(B.phase3), ncol(B.phase3), 
                                           nrow(B.phase3)))
    for (i in 1:nrow(B.phase3)){
      BB.estweighted[,, i] <- tcrossprod(B.phase3[i,], B.estweighted[i,])
    }
    drond.f.gamma0    <- apply(X = BB.estweighted, MARGIN = c(1,2), FUN = sum)
    gamma0            <- gamma0 - c(f %*% solve(drond.f.gamma0))
  }
  
  estimated.weights   <- as.numeric(exp(B.phase3 %*% gamma0))
  B.estweighted       <- sweep(B.phase3, 1, estimated.weights, "*")
  
  return(list(gamma.hat = gamma0, estimated.weights = estimated.weights, 
              estimated.total = colSums(B.estweighted)))
}


## -----------------------------------------------------------------------------
## Function: scenarios()
## -----------------------------------------------------------------------------
## Description:   This function creates a dataframe with the different scenarios
##                to be investigated over simulations
## -----------------------------------------------------------------------------
## Arguments:
##
##  n                     vector with the cohort sizes
##
##  prob.y                vector with the outcome probabilities (at the end of 
##                        study duration)
##
##  noncases.per.case     vector with the desired numbers of non-cases to be 
##                        sampled for each case
##
##  part                  vector with the number of parts to break down the 
##                        replications. By default there is no breaking down
## -----------------------------------------------------------------------------

scenarios <- function (n = c(5000, 10000), prob.y = c(0.02, 0.05, 0.1),
                       noncases.per.case = c(2, 4), part = NULL) {
  
  parameters <- NULL 
  
  # ----------------------------------------------------------------------------
  # Different combinations of parameters ---------------------------------------
  
  if(is.null(part)){
    
    for (l in 1:length(prob.y)) {
      
      for (m in 1:length(n)) {
        
        for (i in 1:length(noncases.per.case)) {
          
          parameters <- rbind(parameters, c(prob.y[l], n[m], 
                                            noncases.per.case[i]))
            
        }
      }
    }
    colnames(parameters) = c("prob.y", "n", "noncases.per.case")
  } else {
    
    for (l in 1:length(prob.y)) {
      
      for (m in 1:length(n)) {
        
        for (i in 1:length(noncases.per.case)) {
          
          for (j in 1:part) {
            
            parameters <- rbind(parameters, c(prob.y[l], n[m], 
                                              noncases.per.case[i], j))
            
          }
        }
      }
    }
    colnames(parameters) = c("prob.y", "n", "noncases.per.case", "part")
  }


  
  return(parameters)
}


## -----------------------------------------------------------------------------
## Function: X.generation()
## -----------------------------------------------------------------------------
## Description:   This function generates a gaussian covariate value, with unit
##                variance and given mean
## -----------------------------------------------------------------------------
## Arguments:
##
##  mu                    mean of the covariate
## -----------------------------------------------------------------------------

X.generation <- function (mu) { 
  return(rnorm(n = 1, mean = mu, sd = 1))
}


## -----------------------------------------------------------------------------
## Functions: densityX1X3() and densityX1()
## -----------------------------------------------------------------------------
## Description:   Functions needed to compute the time-fixed baseline hazard in
##                the design used in the simulation study in Section 5 and Web 
##                Appendix E
## -----------------------------------------------------------------------------

densityX1X3 <- function (x) {
  return(1 / sqrt(2 * pi) * exp(-1/2 * (x - (beta1 + beta3 * alpha1))^2))
}

densityX1 <- function (x) {
  return(1 / sqrt(2 * pi) * exp(-1/2 * (x)^2))
}

## -----------------------------------------------------------------------------
## Functions: E.RH()
## -----------------------------------------------------------------------------
## Description:   This function computes the quantity 
##                E{exp(beta1 * X1 + beta2 * X2 + beta3 * X3)}, needed to
##                compute the time-fixed baseline hazard and the sampling
##                fraction to be used for an un-stratified sampling of the 
##                case-cohort in the design used in the simulation study in
##                Section 5 and Web Appendix E
## -----------------------------------------------------------------------------
## Arguments:
##
##  alpha1                effect of covariate X1 on X3
##
##  alpha2                effect of covariate X2 on X3
##
##  beta1                 effect of covariate X1 on the outcome
##
##  beta2                 effect of covariate X2 on the outcome
##
##  beta3                 effect of covariate X1=3 on the outcome
##
##  Prob.x2               probabilities for the three possible values of 
##                        covariate X2, conditional on the value of covariate X1
## -----------------------------------------------------------------------------

E.RH <- function (alpha1 = 0.05, alpha2 = -0.35, beta1 = - 0.2, beta2 = 0.25, 
                    beta3 = - 0.3, Prob.x2 = matrix(c(0.7, 0.05, 0.25, 0.45, 
                                                      0.2, 0.35, 0.40, 0.3, 
                                                      0.3), nrow = 3, 
                                                    byrow = TRUE)) {
  
 return(exp(1/2 * beta3^2 + 1/2 * (beta1 + beta3 * alpha1)^2) * 
          (integrate(f = densityX1X3, lower = -10, upper = -2)$value * 
             sum(exp((beta2 + beta3 * alpha2) * 0:2) * Prob.x2[1,]) + 
             integrate(densityX1X3, -2, 1)$value * 
             sum(exp((beta2 + beta3 * alpha2) * 0:2) * Prob.x2[2,]) + 
             integrate(densityX1X3, 1, 10)$value * 
             sum(exp((beta2 + beta3 * alpha2) * 0:2) * Prob.x2[3,]))) 
}



## -----------------------------------------------------------------------------
## Functions: E.RH.w()
## -----------------------------------------------------------------------------
## Description:   This function computes the quantities 
##                E{exp(beta1 * X1 + beta2 * X2 + beta3 * X3) | W = w}, w = {0,
##                1, 2, 3} that are needed to compute the stratum-specific 
##                sampling fractions to be used for a stratified sampling of the
##                case-cohort in the design used in the simulation study in
##                Section 5 and Web Appendix E
## -----------------------------------------------------------------------------
## Arguments:
##
##  alpha1                effect of covariate X1 on X3
##
##  alpha2                effect of covariate X2 on X3
##
##  beta1                 effect of covariate X1 on the outcome
##
##  beta2                 effect of covariate X2 on the outcome
##
##  beta3                 effect of covariate X1=3 on the outcome
##
##  Prob.x2               probabilities for the three possible values of 
##                        covariate X2, conditional on the value of covariate X1
## -----------------------------------------------------------------------------

E.RH.w <- function (alpha1 = 0.05, alpha2 = -0.35, beta1 = - 0.2, beta2 = 0.25, 
                    beta3 = - 0.3, Prob.x2 = matrix(c(0.7, 0.05, 0.25, 0.45, 
                                                    0.2, 0.35, 0.40, 0.3, 
                                                    0.3), nrow = 3, 
                                                    byrow = TRUE)) {
  
  E.exp.w0 <- exp(1/2 * beta3^2 + 1/2 * (beta1 + beta3 * alpha1)^2) * 
    (integrate(f = densityX1X3, lower = 0, upper = 1)$value * 
       sum((exp((beta2 + beta3 * alpha2) * 0:2) * Prob.x2[2,])[1]) / 
       sum(Prob.x2[2,][1]) / 
       integrate(f = densityX1, lower = 0, upper = 10)$value +
       integrate(f = densityX1X3, lower = 1, upper = 10)$value * 
       sum((exp((beta2 + beta3 * alpha2) * 0:2) * Prob.x2[3,])[1]) / 
       sum(Prob.x2[3,][1]) / 
       integrate(f = densityX1, lower = 0, upper = 10)$value) 
  
  E.exp.w1 <- exp(1/2 * beta3^2 + 1/2 * (beta1 + beta3 * alpha1)^2) * 
    (integrate(f = densityX1X3, lower = -10, upper = -2)$value * 
       sum((exp((beta2 + beta3 * alpha2) * 0:2) * Prob.x2[1,])[1:2]) / 
       sum(Prob.x2[1,][1:2]) / 
       integrate(f = densityX1, lower = -10, upper = 0)$value + 
       integrate(f = densityX1X3, lower = -2, upper = 0)$value * 
       sum((exp((beta2 + beta3 * alpha2) * 0:2) * Prob.x2[2,])[1:2]) / 
       sum(Prob.x2[2,][1:2]) / 
       integrate(f = densityX1, lower = -10, upper = 0)$value)
  
  E.exp.w2 <- exp(1/2 * beta3^2 + 1/2 * (beta1 + beta3 * alpha1)^2) * 
    (integrate(f = densityX1X3, lower = 0, upper = 1)$value * 
       sum((exp((beta2 + beta3 * alpha2) * 0:2) * Prob.x2[2,])[2:3]) / 
       sum(Prob.x2[2,][2:3]) / 
       integrate(f = densityX1, lower = 0, upper = 10)$value + 
       integrate(f = densityX1X3, lower = 1, upper = 10)$value * 
       sum((exp((beta2 + beta3 * alpha2) * 0:2) * Prob.x2[3,])[2:3]) / 
       sum(Prob.x2[3,][2:3]) / 
       integrate(f = densityX1, lower = 0, upper = 10)$value)
  
  E.exp.w3 <- exp(1/2 * beta3^2 + 1/2 * (beta1 + beta3 * alpha1)^2) * 
    (integrate(f = densityX1X3, lower = -10, upper = -2)$value * 
       ((exp((beta2 + beta3 * alpha2) * 0:2) * Prob.x2[1,])[3]) / 
       sum(Prob.x2[1,][3]) / 
       integrate(f = densityX1, lower = -10, upper = 0)$value + 
       integrate(f = densityX1X3, lower = -2, upper = 0)$value * 
       ((exp((beta2 + beta3 * alpha2) * 0:2) * Prob.x2[2,])[3]) / 
       sum(Prob.x2[2,][3]) / 
       integrate(f = densityX1, lower = -10, upper = 0)$value)

  return(list(E.exp.w0 = E.exp.w0, E.exp.w1 = E.exp.w1, E.exp.w2 = E.exp.w2,
              E.exp.w3 = E.exp.w3)) 
}


## -----------------------------------------------------------------------------
## Function: scenarios.missingdata()
## -----------------------------------------------------------------------------
## Description:   This function creates a dataframe with the different scenarios
##                to be investigated over simulations in Web Appendix I (when 
##                covariate data is missing for individuals in the case cohort)
## -----------------------------------------------------------------------------
## Arguments:
##
##  n                     vector with the cohort sizes
##
##  prob.y                vector with the outcome probabilities (at the end of 
##                        study duration)
##
##  noncases.per.case     vector with the desired numbers of non-cases to be 
##                        sampled for each case
##
##  part                  vector with the number of parts to break down the 
##                        replications. By default there is no breaking down
## -----------------------------------------------------------------------------

scenarios.missingdata <- function (n = c(5000, 10000), 
                                   prob.y = c(0.02, 0.05, 0.1),
                                   noncases.per.case = c(2, 4), 
                                   strata.prob.missingness = rbind(c(0.9, 0.8), 
                                                                   c(0.98, 0.9)),
                                   part = NULL) {
  
  parameters <- NULL 
  
  # ----------------------------------------------------------------------------
  # Different combinations of parameters ---------------------------------------
  
  if(is.null(part)){
    
    for (l in 1:length(prob.y)) {
      
      for (m in 1:length(n)) {
        
        for (i in 1:length(noncases.per.case)) {
          
          for(k in 1:length(strata.prob.missingness[, 1])) {
            
            parameters <- rbind(parameters, c(prob.y[l], n[m], 
                                              noncases.per.case[i],
                                              strata.prob.missingness[k,]))
            
          }
        }
      }
    }
    colnames(parameters) = c("prob.y", "n", "noncases.per.case", 
                             "prob.missingness.noncases", 
                             "prob.missingness.cases")
  } else {
    
    for (l in 1:length(prob.y)) {
      
      for (m in 1:length(n)) {
        
        for (i in 1:length(noncases.per.case)) {
          
          for(k in 1:length(strata.prob.missingness[, 1])) {
            
            for (j in 1:part) {
              
              parameters <- rbind(parameters, c(prob.y[l], n[m], 
                                                noncases.per.case[i], 
                                                strata.prob.missingness[k,], j))
              
              
            }
          }
        }
      }
    }
    colnames(parameters) = c("prob.y", "n", "noncases.per.case", 
                             "prob.missingness.noncases", 
                             "prob.missingness.cases", "part")
  }
  return(parameters)
}

## -----------------------------------------------------------------------------
## Function: robustvariance()
## -----------------------------------------------------------------------------
## Description: This function returns the robust variance estimate, i.e., 
##              computes the sum of the squared influence functions, for
##              parameters such as log-relative hazard, cumulative baseline 
##              hazard or covariate specific pure-risk. The overall influences
##              should be provided. The function works with design weights or 
##              calibrated weights, or with missing covariates in phase two data  
##              (i.e. with influences obtained from the influences.missingdata R 
##              function)
## -----------------------------------------------------------------------------
## Argument:
##
##  infl        overall influences on parameters such as log-relative hazard, 
##              cumulative baseline hazard or covariate specific pure-risk. 
## -----------------------------------------------------------------------------

robustvariance <- function (infl) {
  
  return(robust.var = diag(crossprod(infl)))
  
}


## -----------------------------------------------------------------------------
## Function: variance()
## -----------------------------------------------------------------------------
## Description: This function returns the variance estimate that influence-based
##              and follows the complete variance decomposition, for 
##              parameters such as log-relative hazard, cumulative baseline 
##              hazard or covariate specific pure-risk. The function works with 
##              design weights or calibrated weights
## -----------------------------------------------------------------------------
## Arguments:
##
## n                number of individuals in the whole cohort
##
## casecohort       data frame with id (identification) and weights (design or 
##                  calibrated, if they are not provided in the argument below),
##                  for each individual in the case cohort data. If stratified =
##                  TRUE, it should also contain W (strata), strata.m (number of
##                  sampled individuals in the stratum) and strata.n (stratum 
##                  size). If stratified = FALSE, it should also contain  m 
##                  (number of sampled individuals) and n (cohort size).
##
## weights          vector with design or calibrated weights for each individual 
##                  in the case cohort data
##
## infl             overall influences on a parameter such as log-relative 
##                  hazard, cumulative baseline hazard or covariate specific 
##                  pure-risk 
##
## calibrated       are calibrated weights used for the estimation of the 
##                  parameter? If calibrated = TRUE, the arguments below need to 
##                  be provided. Default is FALSE
##
## infl2            phase-two influences. Should be provided if calibrated 
##                  weights are used
##
## cohort           if stratified = TRUE, data frame with status (case status)
##                  and subcohort (subcohort sampling indicators) for each 
##                  individual in the stratified case cohort data. If stratified
##                  = FALSE, data frame with status (case status) and 
##                  unstrat.subcohort (subcohort unstratified sampling 
##                  indicators). Should be provided if calibrated weights are 
##                  used
##                          
## stratified       was the sampling of the case-cohort stratified on W? Default
##                  is FALSE
##
## variance.phase2  should the phase-two variance component also be returned? 
##                  Default is FALSE
## -----------------------------------------------------------------------------

variance <- function(n, casecohort, weights = NULL, infl, calibrated = NULL, 
                     infl2 = NULL, cohort = NULL, stratified = NULL, 
                     variance.phase2 = NULL) {
  
  if (is.null(calibrated)) {
    calibrated <- FALSE
  }
  if (is.null(stratified)) {
    stratified <- FALSE
  }
  if (is.null(variance.phase2)) {
    variance.phase2 <- FALSE
  }
  if (!is.null(weights)){
    casecohort$weights = weights
    } else {
      if (is.null(casecohort$weights)) {
        print("the individuals weights must be provided")
        break
      }
    }
    
  prod.covar.weight <- product.covar.weight(casecohort, stratified)
  
  if(calibrated == FALSE) {
    superpop.var  <- n / (n - 1) * crossprod((infl / sqrt(casecohort$weights)))
    phase2.var    <- with(casecohort, t(infl[status == 0,]) %*% 
                            prod.covar.weight %*% infl[status == 0,])
    var           <- diag(superpop.var + phase2.var)
  } else {
    if (is.null(infl2) | is.null(cohort)) {
      print("infl2 and cohort must be provided if calibrated weights have been 
            used")
      break
    } else {
      omega         <- rep(1, n)
      omega[casecohort$id]  <- casecohort$weights
      
      superpop.var  <- n / (n - 1) * colSums((infl - infl2) ^ 2 + 2 * 
                                               (infl - infl2) * infl2 + 
                                               (infl2 / sqrt(omega)) ^ 2)
      if (stratified == TRUE) {
        phase2.var    <- with(cohort, 
                              t(infl2[(status == 0) & (subcohort == TRUE),]) %*% 
                                prod.covar.weight %*% 
                                infl2[(status == 0) & (subcohort == TRUE),])
      } else {
        phase2.var    <- with(cohort, 
                              t(infl2[(status == 0) & 
                                        (unstrat.subcohort == TRUE),]) %*% 
                                prod.covar.weight %*% 
                                infl2[(status == 0) & 
                                        (unstrat.subcohort == TRUE),])
      }
      var           <- superpop.var + diag(phase2.var)
    }
  }
  
  if (variance.phase2 == FALSE) {
    return(variance = var)
  } else {
    return(list(variance = var, variance.phase2 = diag(phase2.var)))
  }
}


## -----------------------------------------------------------------------------
## Function: variance.missingdata()
## -----------------------------------------------------------------------------
## Description: This function returns the variance estimate that influence-based
##              and follows the complete variance decomposition, for 
##              parameters such as log-relative hazard, cumulative baseline 
##              hazard or covariate specific pure-risk, when covariate 
##              information is missing for individuals in phase-two sample
## -----------------------------------------------------------------------------
## Arguments:
##
## n                  number of individuals in the whole cohort
##
## casecohort         data frame with id (identification) and weights (design or 
##                    calibrated, if they are not provided in the argument 
##                    below), for each individual in the case cohort data. If 
##                    stratified = TRUE, it should also contain W (strata), 
##                    strata.m (number of sampled individuals in the stratum) 
##                    and strata.n (stratum size). If stratified = FALSE, it 
##                    should also contain  m (number of sampled individuals) and 
##                    n (cohort size)
##
## casecohort.phase2  if stratified = TRUE, data frame with W (strata), strata.m
##                    (number of sampled individuals in the stratum) and 
##                    strata.n (stratum size), for each individual in the phase 
##                    two sample. If stratified = FALSE, data frame with m 
##                    (number of sampled individuals) and n (cohort size), for 
##                    each individual in the case cohort data
##
## weights            vector with design weights for each individual in the case
##                    cohort data
##
## weights.phase2     vector with design weights for each individual in the 
##                    phase two sample
##
## weights.p2.phase2  vector with phase-two design weights for each individual 
##                    in the phase two sample
##
## infl2              phase-two influences on a parameter such as log-relative 
##                    hazard, cumulative baseline hazard or covariate specific 
##                    pure-risk  
##
## infl3              phase-three influences on a parameter such as log-relative 
##                    hazard, cumulative baseline hazard or covariate specific 
##                    pure-risk  
##
## stratified.p2      was the second phase of sampling stratified on W? Default 
##                    is FALSE
##
## estimated.weights  were the phase-three weights estimated? Default is FALSE
## -----------------------------------------------------------------------------


variance.missingdata <- function (n, casecohort, casecohort.phase2, 
                                  weights, weights.phase2, weights.p2.phase2, 
                                  infl2, infl3, stratified.p2 = NULL,
                                  estimated.weights = NULL) {
  
  if (is.null(stratified.p2)) {
    stratified.p2 <- FALSE
  }
  if (is.null(estimated.weights)) {
    estimated.weights <- FALSE
  }
  
  if (estimated.weights == FALSE){
    if (stratified.p2 == TRUE){
      casecohort.W    <- casecohort$W
      casecohort.W3   <- casecohort$W3
      casecohort.n    <- casecohort$strata.n
      casecohort.m    <- casecohort$strata.m
      casecohort.m.n  <- casecohort.m / casecohort.n
      casecohort.n.m  <- casecohort.n / casecohort.m
      
      product.covar.weight <- function (indiv) { 
        return((casecohort.W == casecohort.W[indiv]) * 
                 (casecohort.W3 == casecohort.W3[indiv]) * 
                 (casecohort.W3[indiv] == 0) * 
                 (1 - casecohort.m.n * (casecohort.n - 1) / (casecohort.m - 1) ))
      }
      prod.covar.weight       <- do.call(rbind, lapply(1:nrow(casecohort), 
                                                       product.covar.weight))
      diag(prod.covar.weight) <- (1 - 1 / weights)
      
      superpop.var  <- n / (n-1) * (crossprod(infl2 / sqrt(weights.p2.phase2)) +
                                      2 * crossprod(infl2 /  weights.p2.phase2,
                                                    infl3) + 
                                      crossprod(infl3 / sqrt(weights.phase2)))
      phase23.var   <- with(casecohort.phase2, 
                            t(infl3[phase3 == 1,]) %*% prod.covar.weight %*% 
                              infl3[phase3 == 1,])
      var           <- diag(superpop.var + phase23.var)
      
    } else {
      casecohort.W3   <- casecohort$W3
      casecohort.n    <- casecohort$n
      casecohort.m    <- casecohort$m
      casecohort.m.n  <- casecohort.m / casecohort.n
      casecohort.n.m  <- casecohort.n / casecohort.m
      
      product.covar.weight <- function (indiv) { 
        return((casecohort.W3 == casecohort.W3[indiv]) * 
                 (casecohort.W3[indiv] == 0) * 
                 (1 - casecohort.m.n * (casecohort.n - 1) / (casecohort.m - 1) ))
      }
      prod.covar.weight       <- do.call(rbind, lapply(1:nrow(casecohort), 
                                                       product.covar.weight))
      diag(prod.covar.weight) <- (1 - 1 / weights)
      
      superpop.var  <- n / (n-1) * (crossprod(infl2 / sqrt(weights.p2.phase2)) +
                                      2 * crossprod(infl2 /  weights.p2.phase2,
                                                    infl3) +
                                      crossprod(infl3 / sqrt(weights.phase2)))
      phase23.var   <- with(casecohort.phase2, 
                            t(infl3[unstrat.phase3 == 1,]) %*% 
                              prod.covar.weight %*% infl3[unstrat.phase3 == 1,])
      var           <- diag(superpop.var + phase23.var)
    }
  } else {
    
    weights.p3.phase2 <- weights.phase2 / weights.p2.phase2
    infl <- infl2 + infl3
    
    if (stratified.p2 == TRUE) {
      phase2.W   <- casecohort.phase2$W
      phase2.W3   <- casecohort.phase2$W3
      phase2.n    <- casecohort.phase2$strata.n
      phase2.m    <- casecohort.phase2$strata.m
      phase2.m.n  <- phase2.m / phase2.n
      phase2.n.m  <- phase2.n / phase2.m
      product.covar.weight <- function (indiv) {
        return((phase2.W == phase2.W[indiv]) * 
                 (phase2.W3 == phase2.W3[indiv]) * 
                 (phase2.W3[indiv] == 0) * 
                 (1 - phase2.m.n * (phase2.n - 1) / (phase2.m - 1) ))
      }
      prod.covar.weights        <- do.call(rbind, 
                                           lapply(1:nrow(casecohort.phase2), 
                                                  product.covar.weight))
      diag(prod.covar.weights)  <- (1 - 1 / weights.phase2)
      
      superpop.var  <- n / (n-1) * (crossprod(infl2 / sqrt(weights.p2.phase2)) +
                                      2 * crossprod(infl2 / weights.p2.phase2, 
                                                    infl3) + 
                                      crossprod(infl3 / sqrt(weights.phase2)))
      phase23.var   <- t(infl) %*% prod.covar.weights %*% infl - t(infl2) %*% 
        diag(c(1 / weights.p2.phase2 * (1 - 1 / weights.p3.phase2))) %*% 
        (infl2 + 2 * infl3)
      var           <- diag(superpop.var + phase23.var)
    } else {
      phase2.W3   <- casecohort.phase2$W3
      phase2.n    <- casecohort.phase2$n
      phase2.m    <- casecohort.phase2$m
      phase2.m.n  <- phase2.m / phase2.n
      phase2.n.m  <- phase2.n / phase2.m
      product.covar.weight <- function (indiv) {
        return((phase2.W3 == phase2.W3[indiv]) * 
                 (phase2.W3[indiv] == 0) * 
                 (1 - phase2.m.n * (phase2.n - 1) / (phase2.m - 1) ))
      }
      prod.covar.weights        <- do.call(rbind, 
                                           lapply(1:nrow(casecohort.phase2), 
                                                  product.covar.weight))
      diag(prod.covar.weights)  <- (1 - 1 / weights.phase2)
      
      superpop.var  <- n / (n-1) * (crossprod(infl2 / sqrt(weights.p2.phase2)) +
                                      2 * crossprod(infl2 / weights.p2.phase2, 
                                                    infl3) + 
                                      crossprod(infl3 / sqrt(weights.phase2)))
      phase23.var   <- t(infl) %*% prod.covar.weights %*% infl - t(infl2) %*% 
        diag(c(1 / weights.p2.phase2 * (1 - 1 / weights.p3.phase2))) %*% 
        (infl2 + 2 * infl3)
      var           <- diag(superpop.var + phase23.var)
    }
  }
  return(variance = var)
}

