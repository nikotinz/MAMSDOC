################################################################################
################################ Check input ###################################
################################################################################

mams.check.dtl <- function(obj) {

  m <- match(c(
    "K", "J", "alpha", "power", "r", "r0", "p", "p0",
    "delta", "delta0", "sd", "ushape", "lshape", "ufix", "lfix",
    "nstart", "nstop", "sample.size", "Q", "type", "method",
    "parallel", "print", "nsim", "H0", "obj", "par", "sim"
  ), names(obj), 0)

  mc <- obj[c(m)]

if (mc[["J"]] <= 1) {
stop("For method = 'dtl', the number of stages 'J' should be greater than 1.")
}
if (length(mc[["K"]]) != mc[["J"]]) {
stop("`K` need to be defined as vector length of `J`")
}

obj$Kv <- mc[["K"]]

    if (any(obj$Kv%%1 != 0, !is.numeric(obj$Kv), obj$Kv < 1,
          is.infinite(obj$Kv),
          length(obj$Kv) < 2,
          obj$Kv[2:length(obj$Kv)] >= obj$Kv[1:(length(obj$Kv) - 1)],
          obj$Kv[length(obj$Kv)] != 1)) {
    stop("K must be a monotonically decreasing vector, of length at least 2, ",
          "containing only integers, with final element equal to 1.")
  }

  if (mc[["alpha"]] < 0 | mc[["alpha"]] > 1 | mc[["power"]] < 0 | 
  mc[["power"]] > 1) {
    stop("Error rate or power not between 0 and 1.")
  }
  if (length(mc[["r"]]) != length(mc[["r0"]])) {
    stop("Different length of allocation ratios on control and experimental ",
          "treatments.")
  }
  if (length(mc[["r"]]) != mc[["J"]]) {
    stop("Length of allocation ratios does not match number of stages.")
  }

  if (is.numeric(mc[["p"]]) & is.numeric(mc[["p0"]]) & 
      is.numeric(mc[["delta"]]) &
      is.numeric(mc[["delta0"]]) &
      is.numeric(mc[["sd"]])) {
    stop("Specify the effect sizes either via p or via (delta, sd) and set the",
          " other parameter(s) to NULL.")
  }
  if (is.numeric(mc[["p"]]) & is.numeric(mc[["p0"]])) {
    if (mc[["p"]] < 0 | mc[["p"]] > 1) {
      stop("Treatment effect parameter not within 0 and 1.")
    }
  } else {
    if (is.numeric(obj[["delta"]]) ||
    (!is.null(obj$par) && !is.null(obj$par[["delta"]]) 
    && is.numeric(obj$par[["delta"]])) & 
    is.numeric(mc[["sd"]])) {
      if (mc[["sd"]] <= 0) {
        stop("Standard deviation must be positive.")
      }
    } else {
      stop("Specify the effect sizes either via p or via (delta, sd) and set ",
            "the other parameter(s) to NULL.")
    }
  }
  if (length(mc[["r0"]]) != mc[["J"]]) {
    stop("`r0` needs to be defined as a vector of length `J`.")
  }
  if (any(diff(mc[["r0"]]) <= 0)) {
    stop("`r0` must be a monotonically increasing vector.")
  }
  if (length(mc[["r"]]) != mc[["J"]]) {
    stop("`r` needs to be defined as a vector of length `J`.")
  }
  if (any(diff(mc[["r"]]) <= 0)) {
    stop("`r` must be a monotonically increasing vector.")
  }
    if (mc[["r0"]][1] < 1) {
    stop("`r0[1]` must be >= 1.")
  }
  if (any(mc[["r"]] <= 0)) {
    stop("`r` values must be >= 1.")
  }
  if (mc[["r0"]][1] %% 1 != 0) {
    stop("First element of `r0` must be integers.")
  }

  return(obj)
}

###############################################################################
###################### Fit function ###########################################
###############################################################################

mams.fit.dtl <- function(obj) {
##### Initialise internal functions ##########################################
# Used to find the critical boundary e
  dtl_find_e <- function(e, Kv, alpha, J, outcomes, lowers, uppers, means_HG,
                          Lambda, print) {
    if (print) {
    message(".", appendLF = FALSE)
    }
    # Update the lower and upper boundaries of the integrals for the current
    # value of e, and evaluate the probability of each of the possible outcomes
    for (i in 1:Kv[J]) {
      for (k in 1:Kv[J]) {
        if (outcomes[i, Kv[1] + k] == 1L) {
          lowers[[i]][nrow(Lambda) - Kv[J] + k] <- e
        } else {
          lowers[[i]][nrow(Lambda) - Kv[J] + k] <- -Inf
          uppers[[i]][nrow(Lambda) - Kv[J] + k] <- e
        }
      }
      outcomes[i, Kv[1] + Kv[J] + 1L]           <-
        mvtnorm::pmvnorm(lowers[[i]], uppers[[i]], means_HG, sigma = Lambda)[1]
    }
    # Return the difference between the computed FWER and the nominal alpha
    return(factorial(Kv[1])*sum(outcomes[, Kv[1] + Kv[J] + 1L]) - alpha)
  }

  # Used to find the required sample size n
  dtl_find_n <- function(n, Kv, power, J, outcomes, lowers, uppers, means_LFC,
                          Lambda) {
    # Evaluate the probability of each of the possible outcomes
    for (i in 1:Kv[J]) {
      outcomes[i, Kv[1] + Kv[J] + 2L] <-
        mvtnorm::pmvnorm(lowers[[i]], uppers[[i]], sqrt(n)*means_LFC,
                          sigma = Lambda)[1]
    }
    # Return the difference between the computed power and the nominal power
    return(factorial(Kv[1] - 1)*sum(outcomes[, Kv[1] + Kv[J] + 2L]) - power)
  }

  ##### Convert treatment effects ##############################################
  if (is.numeric(obj$p) & is.numeric(obj$p0)) {
    delta <- sqrt(2) * qnorm(obj$p)
    delta0 <- sqrt(2) * qnorm(obj$p0)
    sig <- 1
    p0  <- obj$p0
  } else {
    delta <- ifelse(is.null(obj[["delta"]]), obj$par[["delta"]],
                                              obj[["delta"]]) 
    delta0 <- ifelse(is.null(obj[["delta0"]]), obj$par[["delta0"]],
                                              obj[["delta0"]]) 
    # for subsequent if(J==1 & p0==0.5)
    p0 <- pnorm(delta0 / sqrt(2 * obj$sd^2))
    sig <- obj$sd
  }

  ##### Create no-drop covariance matrix #######################################

  Lambda                         <- matrix(0, obj$J*obj$Kv[1], obj$J*obj$Kv[1])
  add                            <- obj$r / (obj$r + obj$r0)
  minus                          <- 1 - add
  for (j in 1:obj$J) {
    range                        <- (1 + (j - 1)*obj$Kv[1]):(j*obj$Kv[1])
    Lambda[range, range]         <- matrix(add[j], obj$Kv[1], obj$Kv[1]) +
      diag(minus[j], obj$Kv[1], obj$Kv[1])
  }
  for (j1 in 2:obj$J) {
    for (j2 in 1:(j1 - 1)) {
      range_j1                   <- (1 + (j1 - 1)*obj$Kv[1]):(j1*obj$Kv[1])
      range_j2                   <- (1 + (j2 - 1)*obj$Kv[1]):(j2*obj$Kv[1])
      cov_factor                 <- sqrt(obj$r[j1]*obj$r0[j1] /
                                        (obj$r[j1] + obj$r0[j1]))*
                                        (1/obj$r[j1] + 1/obj$r0[j1])*
                          sqrt(obj$r[j2]*obj$r0[j2] / (obj$r[j2] + obj$r0[j2]))
      cov_Zj1Zj2                 <- matrix(add[j2]*cov_factor,
                                          obj$Kv[1], obj$Kv[1]) +
        diag(minus[j2]*cov_factor, obj$Kv[1], obj$Kv[1])
      Lambda[range_j2, range_j1] <- Lambda[range_j1, range_j2] <- cov_Zj1Zj2
    }
  }

  ##### Create matrix of unique rejection outcomes and associated MVN ##########
  ##### objects ################################################################

  rej                           <- matrix(0L, obj$Kv[obj$J], obj$Kv[obj$J])
  rej[lower.tri(rej, diag = TRUE)] <- 1L
  outcomes                      <- cbind(matrix(1:obj$Kv[1], obj$Kv[obj$J],
                                                obj$Kv[1],
                                                byrow = TRUE),
                                          rej, matrix(0L, obj$Kv[obj$J], 2))
  conditions                    <- sum(obj$Kv[1:(obj$J - 1)]) - (obj$J - 1) +
                                      obj$Kv[obj$J] +
                             as.numeric(obj$Kv[obj$J] > 1) * (obj$Kv[obj$J] - 1)
  A                             <- matrix(0L, conditions, obj$J*obj$Kv[1])
  counter                       <- 1
  for (j in 1:(obj$J - 1)) {
    dropped                     <- obj$Kv[j]:(obj$Kv[j + 1] + 1)
    if (length(dropped) > 1) {
      for (cond in 1:(obj$Kv[j] - obj$Kv[j + 1] - 1)) {
        A[counter,
          obj$Kv[1] * (j - 1) +
            which(outcomes[1, 1:obj$Kv[1]] == dropped[cond + 1])]        <- 1L
        A[counter,
          obj$Kv[1] * (j - 1) + which(outcomes[1,
                                        1:obj$Kv[1]] == dropped[cond])]  <- -1L
        counter                 <- counter + 1
      }
    }
    continue                    <- obj$Kv[j + 1]:1
    for (cond in 1:obj$Kv[j + 1]) {
      A[counter,
        obj$Kv[1] * (j - 1) + which(outcomes[1,
                                    1:obj$Kv[1]] == continue[cond])]   <- 1L
      A[counter,
        obj$Kv[1] * (j - 1) +
          which(outcomes[1, 1:obj$Kv[1]] == dropped[length(dropped)])] <- -1L
      counter                   <- counter + 1
    }
  }

  if (obj$Kv[obj$J] > 1) {
    condition                  <- matrix(0L, obj$Kv[obj$J] - 1, obj$J*obj$Kv[1])
    for (cond in 1:(obj$Kv[obj$J] - 1)) {
      A[counter,
        obj$Kv[1] * (obj$J - 1) + which(outcomes[1,
                                1:obj$Kv[1]] == obj$Kv[obj$J] - cond)]     <- 1L
      A[counter,
        obj$Kv[1] * (obj$J - 1) + which(outcomes[1,
                              1:obj$Kv[1]] == obj$Kv[obj$J] - cond + 1)] <- -1L
      counter                   <- counter + 1
    }
  }
  for (k in 1:obj$Kv[obj$J]) {
    A[counter, (obj$J - 1)*obj$Kv[1] + which(outcomes[1,
                                            1:obj$Kv[1]] == k)]         <- 1L
    counter                     <- counter + 1
  }
  means_HG                      <- numeric(conditions)
  means_LFC                     <-
    as.numeric(A %*% (rep(c(delta, rep(delta0, obj$Kv[1] - 1)), obj$J)*
                      sqrt(rep((1 / (sig^2/obj$r0 +
                                    sig^2/obj$r)),
                                    each = obj$Kv[1]))))
  Lambda                        <- A%*%Lambda%*%t(A)
  lowers_i                      <- numeric(conditions)
  uppers_i                      <- rep(Inf, conditions)
  lowers                        <- uppers <- list()
  for (i in 1:obj$Kv[obj$J]) {
    lowers[[i]]                 <- lowers_i
    uppers[[i]]                 <- uppers_i
  }
  ##### Compute critical boundary ##############################################
  # Search using uniroot
    if (obj$print) {
    message("   i) find new lower and upper boundaries\n      ", 
            appendLF = FALSE)
  }
  e                                           <-
    stats::uniroot(f        = dtl_find_e, interval = c(0, 5),
                    Kv       = obj$Kv,
                    alpha    = obj$alpha,
                    J        = obj$J,
                    outcomes = outcomes,
                    lowers   = lowers,
                    uppers   = uppers,
                    means_HG = means_HG,
                    Lambda   = Lambda,
                    print    = obj$print)$root
  # Update the lower and upper boundaries of the integrals for the determined
  # value of e
  for (i in 1:obj$Kv[obj$J]) {
    for (k in 1:obj$Kv[obj$J]) {
      if (outcomes[i, obj$Kv[1] + k] == 1L) {
        lowers[[i]][nrow(Lambda) - obj$Kv[obj$J] + k] <- e
      } else {
        lowers[[i]][nrow(Lambda) - obj$Kv[obj$J] + k] <- -Inf
        uppers[[i]][nrow(Lambda) - obj$Kv[obj$J] + k] <- e
      }
    }
  }

  ##### Compute required sample size ###########################################
  if (obj$sample.size) {
  if (obj$print) {
        message("\n   ii) perform sample size calculation\n", appendLF = FALSE)
      }
    # If nstop == NULL, then set nstop as 3 x the sample size required by a
    # corresponding single-stage design
    if (is.null(obj$nstop)) {
      rho       <- obj$r[1] / (obj$r[1] + 1)
      corr      <- matrix(rho, obj$Kv[1], obj$Kv[1]) + diag(1 - rho, obj$Kv[1])
      quan      <- mvtnorm::qmvnorm(1 - obj$alpha, corr = corr)$quantile

      obj$nstop     <- 3*ceiling(((quan*sig[1]*sqrt(1 + 1/obj$r[1]) +
                                stats::qnorm(obj$power)*
                  sqrt(sig[1]^2 + sig[1]^2/obj$r[1]))/delta)^2)
    }

      if (obj$r[1] > obj$r0[1]) {
      obj$r <- obj$r / obj$r0[1]
      obj$r0 <- obj$r0 / obj$r0[1]
      message("Allocation ratio for control arm at first stage greater than for 
      treatment arm(s), using normalisation by r0[1] \n")
    }
    # Loop until a value of n providing the desired power is found
        n             <- obj$nstart
    power_check   <- FALSE
    while (all(!power_check, n <= obj$nstop)) {
      power_check <- (dtl_find_n(n, obj$Kv, obj$power, obj$J, outcomes, 
                                  lowers, uppers, means_LFC, Lambda) >= 0)
      n           <- n + 1L
    }
    n             <- n - 1L
    n <- n * obj$r0[1]

    if (n == obj$nstop) {
      warning("The sample size was limited by nstop.")
    }
  } else {
    n             <- NULL
  }

  
  ##### Output #################################################################
  Kv_diff     <- c(obj$Kv[1:(obj$J - 1)] - obj$Kv[2:obj$J], obj$Kv[obj$J])
  

  res         <- list(K           = obj$K,
                      l           = c(rep(NA, obj$J - 1), e), # May want to
                      u           = c(rep(NA, obj$J - 1), e), # change this
                      n           = n,
                      r           = obj$r,
                      r0          = obj$r0,
                      Q           = obj$Q,
                      Kv          = obj$Kv,
                      J           = obj$J,
                      p           = obj$p,
                      p0          = obj$p0,
                      delta       = obj[["delta"]],
                      delta0      = obj[["delta0"]],
                      alpha       = obj$alpha,
                      alpha.star  = c(rep(0, obj$J - 1), obj$alpha),
                      lshape      = ifelse(is.function(obj$lshape),
                                      "self-defined",obj$lshape),
                      ushape      = ifelse(is.function(obj$ushape),
                                      "self-defined",obj$ushape),
                      ufix        = obj$ufix,
                      lfix        = obj$lfix,
                      sample.size = obj$sample.size,
                      print       = obj$print,
                      nstart      = obj$nstart,
                      H0          = obj$H0)

  h <- min(obj$r0) # check that here we are not using r0[1]
  r_norm <- obj$r / h
  r0_norm <- obj$r0 / h
  res$rMat  <- rbind(r0_norm, matrix(r_norm, obj$Kv[1], obj$J, byrow = TRUE))

  res$N           = ceiling(n*r0_norm[obj$J]) + sum(ceiling(n*r_norm*Kv_diff))
  dimnames(res$rMat) <- list(
  c("Control", paste0("T", 1:obj$K[1])),
  paste("Stage", 1:obj$J)
)


  if (obj$sample.size) {
    res$power <- obj$power
  } else {
    res$power <- NA
  }
  res$type    <- obj$type

  res$par = list(p=obj$p, p0=p0, delta=delta, delta0=delta0,
              sigma=sig,
              ushape=ifelse(is.function(obj$ushape),"self-defined",obj$ushape),
              lshape=ifelse(is.function(obj$lshape),"self-defined",obj$lshape))

  class(res)<-"MAMS"
  attr(res, "mc") <- attr(obj, "mc")
  attr(res, "method") <- obj$method
  #############################################################
  ##  simulation to define operating characteristics
  #############################################################

  res$nsim  = obj$nsim
  res$ptest = 1
  if (obj$sample.size) {
      if (is.numeric(obj$nsim)) {
        if (obj$print) {
          message("   iii) run simulation \n",appendLF=FALSE)
        }
        sim  = mams.sim.dtl(obj=res,nsim=obj$nsim,ptest=1,H0=obj$H0, 
                                  K = obj$K)

        sim  <- unpack_object(sim)
        res$sim  <- sim$sim
      } else {
        res$sim = NULL
      }
  } else {
      res$sim = NULL
  }
  class(res)<-"MAMS"
  attr(res, "mc") <- attr(obj, "mc")
  attr(res, "method") <- obj$method

  #############################################################
  ##  out
  #############################################################
  
  return(pack_object(res))
}

###############################################################################
###################### simulation function ####################################
###############################################################################

# 'dtl.sim' evaluates 'dtl_sim' x times and computes average
mams.sim.dtl <- function(obj=NULL, nsim = NULL, K = NULL, nMat = NULL,
                              u = NULL, l = NULL, pv = NULL,
                              deltav = NULL, sd = NULL, ptest = NULL, H0=NULL) {

  if (!is.null(obj) & !is.null(obj$input)) {
  obj  <- unpack_object(obj)
  }
  
res  <- list()
if (!is.null(deltav) | !is.null(pv)) {
  attr(res, "altered") <- "mams.sim"
}
  defaults <- list(
    K = c(4, 1), 
    nsim = 50000,
    nMat = matrix(c(13, 26), 2, 5),
    u = c(NA, 2.169),
    l = c(NA, 2.169),
    pv = NULL,
    deltav = NULL,
    sd = NULL,
    ptest = 1,
    H0 = TRUE
  )

  user_defined <- list(
    K = K,
    nsim = nsim,
    nMat = nMat,
    u = u,
    l = l,
    pv = pv,
    deltav = deltav,
    sd = sd,
    ptest = ptest,
    H0 = H0
  )

  if (is.numeric(pv) & is.numeric(deltav)) {
    stop("Specify the effect sizes either via 'pv' or via 'deltav' and 'sd', and
          set the other argument(s) to NULL.")
  }

  user_defined <- Filter(Negate(is.null), user_defined)

  # Initialize par list
  par <- list()


  if (is.null(obj)) {
    Kv  <- user_defined$K
    K  <- Kv[1]
    final_params <- modifyList(defaults, user_defined)
    K <- ncol(final_params$nMat) - 1
        if (is.null(final_params$pv) & is.null(final_params[["deltav"]])) {
      final_params$pv <- c(0.75, rep(0.5, ifelse(K > 1, K - 1, K)))
    }

    J <- ncol(t(final_params$nMat))
    par <- final_params
  }
  # If MAMS object is provided
  if (!is.null(obj)) {
    if (!inherits(obj, "MAMS")) {
      stop("Object specified under 'obj' has to be of class 'MAMS'")
    }

    par <- obj$par
    J <- obj$J
    K <- obj$K[1]
    Kv  <- obj$K
    obj_params <- list(
      nsim = obj$nsim,
      ptest = obj$ptest,
      nMat = t(obj$rMat * obj$n),
      u = obj$u,
      l = obj$l,
      H0 = ifelse(!is.null(obj$H0),obj$H0, obj$par$H0),
      sd = obj$par$sig,
      pv = if (is.null(pv) & is.null(deltav)) {
        if (is.numeric(par$p) & is.numeric(par$p0)) {
          c(par$p, rep(par$p0, K - 1))
        } else if (is.numeric(par$pv)) {
          par$pv
        } 
      } else {
        pv
      },
      deltav = if (is.null(pv) & is.null(deltav)) {
        if (is.numeric(par[["delta"]]) & is.numeric(par[["delta0"]])) {
          c(par[["delta"]], rep(par[["delta0"]], K - 1))
        } else if (is.numeric(par[["deltav"]])) {
          par[["deltav"]]
        } else {
              stop("Please provide either pv or deltav or delta, delta0 or 
              p, p0 parameters")
        }
      } else {
        deltav
      }
    )

    # Merge object parameters with user-defined values
    final_params <- modifyList(obj_params, user_defined)

    par <- final_params
  }
nMat  <- if (length(par$nMat) == 0) {
  NULL
} else {
  par$nMat
}
  # sample sizes
  if (is.null(nMat)) {
    if (!is.null(obj)) {
      if (is.null(obj$n)) {
        stop("Either provide an entry to the argument 'nMat' or generate a MAMS
              object with argument 'sample.size = TRUE' ")
      } else {
        nMat <- t(obj$rMat * obj$n)
      }
    } else {
      stop("'nMat' and 'obj' can't both be set to NULL")
    }
  } else {
    if (!is.null(obj)) {
      if (nrow(nMat) != obj$J) {
        stop("number of rows of 'nMat' should match the number of stages
        considered when generating the MAMS object indicated under 'obj'.")
      }
      if (ncol(nMat) != (obj$K[1] + 1)) {
      stop("number of columns of 'nMat' should match the number of groups (K+1)
            considered when generating the MAMS object indicated under 'obj'.")
      }
    }
  }

  # effect sizes
  
  sd <- par$sd[1]
  if (!is.numeric(sd)) {
    if (!is.null(obj)) {
      sd <- par$sigma <- obj$par$sigma
    } else {
      sd <- par$sigma <- 1
      warning("Standard deviation set to 1")
    }
  } else {
    if (sd <= 0) {
      stop("Standard deviation must be positive.")
    }
    par$sigma <- sd
  }
  
  
    # pv
    pv <- par$pv
    deltav <- par[["deltav"]]

    if ((!is.numeric(pv) & !is.null(pv)) |
        (!is.numeric(deltav) & !is.null(deltav))) {
        stop("The parameter 'pv' or 'deltav' should be a numeric vector.")
        }

    if (is.numeric(pv)) {
      if (any(pv < 0) | any(pv > 1)) {
        stop("Treatment effect parameter not within
                                  0 and 1.")
      }
      if (length(pv) != (ncol(nMat) - 1)) {
      stop("Length of pv is not equal to K.")
      }
      deltav <- sqrt(2) * qnorm(pv)
      par$p <- pv[1]
      par$p0 <- pv[2]
      par[["delta"]] <- deltav[1]
      par[["delta0"]] <- deltav[2]
      # deltav
    } 
    if (is.numeric(deltav)) {
      if (length(deltav) != (ncol(nMat) - 1)) stop("Length of deltav is
                                              not equal to K.")
      pv <- pnorm(deltav / (sqrt(2) * sd))
      par[["delta"]] <- deltav[1]
      par[["delta0"]] <- deltav[2]
      par$p <- pv[1]
      par$p0 <- pv[2]
    }
  # limits
  if (is.null(u)) {
    if (!is.null(obj)) {
      u <- obj$u
      par$ushape <- obj$par$ushape
    } else if (!is.null(par$u)) {
    u  <- par$u
    } else {
      stop("'u' and 'obj' can't both be set to NULL")
    }
  } else {
    if (!is.null(obj)) {
      if (length(u) != obj$J) {
        stop("the length of 'u' should match the number of stages considered
          when generating the MAMS object indicated under 'obj'.")
      }
    }
    par$ushape <- "provided"
  }

  if (is.null(l)) {
    if (!is.null(obj)) {
      l <- obj$l
      par$lshape <- obj$par$lshape
    } else if (!is.null(par$l)) {
    l  <- par$l
    } else {
      stop("'l' and 'obj' can't both be set to NULL")
    }
  } else {
    if (!is.null(obj)) {
      if (length(l) != obj$J) {
        stop("the length of 'l' should match the number of stages considered 
        when generating the MAMS object indicated under 'obj'.")
      }
    }
    par$lshape <- "provided"
  }

  # may want to change the definition of l and u above

  # 'dtl_sim' simulates the trial once. For general number of patients per arm
  # per stage - allocation given by the matrix R for active treatment arms.
  # R[i,j] = allocation to stage i treatment j and vector r0 for control arm.
  # Treatment effect is specified by delta, delta0 and sig. Output is:
  # (1) rejection any hypothesis yes/no, (2) rejection first hypothesis yes/no,
  # (3) total sample size, (4) matrix of rejected hypotheses
  dtl_sim <- function(n, Kv, l, u, R, r0, delta, sig) {
    J                            <- dim(R)[1]
    K                            <- dim(R)[2]
    Rdiff                        <- R - rbind(0, R[-J, ])
    r0diff                       <- r0 - c(0, r0[-J])
    all.remaining <- list()
    # Create test statistics using independent normal increments in sample
    # means

    mukhats <- apply(matrix(rnorm(J*K), J, K)*sig*sqrt(Rdiff*n) +
                      Rdiff*n*matrix(delta, J, K, byrow = TRUE), 2,
                                                                cumsum) / (R*n)
    mu0hats <- cumsum(rnorm(J, 0, sig*sqrt(r0diff*n))) / (r0*n)

    dmat <- mukhats - mu0hats
    zks  <- (mukhats - mu0hats) / (sig*sqrt((R + r0) / (R*r0*n)))

    # At each interim (j) determine which treatment are dropped, and at the
    # final analysis which are better than control
    j    <- 1
    # matrix of rejected hypotheses
    hmat <- matrix(0, J, K)
    # matrix of futility and efficacy
    fmat = emat = matrix(NA,nrow=J,ncol=K)

    remaining                    <- rep(TRUE, K)
    # all.remaining[[1]] <- remaining
    ss                           <- 0
    for (j in 1:(J - 1)) {
      ss <- sum((n*Rdiff[j, ])[remaining]) + n*r0diff[j] + ss

      emat[j,remaining] <- ((zks[j,remaining]) > u[j])
      fmat[j,remaining] <- ((zks[j,remaining]) < l[j])

      remaining[order(zks[j, ], decreasing = TRUE)[-(1:Kv[j + 1])]] <- FALSE

      zks[(j + 1):J, !remaining] <- -Inf
    }

    ss <- sum((n*Rdiff[J, ])[remaining]) + n*r0diff[J] + ss
    for (k in (1:K)[remaining]) {
      if (zks[J, k] > u[J]) {
        hmat[J, k] <- 1
      }
    }

    emat[J,remaining] <- ((zks[J,remaining]) > u[J])
    fmat[J,remaining] <- ((zks[J,remaining]) < l[J])

    old.rej <- any(hmat[J, ] == 1)
    pow <- (hmat[J, 1] == 1)
    # any arm > control?
    rej=ifelse(any(emat[J,],na.rm=TRUE),J,0)
# if yes, is T1 also the arm with the largest test statistics
# among remaing arms?
    first=ifelse(rej>0,ifelse(!is.na(emat[J,1])&emat[J,1],
                      zks[J,1]==max(zks[J,remaining]),
                              FALSE),FALSE)

    all.remaining <- cbind(matrix(rep(TRUE, J), nrow = J), zks != -Inf)
    # out
    return(list(stage = rej, rej = old.rej, first = first, pow = pow, ess = ss,
                hmat = hmat, zks = zks,
                remaining = matrix(unlist(all.remaining), nrow = J, ncol = K+1),
                efficacy = emat, futility = fmat))
  }

  ##### Perform the simulation study ###########################################
  if (!is.null(obj)) nMat <- t(obj$rMat*obj$n)
  r0       <- nMat[, 1]/nMat[1, 1]
  if (ncol(nMat) == 2) {
    R      <- t(t(nMat[, -1]/nMat[1, 1]))
  } else {
    R      <- nMat[, -1]/nMat[1, 1]
  }
  if (!is.matrix(R) && is.vector(R)) {
    R      <- t(as.matrix(nMat[, -1]/nMat[1, 1]))
  }
  n        <- nMat[1, 1]
  if (is.numeric(pv) && is.null(deltav)) {
    deltas <- sqrt(2)*stats::qnorm(pv) 
    sig    <- 1
  } else {
    deltas <- deltav
    sig    <- sd
  }

# Useful information
  n <- nMat[1, 1]
  sig <- sd
  J <- dim(R)[1]
  # K <- dim(R)[2]
  Rdiff <- R - rbind(0, R[-J, , drop = FALSE])
  r0diff <- r0 - c(0, r0[-J])
  nsim = obj$nsim
  H0  <- obj$H0
  Kv  <- obj$K

  ## H1
  if (!all(pv==0.5)) {
    # sim
    H1 = list()

    H1$full<-sapply(rep(n, nsim), dtl_sim, Kv, l, u, R, r0, deltas, sig)
    # main results
    H1$main = list()
    # sample size
    tmp <- lapply(H1$full["remaining",], function(x) x*n*cbind(r0, R))
    # tmp <- sapply(tmp,function(x) apply(x,2,sum))
    tmp <- sapply(tmp, function(x) apply(x, 2, max)) 
    H1$main$ess = data.frame(ess  = apply(tmp,1,mean),
                            sd   = sqrt(apply(tmp,1,var)),
                            low  = apply(tmp,1,quantile,prob=0.025),
                            high = apply(tmp,1,quantile,prob=0.975)
    )
    # futility
    tmp0 = array(unlist(H1$full["futility",]),dim=c(J,K,nsim))
    tmp1 = t(apply(apply(tmp0,2:1,sum,na.rm=TRUE)/nsim,1,cumsum))
    if (J>1) {
      tmp2 = apply(apply(tmp0,c(3,1),sum,na.rm=TRUE),1,cumsum)
      tmp3 = apply(tmp2>0,1,mean)
      tmp4 = apply(tmp2==K,1,mean)
    } else {
      tmp1 = t(tmp1)
      tmp2 = apply(tmp0,c(3,1),sum,na.rm=TRUE)
      tmp3 = mean(tmp2>0)
      tmp4 = mean(tmp2==K)
    }
    H1$main$futility = as.data.frame(rbind(tmp1,tmp3,tmp4))
    dimnames(H1$main$futility)=list(
    c(paste0("T", 1:K, "  rejected"), "Any rejected", "All rejected"),
      paste("Stage",1:J))
    # efficacy
    tmp0 = array(unlist(H1$full["efficacy",]),dim=c(J,K,nsim))
    tmp1 = t(apply(apply(tmp0,2:1,sum,na.rm=TRUE)/nsim,1,cumsum))
    tmp5 = tapply(unlist(H1$full["first",]),unlist(H1$full["stage",]),sum)
    tmp6 = rep(0,J); names(tmp6) = 1:J
    tmp6[names(tmp5)[names(tmp5)!="0"]] = tmp5[names(tmp5)!="0"]
    if (J>1) {
      tmp2 = apply(apply(tmp0,c(3,1),sum,na.rm=TRUE),1,cumsum)
      tmp3 = apply(tmp2>0,1,mean)
      tmp4 = apply(tmp2==K,1,mean)
    } else {
      tmp1 = t(tmp1)
      tmp2 = apply(tmp0,c(3,1),sum,na.rm=TRUE)
      tmp3 = mean(tmp2>0)
      tmp4 = mean(tmp2==K)
    }
    H1$main$efficacy  = as.data.frame(rbind(tmp1,tmp3,cumsum(tmp6)/nsim,tmp4))
    dimnames(H1$main$efficacy)=list(
  c(paste0("T", 1:K, "  rejected"), "Any rejected", "T1  is best", 
            "All rejected"),
      paste("Stage",1:J))
    if (length(ptest)>1) {
      if (J>1) {
        tmp7 = apply(apply(tmp0[,ptest,,drop=FALSE],c(3,1),sum,na.rm=TRUE),1,
                                                                      cumsum)
        tmp8 = apply(tmp7>0,1,mean)
      } else {
        tmp7 = apply(tmp0[,ptest,,drop=FALSE],c(3,1),sum,na.rm=TRUE)
        tmp8 = mean(tmp7>0)
      }
      H1$main$efficacy = rbind(H1$main$efficacy,tmp8)
      rownames(H1$main$efficacy)[nrow(H1$main$efficacy)] =
        paste(paste0("T",ptest),collapse=" AND/OR ")
    }

  } else {
    H1<-NULL
  }

  ## H0
K  <- Kv[1]
  if (all(pv==0.5)|H0) {
    # sim
    H0 = list()
      H0$full<-sapply(rep(n, nsim), dtl_sim, Kv, l, u, R, r0, rep(0,K), sig)

    # main results
    H0$main = list()
    # sample size
    tmp <- lapply(H0$full["remaining",], function(x) x*n)
    tmp <- sapply(tmp,function(x) apply(x,2,sum))
    H0$main$ess = data.frame(ess  = apply(tmp,1,mean),
                            sd   = sqrt(apply(tmp,1,var)),
                            low  = apply(tmp,1,quantile,prob=0.025),
                            high = apply(tmp,1,quantile,prob=0.975)
    )
    # futility
    tmp0 = array(unlist(H0$full["futility",]),dim=c(J,K,nsim))
    tmp1 = t(apply(apply(tmp0,2:1,sum,na.rm=TRUE)/nsim,1,cumsum))
    if (J>1) {
      tmp2 = apply(apply(tmp0,c(3,1),sum,na.rm=TRUE),1,cumsum)
      tmp3 = apply(tmp2>0,1,mean)
      tmp4 = apply(tmp2==K,1,mean)
    } else {
      tmp1 = t(tmp1)
      tmp2 = apply(tmp0,c(3,1),sum,na.rm=TRUE)
      tmp3 = mean(tmp2>0)
      tmp4 = mean(tmp2==K)
    }
    H0$main$futility = as.data.frame(rbind(tmp1,tmp3,tmp4))
    dimnames(H0$main$futility)=list(
    c(paste0("T", 1:K, "  rejected"), "Any rejected", "All rejected"),
    paste("Stage",1:J))
    # efficacy
    tmp0 = array(unlist(H0$full["efficacy",]),dim=c(J,K,nsim))
    tmp1 = t(apply(apply(tmp0,2:1,sum,na.rm=TRUE)/nsim,1,cumsum))
    tmp5 = tapply(unlist(H0$full["first",]),unlist(H0$full["stage",]),sum)
    tmp6 = rep(0,J); names(tmp6) = 1:J
    tmp6[names(tmp5)[names(tmp5)!="0"]] = tmp5[names(tmp5)!="0"]
    if (J>1) {
      tmp2 = apply(apply(tmp0,c(3,1),sum,na.rm=TRUE),1,cumsum)
      tmp3 = apply(tmp2>0,1,mean)
      tmp4 = apply(tmp2==K,1,mean)
    } else {
      tmp1 = t(tmp1)
      tmp2 = apply(tmp0,c(3,1),sum,na.rm=TRUE)
      tmp3 = mean(tmp2>0)
      tmp4 = mean(tmp2==K)
    }
    H0$main$efficacy  = as.data.frame(rbind(tmp1,tmp3,cumsum(tmp6)/nsim,tmp4))
    dimnames(H0$main$efficacy)=list(
  c(paste0("T", 1:K, "  rejected"), "Any rejected", "T1  is best", 
          "All rejected"),
                                    paste("Stage",1:J))
    if (length(ptest)>1) {
      if (J>1) {
        tmp7 = apply(apply(tmp0[,ptest,,drop=FALSE],c(3,1),sum,na.rm=TRUE),1,
                                                                        cumsum)
        tmp8 = apply(tmp7>0,1,mean)
      } else {
        tmp7 = apply(tmp0[,ptest,,drop=FALSE],c(3,1),sum,na.rm=TRUE)
        tmp8 = mean(tmp7>0)
      }
      H0$main$efficacy = rbind(H0$main$efficacy,tmp8)
      rownames(H0$main$efficacy)[nrow(H0$main$efficacy)] =
        paste(paste0("T",ptest),collapse=" AND/OR ")
    }
  } else {
    H0<-NULL
  }
  ##### Output #################################################################

  Kv_diff <- c(Kv[1:(J - 1)] - Kv[2:J], Kv[J])
  res$N <- ceiling(n * r0[obj$J]) + sum(ceiling((n * R)[,1] * Kv_diff))

  res$l <- l
  res$u <- u
  res$n <- n
  res$rMat <- rbind(r0, t(R))
  res$K <- obj$K
  res$J <- dim(R)[1]
  res$alpha <- ifelse(is.null(obj), 0.05, obj$alpha)
  res$alpha.star <- NULL
  res$type <- "normal"
  res$par <- par
  res$nsim <- nsim
  res$ptest <- par$ptest
  res$sim <- list(H0 = H0, H1 = H1)
  res$H0  <- obj$H0

  if (!is.null(H0)) {
    res$typeI <- mean(as.numeric(unlist(H0$full["rej",])), na.rm = TRUE)
  } else {
    res$typeI <- NULL
  }

  if (!is.null(H1)) {
    res$power <- mean(as.numeric(unlist(H1$full["rej",])), na.rm = TRUE)
  } else {
    res$power <- NULL
  }

  res$prop.rej <- ifelse(!is.null(H0), sum(unlist(H0$full["rej",])),
                                      sum(unlist(H1$full["rej",])))/nsim

  res$exss <- ifelse(!is.null(H0), mean(unlist(H0$full["ess", ])),
                                  mean(unlist(H1$full["ess", ])))

  res$sim$H0$full  <-  NULL
  res$sim$H1$full  <-  NULL
  attr(res, "method") <- "dtl"
  class(res) <- "MAMS"
  attr(res, "mc") <- attr(obj, "mc")

  return(pack_object(res))
}

###############################################################################
###################### print function #########################################
###############################################################################
mams.print.dtl  <- function(x,
                                digits = max(3, getOption("digits") - 4),
                                ...) {

  x  <- unpack_object(x)
  cat(paste("Design parameters for a ", x$J, " stage drop-the-losers trial ",
            "with Kv = (", paste(x$Kv, collapse = ", "), ") \n\n", sep = ""))
  if (!is.na(x$power)) {
    res             <- matrix(NA, 2, x$J)
    colnames(res)   <- paste("Stage", 1:x$J)
    if (x$type == "tite") {
      rownames(res) <- c("Cumulative number of events per stage (control):",
                        "Cumulative number of events per stage (active):")
    } else {
      rownames(res) <- c("Cumulative sample size per stage (control):",
                        "Cumulative sample size per stage (active):")
    }
    res[1,]         <- ceiling(x$n*x$rMat[1, ])
    res[2,]         <- ceiling(x$n*x$rMat[2, ])
    print(res)
    if (x$type == "tite") {
      cat(paste("\nRequired total number of events: ", x$N, "\n\n"))
    } else {
      cat(paste("\nRequired total sample size: ", x$N, "\n\n"))
    }
  }
  res               <- matrix(NA, 2, x$J)
  colnames(res)     <- paste("Stage", 1:x$J)
  rownames(res)     <- c("Upper bound:", "Lower bound:")
  res[1, ]          <- round(x$u, digits)
  res[2, ]          <- round(x$l, digits)
  print(res)

  if (!is.null(x$sim)) {
    cat(paste("\n\nSimulated error rates based on ", as.integer(x$nsim),
              " simulations:\n",sep=""))

    res <- matrix(NA,nrow=4,ncol=1)
    hyp   = ifelse(is.null(x$sim$H1),"H0","H1")
    K     = x$K
    ptest = ifelse(length(x$ptest)==1,paste0("T", x$ptest, "  rejected"),
                  paste(paste0("T",x$ptest),collapse=" AND/OR "))

    res[1,1] <- round(x$sim[[hyp]]$main$efficacy["Any rejected",x$J],digits)
    res[2,1] <- round(x$sim[[hyp]]$main$efficacy["T1  is best",x$J],digits)
    res[3,1] <- round(x$sim[[hyp]]$main$efficacy[ptest,x$J],digits)
    res[4,1] <- round(sum(x$sim[[hyp]]$main$ess[,"ess"]),digits)


    if (length(x$ptest)==1) {
      rownames(res) <- c("Prop. rejecting at least 1 hypothesis:",
                        "Prop. rejecting first hypothesis (Z_1>Z_2,...,Z_K)",
                        paste("Prop. rejecting hypothesis ",x$ptest,":",
                              sep=""),"Expected sample size:")
    } else {
      rownames(res) <- c("Prop. rejecting at least 1 hypothesis:",
                        "Prop. rejecting first hypothesis (Z_1>Z_2,...,Z_K)",
                        paste("Prop. rejecting hypotheses ",
                              paste(as.character(x$ptest),collapse=" or "),":",
                              sep=""),"Expected sample size:")
    }
    colnames(res)<-""

    print(res)
  }
  cat("\n")
}

###############################################################################
###################### summary function #######################################
###############################################################################
mams.summary.dtl <- function(object, digits, extended=FALSE, ...) {
  

  object  <- unpack_object(object) 
    
  if (is.null(object$sim)) {
      stop("No simulation data provided")  
}
  object$Kv =   object$K
  object$K <-  object$K[1]

  cli_h1(col_red("MAMS design"))
  ## normal
  if (object$type=="normal") {

    # main
    cli_h2(col_blue("Design characteristics"))
    ulid1 <- cli_ul()
    cli_li("Normally distributed endpoint")
    cli_li("Drop the losers")
    cli_li(paste0(col_red(object$J),ifelse(object$J>1," stages","stage")))
    cli_li(paste0(col_red(object$K),ifelse(object$K>1," treatment arms",
                                              " treatment arm")))
    if (!is.na(object$alpha)) {
      cli_li(paste0(col_red(round(object$alpha*100,2)), col_red("%"),
                                                      " overall type I error"))
    }
    if (!is.na(object$power)) {
      cli_li(paste0(col_red(round(object$power*100,2)), col_red("%"),
      " power of detecting Treatment 1 as the best arm"))
    }
    cli_li("Assumed effect sizes per treatment arm:")
    cli_end(ulid1)
    cat("\n")
    out = data.frame(abbr = paste0("T",1:object$K),
                    row.names = paste0("  Treatment ",1:object$K))
    hyp = NULL

    if (!isTRUE(object$H0)) {
    if (is.null(object$par$pv)&!is.null(object$par[["delta"]])) {
      object$par$pv = pnorm(c(object$par[["delta"]], rep(object$par[["delta0"]],
                                object$Kv[1] - 1)) / (sqrt(2)*object$par$sigma))
    }
    } else {
      if (is.null(object$par$pv)&!is.null(object$par[["delta"]])) {
        object$par$pv = pnorm(rep(object$par[["delta0"]],
                              object$Kv[1]) / (sqrt(2)*object$par$sigma))
      }
    }
    if (any(object$par$pv!=0.5)|!is.null(object$sim$H1)) {
      out = cbind(out, "|" = "|",
                  cohen.d = round(c(object$par[["delta"]] / object$par$sig, 
                                  rep(object$par[["delta0"]] / object$par$sig,
                                  object$Kv[1] - 1)),digits),
                  prob.scale = round(pnorm(c(object$par[["delta"]], 
                              rep(object$par[["delta0"]],
                              object$Kv[1] - 1)) / (sqrt(2)*object$par$sigma)),
                              digits = digits))
      hyp = c(hyp,"H1")
    }
    if (!is.null(object$sim$H0) | (all(object$par$pv==0.5))) {
      out = cbind(out, "|" = "|",
                  cohen.d = round(rep(0,object$K),digits),
                  prob.scale = round(rep(0.5,object$K),digits))
      hyp = c(hyp,"H0")
    }

    space = cumsum(apply(rbind(out,colnames(out)),2,
                                function(x) max(nchar(x)))+1)+12
    main  = paste0(paste0(rep(" ",space[2]),collapse=""),"| ",
                              paste0("Under ",hyp[1]))
    if (length(hyp)==2) {
      main = paste0(main,paste0(rep(" ",space[4]-space[1]-10),collapse=""),"| ",
                    paste0("Under ",hyp[2]))
    }
    cat(main,"\n")
    print(out)
    cat("\n")
# Design table
    cli_h2(col_blue("Arms allocation per stage"))
    header_entries <- c("", paste("Stage", 1:length(object$Kv)))
    treatment_entries <- c("Treatment", as.character(object$Kv))
    control_entries <- c("Control", rep(1, object$J))

    # Determine the maximum widths
    max_widths <- sapply(1:length(header_entries), function(i) {
      max(nchar(header_entries[i]), nchar(treatment_entries[i]),
                                    nchar(control_entries[i]))
    })

    # Create a helper function to center text
    center_text <- function(text, width) {
      pad_length <- (width - nchar(text)) %/% 2
      paste0(strrep(" ", pad_length), text, strrep(" ", pad_length +
                                            (width - nchar(text)) %% 2))
    }

    # Create the header row with dynamic-width columns
    header <- paste0("",
                      sprintf("%-*s", max_widths[1], header_entries[1]),
                      "  ", paste(sapply(2:length(header_entries), function(i) {
                        sprintf("%-*s", max_widths[i], header_entries[i])
                      }), collapse = " | "))

    # Create the treatment arms row with centered numbers and
    # dynamic-width columns
    treatment_row <- paste0(sprintf("%-*s", max_widths[1],
                            treatment_entries[1]),
                            "  ", paste(sapply(2:length(treatment_entries),
                            function(i) {
                              center_text(treatment_entries[i], max_widths[i])
                            }), collapse = " | "))

    # Create the control row with centered text and dynamic-width columns
    control_row <- paste0(sprintf("%-*s", max_widths[1], control_entries[1]),
                          "  ", paste(sapply(2:length(control_entries),
                          function(i) {
                            center_text(control_entries[i], max_widths[i])
                          }), collapse = " | "))

    # Print the table
    cat(header, "\n")
    cat(control_row, "\n")
    cat(treatment_row, "\n", "\n")
        # limits
    cli_h2(col_blue("Limits"))
    out = as.data.frame(matrix(round(c(object$u,object$l),digits),nrow=2,
                                byrow=TRUE,
                                dimnames=list(c("Upper bounds", "Lower bounds"),
                                              paste("Stage",1:object$J))))

      out$shape = c("dtl", "dtl")

    print(out)
    cat("\n")
    # sample sizes
    if (!is.null(object$n)) {
      cli_h2(col_blue("Sample sizes"))
      out = as.data.frame(object$rMat*object$n)
      dimnames(out)=list(c("Control", paste("Treatment",1:object$K)),
                          if (object$J==1) {
                          "  Stage 1"
                          } else {
                            paste("Stage",1:object$J)})
      dimnames(out)[[2]][object$J] <- paste0(dimnames(out)[[2]][object$J],
                                                                "\u2020")
      shift = 12
      if (!is.null(object$sim)) {
        if (!is.null(object$sim$H1)) {
          tmp = cbind(NA,round(object$sim$H1$main$ess[,c("low","ess","high")],
                                digits))
          colnames(tmp) = c("|","low","mid","high")
          out = cbind(out,tmp)
        }
        if (!is.null(object$sim$H0)) {
          tmp = cbind(NA,round(object$sim$H0$main$ess[,c("low","ess","high")],
                                digits))
          colnames(tmp) = c("|","low","mid","high")
          out = cbind(out,tmp)
        }
        out = as.data.frame(apply(out,2,function(x) format(round(x,digits))))
        # !FIXME
        # total <- cumsum(object$n*object$Kv + rep(1, length(object$Kv))*object$n)
                total <- cumsum(object$n/object[["rMat"]][1,1] * 
                                object[["rMat"]][2,1]*object$Kv + 
                        rep(1, length(object$Kv))*object$n)
        total[length(total)] <- object$N
               
        if (object$H0) {
            total <- c(total, NA, " ", format(object$N, nsmall = digits), " ",
                              NA, " ", format(object$N, nsmall = digits), " ")
        } else {
          total <- c(total, NA, " ", format(object$N, nsmall = digits), " ")
        }
        out <- rbind(out,"TOTAL\u2021" = total)


        out[,colnames(out)=="|"] = "|"
        space = cumsum(apply(rbind(out,colnames(out)),2,
                              function(x) max(nchar(x)))+1)
        bar   = which(names(space)=="|")
        if (!is.null(object$sim$H1)&!is.null(object$sim$H0)) {
          cat(paste0(rep(" ",space[bar[1]-1]+shift),collapse=""),
              "| Expected (\u00A7)\n", sep="")
          cat(paste0(rep(" ",shift),collapse=""),"Cumulative",
              paste0(rep(" ",space[bar[1]]-12),sep=""),"| Under H1",
              paste0(rep(" ",space[bar[2]-1]-space[bar[1]-1]-10),collapse=""),
              "| Under H0\n",sep="")
        } else {
          cat(paste0(rep(" ",shift),collapse=""),"Cumulative",
              paste0(rep(" ",space[bar[1]]-12),collapse=""),
              "| Expected (\u00A7)\n",sep="")
        }
      } else {
        out = rbind(out,"TOTAL\u2021" = apply(out,2,sum))
        cat(paste0(rep(" ",shift),collapse=""),"Cumulated\n",sep="")
      }
      print(out)
      cat("\u2020", "Max cumulative size per arm", "\n")
      cat("\u2021", "Based on arms allocation at each stage", "\n")
      cat("\n")
    } else {
      cli_h2("Allocation ratios")

      dimnames(out)=list(c("Control", paste("Treatment",1:object$K)),
                          paste("Stage",1:object$J))
      cat(paste0(rep(" ",shift),collapse=""),"Cumulated\n",sep="")
      print(out)
      cat("\n")
    }
    # Futility
    if (!is.null(object$sim)) {
      cli_h2(col_blue("Futility cumulated probabilities (\u00A7)"))

      if (!is.null(object$sim$H1)) {
            out = round(object$sim$H1$main$futility,digits)
            }
      if (!is.null(object$sim$H0)) {
          out = cbind(out,"|"="|", round(object$sim$H0$main$futility,digits))
      }
      shift = 13
      space = cumsum(apply(rbind(out,colnames(out)),2,
                            function(x) max(nchar(x)))+1)
      bar   = which(names(space)=="|")
      if (!is.null(object$sim$H1)&!is.null(object$sim$H0)) {
        cat(paste0(rep(" ",shift),collapse=""),paste0("Under ",hyp[1]),
            paste0(rep(" ",space[bar[1]-1]-8),sep=""),"| ",
            paste0("Under ", hyp[2]),"\n",sep="")
      }

      print(out)
      cat("\n")
    }
    # Efficacy
    if (!is.null(object$sim)) {
      cli_h2(col_blue("Efficacy cumulated probabilities (\u00A7)"))

      if (!is.null(object$sim$H1)) {
        out = round(object$sim$H1$main$efficacy,digits)
        }
      if (!is.null(object$sim$H0)) {
        out = cbind(out,"|"="|", round(object$sim$H0$main$efficacy,digits))
      }


      shift = max(nchar(rownames(out)))+1
      space = cumsum(apply(rbind(out,colnames(out)),2,
                            function(x) max(nchar(x)))+1)
      bar   = which(names(space)=="|")
      if (!is.null(object$sim$H1)&!is.null(object$sim$H0)) {
        cat(paste0(rep(" ",shift),collapse=""),paste0("Under ",hyp[1]),
            paste0(rep(" ",space[bar[1]-1]-8),sep=""),"| ",
            paste0("Under ", hyp[2]),"\n",sep="")
      }

      print(out)
      cat("\n")

      # estimated power and overall type I error
      ulid1 <- cli_ul()
      if (any(hyp=="H1")) {
        prob = object$sim$H1$main$efficacy["T1  is best",object$J]
        text = paste0("Estimated prob. T1  is best (\u00A7) = ", 
                      round(prob*100,digits),
                      "%, [", paste0(round(qbinom(c(0.025,.975),object$nsim,
                                          prob)/object$nsim*100,digits),
                              collapse=", "),"] 95% CI")
        cli_li(text)
      }
      if (any(hyp=="H0")) {
        prob = object$sim$H0$main$efficacy["Any rejected",object$J]
        text = paste0("Estimated overall type I error (\u00A7) = ",
                      round(prob*100,digits),"%, [",
                      paste0(round(qbinom(c(0.025,.975),object$nsim,
                                          prob)/object$nsim*100,digits),
                              collapse=", "),"] 95% CI")
        cli_li(text)
      }
      cli_end(ulid1)
      #cat("\n")
    }

    # biases
    if (!is.null(object$sim)&extended) {
      cli_h2("Delta expected values (\u00A7)")
      # futility
      cli_h3("After futility stop")
      cat("\n")
      out = data.frame(assumed = object$par[["deltav"]],
                        row.names = paste0("  Treatment ",1:object$K))
      hyp = NULL
      if (any(object$par$pv!=0.5)) {
        out = cbind(out, "|" = "|",
                    round(object$sim$H1$main$bias$futility,digits))
        hyp = c(hyp,"H1")
      }
      if (!is.null(object$sim$H0) | (all(object$par$pv==0.5))) {
        out = cbind(out, "|" = "|",
                    round(object$sim$H0$main$bias$futility,digits))
        hyp = c(hyp,"H0")
      }
      space = cumsum(apply(rbind(out,colnames(out)),2,
                            function(x) max(nchar(x)))+1)+12
      main  = paste0(paste0(rep(" ",space[2]),collapse=""),"| ",paste0("Under ",
                                                                      hyp[1]))
      if (length(hyp)==2) {
        main = paste0(main,paste0(rep(" ",space[4]-space[1]-10),
                      collapse=""),"| ",
                      paste0("Under ",hyp[2]))
      }
      cat(main,"\n")
      print(out)
      # efficacy
      cli_h3("After efficacy stop")
      cat("\n")
      out = data.frame(assumed = object$par[["deltav"]],
                        row.names = paste0("  Treatment ",1:object$K))
      hyp = NULL
      if (any(object$par$pv!=0.5)) {
        out = cbind(out, "|" = "|",
                    round(object$sim$H1$main$bias$efficacy,digits))
        hyp = c(hyp,"H1")
      }
      if (!is.null(object$sim$H0) | (all(object$par$pv==0.5))) {
        out = cbind(out, "|" = "|",
                    round(object$sim$H0$main$bias$efficacy,digits))
        hyp = c(hyp,"H0")
      }
      space = cumsum(apply(rbind(out,colnames(out)),2,
                            function(x) max(nchar(x)))+1)+12
      main  = paste0(paste0(rep(" ",space[2]),collapse=""),"| ",
                      paste0("Under ",hyp[1]))
      if (length(hyp)==2) {
        main = paste0(main,paste0(rep(" ",space[4]-space[1]-10),
                      collapse=""),"| ",
                      paste0("Under ",hyp[2]))
      }
      cat(main,"\n")
      print(out)
    }

    if (!is.null(object$sim$TIME)&extended) {
    cli_h2("Estimated study duration and number of enrolled participants (**)")
      # futility
      cli_h3("Study duration")
      cat("\n")
      tmp = round(object$sim$TIME$time,digits)
      out = as.data.frame(tmp[,1:2])
      for (jw in 2:object$J) {
        out = cbind(out,"|"="|",tmp[, (jw-1)*2+1:2])
      }
      shift = 12
      space = cumsum(apply(rbind(out,colnames(out)),2,
      function(x) max(nchar(x)))+1)
      bar   = which(names(space)=="|")
      cat(paste0(rep(" ",shift),collapse=""),paste0("Stage 1     | Stage 2\n"))
      print(out)
      cat("\n")

      # efficacy
      cli_h3("Number of enrolled participants at end of each stage")
      cat("\n")
      tmp = round(object$sim$TIME$enrolled,digits)
      out = as.data.frame(tmp[,1:2])
      for (jw in 2:object$J) {
        out = cbind(out,"|"="|",tmp[, (jw-1)*2+1:2])
      }
      shift = 12
      space = cumsum(apply(rbind(out,colnames(out)),2,
      function(x) max(nchar(x)))+1)
      bar   = which(names(space)=="|")
      cat(paste0(rep(" ",shift),collapse=""),paste0("Stage 1      | Stage 2\n"))
      print(out)
      cat("\n")
    }

    # simulation
    if (!is.null(object$sim)) {
      cat("\n(\u00A7) Operating characteristics estimated by a simulation\n",
          "   considering",as.integer(object$nsim),"Monte Carlo samples\n")
      if (!is.null(object$sim$TIME)&extended) {
        cat("\n(**) Operating characteristics estimated by a simulation\n",
            "   considering 1000 Monte Carlo samples\n")
      }
    }

    # other types
  } else {
    cat(paste("Design parameters for a ", object$J, " stage trial with ",
              object$K, " treatments\n\n",sep=""))

    if (object$type!="new.bounds") {
      if (!is.null(object$n)) {
        res <- matrix(NA,nrow=2,ncol=object$J)
        colnames(res)<-paste("Stage",1:object$J)
        if (object$type=="tite") {
          rownames(res) <- c("Cumulative number of events per stage (control):",
                            "Cumulative number of events per stage (active):")
        } else {
          rownames(res) <- c("Cumulative sample size per stage (control):",
                            "Cumulative sample size per stage (active):")
        }

        res[1,] <- ceiling(object$n*object$rMat[1,])
        res[2,] <- ceiling(object$n*object$rMat[2,])

        print(res)

        if (object$type=="tite") {
          cat(paste("\nMaximum total number of events: ", object$N,"\n\n"))
        } else {
          cat(paste("\nMaximum total sample size: ", object$N,"\n\n"))
        }

      }
    }

    res <- matrix(NA,nrow=2,ncol=object$J)
    colnames(res)<-paste("Stage",1:object$J)
    rownames(res) <- c("Upper bound:", "Lower bound:")
    res[1,] <- round(object$u,digits)
    res[2,] <- round(object$l,digits)

    print(res)
  }
  cli_rule()
}
