################################################################################
################################ Check input ###################################
################################################################################
mams.check.sep <- function(obj) {

  m <- match(c(
    "K", "J", "alpha", "power", "r", "r0", "p", "p0",
    "delta", "delta0", "sd", "ushape", "lshape", "ufix", "lfix",
    "nstart", "nstop", "sample.size", "Q", "type", "method",
    "parallel", "print", "nsim", "H0", "obj", "par", "sim"
  ), names(obj), 0)
  mc <- obj[c(m)]

  if (mc[["K"]] %% 1 > 0 | mc[["J"]] %% 1 > 0) {
    stop("K and J need to be integers.")
  }
  if (mc[["K"]] < 1 | mc[["J"]] < 1) {
    stop("The number of stages and treatments must be at least 1.")
  }
  if (mc[["Q"]] <= 3) {
    stop("Number of points for integration by quadrature to small or negative.")
  }
  if (mc[["Q"]] > 3 & mc[["Q"]] <= 10) {
    warning(
      "Number of points for integration by quadrature is small which may",
      " result in inaccurate solutions."
    )
  }
  if (mc[["alpha"]] < 0 | mc[["alpha"]] > 1 | mc[["alpha"]] < 0 | 
      mc[["alpha"]] > 1) {
    stop("Error rate or power not between 0 and 1.")
  }
  if (length(mc[["r"]]) != length(mc[["r0"]])) {
    stop(
      "Different length of allocation ratios on control and experimental ",
      "treatments."
    )
  }
  if (length(mc[["r"]]) != mc[["J"]]) {
    stop("Length of allocation ratios does not match number of stages.")
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

  if (is.numeric(mc[["p"]]) & is.numeric(mc[["delta"]]) & 
      is.numeric(mc[["sd"]])) {
    stop(
      "Specify the effect sizes either via p or via (delta, sd) and set the",
      " other parameter(s) to NULL."
    )
  }
  if (is.numeric(mc[["p"]])) {
    if (mc[["p"]] < 0 | mc[["p"]]> 1) {
      stop("Treatment effect parameter not within 0 and 1.")
    }
  } else {
    if (is.numeric(mc[["delta"]]) ||
    (!is.null(mc[["par"]]) && !is.null(obj$par[["delta"]]) 
    && is.numeric(obj$par[["delta"]])) & 
    is.numeric(mc[["sd"]])) {
      if (obj$sd <= 0) {
        stop("Standard deviation must be positive.")
      }
    } else {
      stop(
        "Specify the effect sizes either via p or via (delta, sd) and set ",
        "the other parameter(s) to NULL."
      )
    }
  }
  if (is.function(mc[["ushape"]]) & is.function(mc[["lshape"]])) {
    warning(
      "You have specified your own functions for both the lower and ",
      "upper boundary. Please check carefully whether the resulting ",
      "boundaries are sensible."
    )
  }
  if (!is.function(mc[["ushape"]])) {
    if (!mc[["ushape"]] %in% c("pocock", "obf", "triangular", "fixed")) {
      stop("Upper boundary does not match the available options.")
    }
    if (mc[["ushape"]] == "fixed" & is.null(mc[["ufix"]])) {
      stop("ufix required when using a fixed upper boundary shape.")
    }
  } else {
      ushape  <- mc[["ushape"]]
    b <- ushape(mc[["J"]])
    if (!all(sort(b, decreasing = TRUE) == b)) {
      stop("Upper boundary shape is increasing.")
    }
  }
  if (!is.function(mc[["lshape"]])) {
    if (!mc[["lshape"]] %in% c("pocock", "obf", "triangular", "fixed")) {
      stop("Lower boundary does not match the available options.")
    }
    if (mc[["lshape"]] == "fixed" & is.null(mc[["lfix"]])) {
      stop("lfix required when using a fixed lower boundary shape.")
    }
  } else {
    lshape  <- mc[["lshape"]]
    b <- lshape(mc[["J"]])
    if (!all(sort(b, decreasing = FALSE) == b)) {
      stop("Lower boundary shape is decreasing.")
    }
  }

  if (obj$nsim<1000) {

      stop("The number of simulations should be equal to or greater than 1000.")
  }

  return(obj)
}
###############################################################################
###################### Fit function ###########################################
###############################################################################
mams.fit.sep <- function(obj) {

  mc <- attr(obj, "mc")

  K <- obj$K
  J <- obj$J
  alpha <- obj$alpha
  power <- obj$power
  r <- obj$r
  r0 <- obj$r0
  p <- obj$p
  delta <- obj[["delta"]]
  delta0 <- obj[["delta0"]]
  sd <- obj$sd
  ushape <- obj$ushape
  lshape <- obj$lshape
  ufix <- obj$ufix
  lfix <- obj$lfix
  nstart <- obj$nstart
  nstop <- obj$nstop
  sample.size <- obj$sample.size
  Q <- obj$Q
  type <- obj$type
  H0 <- obj$H0
  print  <-  obj$print
  ##### Initialize internal functions ##########################################

  # 'mesh' creates the points and respective weights to use in the outer
  # quadrature integrals you provide one-dimensional set of points (x) and
  # weights (w) (e.g. midpoint rule, gaussian quadrature, etc) and the number of
  # dimensions (d) and it gives d-dimensional collection of points and weights
  mesh <- function(x, d, w = 1 / length(x) + x * 0) {
    n <- length(x)
    W <- X <- matrix(0, n^d, d)
    for (i in 1:d) {
      X[, i] <- x
      W[, i] <- w
      x <- rep(x, rep(n, length(x)))
      w <- rep(w, rep(n, length(w)))
    }
    w <- exp(rowSums(log(W)))
    list(X = X, w = w)
  }

  prodsum <- function(x, l, u, r, r0, r0diff, J, K, Sigma) {
    # x is vector of dummy variables ( the t_j in gdt paper ), l and u are
    # boundary vectors
    int <- prod(sapply(x, dnorm))
    L <- sqrt(r[1]) / r0[1] * sqrt(r0diff[1]) * x[1] +
      l[1] * sqrt(1 + r[1] / r0[1])
    U <- sqrt(r[1]) / r0[1] * sqrt(r0diff[1]) * x[1] +
      u[1] * sqrt(1 + r[1] / r0[1])
    insum <- stats::pnorm(L)
    if (J > 1) {
      for (j in 2:J) {
        L[j] <- sqrt(r[j]) / r0[j] * (sqrt(r0diff[1:j]) %*% x[1:j]) +
          l[j] * sqrt(1 + r[j] / r0[j])
        U[j] <- sqrt(r[j]) / r0[j] * (sqrt(r0diff[1:j]) %*% x[1:j]) +
          u[j] * sqrt(1 + r[j] / r0[j])
        insum <- insum + mvtnorm::pmvnorm(
          lower = c(L[1:j - 1], -Inf),
          upper = c(U[1:j - 1], L[j]),
          sigma = Sigma[1:j, 1:j]
        )[1]
      }
    }
    int <- int * insum^K
    return(int)
  }

  # 'typeI' performs the outer quadrature integral in type I error equation
  # using 'mesh' and 'prodsum' and calculates the difference with the nominal
  # alpha. The accuracy of the quadrature integral will depend on the choice of
  # points and weights. Here, the number of points in each dimension is an
  # input, Q. The midpoint rule is used with a range of -6 to 6 in each
  # dimension.
  typeI <- function(C, alpha, Q, r, r0, r0diff, J, K, Sigma, ushape, lshape,
                    lfix = NULL, ufix = NULL) {
    # the form of the boundary constraints are determined as functions of C.
    if (!is.function(ushape)) {
      if (ushape == "obf") {
        u <- C * sqrt(r[J] / r)
      } else if (ushape == "pocock") {
        u <- rep(C, J)
      } else if (ushape == "fixed") {
        u <- c(rep(ufix, J - 1), C)
      } else if (ushape == "triangular") {
        u <- C * (1 + r / r[J]) / sqrt(r)
      }
    } else {
      u <- C * ushape(J)
    }

    if (!is.function(lshape)) {
      if (lshape == "obf") {
        l <- c(-C * sqrt(r[J] / r[1:(J - 1)]), u[J])
      } else if (lshape == "pocock") {
        l <- c(rep(-C, J - 1), u[J])
      } else if (lshape == "fixed") {
        l <- c(rep(lfix, J - 1), u[J])
      } else if (lshape == "triangular") {
        if (ushape == "triangular") {
          l <- -C * (1 - 3 * r / r[J]) / sqrt(r)
        } else {
          l <- -C * (1 - 3 * r / r[J]) / sqrt(r) / (-1 * (1 - 3) / sqrt(J))
        }
      }
    } else {
      l <- c(C * lshape(J)[1:(J - 1)], u[J])
    }
    mmp <- mesh((1:Q - 0.5) / Q * 12 - 6, J, rep(12 / Q, Q))
    evs <- apply(mmp$X, 1, prodsum,
      l = l, u = u, r = r, r0 = r0,
      r0diff = r0diff, J = J, K = K, Sigma = Sigma
    )
    truealpha <- 1 - mmp$w %*% evs
      if (print) {
    message(".", appendLF = FALSE)
  }
    return(truealpha - alpha)
  }

  # 'typeII' evaluates the power of a design for a particular given group size
  # n, and returns the difference from the nominal power. It achieves this using
  # mvtnorm::pmvnorm() and the fact that for a separate stopping rule power can
  # be evaluated using methods for a conventional two-arm group-sequential
  # design.
  typeII <- function(n, beta, l, u, r, r0, J, delta, sig, Sigma) {
    delta_sqrt_I <- delta * sqrt(1 / (sig^2 / (r0 * n) + sig^2 / (r * n)))
    
    pi <- stats::pnorm(u[1], mean = delta_sqrt_I[1], lower.tail = FALSE)
    if (J > 1) {
      for (j in 2:J) {
        pi <- pi + mvtnorm::pmvnorm(
          lower = c(l[1:(j - 1)], u[j]),
          upper = c(u[1:(j - 1)], Inf),
          mean = delta_sqrt_I[1:j],
          sigma = Sigma[1:j, 1:j]
        )[1]
      }
    }
    return(1 - beta - pi)
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


  ##### Ensure equivalent allocation ratios yield same sample size #############

  # h <- min(c(r, r0))
  # r <- r / h
  # r0 <- r0 / h

  ##### Create the variance covariance matrix from allocation proportions ######

  bottom <- matrix(r, J, J)
  top <- matrix(rep(r, rep(J, J)), J, J)
  top[upper.tri(top)] <- t(top)[upper.tri(top)]
  bottom[upper.tri(bottom)] <- t(bottom)[upper.tri(bottom)]
  Sigma <- sqrt(top / bottom)

  ##### Create r0diff: the proportion of patients allocated to each particular #
  ##### stage ##################################################################

  r0lag1 <- c(0, r0[1:J - 1])
  r0diff <- r0 - r0lag1

  ##### Find boundaries using 'typeI' ##########################################
  if (print) {
    message("   i) find lower and upper boundaries\n      ", 
            appendLF = FALSE)
  }

  # Quick & dirty fix to enable single-stage design with specification
  if (!is.function(lshape)) {
    if (J == 1 & lshape == "obf") {
      lshape <- "pocock"
    }
  }
  uJ <- NULL
  # making sure that lfix is not larger then uJ
  try(uJ <- uniroot(typeI, c(stats::qnorm(1 - alpha) / 2, 5),
    alpha = alpha,
    Q = Q, r = r, r0 = r0, r0diff = r0diff, J = J, K = K,
    Sigma = Sigma, ushape = ushape, lshape = lshape,
    lfix = lfix, ufix = ufix, tol = 0.001
  )$root, silent = TRUE)
  if (is.null(uJ)) {
    stop("No boundaries can be found.")
  }
  if (!is.function(ushape)) {
    if (ushape == "obf") {
      u <- uJ * sqrt(r[J] / r)
    } else if (ushape == "pocock") {
      u <- rep(uJ, J)
    } else if (ushape == "fixed") {
      u <- c(rep(ufix, J - 1), uJ)
    } else if (ushape == "triangular") {
      u <- uJ * (1 + r / r[J]) / sqrt(r)
    }
  } else {
    u <- uJ * ushape(J)
  }
  if (!is.function(lshape)) {
    if (lshape == "obf") {
      l <- c(-uJ * sqrt(r[J] / r[1:(J - 1)]), u[J])
    } else if (lshape == "pocock") {
      l <- c(rep(-uJ, J - 1), u[J])
    } else if (lshape == "fixed") {
      l <- c(rep(lfix, J - 1), u[J])
    } else if (lshape == "triangular") {
      if (ushape == "triangular") {
        l <- -uJ * (1 - 3 * r / r[J]) / sqrt(r)
      } else {
        l <- -uJ * (1 - 3 * r / r[J]) / sqrt(r) / (-1 * (1 - 3) / sqrt(J))
      }
    }
  } else {
    l <- c(uJ * lshape(J)[1:(J - 1)], u[J])
  }

  ##### Now find sample size for arm 1 stage 1 (n)  using 'typeII'. Sample #####
  ##### sizes for all stages are then determined by r*n and r0*n ###############
  if (obj$print & sample.size) {
        message("\n   ii) perform sample size calculation\n", appendLF = FALSE)
      }
  if (J == 1) {
    if (r0 > r) {
      r <- r / r0
      r0 <- r0 / r0
    }
    rho <- r / (r + r0)
    corr <- matrix(rho, K, K) + diag(1 - rho, K)
    if (obj$K == 1) {
      quan <- qmvnorm(1 - obj$alpha, sigma = 1)$quantile
    } else {
      quan <- mvtnorm::qmvnorm(1 - alpha, corr = corr)$quantile   
    }

if (is.null(obj$p)) {
    p  <- pnorm(delta / (sqrt(2) * obj[["sd"]]))
    } else {
    p <- obj[["p"]]
    }
    n <- ceiling(
            ((quan + qnorm(power)) / (qnorm(p) * sqrt(2)))^2 * (1 + 1 / r))


  } else {
    if (r[1] > r0[1]) {
      r <- r / r0[1]
      r0 <- r0 / r0[1]
      message("Allocation ratio for control arm at first stage greater than for 
      treatment arm(s), using normalisation by r0[1] \n")
    }
    n <- nstart
    pow <- 0
    if (sample.size) {
      if (is.null(nstop)) {
        nx <- nstart
        po <- 0
        while (po == 0) {
          nx <- nx + 1
          po <- (typeII(nx, 1 - power, l, u, r, r0, 1, delta, sig, Sigma) < 0)
        }
        nstop <- 3 * nx
      }
      while (pow == 0 & n <= nstop) {
        n <- n + 1
        pow <- (typeII(n, 1 - power, l, u, r, r0, J, delta, sig, Sigma) < 0)
      }

      n <- n * r0[1]
      if (n - 1 == nstop) {
        warning("The sample size was limited by nstop.")
      }
    } else {
      n <- NULL
    }
  }
  ##### Output #################################################################

  res <- NULL
  res  <- list(K=obj$K, J=obj$J, alpha=obj$alpha, power=obj$power, 
                r=obj$r, r0=obj$r0, p=obj$p, p0=obj$p0, 
                delta=obj[["delta"]], delta0=obj[["delta0"]], sd=obj$sd, 
                ushape=obj$ushape, lshape=obj$lshape, 
                ufix=obj$ufix, lfix=obj$lfix,
                nstart=obj$nstart, nstop=obj$nstop, 
                sample.size=obj$sample.size, Q=obj$Q, 
                type=obj$type, parallel=obj$parallel, print=obj$print, 
                nsim=obj$nsim, H0=obj$H0)
  res$l <- l
  res$u <- u
  res$n <- n
  ## allocation ratios


h <- min(obj$r0) # check that here we are not using r0[1]
r_norm <- obj$r / h
r0_norm <- obj$r0 / h
res$rMat <- rbind(r0_norm, matrix(r_norm, ncol = obj$J, nrow = obj$K, 
                                  byrow = TRUE))
  
dimnames(res$rMat) <- list(
  c("Control", paste0("T", 1:obj$K)),
  paste("Stage", 1:obj$J)
)

  # res$rMat <- rbind(r0, matrix(r, ncol = obj$J, nrow = obj$K, byrow = TRUE))
  ## maximum total sample sizeres$Q <- K*r[J]*n+r0[J]*n ## maximum total
  ## sample size
  res$N <- sum(ceiling(res$rMat[, obj$J] * res$n))
  # res$alpha.star <- alpha.star

  res$type <- obj$type
  res$par <- list(
    p = obj$p, p0 = p0, delta = delta, delta0 = delta0,
    sigma = sig,
    ushape = ifelse(is.function(obj$ushape), "self-defined", obj$ushape),
    lshape = ifelse(is.function(obj$lshape), "self-defined", obj$lshape)
  )

  class(res) <- "MAMS"
  attr(res, "mc") <- attr(obj, "mc")
  attr(res, "method") <- "sep"


  #############################################################
  ##  simulation to define operating characteristics
  #############################################################
  res$ptest <- 1
  if (obj$sample.size) {
    if (is.numeric(res$nsim)) {
        if (obj$print) {
          message("   iii) run simulation \n", appendLF = FALSE)
        }
        sim <- mams.sim.sep(
          obj = res, nsim = res$nsim, ptest = res$ptest,
          H0 = obj$H0
        )

        sim  <- unpack_object(sim)
        res$sim  <- sim$sim
        res$par <- sim$par 
      } else {
        res$sim <- NULL
      }
  } else {
      res$sim = NULL
  }

  #############################################################
  ##  out
  #############################################################
  return(pack_object(res))
}

###############################################################################
###################### simulation function ####################################
###############################################################################

# 'mams_sep.sim' evaluates 'sep_sim' x times and computes average
mams.sim.sep <- function(obj = NULL, nsim = NULL, nMat = NULL,
                              u = NULL, l = NULL, pv = NULL, deltav = NULL,
                              sd = NULL, ptest = NULL, H0 = NULL) {

  if (!is.null(obj) & !is.null(obj$input)) {
  obj  <- unpack_object(obj)
  }

res  <- list()
if (!is.null(deltav) | !is.null(pv)) {
  attr(res, "altered") <- "mams.sim"
}

  defaults <- list(
    nsim = 50000,
    nMat = matrix(c(14, 28), 2, 5),
    u = c(3.068, 2.169),
    l = c(0.000, 2.169),
    pv = NULL,
    deltav = NULL,
    sd = NULL,
    ptest = 1,
    H0 = TRUE
  )
  user_defined <- list(
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
    K <- obj$K
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
        } else if (is.numeric(par$pv) & is.null(pv)) {
          par$pv
        } 
      } else {
        pv
      },
      deltav = if (is.null(pv) & is.null(deltav)) {
        if (is.numeric(par[["delta"]]) & is.numeric(par[["delta0"]])) {
          c(par[["delta"]], rep(par[["delta0"]], K - 1))
        } else if (is.numeric(par[["deltav"]]) & is.null(deltav)) {
          par[["deltav"]]
        } else {
              deltav
        }
      } 
    )

    # Merge object parameters with user-defined values
    final_params <- modifyList(obj_params, user_defined)

    if (is.null(final_params$pv) & is.null(final_params[["deltav"]])) {
      stop("Please provide either pv or deltav or delta, delta0 or 
              p, p0 parameters")
      }

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
      if (ncol(nMat) != (obj$K + 1)) {
      stop("number of columns of 'nMat' should match the number of groups (K+1)
            considered when generating the MAMS object indicated under 'obj'.")
      }
    }
  }

  # effect sizes
  sd <- par$sd
  if (!is.numeric(sd)) {
    if (!is.null(obj)) {
      sd <- sig <-  par$sigma <- obj$par$sigma
    } else {
      sd <-  sig <- par$sigma <- 1
      warning("Standard deviation set to 1")
    }
  } else {
    if (sd <= 0) {
      stop("Standard deviation must be positive.")
    }
    par$sigma <-  sig <- sd
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
    } else {
      stop("'l' and 'obj' can't both be set to NULL")
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

# 'sep_sim' simulates the trial once. For general number of patients per arm
# per stage - allocation given by the matrix R for active treatment arms.
# R[i,j] = allocation to stage i treatment j and vector r0 for control arm.
# Treatment effect is specified by delta, delta0 and sig. Output is:
# (1) rejection any hypothesis yes/no, (2) rejection first hypothesis yes/no,
    # (3) total sample size, (4) matrix of rejected hypotheses
    # ^ORIGINAL
    sep_sim <- function(n, l, u, R, r0, delta, sig) {
      J <- dim(R)[1]
      K <- dim(R)[2]
      Rdiff <- R - rbind(0, R[-J, ])
      r0diff <- r0 - c(0, r0[-J])
      # Create test statistics using independent normal increments in sample
      # means
      mukhats <- apply(
        matrix(rnorm(J * K), J, K) * sig * sqrt(Rdiff * n) +
          Rdiff * n * matrix(delta, J, K, byrow = TRUE), 2,
        cumsum
      ) / (R * n)
      mu0hats <- cumsum(rnorm(J, 0, sig * sqrt(r0diff * n))) / (r0 * n)
      zks <- (mukhats - mu0hats) / (sig * sqrt((R + r0) / (R * r0 * n)))
      # At each interim (j) determine whether or not each treatment has been
      # declared significantly better than control
      eff <- 0
      fut <- 0
      # ss <- 0
      j <- 1
      all.remaining <- matrix(FALSE, J, K)
      control <- matrix(FALSE, J, 1)

      # matrix of rejected hypotheses
      emat <- matrix(0, J, K)
      fmat <- matrix(0, J, K)
      remaining <- rep(TRUE, K)
      while ((eff == 0) && (fut == 0) && (j <= J)) {
        # ss        <- sum((n*Rdiff[j, ])[remaining]) + n*r0diff[j] + ss

        eff <- (min(zks[j, remaining]) > u[j])
        fut <- (max(zks[j, remaining]) < l[j])
        if (any(zks[j, remaining] > u[j])) {
          emat[j, which((zks[j, ] > u[j]) & remaining)] <- 1
        }

        if (any(zks[j, remaining] < l[j])) {
          fmat[j, which((zks[j, ] < l[j]) & remaining)] <- 1
        }


        all.remaining[j, ] <- remaining
        control[j] <- TRUE

        remaining <- ((zks[j, ] > l[j]) & (zks[j, ] < u[j]) & remaining)
        if (all(remaining == FALSE)) {
          break
        }
        j <- ifelse(j < J & (eff == 0 & fut == 0), j + 1, break)
      }


      rej <- ifelse(any(emat == 1, na.rm = TRUE), j, 0)

        first <- apply(zks*emat, 2, function(x) {
        x <- x[!is.na(x) & x != 0]
        if (length(x) == 0) 0 else x[length(x)]
        }
        )
        if (K == 1) {
        first <- any(emat[, 1] == 1)
        } else {
        first= all(first[1]>first[2:length(first)])

        }

      all.remaining <- cbind(control, all.remaining)

      return(list(
        stage = rej, remaining = all.remaining,
        futility = fmat, efficacy = emat, first = first
      ))
    }
  #### Perform the simulation study ###########################################
  # nMat[1,1] (1st stage, 1st group = control) is the reference for r0 and R
  # r0: attribution rate per stage of control
  # R: attribution rate per stage for all groups
  r0 <- nMat[, 1] / nMat[1, 1]
  if (ncol(nMat) == 2) {
    R <- t(t(nMat[, -1] / nMat[1, 1]))
  } else {
    R <- nMat[, -1] / nMat[1, 1]
  }
  if (!is.matrix(R) && is.vector(R)) R <- t(as.matrix(nMat[, -1] / nMat[1, 1]))

  # useful information
  n <- nMat[1, 1]
  deltas <- deltav
  sig <- sd
  J <- dim(R)[1]
  K <- dim(R)[2]
  Rdiff <- R - rbind(0, R[-J, , drop = FALSE])
  r0diff <- r0 - c(0, r0[-J])
  nsim = par$nsim
H0  <- par$H0


  ## H1
  if (!all(pv == 0.5)) {
    # sim
    H1 <- list()

    H1$full <- sapply(rep(n, nsim), sep_sim, l, u, R, r0, deltav, sig)

    # main results
    H1$main <- list()
    # sample size
    tmp <- lapply(H1$full["remaining", ], function(x) x * n)
    tmp <- sapply(tmp, function(x) apply(x, 2, sum))
    H1$main$ess <- data.frame(
      ess = apply(tmp, 1, mean),
      sd = sqrt(apply(tmp, 1, var)),
      low = apply(tmp, 1, quantile, prob = 0.025),
      high = apply(tmp, 1, quantile, prob = 0.975)
    )
    # futility
    tmp0 <- array(unlist(H1$full["futility", ]), dim = c(J, K, nsim))
    tmp1 <- t(apply(apply(tmp0, 2:1, sum, na.rm = TRUE) / nsim, 1, cumsum))
    if (J > 1) {
      tmp2 <- apply(apply(tmp0, c(3, 1), sum, na.rm = TRUE), 1, cumsum)
      tmp3 <- apply(tmp2 > 0, 1, mean)
      tmp4 <- apply(tmp2 == K, 1, mean)
    } else {
      tmp1 <- t(tmp1)
      tmp2 <- apply(tmp0, c(3, 1), sum, na.rm = TRUE)
      tmp3 <- mean(tmp2 > 0)
      tmp4 <- mean(tmp2 == K)
    }
    H1$main$futility <- as.data.frame(rbind(tmp1, tmp3, tmp4))
    dimnames(H1$main$futility) <- list(
      c(paste0("T", 1:K, "  rejected"), "Any rejected", "All rejected"),
      paste("Stage", 1:J)
    )
    # efficacy
    tmp0 <- array(unlist(H1$full["efficacy", ]), dim = c(J, K, nsim))
    tmp1 <- t(apply(apply(tmp0, 2:1, sum, na.rm = TRUE) / nsim, 1, cumsum))
    # tmp1 <- apply(tmp0, 2:1, sum, na.rm = TRUE) / nsim
    tmp5 <- tapply(unlist(H1$full["first", ]), unlist(H1$full["stage", ]), sum)
    tmp6 <- rep(0, J)
    names(tmp6) <- 1:J
    tmp6[names(tmp5)[names(tmp5) != "0"]] <-
      tmp5[names(tmp5)[names(tmp5) != "0"]]
    if (J > 1) {
      tmp2 <- apply(apply(tmp0, c(3, 1), sum, na.rm = TRUE), 1, cumsum)
      tmp3 <- apply(tmp2 > 0, 1, mean)
      tmp4 <- apply(tmp2 == K, 1, mean)
    } else {
      tmp1 <- t(tmp1)
      tmp2 <- apply(tmp0, c(3, 1), sum, na.rm = TRUE)
      tmp3 <- mean(tmp2 > 0)
      tmp4 <- mean(tmp2 == K)
    }

H1$main$efficacy <- as.data.frame(rbind(
      tmp1, tmp3, cumsum(tmp6) / nsim,
      tmp4
    ))
    dimnames(H1$main$efficacy) <- list(
      c(
      paste0("T", 1:K, "  rejected"), "Any rejected", "T1  is best", 
      "All rejected"
      ),
      paste("Stage", 1:J)
    )
    if (length(ptest) > 1) {
      if (J > 1) {
        tmp7 <- apply(apply(tmp0[, ptest, , drop = FALSE], c(3, 1), sum,
          na.rm = TRUE
        ), 1, cumsum)
        tmp8 <- apply(tmp7 > 0, 1, mean)
      } else {
        tmp7 <- apply(tmp0[, ptest, , drop = FALSE], c(3, 1), sum, na.rm = TRUE)
        tmp8 <- mean(tmp7 > 0)
      }
      H1$main$efficacy <- rbind(H1$main$efficacy, tmp8)
      rownames(H1$main$efficacy)[nrow(H1$main$efficacy)] <-
        paste(paste0("T", ptest), collapse = " AND/OR ")
    }
    H1$full  <- NULL
  } else {
    H1 <- NULL
  }

  # H0
  if (all(pv == 0.5) | H0) {
    # sim
    H0 <- list()
    H0$full <- sapply(rep(n, nsim), sep_sim, l, u, R, r0, 0, sig)
    # main results
    H0$main <- list()
    # sample size
    tmp <- lapply(H0$full["remaining", ], function(x) x * n)
    tmp <- sapply(tmp, function(x) apply(x, 2, sum))
    H0$main$ess <- data.frame(
      ess = apply(tmp, 1, mean),
      sd = sqrt(apply(tmp, 1, var)),
      low = apply(tmp, 1, quantile, prob = 0.025),
      high = apply(tmp, 1, quantile, prob = 0.975)
    )
    # futility
    tmp0 <- array(unlist(H0$full["futility", ]), dim = c(J, K, nsim))
    tmp1 <- t(apply(apply(tmp0, 2:1, sum, na.rm = TRUE) / nsim, 1, cumsum))

    if (J > 1) {
      tmp2 <- apply(apply(tmp0, c(3, 1), sum, na.rm = TRUE), 1, cumsum)
      tmp3 <- apply(tmp2 > 0, 1, mean)
      tmp4 <- apply(tmp2 == K, 1, mean)
    } else {
      tmp1 <- t(tmp1)
      tmp2 <- apply(tmp0, c(3, 1), sum, na.rm = TRUE)
      tmp3 <- mean(tmp2 > 0)
      tmp4 <- mean(tmp2 == K)
    }
    H0$main$futility <- as.data.frame(rbind(tmp1, tmp3, tmp4))
    dimnames(H0$main$futility) <- list(
      c(paste0("T", 1:K, "  rejected"), "Any rejected", "All rejected"),
      paste("Stage", 1:J)
    )
    # efficacy
    tmp0 <- array(unlist(H0$full["efficacy", ]), dim = c(J, K, nsim))
    tmp1 <- t(apply(apply(tmp0, 2:1, sum, na.rm = TRUE) / nsim, 1, cumsum))
    tmp5 <- tapply(unlist(H0$full["first", ]), unlist(H0$full["stage", ]), sum)
    tmp6 <- rep(0, J)
    names(tmp6) <- 1:J
    tmp6[names(tmp5)[names(tmp5) != "0"]] <- tmp5[names(tmp5) != "0"]

    if (J > 1) {
      tmp2 <- apply(apply(tmp0, c(3, 1), sum, na.rm = TRUE), 1, cumsum)
      tmp3 <- apply(tmp2 > 0, 1, mean)
      tmp4 <- apply(tmp2 == K, 1, mean)
    } else {
      tmp1 <- t(tmp1)
      tmp2 <- apply(tmp0, c(3, 1), sum, na.rm = TRUE)
      tmp3 <- mean(tmp2 > 0)
      tmp4 <- mean(tmp2 == K)
    }
    H0$main$efficacy <- as.data.frame(rbind(
      tmp1, tmp3, cumsum(tmp6) / nsim,
      tmp4
    ))
    dimnames(H0$main$efficacy) <- list(
      c(
      paste0("T", 1:K, "  rejected"), "Any rejected", "T1  is best",
      "All rejected"
      ),
      paste("Stage", 1:J)
    )
    if (length(ptest) > 1) {
      if (J > 1) {
        tmp7 <- apply(apply(tmp0[, ptest, , drop = FALSE],
          c(3, 1), sum,
          na.rm = TRUE
        ), 1, cumsum)
        tmp8 <- apply(tmp7 > 0, 1, mean)
      } else {
        tmp7 <- apply(tmp0[, ptest, , drop = FALSE], c(3, 1), sum, na.rm = TRUE)
        tmp8 <- mean(tmp7 > 0)
      }
      H0$main$efficacy <- rbind(H0$main$efficacy, tmp8)
      rownames(H0$main$efficacy)[nrow(H0$main$efficacy)] <-
        paste(paste0("T", ptest), collapse = " AND/OR ")
    }
    H0$full  <- NULL
  } else {
    H0 <- NULL
  }

  ##### Output #################################################################
  res$l <- l
  res$u <- u
  res$n <- n
  res$rMat <- rbind(r0, t(R))
  res$K <- dim(R)[2]
  res$J <- dim(R)[1]
  res$N <- sum(res$rMat[, res$J] * res$n)
  res$alpha <- ifelse(is.null(obj), 0.05, obj$alpha)
  res$alpha.star <- NULL
  res$power <- ifelse(is.null(obj), 0.9, obj$power)
  res$type <- "normal"
  res$par <- par
  res$nsim <- nsim
  res$ptest <- par$ptest
  res$sim <- list(H0 = H0, H1 = H1)
  res$sample.size  <- TRUE

  class(res) <- "MAMS"
  attr(res, "mc") <- attr(obj, "mc")
  if (!is.null(attr(obj, "method"))) {
    attr(res, "method") <- attr(obj, "method")
  } else {
    attr(res, "method") <- "sep"
  }

  return(pack_object(res))
}

###############################################################################
###################### print function #########################################
###############################################################################

mams.print.sep <- function(x,
                               digits = max(3, getOption("digits") - 4),
                               ...) {

  x  <- unpack_object(x)

  cat(paste("Design parameters for a ", x$J,
    " stage separate stopping rules trial ",
    "with K = ", x$K, "\n\n",
    sep = ""
  ))
  if (!isTRUE(x$sample.size)) {
  res <- matrix(NA, nrow = 2, ncol = x$J)
  colnames(res) <- paste("Stage", 1:x$J)
  rownames(res) <- c("Upper bound:", "Lower bound:")
  res[1, ] <- round(x$u, digits)
  res[2, ] <- round(x$l, digits)

  print(res)
} else {
  if (!is.na(x$power)) {
    res <- matrix(NA, 2, x$J)
    colnames(res) <- paste("Stage", 1:x$J)
    if (x$type == "tite") {
      rownames(res) <- c(
        "Cumulative number of events per stage (control):",
        "Cumulative number of events per stage (active):"
      )
    } else {
      rownames(res) <- c(
        "Cumulative sample size per stage (control):",
        "Cumulative sample size per stage (active):"
      )
    }
    res[1, ] <- ceiling(x$n * x$rMat[1, ])
    res[2, ] <- ceiling(x$n * x$rMat[2, ])
    print(res)
    if (x$type == "tite") {
      cat(paste("\nRequired total number of events: ", x$N, "\n\n"))
    } else {
      cat(paste("\nRequired total sample size: ", x$N, "\n\n"))
    }
  }
  res <- matrix(NA, 2, x$J)
  colnames(res) <- paste("Stage", 1:x$J)
  rownames(res) <- c("Upper bound:", "Lower bound:")
  res[1, ] <- round(x$u, digits)
  res[2, ] <- round(x$l, digits)
  print(res)

  if (!is.null(x$sim)) {
    cat(paste("\n\nSimulated error rates based on ", as.integer(x$nsim),
      " simulations:\n",
      sep = ""
    ))

    res <- matrix(NA, nrow = 4, ncol = 1)
    hyp <- ifelse(is.null(x$sim$H1), "H0", "H1")
    K <- x$K
    ptest <- ifelse(length(x$ptest) == 1, paste0("T", x$ptest, "  rejected"),
      paste(paste0("T", x$ptest), collapse = " AND/OR ")
    )
    res[1, 1] <- round(x$sim[[hyp]]$main$efficacy["Any rejected", x$J], digits)
    res[2, 1] <- round(x$sim[[hyp]]$main$efficacy["T1  is best", x$J], digits)
    res[3, 1] <- round(x$sim[[hyp]]$main$efficacy[ptest, x$J], digits)
    res[4, 1] <- round(sum(x$sim[[hyp]]$main$ess[, "ess"]), digits)

    if (length(x$ptest) == 1) {
      rownames(res) <- c(
        "Prop. rejecting at least 1 hypothesis:",
        "Prop. rejecting first hypothesis (Z_1>Z_2,...,Z_K)",
        paste("Prop. rejecting hypothesis ", x$ptest, ":",
          sep = ""
        ), "Expected sample size:"
      )
    } else {
      rownames(res) <- c(
        "Prop. rejecting at least 1 hypothesis:",
        "Prop. rejecting first hypothesis (Z_1>Z_2,...,Z_K)",
        paste("Prop. rejecting hypotheses ",
          paste(as.character(x$ptest), collapse = " or "), ":",
          sep = ""
        ), "Expected sample size:"
      )
    }
    colnames(res) <- ""

    print(res)
  }
  }
  cat("\n")
}

###############################################################################
###################### summary function #######################################
###############################################################################

mams.summary.sep <- function(object, digits, extended = FALSE, ...) {
  
  object  <- unpack_object(object)     
if (is.null(object$sim)) {
      stop("No simulation data provided")  
}
  cli_h1(col_red("MAMS design"))
  ## normal
  if (object$type == "normal") {
    # main
    cli_h2(col_blue("Design characteristics"))
    ulid1 <- cli_ul()
    cli_li("Normally distributed endpoint")
    cli_li("Separate stopping rules")
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
    out <- data.frame(
      abbr = paste0("T", 1:object$K),
      row.names = paste0("  Treatment ", 1:object$K)
    )
    hyp <- NULL
  if (!is.null(object$sim$H1) | (object$par$p!=0.5)) {
      if (is.null(object$par[["deltav"]])) {
      object$par[["deltav"]] = sqrt(2) * qnorm(object$par$pv)
      }
        out = cbind(out, "|" = "|",
          cohen.d = round(object$par[["deltav"]] / object$par$sig, digits),
          prob.scale = round(pnorm(object$par[["deltav"]] / 
                            (sqrt(2)*object$par$sd)),
                              digits))
            hyp = c(hyp,"H1")
        }

    if (!is.null(object$sim$H0) | (all(object$par$pv == 0.5))) {
      out <- cbind(out,
        "|" = " |",
        cohen.d = round(rep(0, object$K), digits),
        prob.scale = round(rep(0.5, object$K), digits)
      )
      hyp <- c(hyp, "H0")
    }
    space <- cumsum(apply(
      rbind(out, colnames(out)), 2,
      function(x) max(nchar(x))
    ) + 1) + 12
    main <- paste0(
      paste0(rep(" ", space[2]), collapse = ""), "| ",
      paste0("Under ", hyp[1])
    )
    if (length(hyp) == 2) {
      main <- paste0(
        main, paste0(rep(" ", space[4] - space[1] - 9), collapse = ""), "| ",
        paste0("Under ", hyp[2])
      )
    }
    cat(main, "\n")
    print(out)

    # limits
    cli_h2(col_blue("Limits"))
    out <- as.data.frame(matrix(round(c(object$u, object$l), digits),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(
        c("Upper bounds", "Lower bounds"),
        paste("Stage", 1:object$J)
      )
    ))
    if (!all(c(object$par$ushape, object$par$lshape) == "provided")) {
      out$shape <- c(
        ifelse(is.function(object$par$ushape),
          "self-designed", object$par$ushape
        ),
        ifelse(is.function(object$par$lshape),
          "self-designed", object$par$lshape
        )
      )
    }
    print(out)
    cat("\n")
    # sample sizes
    if (!is.null(object$n)) {
      cli_h2(col_blue("Sample sizes"))
      out <- as.data.frame(object$rMat * object$n)
      dimnames(out) <- list(
        c("Control", paste("Treatment", 1:object$K)),
        if (object$J == 1) {
          "  Stage 1"
        } else {
          paste("Stage", 1:object$J)
        }
      )
      shift <- 12
      if (!is.null(object$sim)) {
        if (!is.null(object$sim$H1)) {
          tmp <- cbind(NA, round(
            object$sim$H1$main$ess[, c("low", "ess", "high")],
            digits
          ))
          colnames(tmp) <- c("|", "low", "mid", "high")
          out <- cbind(out, tmp)
        }
        if (!is.null(object$sim$H0)) {
          tmp <- cbind(NA, round(
            object$sim$H0$main$ess[, c("low", "ess", "high")],
            digits
          ))
          colnames(tmp) <- c("|", "low", "mid", "high")
          out <- cbind(out, tmp)
        }
        out <- as.data.frame(apply(
          rbind(out, "TOTAL" = apply(out, 2, sum)), 2,
          function(x) format(round(x, digits))
        ))
        out[, colnames(out) == "|"] <- "|"
        space <- cumsum(apply(
          rbind(out, colnames(out)), 2,
          function(x) max(nchar(x))
        ) + 1)
        bar <- which(names(space) == "|")
        if (!is.null(object$sim$H1) & !is.null(object$sim$H0)) {
          cat(paste0(rep(" ", space[bar[1] - 1] + shift), collapse = ""),
            "| Expected (*)\n",
            sep = ""
          )
          cat(paste0(rep(" ", shift), collapse = ""), "Cumulated",
            paste0(rep(" ", space[bar[1]] - 11), sep = ""), "| Under H1",
            paste0(rep(" ", space[bar[2] - 1] - space[bar[1] - 1] - 10), 
                                                          collapse = ""),
            "| Under H0\n",
            sep = ""
          )
        } else {
          cat(paste0(rep(" ", shift), collapse = ""), "Cumulated",
            paste0(rep(" ", space[bar[1]] - 11), collapse = ""),
            "| Expected (*)\n",
            sep = ""
          )
        }
      } else {
        out <- rbind(out, "TOTAL" = apply(out, 2, sum))
        cat(paste0(rep(" ", shift), collapse = ""), "Cumulated\n", sep = "")
      }
      print(out)
      cat("\n")
    } else {
      cli_h2("Allocation ratios")
      out <- object$rMat
      dimnames(out) <- list(
        c("Control", paste("Treatment", 1:object$K)),
        paste("Stage", 1:object$J)
      )
      cat(paste0(rep(" ", shift), collapse = ""), "Cumulated\n", sep = "")
      print(out)
      cat("\n")
    }

    # Futility
    if (!is.null(object$sim)) {
      cli_h2(col_blue("Futility cumulated probabilities (\u00A7)"))
      out <- round(object$sim[[hyp[1]]]$main$futility, digits)
      if (length(hyp) > 1) {
        out <- cbind(out,
          "|" = "|",
          round(object$sim[[hyp[2]]]$main$futility, digits)
        )
      }
      shift <- max(nchar(rownames(out))) + 1
      space <- cumsum(apply(
        rbind(out, colnames(out)), 2,
        function(x) max(nchar(x))
      ) + 1)
      bar <- which(names(space) == "|")
      if (length(hyp) == 2) {
        cat(paste0(rep(" ", shift), collapse = ""), paste0("Under ", hyp[1]),
          paste0(rep(" ", space[bar[1] - 1] - 8), sep = ""), "| ",
          paste0("Under ", hyp[2]), "\n",
          sep = ""
        )
      }
      print(out)
      cat("\n")
    }

    # Efficacy
    if (!is.null(object$sim)) {
      cli_h2(col_blue("Efficacy cumulated probabilities (\u00A7)"))
      out <- round(object$sim[[hyp[1]]]$main$efficacy, digits)
      if (length(hyp) > 1) {
        out <- cbind(out,
          "|" = "|",
          round(object$sim[[hyp[2]]]$main$efficacy, digits)
        )
      }
      shift <- max(nchar(rownames(out))) + 1
      space <- cumsum(apply(
        rbind(out, colnames(out)), 2,
        function(x) max(nchar(x))
      ) + 1)
      bar <- which(names(space) == "|")
      if (length(hyp) == 2) {
        cat(paste0(rep(" ", shift), collapse = ""), paste0("Under ", hyp[1]),
          paste0(rep(" ", space[bar[1] - 1] - 8), sep = ""), "| ",
          paste0("Under ", hyp[2]), "\n",
          sep = ""
        )
      }
      print(out)
      cat("\n")

      # estimated power and overall type I error
      ulid1 <- cli_ul()
      if (any(hyp == "H1")) {
        prob <- object$sim$H1$main$efficacy["T1  is best", object$J]
        text <- paste0(
          "Estimated prob. T1  is best (\u00A7) = ", 
                                              round(prob * 100, digits), 
          "%, [",
          paste0(
            round(qbinom(
              c(0.025, .975), object$nsim,
              prob
            ) / object$nsim * 100, digits),
            collapse = ", "
          ), "] 95% CI"
        )
        cli_li(text)
      }
      if (any(hyp == "H0")) {
        prob <- object$sim$H0$main$efficacy["Any rejected", object$J]
        text <- paste0(
          "Estimated overall type I error (\u00A7) = ",
          round(prob * 100, digits), "%, [",
          paste0(
            round(qbinom(
              c(0.025, .975), object$nsim,
              prob
            ) / object$nsim * 100, digits),
            collapse = ", "
          ), "] 95% CI"
        )
        cli_li(text)
      }
      cli_end(ulid1)
      # cat("\n")
    }

    # biases
    if (!is.null(object$sim) & extended) {
      cli_h2("Delta expected values (*)")
      # futility
      cli_h3("After futility stop")
      cat("\n")
      out <- data.frame(
        assumed = object$par[["deltav"]],
        row.names = paste0("  Treatment ", 1:object$K)
      )
      hyp <- NULL
      if (any(object$par$pv != 0.5)) {
        out <- cbind(out,
          "|" = "|",
          round(object$sim$H1$main$bias$futility, digits)
        )
        hyp <- c(hyp, "H1")
      }
      if (!is.null(object$sim$H0) | (all(object$par$pv == 0.5))) {
        out <- cbind(out,
          "|" = "|",
          round(object$sim$H0$main$bias$futility, digits)
        )
        hyp <- c(hyp, "H0")
      }
      space <- cumsum(apply(
        rbind(out, colnames(out)), 2,
        function(x) max(nchar(x))
      ) + 1) + 12
      main <- paste0(paste0(rep(" ", space[2]), collapse = ""), "| ", paste0(
        "Under ",
        hyp[1]
      ))
      if (length(hyp) == 2) {
        main <- paste0(
          main, paste0(rep(" ", space[4] - space[1] - 10),
            collapse = ""
          ), "| ",
          paste0("Under ", hyp[2])
        )
      }
      cat(main, "\n")
      print(out)
      # efficacy
      cli_h3("After efficacy stop")
      cat("\n")
      out <- data.frame(
        assumed = object$par[["deltav"]],
        row.names = paste0("  Treatment ", 1:object$K)
      )
      hyp <- NULL
      if (any(object$par$pv != 0.5)) {
        out <- cbind(out,
          "|" = "|",
          round(object$sim$H1$main$bias$efficacy, digits)
        )
        hyp <- c(hyp, "H1")
      }
      if (!is.null(object$sim$H0) | (all(object$par$pv == 0.5))) {
        out <- cbind(out,
          "|" = "|",
          round(object$sim$H0$main$bias$efficacy, digits)
        )
        hyp <- c(hyp, "H0")
      }
      space <- cumsum(apply(
        rbind(out, colnames(out)), 2,
        function(x) max(nchar(x))
      ) + 1) + 12
      main <- paste0(
        paste0(rep(" ", space[2]), collapse = ""), "| ",
        paste0("Under ", hyp[1])
      )
      if (length(hyp) == 2) {
        main <- paste0(
          main, paste0(rep(" ", space[4] - space[1] - 10),
            collapse = ""
          ), "| ",
          paste0("Under ", hyp[2])
        )
      }
      cat(main, "\n")
      print(out)
    }

    if (!is.null(object$sim$TIME) & extended) {
    cli_h2("Estimated study duration and number of enrolled participants (**)")
      # futility
      cli_h3("Study duration")
      cat("\n")
      tmp <- round(object$sim$TIME$time, digits)
      out <- as.data.frame(tmp[, 1:2])
      for (jw in 2:object$J) {
        out <- cbind(out, "|" = "|", tmp[, (jw - 1) * 2 + 1:2])
      }
      shift <- 12
      space <- cumsum(apply(
        rbind(out, colnames(out)), 2,
        function(x) max(nchar(x))
      ) + 1)
      bar <- which(names(space) == "|")
      cat(paste0(rep(" ", shift), collapse = ""), 
                                  paste0("Stage 1     | Stage 2\n"))
      print(out)
      cat("\n")

      # efficacy
      cli_h3("Number of enrolled participants at end of each stage")
      cat("\n")
      tmp <- round(object$sim$TIME$enrolled, digits)
      out <- as.data.frame(tmp[, 1:2])
      for (jw in 2:object$J) {
        out <- cbind(out, "|" = "|", tmp[, (jw - 1) * 2 + 1:2])
      }
      shift <- 12
      space <- cumsum(apply(
        rbind(out, colnames(out)), 2,
        function(x) max(nchar(x))
      ) + 1)
      bar <- which(names(space) == "|")
      cat(paste0(rep(" ", shift), collapse = ""), 
                                  paste0("Stage 1      | Stage 2\n"))
      print(out)
      cat("\n")
    }

    # simulation
    if (!is.null(object$sim)) {
      cat(
        "\n(\u00A7) Operating characteristics estimated by a simulation\n",
        "   considering", as.integer(object$nsim), "Monte Carlo samples\n"
      )
      if (!is.null(object$sim$TIME) & extended) {
        cat(
          "\n(**) Operating characteristics estimated by a simulation\n",
          "   considering 1000 Monte Carlo samples\n"
        )
      }
    }

    # other types
  } else {
    cat(paste("Design parameters for a ", object$J, " stage trial with ",
      object$K, " treatments\n\n",
      sep = ""
    ))

    if (object$type != "new.bounds") {
      if (!is.null(object$n)) {
        res <- matrix(NA, nrow = 2, ncol = object$J)
        colnames(res) <- paste("Stage", 1:object$J)
        if (object$type == "tite") {
          rownames(res) <- c(
            "Cumulative number of events per stage (control):",
            "Cumulative number of events per stage (active):"
          )
        } else {
          rownames(res) <- c(
            "Cumulative sample size per stage (control):",
            "Cumulative sample size per stage (active):"
          )
        }

        res[1, ] <- ceiling(object$n * object$rMat[1, ])
        res[2, ] <- ceiling(object$n * object$rMat[2, ])

        print(res)

        if (object$type == "tite") {
          cat(paste("\nMaximum total number of events: ", object$N, "\n\n"))
        } else {
          cat(paste("\nMaximum total sample size: ", object$N, "\n\n"))
        }
      }
    }

    res <- matrix(NA, nrow = 2, ncol = object$J)
    colnames(res) <- paste("Stage", 1:object$J)
    rownames(res) <- c("Upper bound:", "Lower bound:")
    res[1, ] <- round(object$u, digits)
    res[2, ] <- round(object$l, digits)

    print(res)
  }
  cli_rule()
}

###############################################################################
###################### plot function ##########################################
###############################################################################

mams.plot.sep <- function(x, col = NULL, pch = NULL, lty = NULL,
                              main = NULL, xlab = "Analysis",
                              ylab = "Test statistic", ylim = NULL,
                              type = NULL, las = 1, ...) {

  x  <- unpack_object(x)                               
  if (is.null(type)) type <- "p"
  if (is.null(pch)) pch <- 1
  if (is.null(col)) col <- 1
  if (is.null(lty)) lty <- 2
  if (is.null(las)) las <- 1
  if (is.null(ylim)) {
    r <- range(x$l, x$u)
    ylim <- c(r[1] - diff(r) / 6, r[2] + diff(r) / 6)
  }
  matplot(1:x$J, cbind(x$l, x$u),
    type = type, pch = pch, col = col,
    ylab = ylab, xlab = xlab, ylim = ylim, main = main, axes = FALSE,
    las = las, ...
  )
  mtext(1:x$J, side = 1, at = 1:x$J)
  axis(side = 2, at = seq(-10, 10, 1), las = las)
  lines(x$u, lty = lty)
  lines(x$l[1:(x$J)], lty = lty)
}
