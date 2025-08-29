################################################################################
################################ Check input ###################################
################################################################################
mams.check.simultaneous <- function(obj) {
  # re-order
  m <- match(c(
    "K", "J", "alpha", "power", "r", "r0", "p", "p0",
    "delta", "delta0", "sd", "ushape", "lshape", "ufix", "lfix",
    "nstart", "nstop", "sample.size", "Q", "type", "method",
    "parallel", "print", "nsim", "H0", "obj", "par", "sim"
  ), names(obj), 0)
  mc <- obj[c(m)]
  # sizes
  if (!is.numeric(mc[["K"]]) | !is.numeric(mc[["J"]])) {
    stop("K and J need to be integers.")
  }
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
    warning("Number of points for integration by quadrature is small which may
            result in inaccurate solutions.")
  }
if (is.null(mc[["r"]])) {
    mc[["r"]] <- 1:mc[["J"]]
  }
  if (is.null(mc[["r0"]])) {
    mc[["r0"]] <- 1:mc[["J"]]
  }
  if (length(eval(mc[["r"]])) != length(eval(mc[["r0"]]))) {
    stop("Different length of allocation ratios on control and experimental
      treatments.")
  }
  if (length(eval(mc[["r"]])) != mc[["J"]]) {
    stop("Length of allocation ratios does not match number of stages.")
  }

  if (any(diff(mc[["r0"]]) <= 0)) {
    stop("`r0` must be a monotonically increasing vector.")
  }

  if (any(diff(mc[["r"]]) <= 0)) {
    stop("`r` must be a monotonically increasing vector.")
  }
  if (mc[["r0"]][1] < 1) {
    stop("`r0[1]` must be >= 1.")
  }
  if (any(mc[["r"]] <= 0)) {
    stop("`r` values must be > 0.")
  }
  if (mc[["r0"]][1] %% 1 != 0) {
    stop("First element of `r0` must be integers.")
  }

  # alpha and power
  if (mc[["alpha"]] < 0 | mc[["alpha"]] > 1 |
    mc[["power"]] < 0 | mc[["power"]] > 1) {
    stop("Error rate or power not between 0 and 1.")
  }
  # effect sizes
  if (is.numeric(mc[["p"]]) & is.numeric(mc[["delta"]])) {
    stop("Specify the effect sizes either via (p, p0) or via (delta, delta0, sd)
          and set the other parameters to NULL.")
  }
  #
  if (is.numeric(mc[["p"]])) {
    # p
    if (mc[["p"]] < 0 | mc[["p"]] > 1) {
      stop("Treatment effect parameter 'p0' not within 0.5 and 1.")
    }
    if (mc[["p"]] < 0.5) {
      stop("Interesting treatment effect less than 0.5 which implies that
            reductions in effect over placebo are interesting.")
    }
    # p0
    if (!is.numeric(mc[["p0"]])) {
      mc[["p0"]] <- 0.5
    } else {
      if (mc[["p0"]] < 0 | mc[["p0"]] > 1) {
        stop("Treatment effect parameter 'p0' not within 0 and 1.")
      }
    }
    if (mc[["p"]] <= mc[["p0"]]) {
      stop("Interesting treatment effect must be larger than uninteresting
            effect.")
    }
    if (mc[["p0"]] < 0.5) {
      warning("Uninteresting treatment effect less than 0.5 which implies that
              reductions in effect over placebo are interesting.")
    }
  }
  #
  if (is.numeric(mc[["delta"]])) {
      if (mc[["delta"]] < 0) {
      stop("Negative treatment effect parameter 'delta'.")
    }
    # delta0
    if (!is.numeric(mc[["delta0"]])) {
      mc[["delta0"]] <- 0
      warning("assigned value '0' to argument 'delta0'")
    }
    if (mc[["delta"]] <= mc[["delta0"]]) {
      stop("Interesting treatment effect must be larger than uninteresting
            effect.")
    }
    if (mc[["delta0"]] < 0) {
      warning("Negative uninteresting treatment effect which implies that
              reductions in effect over placebo are interesting.")
    }
    # sd
    if (!is.numeric(mc[["sd"]])) {
      mc[["sd"]] <- 1
      warning("assigned value '1' to argument 'sd'")
    } else {
      if (mc[["sd"]] <= 0) {
        stop("Standard deviation must be positive.")
      }
    }
  }
  # selfdefined ushape and lshape
  if (is.function(mc[["ushape"]]) & is.function(mc[["lshape"]])) {
    warning("You have specified your own functions for both the lower and upper
            boundary. Please check carefully whether the resulting
            boundaries are sensible.")
  }
  # ushape and ufix
  if (!is.function(mc[["ushape"]])) {
    if (!mc[["ushape"]] %in% c("pocock", "obf", "triangular", "fixed")) {
      stop("Upper boundary does not match the available options.")
    }
    if (mc[["ushape"]] == "fixed" & is.null(mc[["ufix"]])) {
      stop("ufix required when using a fixed upper boundary shape.")
    }
  } else {
    b <- do.call(mc[["ushape"]], list(mc[["J"]]))
    if (!all(sort(b, decreasing = TRUE) == b)) {
      stop("Upper boundary shape is increasing.")
    }
  }
  # lshape and lfix
  if (!is.function(mc[["lshape"]])) {
    if (!mc[["lshape"]] %in% c("pocock", "obf", "triangular", "fixed")) {
      stop("Lower boundary does not match the available options.")
    }
    if (mc[["lshape"]] == "fixed" & is.null(mc[["lfix"]])) {
      stop("lfix required when using a fixed lower boundary shape.")
    }
  } else {
    b <- do.call(mc[["lshape"]], list(mc[["J"]]))
    if (!all(sort(b, decreasing = FALSE) == b)) {
      stop("Lower boundary shape is decreasing.")
    }
  }
  # numberenrolled and timetodata
  if (!is.null(mc[["numberenrolled"]])) {
    if (!is.numeric(MAMS.eval(mc[["numberenrolled"]]))) {
      stop("when evaluated, argument 'numberenrolled' should return
          a numeric/integer value")
    }
  }
  if (!is.null(mc[["timetodata"]])) {
    if (!is.numeric(MAMS.eval(mc[["timetodata"]]))) {
      stop("when evaluated, argument 'timetodata' should return
            a numeric/integer value")
    }
  }

  if (mc[["nsim"]]<1000) {

      stop("The number of simulations should be equal to or greater than 1000.")
  }
  # out
  class(mc) <- "MAMS"
  attr(mc, "mc") <- attr(obj, "mc")
  attr(mc, "method") <- attr(obj, "method")
  mc
}

MAMS.eval <- function(expr, ...) {
  eval(parse(text = noquote(expr)), ...)
}

################################################################################
################################ Internal functions ############################
################################################################################

## 'MAMS.mesh' creates the points and respective weights
## to use in the outer quadrature integrals. you provide:
## * one-dimensional set of points (x) and weights (w)
#    (e.g. midpoint rule, gaussian quadrature, etc)
## * the number of dimensions (d)
## d-dimensional collection of points and weights are returned
################################################################################
mams.mesh.simultaneous <- function(x, d, w = 1 / length(x) + x * 0) {
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

################################################################################
## x is vector of dummy variables ( the t_j in gdt paper ),
## l and u are boundary vectors
################################################################################
mams.prodsum1.simultaneous <- function(x, l, u, r, r0, r0diff, J, K,
                                                                        Sigma) {
  .Call("C_prodsum1",
    x2 = as.double(x), l2 = as.double(l), u2 = as.double(u),
    r2 = as.double(r), r02 = as.double(r0), r0diff2 = as.double(r0diff),
    J2 = as.integer(J), K2 = as.integer(K), Sigma2 = Sigma,
    maxpts2 = as.integer(25000), releps2 = as.double(0),
    abseps2 = as.double(0.001), tol2 = as.double(0.0001)
  )
}

################################################################################
## 'MAMS.prodsum2' evaluates the integrand of Pi_1 according to  gdt paper:
################################################################################

mams.prodsum2.simultaneous <- function(x, r, r0, u, K, delta, delta0, n,
                                                                          sig) {
  .Call("C_prodsum2",
    x2 = as.double(x), r2 = as.double(r), r02 = as.double(r0),
    K2 = as.double(K), u2 = as.double(u),
    delta2 = as.double(delta), delta02 = as.double(delta0),
    n2 = as.double(n), sig2 = as.double(sig)
  )
}

################################################################################
## 'MAMS.prodsum3' evaluates the integrand of Pi_j for j>1.
################################################################################

mams.prodsum3.simultaneous <- function(x, l, u, r, r0, r0diff, J, K, delta,
                                            delta0, n, sig, Sigma, SigmaJ) {
  .Call("C_prodsum3",
    x2 = as.double(x), l2 = as.double(l), u2 = as.double(u),
    r2 = as.double(r), r02 = as.double(r0), r0diff2 = as.double(r0diff),
    Jfull2 = as.integer(J), K2 = as.integer(K),
    delta2 = as.double(delta), delta02 = as.double(delta0),
    n2 = as.double(n), sig2 = as.double(sig),
    Sigma2 = Sigma, SigmaJ2 = SigmaJ,
    maxpts2 = as.integer(25000), releps2 = as.double(0),
    abseps2 = as.double(0.001), tol2 = as.double(0.001)
  )
}

################################################################################
##  'typeI' performs the outer quadrature integral in type I error equation
## using 'mesh' and 'prodsum'and calculates the difference with the nominal
## alpha.
##  The accuracy of the quadrature integral will depend on the choice of points
## and weights. Here, the number of points in each dimension is an input, Q.
##  The midpoint rule is used with a range of -6 to 6 in each dimension.
################################################################################

mams.typeI.simultaneous <- function(C, alpha, r, r0, r0diff, J, K, Sigma,
                                        mmp, ushape, lshape, lfix = NULL,
                                        ufix = NULL, parallel = parallel,
                                        print = print) {
  ########################################################################
  ## the form of the boundary constraints are determined as functions of C.
  ########################################################################
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

  if (parallel) {

    evs <- future.apply::future_apply(mmp$X, 1, mams.prodsum1.simultaneous,
      l = l, u = u, r = r, r0 = r0, r0diff = r0diff,
      J = J, K = K, Sigma = Sigma, future.seed = TRUE,
      future.packages = "MAMS"
    )
  } else {
    evs <- apply(mmp$X, 1, mams.prodsum1.simultaneous,
      l = l, u = u, r = r, r0 = r0, r0diff = r0diff, J = J, K = K, Sigma = Sigma
    )
  }
  if (print) {
    message(".", appendLF = FALSE)
  }
  truealpha <- 1 - mmp$w %*% evs
  return(truealpha - alpha)
}

################################################################################
##  'typeII' performs the outer quadrature integrals of Pi_j j=1,...,J  using
##  'mesh' (stored in mmp_j),
##  'prodsum2' and 'prodsum3' as well as summing the Pi_1,...,Pi_J and
##  calculates the difference
##  with the nominal power.
##  The accuracy of the quadrature integral again depends on the choice of
##  points and weights.
##  Here, the number of points in each dimension is an input, N.
##  The midpoint rule is used with a range of -6 to 6 in each dimension.
################################################################################

mams.typeII.simultaneous <- function(n, beta, l, u, r, r0, r0diff, J, K,
                                          delta, delta0,
                                          sig, Sigma, mmp_j,
                                          parallel = parallel) {


  evs <- apply(mmp_j[[1]]$X, 1, mams.prodsum2.simultaneous,
    r = r, r0 = r0, K = K, u = u, delta = delta, delta0 = delta0, n = n,
                                                                  sig = sig
  )
  pi <- mmp_j[[1]]$w %*% evs

  if (J > 1) {
    for (j in 2:J) {
      A <- diag(sqrt(r[j] / (r[j] - r[1:(j - 1)])), ncol = j - 1)
      SigmaJ <- A %*% (Sigma[1:(j - 1), 1:(j - 1)] - Sigma[1:(j - 1), j] %*%
        t(Sigma[1:(j - 1), j])) %*% A

      if (parallel) {
        evs <- future.apply::future_apply(mmp_j[[j]]$X, 1,
          mams.prodsum3.simultaneous,
          l = l, u = u, r = r, r0 = r0, r0diff = r0diff, J = j,
          K = K, delta = delta, delta0 = delta0,
          n = n, sig = sig, Sigma = Sigma[1:j, 1:j],
          SigmaJ = SigmaJ,
          future.seed = TRUE,
          future.packages = "MAMS"
        )
      } else {
        evs <- apply(mmp_j[[j]]$X, 1, mams.prodsum3.simultaneous,
          l = l, u = u, r = r, r0 = r0, r0diff = r0diff, J = j, K = K,
          delta = delta, delta0 = delta0, n = n, sig = sig,
          Sigma = Sigma[1:j, 1:j],
          SigmaJ = SigmaJ
        )
      }
      pi <- pi + mmp_j[[j]]$w %*% evs
    }
  }
  return(1 - beta - pi)
}

###############################################################################
###################### Fit function ###########################################
###############################################################################

mams.fit.simultaneous <- function(obj) {
  mc <- attr(obj, "mc")

  if (obj$print) {
    message("\n", appendLF = FALSE)
  }

  parallel  <- ifelse(is.null(obj$parallel), FALSE, obj$parallel)

  ##############################################################################
  ## Convert treatment effects into absolute effects with standard deviation 1:
  ##############################################################################

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

  ##############################################################################
  ## gaussian quadrature's grid and weights for stages 1:J
  ##############################################################################

  mmp_j <- as.list(rep(NA, obj$J))
  for (j in 1:obj$J) {
    mmp_j[[j]] <- mams.mesh.simultaneous(
      x = (1:obj$Q - .5) / obj$Q * 12 - 6, j,
      w = rep(12 / obj$Q, obj$Q)
    )
  }

mmp_1 <- as.list(rep(NA, 1))
for (j in 1:1) {
  mmp_1[[j]] <- mams.mesh.simultaneous(
    x = (1:obj$Q - .5) / obj$Q * 12 - 6, j,
    w = rep(12 / obj$Q, obj$Q)
  )
}



  ##############################################################################
  ## Ensure equivalent allocation ratios yield same sample size
  ##############################################################################

  r <- obj$r
  r0 <- obj$r0

  ##############################################################################
  ## Create the variance covariance matrix from allocation proportions:
  ##############################################################################
  bottom <- matrix(r, obj$J, obj$J)
  top <- matrix(rep(r, rep(obj$J, obj$J)), obj$J, obj$J)
  top[upper.tri(top)] <- t(top)[upper.tri(top)]
  bottom[upper.tri(bottom)] <- t(bottom)[upper.tri(bottom)]
  Sigma <- sqrt(top / bottom)

  ##############################################################################
  # Create r0diff: the proportion of patients allocated to each particular stage
  ##############################################################################
  r0lag1 <- c(0, r0[1:obj$J - 1])
  r0diff <- r0 - r0lag1
  ##############################################################################
  ## Find boundaries using 'typeI'
  ##############################################################################
  if (obj$print) {
    message("   i) find lower and upper boundaries\n      ", appendLF = FALSE)
  }
  # Quick & dirty fix to enable single-stage design with specification
  # lshape="obf"
  if (!is.function(obj$lshape)) {
    if (obj$J == 1 & obj$lshape == "obf") {
      obj$lshape <- "pocock"
    }
  }

  uJ <- NULL
  ## making sure that lfix is not larger then uJ
  try(uJ <- uniroot(mams.typeI.simultaneous, c(qnorm(1 - obj$alpha) / 2, 5),
    alpha = obj$alpha, r = r, r0 = r0, r0diff = r0diff,
    J = obj$J, K = obj$K, Sigma = Sigma, mmp = mmp_j[[obj$J]],
    ushape = obj$ushape, lshape = obj$lshape, lfix = obj$lfix,
    ufix = obj$ufix,
    parallel = parallel,
    print = obj$print,
    tol = 0.001
  )$root, silent = TRUE)

  if (is.null(uJ)) {
    stop("No boundaries can be found.")
  }

  if (!is.function(obj$ushape)) {
    if (obj$ushape == "obf") {
      u <- uJ * sqrt(r[obj$J] / r)
    } else if (obj$ushape == "pocock") {
      u <- rep(uJ, obj$J)
    } else if (obj$ushape == "fixed") {
      u <- c(rep(obj$ufix, obj$J - 1), uJ)
    } else if (obj$ushape == "triangular") {
      u <- uJ * (1 + r / r[obj$J]) / sqrt(r)
    }
  } else {
    u <- uJ * obj$ushape(obj$J)
  }

  if (!is.function(obj$lshape)) {
    if (obj$lshape == "obf") {
      l <- c(-uJ * sqrt(r[obj$J] / r[1:(obj$J - 1)]), u[obj$J])
    } else if (obj$lshape == "pocock") {
      l <- c(rep(-uJ, obj$J - 1), u[obj$J])
    } else if (obj$lshape == "fixed") {
      l <- c(rep(obj$lfix, obj$J - 1), u[obj$J])
    } else if (obj$lshape == "triangular") {
      if (obj$ushape == "triangular") {
        l <- -uJ * (1 - 3 * r / r[obj$J]) / sqrt(r)
      } else {
        l <- -uJ * (1 - 3 * r / r[obj$J]) / sqrt(r) /
                                            (-1 * (1 - 3) / sqrt(obj$J))
      }
    }
  } else {
    l <- c(uJ * obj$lshape(obj$J)[1:(obj$J - 1)], u[obj$J])
  }

  ##############################################################################
  ## Find alpha.star
  ##############################################################################

  if (obj$print) {
    message("\n  ii) define alpha star\n", appendLF = FALSE)
  }

  alpha.star <- numeric(obj$J)
  alpha.star[1] <- mams.typeI.simultaneous(u[1],
    alpha = 0, r = r[1], r0 = r0[1],
    r0diff = r0diff[1], J = 1, K = obj$K, Sigma = Sigma,
    mmp = mmp_j[[1]], ushape = "fixed", lshape = "fixed",
    lfix = NULL, ufix = NULL, parallel = parallel,
    print = FALSE
  )
  if (obj$J > 1) {
    for (j in 2:obj$J) {
      alpha.star[j] <- mams.typeI.simultaneous(u[j],
        alpha = 0, r = r[1:j], r0 = r0[1:j],
        r0diff = r0diff[1:j], J = j, K = obj$K,
        Sigma = Sigma, mmp = mmp_j[[j]], ushape = "fixed",
        lshape = "fixed", lfix = l[1:(j - 1)],
        ufix = u[1:(j - 1)], parallel = parallel,
        print = FALSE
      )
    }
  }
  ##############################################################################
  ##  Now find samplesize for arm 1 stage 1 (n)  using 'typeII'.
  ##  Sample sizes for all stages are then determined by
  ##  r*n and r0*n.
  ##############################################################################
  if (obj$J == 1 & p0 == 0.5) {
  if (obj$print) {
        message(" iii) perform sample size calculation\n", appendLF = FALSE)
      }
      if (r0 > r) {
      r <- r / r0
      r0 <- r0 / r0
      message("Allocation ratio for control arm at first stage greater than for
      treatment arm(s), using normalisation by r0[1] \n")
    }
    rho <- r / (r + r0)

    corr <- matrix(rho, obj$K, obj$K) + diag(1 - rho, obj$K)
    if (obj$K == 1) {
      quan <- qmvnorm(1 - obj$alpha, mean = rep(0, obj$K), sigma = 1)$quantile
    } else {
      quan <- qmvnorm(1 - obj$alpha, mean = rep(0, obj$K), corr = corr)$quantile
    }

    if (is.null(obj$p)) {
    p  <- pnorm(delta / (sqrt(2) * obj[["sd"]]))
    } else p <- obj[["p"]]
    n <- ceiling(
            ((quan + qnorm(obj$power)) / (qnorm(p) * sqrt(2)))^2 * (1 + 1 / r)) 

  } else {

    if (r[1] > r0[1]) {
      r <- r / r0[1]
      r0 <- r0 / r0[1]
      message("Allocation ratio for control arm at first stage greater than for
      treatment arm(s), using normalisation by r0[1] \n")
    }
    n <- obj$nstart
    ############################################################################
    ## program could be very slow starting at n=0, may want to start at more
    ## sensible lower bound on n unlike type I error, power equation does not
    ## neccessarily have unique solution n therefore search for smallest
    # solution:
    ############################################################################
    pow <- 0
    if (obj$sample.size) {
      if (obj$print) {
        message(" iii) perform sample size calculation\n", appendLF = FALSE)
      }
      if (is.null(obj$nstop)) {
        nx <- obj$nstart
        po <- 0

        while (po == 0) {
          nx <- nx + 1
          po <- (mams.typeII.simultaneous(nx,
            beta = 1 - obj$power, l = l[obj$J], u = u[obj$J],
            r = r[obj$J], r0 = r0[obj$J], r0diff = r0diff[obj$J], J = 1,
            K = obj$K, delta = delta, delta0 = delta0,
            sig = sig, Sigma = matrix(1, 1, 1), 
            mmp_j =  mmp_1,
            parallel = parallel
          ) < 0)
          }

        iterations <- 3 * ceiling(nx)
}
nstop <- if (is.null(obj$nstop)) iterations else obj$nstop
      if (obj$print) {
        message(paste0(
          "      (maximum iteration number = ",
          nstop - obj$nstart + 1, ")\n      "
        ), appendLF = FALSE)
      }
      while (pow == 0 & n <= nstop) {
        n <- n + 1
        pow <- (mams.typeII.simultaneous(n,
          beta = 1 - obj$power, l = l, u = u,
          r = r, r0 = r0, r0diff = r0diff,
          J = obj$J, K = obj$K, delta = delta,
          delta0 = delta0, sig = sig,
          Sigma = Sigma, mmp_j = mmp_j,
          parallel = parallel
        ) < 0)
        if (obj$print) {
          if (any(seq(0, nstop, 50) == n)) {
            message(n, "\n      ",
              appendLF = FALSE
            )
          } else {
            message(".", appendLF = FALSE)
          }
        }
      }

      n <- n * r0[1]

      if (obj$print) {
        message("\n", appendLF = FALSE)
      }
      if ((n - 1) == nstop) {
        warning(paste("The sample size search was stopped because the maximum",
        "sample size (nstop, default: 3 times the sample size of a fixed",
        "sample design) was reached.\n"))
      }
    } else {
      n <- NULL
    }
  }
  #############################################################
  ##  output
  #############################################################

  res <- NULL
  res  <- list(K=obj$K, J=obj$J, alpha=obj$alpha, power=obj$power,
                r=obj$r, r0=obj$r0, p=obj$p, p0=obj$p0,
                delta=obj[["delta"]], delta0=obj[["delta0"]], sd=obj$sd,
                ushape=obj$ushape, lshape=obj$lshape,
                ufix=obj$ufix, lfix=obj$lfix,
                nstart=obj$nstart, nstop=obj$nstop,
                sample.size=obj$sample.size, Q=obj$Q,
                type=obj$type, parallel=parallel, print=obj$print,
                nsim=obj$nsim, H0=obj$H0)
  res$l <- l
  res$u <- u
  res$n <- n
  ## allocation ratios

  h <- min(obj$r0) # check that here we are not using r0[1]
  r_norm <- obj$r / h
  r0_norm <- obj$r0 / h
  # res$rMat <- rbind(r0, matrix(r, ncol = obj$J, nrow = obj$K, byrow = TRUE))
  res$rMat <- rbind(r0_norm, matrix(r_norm, ncol = obj$J, nrow = obj$K,
                                    byrow = TRUE))

  dimnames(res$rMat) <- list(
      c("Control", paste0("T", 1:obj$K)),
      paste("Stage", 1:obj$J)
    )


  ## maximum total sample size
  ## sample size

  res$N <- sum(round(res$rMat[, obj$J] * res$n, digits = 0))


  n_temp <- n / r0[1]
  r0diff <- c(r0[1], diff(r0))
  rdiff <- c(r[1], diff(r))
  rMat <- rbind(r0diff, matrix(rdiff, ncol = obj$J, nrow = obj$K, byrow = TRUE))
  new_rMat <-rMat * n_temp




  res$alpha.star <- alpha.star

  res$type <- obj$type
  res$par <- list(
    p = obj$p, p0 = p0, delta = delta, delta0 = delta0,
    sig = sig,
    ushape = ifelse(is.function(obj$ushape), "self-defined", obj$ushape),
    lshape = ifelse(is.function(obj$lshape), "self-defined", obj$lshape)
  )

  class(res) <- "MAMS"
  attr(res, "mc") <- attr(obj, "mc")
  attr(res, "method") <- "simultaneous"
  #############################################################
  ##  simulation to define operating characteristics
  #############################################################

  # res$nsim <- obj$nsim
  res$ptest <- 1
  if (obj$sample.size) {
  if (is.numeric(res$nsim)) {
        if (obj$print) {
          message("  iv) run simulation \n", appendLF = FALSE)
        }
        sim <- mams.sim.simultaneous(
          obj = res, nsim = res$nsim, ptest = res$ptest,
          parallel = parallel, H0 = obj$H0
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

mams.sim.simultaneous <- function(obj = NULL, nsim = NULL, nMat = NULL,
                                      u = NULL, l = NULL, pv = NULL,
                                      deltav = NULL, sd = NULL, ptest = NULL,
                                      parallel = NULL, H0 = NULL) {

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
    parallel = TRUE,
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
    parallel = parallel,
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
    u   <- par$u
    l   <- par$l
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
      parallel = par$parallel,
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
      sd <- obj$par$sd <- par$sig <- obj$par$sig
    } else {
      sd <- par$sd  <- par$sig <- 1
      warning("Standard deviation set to 1")
    }
  } else {
    if (sd <= 0) {
      stop("Standard deviation must be positive.")
    }
    par$sig <- sd
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
      stop("u' and 'obj' can't both be set to NULL")
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


  ##############################################################################
  ## 'sim' simulates the trial once. For general number of patients per arm per
  ## stage - for active treatment arms, allocation given by the matrix R .
  ##    R[i,j]= allocation to stage i treatment j
  ## - for control arm, allocation given by vector r0
  ## - treatment effect is specified by delta, delta0 and sig
  ##############################################################################


  sim <- function(n, l, u, R, r0, deltas, sig, J, K, Rdiff, r0diff) {
  ##############################################################################
  # Create test statistics using independent normal increments in sample means:
  ##############################################################################
    # directly simulate means per group and time points:
    mukhats <-
      apply(matrix(rnorm(J * K), J, K) * sig * sqrt(Rdiff * n) + Rdiff * n *
        matrix(deltas, nrow = J, ncol = K, byrow = TRUE), 2, cumsum) / (R * n)
    mu0hats <- cumsum(rnorm(J, 0, sig * sqrt(r0diff * n))) / (r0 * n)
    # differences
    dmat <- mukhats - mu0hats
    # test staistics
    zks <- (mukhats - mu0hats) / (sig * sqrt((R + r0) / (R * r0 * n)))
    ## run trial
    fmat <- emat <- matrix(NA, nrow = J, ncol = K)
    nmat <- matrix(0, nrow = J, ncol = K + 1)
    remaining <- rep(TRUE, K)
    for (j in 1:J) {
      nmat[j, c(TRUE, remaining)] <- c(n * r0diff[j], n * Rdiff[j, remaining])
      emat[j, remaining] <- ((zks[j, remaining]) > u[j])
      fmat[j, remaining] <- ((zks[j, remaining]) < l[j])
      remaining <- (zks[j, ] > l[j]) & remaining
      if (any(emat[j, ], na.rm = TRUE) | all(fmat[j, ], na.rm = TRUE)) {
        break
      }
    }


    # any arm > control?
    rej <- ifelse(any(emat[j, ], na.rm = TRUE), j, 0)
    # if yes, is T1 also the arm with the largest test statistics
    # among remaing arms?
    first <- ifelse(rej > 0, ifelse(!is.na(emat[j, 1]) & emat[j, 1],
      zks[j, 1] == max(zks[j, remaining]),
      FALSE
    ), FALSE)
    # out
    return(list(
      any = rej > 0, stage = rej, first = first, efficacy = emat,
      futility = fmat, difference = dmat, statistic = zks, samplesize = nmat
    ))
  }

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
  parallel  <- par$parallel
  parallel  <- ifelse(is.null(obj$parallel) & is.null(parallel),
                                            FALSE, par$parallel)
  ##
  ## simulation
  ##
  ## H1
  if (!all(pv == 0.5)) {
    # sim
    H1 <- list()
    if (parallel) {
      H1$full <- future.apply::future_sapply(rep(n, nsim), sim, l, u, R, r0,
      deltas, sig, J, K, Rdiff, r0diff,
        future.seed = TRUE, future.packages="MAMS"
      )
    } else {
      H1$full <- sapply(rep(n, nsim), sim, l, u, R, r0, deltas, sig, J,
      K, Rdiff, r0diff
      )
    }
    # main results
    H1$main <- list()
    # sample size
    tmp <- sapply(H1$full["samplesize", ], function(x) apply(x, 2, sum))
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
    tmp5 <- tapply(unlist(H1$full["first", ]), unlist(H1$full["stage", ]), sum)
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
    H1$main$efficacy <- as.data.frame(rbind(tmp1, tmp3, cumsum(tmp6) / nsim,
                                                                      tmp4))
    dimnames(H1$main$efficacy) <- list(
      c(
      paste0("T", 1:K, "  rejected"), "Any rejected", "T1  is best",
              "All rejected"
      ),
      paste("Stage", 1:J)
    )
    if (length(ptest) > 1) {
      if (J > 1) {
        tmp7 <- apply(
          apply(tmp0[, ptest, , drop = FALSE], c(3, 1), sum, na.rm = TRUE), 1,
          cumsum
        )
        tmp8 <- apply(tmp7 > 0, 1, mean)
      } else {
        tmp7 <- apply(tmp0[, ptest, , drop = FALSE], c(3, 1), sum, na.rm = TRUE)
        tmp8 <- mean(tmp7 > 0)
      }
      H1$main$efficacy <- rbind(H1$main$efficacy, tmp8)
      rownames(H1$main$efficacy)[nrow(H1$main$efficacy)] <-
        paste(paste0("T", ptest), collapse = " AND/OR ")
    }
    H1$full <- NULL
  } else {
    H1 <- NULL
  }
  ## H0
  if (all(pv == 0.5) | H0) {
    # sim
    H0 <- list()
    if (parallel) {
      H0$full <- future.apply::future_sapply(rep(n, nsim), sim, l, u, R, r0,
                                            rep(0, K), sig, J, K, Rdiff, r0diff,
                                            future.seed = TRUE,
                                            future.packages="MAMS"
      )
  
    } else {
      H0$full <- sapply(rep(n, nsim), sim, l, u, R, r0,
                                            rep(0, K), sig, J, K, Rdiff, r0diff
      )
    }
    # main results
    H0$main <- list()
    # sample size
    tmp <- sapply(H0$full["samplesize", ], function(x) apply(x, 2, sum))
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
    H0$main$efficacy <- as.data.frame(rbind(tmp1, tmp3, cumsum(tmp6) / nsim,
                                                                      tmp4))
    dimnames(H0$main$efficacy) <- list(
      c(
      paste0("T", 1:K, "  rejected"), "Any rejected", "T1  is best",
        "All rejected"
      ),
      paste("Stage", 1:J)
    )
    if (length(ptest) > 1) {
      if (J > 1) {
        tmp7 <- apply(
          apply(tmp0[, ptest, , drop = FALSE], c(3, 1), sum, na.rm = TRUE), 1,
          cumsum
        )
        tmp8 <- apply(tmp7 > 0, 1, mean)
      } else {
        tmp7 <- apply(tmp0[, ptest, , drop = FALSE], c(3, 1), sum, na.rm = TRUE)
        tmp8 <- mean(tmp7 > 0)
      }
      H0$main$efficacy <- rbind(H0$main$efficacy, tmp8)
      rownames(H0$main$efficacy)[nrow(H0$main$efficacy)] <-
        paste(paste0("T", ptest), collapse = " AND/OR ")
    }
    H0$full <- NULL
  } else {
    H0 <- NULL
  }
  ##
  ## output
  ##

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
  res$sample.size  <- TRUE
  res$sim <- list(H0 = H0, H1 = H1)

  class(res) <- "MAMS"
  attr(res, "mc") <- attr(obj, "mc")
  if (!is.null(attr(obj, "method"))) {
    attr(res, "method") <- attr(obj, "method")
  } else {
    attr(res, "method") <- "simultaneous"
  }
  return(pack_object(res))
}

###############################################################################
###################### print function #########################################
###############################################################################

mams.print.simultaneous <- function(x, digits, ...) {
  x  <- unpack_object(x)
  cat(paste("\nDesign parameters for a ", x$J, " stage trial with ", x$K,
    " treatments:\n\n",
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
    res <- matrix(NA, nrow = 2, ncol = x$J)
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
      cat(paste("\nMaximum total number of events: ", x$N, "\n\n"))
    } else {
      cat(paste("\nMaximum total sample size: ", x$N, "\n\n"))
    }
  }

  res <- matrix(NA, nrow = 2, ncol = x$J)
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
mams.summary.simultaneous <- function(object, digits,
                                                  extended = FALSE, ...) {
  object  <- unpack_object(object)

if (object$type %in% c("ordinal", "tite")) {
print(object)
return(invisible(NULL))
}

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
    cli_li("Simultaneous stopping rules")
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
          cohen.d = round(object$par[["deltav"]] / object$par$sig ,digits),
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
###################### plot function #########################################
###############################################################################
mams.plot.simultaneous <- function(x, ask = TRUE, which = 1:2, new = TRUE,
                                        col = NULL, pch = NULL, lty = NULL,
                                        main = NULL, xlab = "Analysis",
                                        ylab = "Test statistic", ylim = NULL,
                                        type = NULL, las = 1, ...) {

  x  <- unpack_object(x)
  # checks:
  which.plot <- c(TRUE, TRUE)
  if (!is.null(which)) {
    which.plot[-which] <- FALSE
  }
  if (is.null(x$sim)) {
    which.plot <- c(TRUE, FALSE)
  }
  if (sum(which.plot) == 1) {
    ask <- FALSE
  }
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  # useful
  col1 <- c("#008B00", "#8B8970", "#8B3A3A")
  col2 <- c("#8B8970", "#008B00", "#8B3A3A")

  ## plot 1
  if (which.plot[1]) {
    # if(!ask){dev.new()}
    if (is.null(ylim)) {
      r <- range(x$l, x$u)
      ylim <- c(r[1] - diff(r) / 6, r[2] + diff(r) / 4)
    }
    if (!new) {
      if (is.null(type)) type <- "p"
      if (is.null(pch)) pch <- 1
      if (is.null(col)) col <- 1
      if (is.null(lty)) lty <- 2
      if (is.null(las)) las <- 1

      matplot(1:x$J, cbind(x$l, x$u),
        type = type, pch = pch, col = col, ylab = ylab,
        xlab = xlab, ylim = ylim, main = main, axes = FALSE, las = las
      ) # , ...)
      mtext(1:x$J, side = 1, at = 1:x$J)
      axis(side = 2, at = seq(-10, 10, 1), las = las)
      lines(x$u, lty = lty)
      lines(x$l[1:(x$J)], lty = lty)
    } else {
      par(mfrow = c(1, 1), mar = c(3, 4, 3, 1))
      plot(1, 1,
        pch = "", ylim = ylim, xlim = c(0.75, x$J + .25), axes = FALSE,
        xlab = "", ylab = ylab
      )
      axis(1, 1:x$J, paste0("Stage ", 1:x$J), tick = FALSE)
      axis(2, seq(-10, 10, .5), las = 2)
      abline(h = seq(-10, 10, .5), col = gray(.95))
      abline(h = 0, col = gray(.75), lwd = 2)
      abline(h = qnorm(.975), col = gray(.75), lwd = 1, lty = 3)
      title(ifelse(is.null(main), "Efficacy and Futility limits", main))
      # legend
      eff <- ifelse(x$par$ushape == "obf", "O'Brien & Fleming",
        ifelse(x$par$ushape == "pocock", "Pocock ",
          ifelse(x$par$ushape == "triangular", "Triangular", "")
        )
      )
      eff <- ifelse(eff == "", "Efficacy limit",
                                                paste0(eff, " efficacy limits"))
      fut <- ifelse(x$par$lshape == "obf", "O'Brien & Fleming",
        ifelse(x$par$lshape == "pocock", "Pocock",
          ifelse(x$par$lshape == "triangular", "Triangular", "")
        )
      )
      fut <- ifelse(fut == "", "Futility limit",
                                                paste0(fut, " futility limits"))
      legend("top",
        ncol = 1, legend = c(eff, fut), col = col1[c(1, 3)], pch = 15, lwd = 2,
        lty = NA, bg = "light gray", box.lwd = NA, pt.cex = 1.5
      )
      # limits
      lines(x$u, col = col1[1], lwd = 2)
      points(x$u, col = rep(c(col1[1], 1), c(x$J - 1, 1)), lwd = 2)
      text(1:x$J, x$u, format(round(x$u, 3)),
        col = rep(c(col1[1], 1), c(x$J - 1, 1)),
        pos = rep(c(3, 4), c(x$J - 1, 1))
      )
      lines(x$l, col = col1[3], lwd = 2)
      points(x$l[1:(x$J - 1)], col = col1[3], lwd = 2)
      text(1:(x$J - 1), x$l[-x$J], format(round(x$l[-x$J], 3)), col = col1[3],
                                                                      pos = 1)
    }
  }
  ## plot 2
  if (which.plot[2]) {
    # if(!ask){dev.new()}
    hyp <- !sapply(x$sim, is.null)
    hyp <- names(hyp)[hyp]
    if (length(hyp) > 1) {
      layout(matrix(c(1, 2, 3, 3), ncol = 2, byrow = TRUE), widths = 1,
                                                            heights = c(1, .15))
    } else {
      layout(matrix(c(1, 2), ncol = 2, byrow = TRUE), widths = c(1, .3),
                                                      heights = 1)
    }
    rem <- "All rejected"
    # if(!any(hyp=="H0")){rem = c(rem, "ANY")}
    col.G <- c(
      "Treatment 1" = "red", "Treatment 2" = "black",
      "Any rejected" = "violet", "T1  is best" = "blue"
    )
    name.G <- c(
      "Treatment 1" = "Treatment 1",
      "Treatment 2" = "Other treatment(s)",
      "Any rejected" = "Any treatment", "T1  is best" = "Treatment 1 is 'best'"
    )
    for (hw in 1:length(hyp)) {
      # inf
      eff <- x$sim[[hyp[hw]]]$main$efficacy
      fut <- x$sim[[hyp[hw]]]$main$futility
      if (x$K > 2) {
        rem <- c(rem, paste0("Treatment ", 3:x$K))
      }
      eff <- eff[is.na(match(rownames(eff), c(rem))), ]
      fut <- fut[is.na(match(rownames(fut), c(rem, "Any rejected"))), ]
      grp <- unique(c(rownames(fut), rownames(eff)))
      n.grp <- length(grp)
      col.g <- col.G[grp]
      #
      par(mar = c(3, 4, 1.5, 4))
      plot(1, 1,
        pch = "", ylim = c(0, 1), xlim = c(0.75, x$J + .25), axes = FALSE,
        xlab = "", ylab = ""
      )
      axis(1, 1:x$J, paste0("Stage ", 1:x$J), tick = FALSE)
      axis(2, seq(0, 1, .1), las = 2)
      axis(2, .5, "Efficacy", padj = -3, tick = FALSE)
      axis(4, seq(0, 1, .1), format(seq(1, 0, -.1)), las = 2, lty = 3)
      axis(4, .5, "Futility", padj = 3, tick = FALSE)
      abline(h = seq(0, 1, .1), col = gray(.95))
      if (length(hyp) > 1) {
        title(paste0("Under ", hyp[hw]))
      }
      for (gw in 1:n.grp) {
        y <- ifelse(grp[gw] == "Treatment 2", 0, 0)
        if (any(rownames(eff) == grp[gw])) {
          lines(unlist(eff[grp[gw], ]) + y, lty = 1, col = col.g[gw])
          points(unlist(eff[grp[gw], ]) + y, col = col.g[gw], cex = .25)
        }
        if (any(rownames(fut) == grp[gw])) {
          lines(1 - unlist(fut[grp[gw], ]) + y, lty = 2, col = col.g[gw])
          points(1 - unlist(fut[grp[gw], ]) + y, col = col.g[gw], cex = .25)
        }
        if (grp[gw] == "T1  is best" & hyp[hw] == "H1") {
          points(x$J, eff[grp[gw], x$J], col = col.g[gw], cex = 1)
          text(x$J, eff[grp[gw], x$J], paste(round(eff[grp[gw], x$J] * 100, 2),
                                                                          "%"),
            col = col.g[gw], pos = 1, cex = .75
          )
        }
        if (grp[gw] == "Any rejected" & hyp[hw] == "H0") {
          points(x$J, eff[grp[gw], x$J], col = col.g[gw], cex = 1)
          text(x$J, eff[grp[gw], x$J], paste(round(eff[grp[gw], x$J] * 100, 2),
                                                                          "%"),
            col = col.g[gw], pos = 3, cex = .75
          )
        }
      }
    }
    # legend
    if (length(hyp) > 1) {
      par(mar = c(0, 4, 0, 4))
    } else {
      par(mar = c(3, 0, 1.5, 0.25))
    }
    plot(1, 1, pch = "", ylim = c(0, 1), xlim = c(0, 1), axes = FALSE,
                                                          xlab = "", ylab = "")
    # rect(-10,-10,10,10,col="light gray",border=NULL)
    legend("top",
      ncol = ifelse(length(hyp) > 1, 4, 1),
      legend = name.G[grp], col = col.g, pch = 15, lwd = 2, lty = NA,
      bg = NA, box.lwd = NA, pt.cex = 1.5, cex = .9
    )
    legend("bottom",
      ncol = ifelse(length(hyp) > 1, 2, 1),
      legend = c("Futility", "Efficacy"), col = 1, pch = NA, lwd = 2,
                lty = c(2, 1), bg = NA, box.lwd = NA, pt.cex = 1.5, cex = .9
    )
  }
}

###############################################################################
###################### new.bounds function ####################################
###############################################################################
#' Function to update boundaries based on observed sample sizes
#' @description The function determines updated boundaries of a multi-arm
#' multi-stage study based on observed number of observations per arm.
#' @param K  Number of experimental treatments (default=3).
#' @param J Number of stages (default=2).
#' @param alpha One-sided familywise error rate (default=0.05).
#' @param nMat Jx(K+1) dimensional matrix of observed/expected sample sizes.
#' Rows correspond to stages and columns to arms. First column is control
#' (default: 2x4 matrix with 10 subjects per stage and arm).
#' @param u Vector of previously used upper boundaries (default=NULL).
#' @param l Vector of previously used upper boundaries (default=NULL).
#' @param ushape Shape of upper boundary. Either a function specifying the
#' shape or one of "pocock", "obf" (the default), "triangular" and "fixed".
#'  See details.
#' @param lshape Shape of lower boundary. Either a function specifying the
#' shape or one of "pocock", "obf", "triangular" and "fixed" (the default).
#' See details.
#' @param ufix 	Fixed upper boundary (default=NULL). Only used if shape="fixed".
#' @param lfix Fixed lower boundary (default=0). Only used if shape="fixed".
#' @param N Number of quadrature points per dimension in the outer integral
#' (default=20).
#' @param parallel 	if TRUE (default), allows parallelisation of the computation
#'  via a user-defined strategy specified by means of the function
#' future::plan(). If not set differently, the default strategy is sequential,
#' which corresponds to a computation without parallelisation.
#' @param print if TRUE (default), indicate at which stage the computation is.
#' @details This function finds the boundaries for a given matrix of sample
#' sizes in multi-arm multi-stage study with K active treatments plus control.
#' The vectors u and l are the boundaries used so far while u.shape and l.shape
#'  specify the shape to the boundaries for the remaining analysis.
#' By specifying u and l as NULL, a design using only the shapes given by
#' ushape and lshape can be found for any sample sizes per stage and arm.
#'
#' The shape of the boundaries (ushape, lshape) are either using the
#' predefined shapes following Pocock (1977), O'Brien & Fleming (1979)
#' or the triangular Test (Whitehead, 1997) using options "pocock", "obf" or
#' "triangular" respectively, are constant (option "fixed") or supplied in as
#' a function. If a function is passed it should require exactly one argument
#' specifying the number of stages and return a vector of the same length.
#' The lower boundary shape is required to be non-decreasing while the upper
#' boundary shape needs to be non-increasing. If a fixed lower boundary is
#' used, lfix must be smaller than \eqn{\Phi^{-1}(1-\alpha)/2}{Phi(1-alpha)/2}
#' to ensure that it is smaller than the upper boundary.
#' @return An object of the class MAMS containing the following components: \cr
#' \item{l}{Lower boundary.}
#' \item{u}{Upper boundary.}
#' \item{n}{Sample size on control in stage 1.}
#' \item{N}{Maximum total sample size.}
#' \item{K}{Number of experimental treatments.}
#' \item{J}{Number of stages in the trial.}
#' \item{alpha}{Familywise error rate.}
#' \item{power}{Power under least favorable configuration.}
#' \item{rMat}{Matrix of allocation ratios. First row corresponds to control
#' and second row to experimental treatments.}
#' @author Thomas Jaki, Dominic Magirr and Dominique-Laurent Couturier
#' @references
#' Jaki T., Pallmann P. and Magirr D. (2019), The R Package MAMS for Designing
#' Multi-Arm Multi-Stage Clinical Trials, Journal of Statistical Software,
#' 88(4), 1-25. Link: doi:10.18637/jss.v088.i04
#'
#' Magirr D., Jaki T. and Whitehead J. (2012), A generalized Dunnett test for
#' multi-arm multi-stage clinical studies with treatment selection, Biometrika,
#'  99(2), 494-501. Link: doi:10.1093/biomet/ass002
#'
#' Magirr D., Stallard N. and Jaki T. (2014), Flexible sequential designs for
#' multi-arm clinical trials, Statistics in Medicine, 33(19), 3269-3279.
#' Link: doi:10.1002/sim.6183
#'
#' Pocock S.J. (1977), Group sequential methods in the design and analysis of
#' clinical trials, Biometrika, 64(2), 191-199.
#'
#' O'Brien P.C., Fleming T.R. (1979), A multiple testing procedure for clinical
#'  trials, Biometrics, 35(3), 549-556.
#'
#' Whitehead J. (1997), The Design and Analysis of Sequential Clinical Trials,
#' Wiley: Chichester, UK.
#' @export
#' @examples
#' \dontrun{
#' # Note that some of these examples may take a few minutes to run
#' # 2-stage design with O'Brien & Fleming efficacy and zero futility boundary
#' with
#' # equal sample size per arm and stage. Results are equivalent to using
#'  mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2, ushape="obf",
#'            lshape="fixed", lfix=0, sample.size=FALSE)
#' new.bounds(K=4, J=2, alpha=0.05, nMat=matrix(c(10, 20), nrow=2, ncol=5),
#' u=NULL, l=NULL,
#'            ushape="obf", lshape="fixed", lfix=0)
#' # A 2-stage design that was designed to use an O'Brien & Fleming efficacy
#' # and zero futility boundary with equal sample size per arm and stage (n=14).
#' # The observed sample size after stage one are 10, 10, 18, 10, 13 for each
#' # arm while the original upper bounds used are (3.068, 2.169) for stage 1.
#' # The updated bounds are (3.068, 2.167).
#' new.bounds(K=4, J=2, alpha=0.05,
#'      nMat=matrix(c(10, 28, 10, 28, 18, 28, 10, 28, 13, 28), nrow=2, ncol=5),
#'      u=3.068, l=0, ushape="obf", lshape="fixed", lfix=0)
#'
#' # same using parallelisation via separate R sessions running in the
#' # background
#' future::plan(multisession)
#' new.bounds(K=4, J=2, alpha=0.05,
#'            nMat=matrix(c(10, 28, 10, 28, 18, 28, 10, 28, 13, 28),
#'            nrow=2, ncol=5),
#'            u=3.068, l=0, ushape="obf", lshape="fixed", lfix=0)
#' future::plan("default")
#' }
new.bounds <- function(K = 3, J = 2, alpha = 0.05,
                           nMat = matrix(c(10, 20), nrow = 2, ncol = 4),
                           u = NULL, l = NULL, ushape = "obf", lshape = "fixed",
                           ufix = NULL, lfix = 0, N = 20, parallel = TRUE,
                                                                print = TRUE) {
  # require(mvtnorm) ## the function pmvnorm is required to evaluate
  # multivariate normal probabilities

  ##############################################################################
  ## 'mesh' creates the points and respective weights to use in the outer
  ## quadrature integrals
  ## you provide one-dimensional set of points (x) and weights (w)
  ## (e.g. midpoint rule, gaussian quadrature, etc)
  ## and the number of dimensions (d) and it gives d-dimensional collection of
  ## points and weights
  ##############################################################################

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

  R_prodsum1_nb <- function(x, l, u, R, r0, r0diff, J, K, Sigma) {
    ############################################################################
    ## x is vector of dummy variables ( the t_j in gdt paper ), l and u are
    ## boundary vectors
    ############################################################################
    .Call("C_prodsum1_nb",
      x2 = as.double(x), l2 = as.double(l), u2 = as.double(u),
      R2 = as.double(R), r02 = as.double(r0), r0diff2 = as.double(r0diff),
      Jfull2 = as.integer(J), K2 = as.integer(K), Sigma2 = Sigma,
      maxpts2 = as.integer(25000), releps2 = as.double(0),
      abseps2 = as.double(0.001), tol2 = as.double(0.0001)
    )
  }

  ##############################################################################
  ##  'typeI' performs the outer quadrature integral in type I error equation
  ##  using 'mesh' and 'prodsum'
  ##  and calculates the difference with the nominal alpha.
  ##  The accuracy of the quadrature integral will depend on the choice of
  ##  points and weights.
  ##  Here, the number of points in each dimension is an input, N.
  ##  The midpoint rule is used with a range of -6 to 6 in each dimension.
  ##############################################################################

  typeI <- function(C, alpha, N, R, r0, r0diff, J, K, Sigma, mmp,
                    u, l, ushape, lshape, lfix = NULL, ufix = NULL,
                    parallel = parallel, print = print) {
    ## number of stages already done
    if (!is.null(l)) {
      j <- length(l)
    } else {
      j <- 0
    }
    ########################################################################
    ## the form of the boundary constraints are determined as functions of C.
    ########################################################################
    Jleft <- J - j
    if (Jleft == 1) {
      ind <- J
    } else {
      ind <- (j + 1):J
    }
    if (!is.function(ushape)) {
      if (ushape == "obf") {
        ub <- c(u, C * sqrt(J / (1:J))[ind])
      } else if (ushape == "pocock") {
        ub <- c(u, rep(C, Jleft))
      } else if (ushape == "fixed") {
        ub <- c(u, rep(ufix, Jleft - 1), C)
      } else if (ushape == "triangular") {
        # ub<-c(u,C*(1+(1:J)/J)/sqrt(1:J)[(j+1):J])
        ub <- c(u, (C * (1 + (1:J) / J) / sqrt(1:J))[ind])
      }
    } else {
      ub <- c(u, C * ushape(J)[(j + 1):J])
    }

    if (!is.function(lshape)) {
      if (lshape == "obf") {
        lb <- c(l, -C * sqrt(J / (1:(J - 1)))[(j + 1):(J - 1)], ub[J])
      } else if (lshape == "pocock") {
        lb <- c(l, rep(-C, Jleft - 1), ub[J])
      } else if (lshape == "fixed") {
        lb <- c(l, rep(lfix, Jleft - 1), ub[J])
      } else if (lshape == "triangular") {
        # lb<-c(l,-C*(1-3*(1:J)/J)/sqrt(1:J)[(j+1):J])
        lb <- c(l, (-C * (1 - 3 * (1:J) / J) / sqrt(1:J))[ind])
      }
    } else {
      lb <- c(l, C * lshape(J)[(j + 1):(J - 1)], ub[J])
    }

    if (parallel) {
      evs <- future.apply::future_apply(mmp$X, 1, R_prodsum1_nb,
        l = lb, u = ub, R = R, r0 = r0, r0diff = r0diff, J = J,
        K = K, Sigma = Sigma, future.seed = TRUE,
        future.packages = "MAMS"
      )
    } else {
      evs <- apply(mmp$X, 1, R_prodsum1_nb,
        l = lb, u = ub, R = R, r0 = r0, r0diff = r0diff, J = J, K = K,
        Sigma = Sigma
      )
    }
    if (print) {
      message(".", appendLF = FALSE)
    }
    truealpha <- 1 - mmp$w %*% evs

    return(truealpha - alpha)
  }

  # check input

  if (K %% 1 > 0 | J %% 1 > 0) {
    stop("K and J need to be integers.")
  }
  if (K < 1 | J < 1) {
    stop("The number of stages and treatments must be at least 1.")
  }
  if (N <= 3) {
    stop("Number of points for integration by quadrature to small or negative.")
  }
  if (N > 3 & N <= 10) {
    warning("Number of points for integration by quadrature is small which may
      result in inaccurate solutions.")
  }
  if (alpha < 0 | alpha > 1) {
    stop("Error rate not between 0 and 1.")
  }

  if (!is.function(ushape)) {
    if (!ushape %in% c("pocock", "obf", "triangular", "fixed")) {
      stop("Upper boundary does not match the available options")
    }
    if (ushape == "fixed" & is.null(ufix)) {
      stop("ufix required when using a fixed upper boundary shape.")
    }
  } else {
    b <- ushape(J)
    if (!all(sort(b, decreasing = TRUE) == b)) {
      stop("Upper boundary shape is increasing")
    }
  }
  if (!is.function(lshape)) {
    if (!lshape %in% c("pocock", "obf", "triangular", "fixed")) {
      stop("Lower boundary does not match the available options")
    }
    if (lshape == "fixed" & is.null(lfix)) {
      stop("lfix required when using a fixed lower boundary shape.")
    }
  } else {
    b <- lshape(J)
    if (!all(sort(b, decreasing = FALSE) == b)) {
      stop("Lower boundary shape is decreasing")
    }
  }

  if (!all(diff(nMat) >= 0)) {
    stop("Total sample size per arm can not decrease between stages")
  }
  if (ncol(nMat) != K + 1) {
    stop("Number of columns in nMat not equal to K+1")
  }
  if (nrow(nMat) != J) {
    stop("Number of rows in nMat not equal to J")
  }
  if (length(l) != length(u)) {
    stop("Length of l must be the same as length of u")
  }
  if (length(l) > (J - 1)) {
    stop("Maximum number of stages, J,
    greater or equal to length of boundary vector")
  }

  r0 <- nMat[, 1] / nMat[1, 1]
  R <- nMat[, -1] / nMat[1, 1]

  ############################################################################
  ## gaussian quadrature's grid and weights for stages 1:J
  ############################################################################

  mmp_j <- as.list(rep(NA, J))
  for (j in 1:J) {
    mmp_j[[j]] <- mesh(x = (1:N - .5) / N * 12 - 6, j, w = rep(12 / N, N))
  }


  ####################################################################
  ## Create the variance covariance matrix from allocation proportions:
  ####################################################################

  bottom <- array(R, dim = c(J, K, J))
  top <- array(rep(t(R), rep(J, K * J)), dim = c(J, K, J))
  for (i in 1:K) {
    top[, i, ][upper.tri(top[, i, ])] <- t(top[, i, ])[upper.tri(top[, i, ])]
    bottom[, i, ][upper.tri(bottom[, i, ])] <- t(bottom[, i,
                                                  ])[upper.tri(bottom[, i, ])]
  }
  tmp <- sqrt(top / bottom)
  Sigma <- array(NA, dim = c(J, J, K))
  for (k in 1:K) {
    Sigma[, , k] <- tmp[, k, ]
  }

  ##############################################################################
  ## Create r0diff: the proportion of patients allocated to each
  ## particular stage
  ##############################################################################

  r0lag1 <- c(0, r0[1:J - 1])
  r0diff <- r0 - r0lag1


  ################################
  ## Find boundaries using 'typeI'
  ################################

  if (print) {
    message("   i) find new lower and upper boundaries\n      ",
            appendLF = FALSE)
  }

  uJ <- NULL
  try(uJ <- uniroot(typeI, c(0, 5),
    alpha = alpha, N = N, R = R, r0 = r0, r0diff = r0diff,
    J = J, K = K, Sigma = Sigma, mmp = mmp_j[[J]],
    u = u, l = l, ushape = ushape, lshape = lshape, lfix = lfix, ufix = ufix,
    parallel = parallel, print = print, tol = 0.001
  )$root, silent = TRUE)

  if (is.null(uJ)) {
    stop("No solution found")
  }


  ## number of stages already done
  if (!is.null(l)) {
    j <- length(l)
  } else {
    j <- 0
  }
  ## number of stages left
  Jleft <- J - j
  if (Jleft == 1) {
    ind <- J
  } else {
    ind <- (j + 1):J
  }

  ########################################################################
  ## the form of the boundary constraints are determined as functions of C.
  ########################################################################

  if (!is.function(ushape)) {
    if (ushape == "obf") {
      ub <- c(u, uJ * sqrt(J / (1:J))[ind])
    } else if (ushape == "pocock") {
      ub <- c(u, rep(uJ, Jleft))
    } else if (ushape == "fixed") {
      ub <- c(u, rep(ufix, Jleft), uJ)
    } else if (ushape == "triangular") {
      ub <- c(u, (uJ * (1 + (1:J) / J) / sqrt(1:J))[ind])
    }
  } else {
    ub <- c(u, uJ * ushape(J)[(j + 1):J])
  }

  if (!is.function(lshape)) {
    if (lshape == "obf") {
      lb <- c(l, -uJ * sqrt(J / (1:(J - 1)))[(j + 1):(J - 1)], ub[J])
    } else if (lshape == "pocock") {
      lb <- c(l, rep(-uJ, Jleft - 1), ub[J])
    } else if (lshape == "fixed") {
      lb <- c(l, rep(lfix, Jleft - 1), ub[J])
    } else if (lshape == "triangular") {
      lb <- c(l, (-uJ * (1 - 3 * (1:J) / J) / sqrt(1:J))[ind])
    }
  } else {
    lb <- c(l, uJ * lshape(J)[(j + 1):(J - 1)], ub[J])
  }

  #########################################################
  ## Find alpha.star
  #########################################################
  if (print) {
    message("\n  ii) define alpha star\n", appendLF = FALSE)
  }

  alpha.star <- numeric(J)
  alpha.star[1] <- typeI(ub[1],
    alpha = 0, N = N, R = t(as.matrix(R[1, ])),
    r0 = r0[1], r0diff = r0diff[1],
    J = 1, K = K, Sigma = Sigma, mmp = mmp_j[[1]],
    u = NULL, l = NULL, ushape = "fixed",
    lshape = "fixed", lfix = NULL, ufix = NULL,
    parallel = parallel, print = FALSE
  )
  if (J > 1) {
    for (j in 2:J) {
      alpha.star[j] <- typeI(ub[j],
        alpha = 0, N = N, R = R[1:j, ],
        r0 = r0[1:j], r0diff = r0diff[1:j],
        J = j, K = K, Sigma = Sigma, mmp = mmp_j[[j]],
        u = NULL, l = NULL, ushape = "fixed",
        lshape = "fixed", lfix = lb[1:(j - 1)],
        ufix = ub[1:(j - 1)],
        parallel = parallel, print = FALSE
      )
    }
  }

  res <- NULL
  res$l <- lb
  res$u <- ub
  res$n <- nMat[1, 1]
  res$rMat <- rbind(r0, t(R)) ## allocation ratios
  res$N <- sum(res$rMat[, J] * res$n) ## maximum total sample size

  res$K <- K
  res$J <- J
  res$alpha <- alpha
  res$alpha.star <- alpha.star
  res$power <- NA
  res$type <- "new.bounds"

  class(res) <- "MAMS"
  attr(res, "method") <- "simultaneous"

  return(res)
}

###############################################################################
###################### ordinal.mams function ##################################
###############################################################################
#' Function to design multi-arm multi-stage studies with ordinal or binary
#'  endpoints
#' @description The function determines (approximately) the boundaries of a
#' multi-arm multi-stage study with ordinal or binary endpoints for a given
#' boundary shape and finds the required number of subjects.
#' @param prob Vector of expected probabilities of falling into each category
#' under control conditions. The elements must sum up to one
#' (default=c(0.35, 0.4, 0.25)).
#' @param or Interesting treatment effect on the scale of odds ratios
#' (default=2).
#' @param or0 Uninteresting treatment effect on the scale of odds ratios
#' (default=1.2).
#' @param K Number of experimental treatments (default=4).
#' @param J Number of stages (default=2).
#' @param alpha One-sided familywise error rate (default=0.05).
#' @param power Desired power (default=0.9).
#' @param r Vector of allocation ratios (default=1:2).
#' @param r0 	Vector ratio on control (default=1:2).
#' @param ushape 	Shape of upper boundary. Either a function specifying the
#' shape or one of "pocock", "obf" (the default), "triangular" and "fixed".
#' @param lshape 	Shape of lower boundary. Either a function specifying the
#'  shape or one of "pocock", "obf", "triangular" and "fixed" (the default).
#' @param ufix Fixed upper boundary (default=NULL). Only used if shape="fixed".
#' @param lfix Fixed lower boundary (default=0). Only used if shape="fixed".
#' @param nstart Starting point for finding the sample size (default=1).
#' @param nstop Stopping point for finding the sample size (default=NULL).
#' @param sample.size Logical if sample size should be found as well
#' (default=TRUE).
#' @param Q Number of quadrature points per dimension in the outer integral
#' (default=20).
#' @param parallel if TRUE (default), allows parallelisation of the computation
#' via a user-defined strategy specified by means of the function
#' future::plan(). If not set differently, the default strategy is sequential,
#'  which corresponds to a computation without parallelisation.
#' @param print if TRUE (default), indicate at which stage the computation is.
#' @details This function finds the (approximate) boundaries and sample size of
#'  a multi-arm multi-stage study with ordinal or binary endpoints with K active
#'  treatments plus control in which all promising treatments are continued at
#'  interim analyses as described in Magirr et al (2012). It is a wrapper around
#'  the basic mams function to facilitate its use with ordinal and binary
#' endpoints, following ideas of Whitehead & Jaki (2009) and Jaki & Magirr
#' (2013). For a binary endpoint the vector prob has only two elements
#' (success/failure, yes/no, etc.). See mams for further details on the basic
#' methodology.
#' @return An object of the class MAMS containing the following components: \cr
#' \item{prob}{Vector of expected probabilities of falling into each category
#' under control conditions. The elements must sum up to one
#' (default=\code{c(0.35, 0.4, 0.25)}).}
#'  \item{or}{Interesting treatment effect on the scale of odds ratios
#' (default=\code{2}).}
#'  \item{or0}{Uninteresting treatment effect on the scale of odds ratios
#' (default=\code{1.2}).}
#'  \item{K}{Number of experimental treatments (default=\code{4}).}
#'  \item{J}{Number of stages (default=\code{2}).}
#'  \item{alpha}{One-sided familywise error rate (default=\code{0.05}).}
#'  \item{power}{Desired power (default=\code{0.9}).}
#'  \item{r}{Vector of allocation ratios (default=\code{1:2}).}
#'  \item{r0}{Vector ratio on control (default=\code{1:2}).}
#'  \item{ushape}{Shape of upper boundary. Either a function specifying the
#' shape or one of \code{"pocock"}, \code{"obf"} (the default),
#' \code{"triangular"} and \code{"fixed"}.}
#'  \item{lshape}{Shape of lower boundary. Either a function specifying the
#' shape or one of \code{"pocock"}, \code{"obf"}, \code{"triangular"} and
#' \code{"fixed"} (the default).}
#'  \item{ufix}{Fixed upper boundary (default=\code{NULL}). Only used if
#' \code{shape="fixed"}.}
#'  \item{lfix}{Fixed lower boundary (default=\code{0}). Only used if
#' \code{shape="fixed"}.}
#'  \item{nstart}{Starting point for finding the sample size
#' (default=\code{1}).}
#'  \item{nstop}{Stopping point for finding the sample size
#' (default=\code{NULL}).}
#'  \item{sample.size}{Logical if sample size should be found as well
#' (default=\code{TRUE}).}
#'  \item{N}{Number of quadrature points per dimension in the outer integral
#' (default=\code{20}).}
#'  \item{parallel}{if \code{TRUE} (default), allows parallelisation of the
#' computation via a user-defined strategy specified by means of the function
#' \code{\link[future:plan]{future::plan()}}. If not set differently,
#' the default strategy is \code{sequential}, which corresponds to a
#' computation without parallelisation.}
#'  \item{print}{if \code{TRUE} (default), indicate at which stage the
#' computation is.}
#' @author Philip Pallmann
#' @references
#' Jaki T., Pallmann P. and Magirr D. (2019), The R Package MAMS for Designing
#' Multi-Arm Multi-Stage Clinical Trials, Journal of Statistical Software,
#' 88(4), 1-25. Link: doi:10.18637/jss.v088.i04
#'
#' Magirr D., Jaki T. and Whitehead J. (2012), A generalized Dunnett test for
#' multi-arm multi-stage clinical studies with treatment selection, Biometrika,
#' 99(2), 494-501. Link: doi:10.1093/biomet/ass002
#'
#' Magirr D., Stallard N. and Jaki T. (2014), Flexible sequential designs for
#' multi-arm clinical trials, Statistics in Medicine, 33(19), 3269-3279. Link:
#' doi:10.1002/sim.6183
#'
#' Pocock S.J. (1977), Group sequential methods in the design and analysis of
#' clinical trials, Biometrika, 64(2), 191-199.
#'
#' O'Brien P.C., Fleming T.R. (1979), A multiple testing procedure for clinical
#' trials, Biometrics, 35(3), 549-556.
#'
#' Whitehead J. (1997), The Design and Analysis of Sequential Clinical Trials,
#' Wiley: Chichester, UK.
#' @export
#' @examples
#' \dontrun{
#' ## An example based on the example in Whitehead & Jaki (2009)
#' # 2-stage design with triangular efficacy and futility boundaries
#' prob <- c(0.075, 0.182, 0.319, 0.243, 0.015, 0.166)
#' ordinal.mams(prob=prob, or=3.06, or0=1.32, K=3, J=2, alpha=0.05,
#'                  power=0.9, r=1:2, r0=1:2, ushape="triangular",
#'                  lshape="triangular")
#' # same example with parallelisation via separate R sessions running in the
#' # background
#' future::plan(multisession)
#' ordinal.mams(prob=prob, or=3.06, or0=1.32, K=3, J=2, alpha=0.05,
#'                  power=0.9, r=1:2, r0=1:2, ushape="triangular",
#'                  lshape="triangular", parallel=TRUE)
#' future::plan("default")
#' }
ordinal.mams <- function(prob = c(0.35, 0.4, 0.25), or = 2, or0 = 1.2,
                              K = 4, J = 2, alpha = 0.05, power = 0.9, r = 1:2,
                              r0 = 1:2, ushape = "obf", lshape = "fixed",
                              ufix = NULL, lfix = 0, nstart = 1, nstop = NULL,
                              sample.size = TRUE, Q = 20, parallel = TRUE,
                              print = TRUE) {
  if (sum(prob) != 1) {
    stop("The elements of prob must sum to one.")
  }
  if (or0 < 1) {
    stop("The uninteresting effect must be 1 or larger.")
  }
  if (or <= or0) {
    stop("The interesting effect must be larger than the uninteresting one.")
  }

  q <- (1 - sum(prob^3)) / 3
  sigma <- 1 / sqrt(q)

  p <- pnorm(log(or) / sqrt(2 * sigma^2))
  p0 <- pnorm(log(or0) / sqrt(2 * sigma^2))

  mams(
    K = K, J = J, alpha = alpha, power = power, r = r, r0 = r0, p = p, p0 = p0,
    delta = NULL, delta0 = NULL, sd = NULL, ushape = ushape, lshape = lshape,
    ufix = ufix, lfix = lfix, nstart = nstart, nstop = nstop,
    sample.size = sample.size,
    Q = Q, type = "ordinal", parallel = parallel, print = print
  )
}

###############################################################################
###################### stepdown.mams function #################################
###############################################################################
#' Function to find stopping boundaries for a 2- or 3-stage (step-down)
#' multiple-comparisons-with-control test.
#' @description The function determines stopping boundaries for all intersection
#'  hypothesis tests in a multi-arm multi-stage study, given the amount of alpha
#'  (familywise error rate) to be spent at each analysis.
#' @usage stepdown.mams(nMat=matrix(c(10, 20), nrow=2, ncol=4),
#'              alpha.star=c(0.01, 0.025), lb=0,
#'              selection="all.promising")
#' @param nMat 	Matrix containing the cumulative sample sizes in each treatment
#' arm columns: control, trt 1, ..., trt K), at each analysis (rows).
#' The number of analyses must be either 2 or 3 (default=matrix(c(10, 20),
#' nrow=2, ncol=4)).
#' @param alpha.star 	Cumulative familywise error rate to be spent at each
#' analysis (default=c(0.01, 0.025)).
#' @param lb Fixed lower boundary (default=0).
#' @param selection How are treatments selected for the next stage? Using the
#' default "all.promising" method, all treatments with a test statistic
#' exceeding the lower boundary are taken forward to the next stage.
#' If "select.best", only the treatment with the largest statistic may be
#' selected for future stages. (default="all.promising").
#' @details The function implements the methods described in Magirr et al (2014)
#'  to find individual boundaries for all intersection hypotheses.
#' @return An object of the class MAMS.stepdown containing the following
#' components: \cr
#' \item{l}{Lower boundaries.}
#'  \item{u}{Upper boundaries.}
#'  \item{nMat}{Cumulative sample sizes on each treatment arm.}
#'  \item{K}{Number of experimental treatments.}
#'  \item{J}{Number of stages in the trial.}
#'  \item{alpha.star}{Cumulative familywise error rate spent at each analysis.}
#'  \item{selection}{Pre-specified method of treatment selection.}
#'  \item{zscores}{A list containing the observed test statistics at analyses
#' so far (at the design stage this is NULL).}
#'  \item{selected.trts}{A list containing the treatments selected for
#' each stage.}
#' @author Dominic Magirr
#' @references
#' Jaki T., Pallmann P. and Magirr D. (2019), The R Package MAMS for Designing
#' Multi-Arm Multi-Stage Clinical Trials, Journal of Statistical Software,
#' 88(4), 1-25. Link: doi:10.18637/jss.v088.i04
#'
#' Magirr D., Jaki T. and Whitehead J. (2012), A generalized Dunnett test for
#' multi-arm multi-stage clinical studies with treatment selection, Biometrika,
#'  99(2), 494-501. Link: doi:10.1093/biomet/ass002
#'
#' Magirr D., Stallard N. and Jaki T. (2014), Flexible sequential designs for
#' multi-arm clinical trials, Statistics in Medicine, 33(19), 3269-3279.
#' Link: doi:10.1002/sim.6183
#'
#' Stallard N. and Todd S. (2003), Sequential designs for phase III clinical
#' trials incorporating treatment selection, Statistics in Medicine, 22(5),
#' 689-703.
#' @export
#' @examples
#' \dontrun{
#' # Note that some of these examples may take a few minutes to run
#' # 2-stage 3-treatments versus control design, all promising treatments
#' # are selected:
#' stepdown.mams(nMat=matrix(c(10, 20), nrow=2, ncol=4),
#'               alpha.star=c(0.01, 0.05), lb=0,
#'               selection="all.promising")
#' # select the best treatment after the first stage:
#' stepdown.mams(nMat=matrix(c(10, 20), nrow=2, ncol=4),
#'               alpha.star=c(0.01, 0.05), lb=0,
#'               selection="select.best")
#' # 3 stages and unequal randomization:
#' stepdown.mams(nMat=matrix(c(20, 40, 60, rep(c(10, 20, 30), 3)),
#'               nrow=3, ncol=4),
#'               alpha.star=c(0.01, 0.025, 0.05), lb=c(0, 0.75),
#'               selection="all.promising")
#' }
stepdown.mams <- function(nMat = matrix(c(10, 20), nrow = 2, ncol = 4),
                              alpha.star = c(0.01, 0.025), lb = 0,
                              selection = "all.promising") {
  # checking input parameters
  if (!all(diff(nMat) >= 0)) {
    stop("total sample size per arm cannot decrease between stages.")
  }
  J <- dim(nMat)[1]
  K <- dim(nMat)[2] - 1
  if ((J != 2) && (J != 3)) {
    stop("number of stages must be 2 or 3")
  }
  if (K < 2) {
    stop("must have at least two experimental treatments")
  }
  if (length(alpha.star) != J) {
    stop("length of error spending vector must be same as number of stages")
  }
  if (!all(diff(alpha.star) >= 0)) {
    stop("cumulative familywise error must increase.")
  }
  if (length(lb) != J - 1) {
 stop("lower boundary must be specified at all analysis points except the last")
  }
  match.arg(selection, c("all.promising", "select.best"))

  # find the nth intersection hypothesis (positions of 1s in binary n)
  get.hyp <- function(n) {
    indlength <- ceiling(log(n) / log(2) + .0000001)
    ind <- rep(0, indlength)
    newn <- n

    for (h in seq(1, indlength)) {
      ind[h] <- (newn / (2^(h - 1))) %% 2
      newn <- newn - ind[h] * 2^(h - 1)
    }
    seq(1, indlength)[ind == 1]
  }
  # for argument c(i,j) this gives covariance between statistics in stage i
  # with statistics in stage j
  create.block <- function(control.ratios = 1:2,
                            active.ratios = matrix(1:2, 2, 3)) {
    K <- dim(active.ratios)[2]
    block <- matrix(NA, K, K)
    for (i in 1:K) {
      block[i, i] <- sqrt(active.ratios[1, i] * control.ratios[1] *
        (active.ratios[2, i] + control.ratios[2]) / (active.ratios[1, i] +
          control.ratios[1]) / active.ratios[2, i] / control.ratios[2])
    }
    for (i in 2:K) {
      for (j in 1:(i - 1)) {
        block[i, j] <- sqrt(active.ratios[1, i] * control.ratios[1] *
          active.ratios[2, j] / (active.ratios[1, i] + control.ratios[1]) /
          (active.ratios[2, j] + control.ratios[2]) / control.ratios[2])
        block[j, i] <- sqrt(active.ratios[1, j] * control.ratios[1] *
          active.ratios[2, i] / (active.ratios[1, j] + control.ratios[1]) /
          (active.ratios[2, i] + control.ratios[2]) / control.ratios[2])
      }
    }
    block
  }

  # create variance-covariance matrix of the test statistics
  create.cov.matrix <- function(control.ratios = 1:2,
                                active.ratios = matrix(1:2, 2, 3)) {
    J <- dim(active.ratios)[1]
    K <- dim(active.ratios)[2]

    cov.matrix <- matrix(NA, J * K, J * K)
    for (i in 1:J) {
      for (j in i:J) {
        cov.matrix[((i - 1) * K + 1):(i * K), ((j - 1) * K + 1):(j * K)] <-
          create.block(control.ratios[c(i, j)], active.ratios[c(i, j), ])

        cov.matrix[((j - 1) * K + 1):(j * K), ((i - 1) * K + 1):(i * K)] <-
          t(cov.matrix[((i - 1) * K + 1):(i * K), ((j - 1) * K + 1):(j * K)])
      }
    }
    cov.matrix
  }
  # find the probability that no test statistic crosses the upper boundary +
  # only treatments in surviving_subsetj reach the jth stage
  get.path.prob <- function(
      surviving.subset1, surviving.subset2 = NULL,
      cut.off, treatments, cov.matrix, lb, upper.boundary, K, stage) {
    treatments2 <- treatments[surviving.subset1]
    if (stage == 2) {
      lower <- c(rep(-Inf, length(treatments)), rep(-Inf, length(treatments2)))
      lower[surviving.subset1] <- lb[1]

      upper <- c(
        rep(lb[1], length(treatments)),
        rep(cut.off, length(treatments2))
      )
      upper[surviving.subset1] <- upper.boundary[1]

      return(pmvnorm(
        lower = lower, upper = upper,
        sigma = cov.matrix[
          c(treatments, K + treatments2),
          c(treatments, K + treatments2)
        ]
      )[1])
    }
    treatments3 <- treatments2[surviving.subset2]

    lower <- c(
      rep(-Inf, length(treatments)), rep(-Inf, length(treatments2)),
      rep(-Inf, length(treatments3))
    )
    lower[surviving.subset1] <- lb[1]
    lower[length(treatments) + surviving.subset2] <- lb[2]

    upper <- c(
      rep(lb[1], length(treatments)), rep(lb[2], length(treatments2)),
      rep(cut.off, length(treatments3))
    )
    upper[surviving.subset1] <- upper.boundary[1]
    upper[length(treatments) + surviving.subset2] <- upper.boundary[2]

    pmvnorm(lower = lower, upper = upper, sigma = cov.matrix[c(
      treatments,
      K + treatments2, 2 * K + treatments3
    ), c(
      treatments, K + treatments2,
      2 * K + treatments3
    )])[1]
  }

  # for the "select.best" method, find the probability that "select.treatment"
  # is selected and subsequently crosses the upper boundary
  rejection.paths <- function(
      selected.treatment, cut.off, treatments,
      cov.matrix, lb, upper.boundary, K, stage) {
    contrast <- diag(-1, K + stage - 1)
    contrast[1:K, selected.treatment] <- 1
    for (i in 1:(stage - 1)) contrast[K + i, K + i] <- 1

    bar.cov.matrix <- contrast %*% cov.matrix[c(1:K, 1:(stage - 1) *
      K + selected.treatment), c(1:K, 1:(stage - 1) * K +
                                                        selected.treatment)] %*%
      t(contrast)

    lower <- c(rep(0, length(treatments)), cut.off)
    if (stage > 2) {
      lower <- c(
        rep(0, length(treatments)), lb[2:(stage - 1)],
        cut.off
      )
    }
    lower[which(treatments == selected.treatment)] <- lb[1]

    upper <- c(rep(Inf, length(treatments)), Inf)
    if (stage > 2) {
      upper <- c(
        rep(Inf, length(treatments)),
        upper.boundary[2:(stage - 1)], Inf
      )
    }
    upper[which(treatments == selected.treatment)] <- upper.boundary[1]

    pmvnorm(lower = lower, upper = upper, sigma = bar.cov.matrix[c(
      treatments,
      K + 1:(stage - 1)
    ), c(treatments, K + 1:(stage - 1))])[1]
  }

  # for "all.promising" rule, this gives the cumulative typeI error for
  # 'stage' stages
  # for "select.best" rule, this gives the Type I error spent at the
  # 'stage'th stage
  excess.alpha <- function(
      cut.off, alpha.star, treatments, cov.matrix, lb,
      upper.boundary, selection, K, stage) {
    if (stage == 1) {
      return(1 - alpha.star[1] - pmvnorm(
        lower =
          rep(-Inf, length(treatments)), upper = rep(
          cut.off,
          length(treatments)
        ), sigma = cov.matrix[treatments, treatments]
      )[1])
    }
    if (selection == "select.best") {
      return(alpha.star[stage] -
        alpha.star[stage - 1] - sum(unlist(lapply(treatments, rejection.paths,
          cut.off = cut.off, treatments = treatments, cov.matrix = cov.matrix,
          lb = lb, upper.boundary = upper.boundary, K = K, stage = stage
        ))))
    }
    # any of 'treatments' could be selected, so we add all these probabilities
    if (stage == 2) {
      # list all possible subsets of surviving treatments after the first stage
      surviving.subsets <- c(list(numeric(0)), lapply(as.list(1:(2^
        length(treatments) - 1)), get.hyp))
      return(1 - alpha.star[2] - sum(unlist(lapply(surviving.subsets,
        get.path.prob,
        cut.off = cut.off, treatments = treatments,
        cov.matrix = cov.matrix, lb = lb, upper.boundary = upper.boundary,
        K = K, stage = stage
      ))))
    }
    # all possible subsets of surviving treatments after the first stage
    surviving.subsets1 <- c(list(numeric(0)), lapply(as.list(1:(2^
      length(treatments) - 1)), get.hyp))
    surviving.subsets2 <- c(
      list(list(numeric(0))),
      lapply(surviving.subsets1[-1], function(x) {
        c(
          list(numeric(0)),
          lapply(as.list(1:(2^length(x) - 1)), get.hyp)
        )
      })
    )

# for each possible subset of survivng subsets after stage 1, list the possible
# subsets still surviving after stage 2
    1 - alpha.star[3] - sum(unlist(Map(function(x, y) {
      sum(unlist(lapply(y, get.path.prob,
        surviving.subset1 = x,
        cut.off = cut.off,
        treatments = treatments,
        cov.matrix = cov.matrix,
        lb = lb,
        upper.boundary = upper.boundary,
        K = K,
        stage = stage
      )))
    }, surviving.subsets1, surviving.subsets2)))
  }

  # get sample size ratios
  R <- nMat[, -1] / nMat[1, 1]
  r0 <- nMat[, 1] / nMat[1, 1]

  cov.matrix <- create.cov.matrix(r0, R)

  l <- u <- as.list(1:(2^K - 1))

  alpha.star <- rep(list(alpha.star), 2^K - 1)

  for (i in 1:(2^K - 1)) {
    names(u)[i] <- paste("U_{", paste(get.hyp(i), collapse = " "), "}",
                                                                sep = "")
    names(l)[i] <- paste("L_{", paste(get.hyp(i), collapse = " "), "}",
                                                                sep = "")
    names(alpha.star)[i] <- paste("alpha.star.{",
      paste(get.hyp(i), collapse = " "), "}",
      sep = ""
    )

    for (j in 1:J) {
      try(
        new.u <- uniroot(excess.alpha, c(0, 10),
          alpha.star = alpha.star[[i]],
          treatments = get.hyp(i), cov.matrix = cov.matrix, lb = lb,
          upper.boundary = u[[i]], selection = selection, K = K, stage = j
        )$root,
        silent = TRUE
      )
      if (is.null(new.u)) {
        stop("upper boundary not between 0 and 10")
      }
      u[[i]][j] <- round(new.u, 2)
    }

    l[[i]] <- c(lb, u[[i]][J])
  }

  res <- NULL
  res$l <- l
  res$u <- u
  res$sample.sizes <- nMat
  res$K <- K
  res$J <- J
  res$alpha.star <- alpha.star
  res$selection <- selection
  res$zscores <- NULL
  res$selected.trts <- list(1:K)

  class(res) <- "MAMS.stepdown"

  return(res)
}


###############################################################################
###################### stepdown.update function ###############################
###############################################################################
#' Update the stopping boundaries of multi-arm multi-stage study at an interim
#' analysis, allowing for unplanned treatment selection and/or sample-size
#' reassessment.
#' @description Function to update a planned multi-arm multi-stage design to
#' account for unplanned adaptations.
#' @param current.mams The planned step-down MAMS design prior to the current
#' interim analysis (=defaultstepdown.mams()).
#' @param nobs Cumulative sample sizes observed on each treatment arm up to and
#' including the current interim analysis.
#' @param zscores Observed vector of test statistics at the current interim
#' analysis.
#' @param selected.trts The set of experimental treatments to be taken forward
#' to the next stage of testing. This argument should be omitted at the final
#' analysis.
#' @param nfuture A matrix of future cumulative sample sizes. The number of rows
#'  must be equal to the originally planned number of stages (2 or 3) minus the
#' number of stages already observed. The number of columns must be equal to
#' the number of treatment arms (default=NULL).
#' @details The function implements the ideas described in Magirr et al. (2014)
#' to update a design according to unplanned design modifications. It takes as
#' input the planned multi-arm multi-stage design prior to the interim analysis,
#'  together with the actually observed cumulative sample sizes and test
#' statistics. Treatments to be included in future stages, as well as future
#' sample sizes, can be chosen without following pre-specified rules. The output
#'  is a new multi-arm multi-stage design for the remaining stages such that the
#'  familywise error remains controlled at the pre-specified level.
#' @author Dominic Magirr
#' @references
#' Jaki T., Pallmann P. and Magirr D. (2019), The R Package MAMS for Designing
#' Multi-Arm Multi-Stage Clinical Trials, Journal of Statistical Software,
#' 88(4), 1-25. Link: doi:10.18637/jss.v088.i04
#'
#' Magirr D., Jaki T. and Whitehead J. (2012), A generalized Dunnett test for
#' multi-arm multi-stage clinical studies with treatment selection, Biometrika,
#' 99(2), 494-501. Link: doi:10.1093/biomet/ass002
#'
#' Magirr D., Stallard N. and Jaki T. (2014), Flexible sequential designs for
#' multi-arm clinical trials, Statistics in Medicine, 33(19), 3269-3279.
#' Link: doi:10.1002/sim.6183
#'
#' Stallard N. and Todd S. (2003), Sequential designs for phase III clinical
#' trials incorporating treatment selection, Statistics in Medicine,
#' 22(5), 689-703.
#' @export
#' @examples
#' \dontrun{
#' # 2-stage 3-treatments versus control design
#' # all promising treatments are selected:
#' orig_mams <- stepdown.mams(nMat=matrix(c(10, 20), nrow=2, ncol=4),
#'                            alpha.star=c(0.01, 0.05), lb=0,
#'                            selection="all.promising")
#'
#' # make adjustment for the observed sample sizes
#' # not being exactly as planned:
#' stepdown.update(orig_mams, nobs=c(9, 8, 13, 11),
#'                  zscores=c(1.1, -0.5, 0.2),
#'                  selected.trts=1:3, nfuture=NULL)
#'
#' # make adjustment for the observed sample sizes
#' # not being exactly as planned. In addition, drop treatment 2:
#' stepdown.update(orig_mams, nobs=c(9, 8, 13, 11),
#'                  zscores=c(1.1, -0.5, 0.2),
#'                  selected.trts=c(1, 3), nfuture=NULL)
#'
#' # make adjustment for the observed sample sizes not being
#' # exactly as planned. In addition, drop treatment 2. In addition,
#' # double the planed cumulative second stage sample sizes:
#' updated_mams <- stepdown.update(orig_mams, nobs=c(9, 8, 13, 11),
#'                                  zscores=c(1.1, -0.5, 0.2),
#'                                  selected.trts=c(1, 3),
#'                                  nfuture=matrix(c(40, 40, 13, 40),
#'                                  nrow=1, ncol=4))
#'
#' # Account for the observed second stage sample sizes:
#' stepdown.update(updated_mams, nobs=c(38, 41, 13, 36),
#'                 zscores=c(1.9, -Inf, 1.2),
#'                 selected.trts=NULL)
#'
#' # 'select.best' design. Account for actually observed sample sizes
#' # in first stage, and drop treatment 2:
#' orig_mams <- stepdown.mams(nMat=matrix(c(10, 20), nrow=2, ncol=4),
#'                            alpha.star=c(0.01, 0.05), lb=0,
#'                            selection="select.best")
#'
#' stepdown.update(orig_mams, nobs=c(9, 8, 13, 11),
#'                  zscores=c(1.1, -0.5, 0.2),
#'                  selected.trts=c(1, 3), nfuture=NULL)
#' }
stepdown.update <- function(current.mams = stepdown.mams(), nobs = NULL,
                                zscores = NULL, selected.trts = NULL,
                                                nfuture = NULL) {
  zscores <- c(current.mams$zscores, list(zscores))

  if (!is.null(selected.trts)) {
    selected.trts <- c(
      current.mams$selected.trts,
      list(selected.trts)
    )
  }

  # checking input parameters
  if (!is(current.mams, "MAMS.stepdown")) {
    stop("current.mams must be a 'MAMS.stepdown' object")
  }
  if (length(nobs) != current.mams$K + 1) {
    stop("must provide observed cumulative sample size for each treatment")
  }

  completed.stages <- length(zscores)
  for (i in 1:completed.stages) {
    if (length(zscores[[i]]) != current.mams$K) {
      stop("vector of statistics is wrong length")
    }
  }

  if (is.null(selected.trts)) {
    if (current.mams$J > completed.stages) {
      stop("must specify treatments selected for next stage")
    }
  }

  for (i in seq_along(selected.trts)) {
    if (length(setdiff(selected.trts[[i]], 1:current.mams$K) > 0)) {
      stop("inappropriate treatment selection")
    }
  }
  if (is.matrix(nfuture)) {
    if (dim(nfuture)[1] != current.mams$J - completed.stages) {
      stop("must provide future sample sizes for all remaining stages")
    }
    if (dim(nfuture)[2] != current.mams$K + 1) {
      stop("must provide future sample sizes for all treatment arms")
    }
  }

  # load all necessary functions
  # find the nth intersection hypothesis (positions of 1s in binary n)
  get.hyp <- function(n) {
    indlength <- ceiling(log(n) / log(2) + .0000001)
    ind <- rep(0, indlength)
    newn <- n

    for (h in seq(1, indlength)) {
      ind[h] <- (newn / (2^(h - 1))) %% 2
      newn <- newn - ind[h] * 2^(h - 1)
    }
    seq(1, indlength)[ind == 1]
  }

  create.block <- function(control.ratios = 1:2,
                            active.ratios = matrix(1:2, 2, 3)) {
    # for argument c(i,j) this gives covariance between statistics in stage i
    # with statistics in stage j
    K <- dim(active.ratios)[2]
    block <- matrix(NA, K, K)
    for (i in 1:K) {
      block[i, i] <- sqrt(active.ratios[1, i] * control.ratios[1] *
        (active.ratios[2, i] + control.ratios[2]) / (active.ratios[1, i] +
          control.ratios[1]) / active.ratios[2, i] / control.ratios[2])
    }
    for (i in 2:K) {
      for (j in 1:(i - 1)) {
        block[i, j] <- sqrt(active.ratios[1, i] * control.ratios[1] *
          active.ratios[2, j] / (active.ratios[1, i] + control.ratios[1]) /
          (active.ratios[2, j] + control.ratios[2]) / control.ratios[2])
        block[j, i] <- sqrt(active.ratios[1, j] * control.ratios[1] *
          active.ratios[2, i] / (active.ratios[1, j] + control.ratios[1]) /
          (active.ratios[2, i] + control.ratios[2]) / control.ratios[2])
      }
    }
    block
  }

  # create variance-covariance matrix of the test statistics
  create.cov.matrix <- function(control.ratios = 1:2,
                                active.ratios = matrix(1:2, 2, 3)) {
    J <- dim(active.ratios)[1]
    K <- dim(active.ratios)[2]

    cov.matrix <- matrix(NA, J * K, J * K)
    for (i in 1:J) {
      for (j in i:J) {
        cov.matrix[((i - 1) * K + 1):(i * K), ((j - 1) * K + 1):(j * K)] <-
          create.block(control.ratios[c(i, j)], active.ratios[c(i, j), ])
        cov.matrix[((j - 1) * K + 1):(j * K), ((i - 1) * K + 1):(i * K)] <-
          t(cov.matrix[((i - 1) * K + 1):(i * K), ((j - 1) * K + 1):(j * K)])
      }
    }
    cov.matrix
  }
  create.cond.cov.matrix <- function(cov.matrix, K, completed.stages) {
  # find the conditional covariance of future test statistics given data so far

    sigma_1_1 <- cov.matrix[((completed.stages - 1) * K + 1):(completed.stages *
      K), ((completed.stages - 1) * K + 1):(completed.stages * K)]
    sigma_1_2 <- cov.matrix[((completed.stages - 1) * K + 1):(completed.stages *
      K), -(1:(completed.stages * K))]
    sigma_2_1 <- t(sigma_1_2)
    sigma_2_2 <- cov.matrix[
      -(1:(completed.stages * K)),
      -(1:(completed.stages * K))
    ]
    sigma_2_2 - sigma_2_1 %*% solve(sigma_1_1) %*% sigma_1_2
  }

  # find the conditional mean of future test statistics given data so far
  create.cond.mean <- function(cov.matrix, K, completed.stages, zscores) {
    sigma_1_1 <- cov.matrix[((completed.stages - 1) * K + 1):(completed.stages *
      K), ((completed.stages - 1) * K + 1):(completed.stages * K)]
    sigma_1_2 <- cov.matrix[((completed.stages - 1) * K + 1):(completed.stages *
      K), -(1:(completed.stages * K))]
    sigma_2_1 <- t(sigma_1_2)
    sigma_2_1 %*% solve(sigma_1_1) %*% zscores
  }

  # find the probability that no test statistic crosses the upper boundary +
  # only treatments in surviving_subsetj reach the jth stage
  get.path.prob <- function(
      surviving.subset1, surviving.subset2 = NULL,
      cut.off, treatments, cov.matrix, lower.boundary, upper.boundary, K, stage,
      z.means) {
    treatments2 <- treatments[surviving.subset1]
    if (stage == 2) {
      lower <- c(rep(-Inf, length(treatments)), rep(-Inf, length(treatments2)))
      lower[surviving.subset1] <- lower.boundary[1]

      upper <- c(rep(lower.boundary[1], length(treatments)), rep(
        cut.off,
        length(treatments2)
      ))
      upper[surviving.subset1] <- upper.boundary[1]

      return(pmvnorm(lower = lower, upper = upper, mean = z.means[c(
        treatments,
        K + treatments2
      )], sigma = cov.matrix[
        c(treatments, K + treatments2),
        c(treatments, K + treatments2)
      ])[1])
    }
    treatments3 <- treatments2[surviving.subset2]

    lower <- c(
      rep(-Inf, length(treatments)), rep(-Inf, length(treatments2)),
      rep(-Inf, length(treatments3))
    )
    lower[surviving.subset1] <- lower.boundary[1]
    lower[length(treatments) + surviving.subset2] <- lower.boundary[2]

    upper <- c(
      rep(lower.boundary[1], length(treatments)),
      rep(lower.boundary[2], length(treatments2)),
      rep(cut.off, length(treatments3))
    )
    upper[surviving.subset1] <- upper.boundary[1]
    upper[length(treatments) + surviving.subset2] <- upper.boundary[2]

    pmvnorm(lower = lower, upper = upper, mean = z.means[c(treatments, K +
      treatments2, 2 * K + treatments3)], sigma = cov.matrix[c(treatments, K +
      treatments2, 2 * K + treatments3), c(treatments, K + treatments2, 2 * K +
      treatments3)])[1]
  }

  # for the "select.best" method, find the probability that "select.treatment"
  # is selected and subsequently crosses the upper boundary
  rejection.paths <- function(
      selected.treatment, cut.off, treatments,
      cov.matrix, lower.boundary, upper.boundary, K, stage, z.means) {
    contrast <- diag(-1, K + stage - 1)
    contrast[1:K, selected.treatment] <- 1
    for (i in 1:(stage - 1)) contrast[K + i, K + i] <- 1

    bar.mean <- contrast %*% z.means[c(1:K, 1:(stage - 1) * K +
      selected.treatment)]

    bar.cov.matrix <- contrast %*% cov.matrix[c(1:K, 1:(stage - 1) * K +
      selected.treatment), c(1:K, 1:(stage - 1) * K + selected.treatment)] %*%
      t(contrast)

    lower <- c(rep(0, length(treatments)), cut.off)
    if (stage > 2) {
      lower <- c(
        rep(0, length(treatments)),
        lower.boundary[2:(stage - 1)], cut.off
      )
    }
    lower[which(treatments == selected.treatment)] <- lower.boundary[1]

    upper <- c(rep(Inf, length(treatments)), Inf)
    if (stage > 2) {
      upper <- c(
        rep(Inf, length(treatments)),
        upper.boundary[2:(stage - 1)], Inf
      )
    }
    upper[which(treatments == selected.treatment)] <- upper.boundary[1]

    pmvnorm(lower = lower, upper = upper, mean = bar.mean[c(
      treatments,
      K + 1:(stage - 1)
    )], sigma = bar.cov.matrix[c(treatments, K +
      1:(stage - 1)), c(treatments, K + 1:(stage - 1))])[1]
  }

  # for "all.promising" rule, this gives the cumulative typeI error for 'stage'
  # stages
  excess.alpha <- function(
      cut.off, alpha.star, treatments, cov.matrix,
      lower.boundary, upper.boundary, selection.method, K, stage, z.means) {
    # for "select.best" rule, this gives the Type I error spent at the
    # 'stage'th stage
    if (stage == 1) {
      return(1 - alpha.star[1] - pmvnorm(
        lower = rep(
          -Inf,
          length(treatments)
        ), upper = rep(cut.off, length(treatments)),
        mean = z.means[treatments], sigma = cov.matrix[treatments, treatments]
      )[1])
    }
    if (selection.method == "select.best") {
      return(sum(unlist(lapply(treatments,
        rejection.paths,
        cut.off = cut.off, treatments = treatments,
        cov.matrix = cov.matrix, lower.boundary = lower.boundary,
        upper.boundary = upper.boundary, K = K, stage = stage,
        z.means = z.means
      ))) - (alpha.star[stage] - alpha.star[stage - 1]))
    }
    # any of 'treatments' could be selected, so we add all these probabilities

    if (stage == 2) {
      surviving.subsets <- c(
        list(numeric(0)),
        lapply(as.list(1:(2^length(treatments) - 1)), get.hyp)
      )
      # list all possible subsets of surviving treatments after the first stage
      return(1 - alpha.star[2] - sum(unlist(lapply(surviving.subsets,
        get.path.prob,
        cut.off = cut.off, treatments = treatments,
        cov.matrix = cov.matrix, lower.boundary = lower.boundary,
        upper.boundary = upper.boundary, K = K, stage = stage,
        z.means = z.means
      ))))
    }
    surviving.subsets1 <- c(
      list(numeric(0)),
      lapply(as.list(1:(2^length(treatments) - 1)), get.hyp)
    )
    # all possible subsets of surviving treatments after the first stage
    surviving.subsets2 <- c(
      list(list(numeric(0))),
      lapply(surviving.subsets1[-1], function(x) {
        c(
          list(numeric(0)),
          lapply(as.list(1:(2^length(x) - 1)), get.hyp)
        )
      })
    )

    # for each possible subset of survivng subsets after stage 1,
    # list the possible subsets still surviving after stage 2
    1 - alpha.star[3] - sum(unlist(Map(function(x, y) {
      sum(unlist(lapply(y, get.path.prob,
        surviving.subset1 = x,
        cut.off = cut.off,
        treatments = treatments,
        cov.matrix = cov.matrix,
        lower.boundary = lower.boundary,
        upper.boundary = upper.boundary,
        K = K,
        stage = stage,
        z.means = z.means
      )))
    }, surviving.subsets1, surviving.subsets2)))
  }

  # give everything the correct name
  alpha.star <- current.mams$alpha.star
  l <- current.mams$l
  u <- current.mams$u
  selection.method <- current.mams$selection
  sample.sizes <- current.mams$sample.sizes
  sample.sizes[completed.stages, ] <- nobs
  # Update given the sample sizes actually observed
  if (!all(diff(sample.sizes) >= 0)) {
    stop("total sample size per arm cannot decrease between stages.")
  }
  J <- dim(sample.sizes)[1]
  K <- dim(sample.sizes)[2] - 1
  R <- sample.sizes[, -1] / sample.sizes[1, 1]
  r0 <- sample.sizes[, 1] / sample.sizes[1, 1]

  # get conditional distributions BEFORE seeing the new z scores
  cond.cov.matrix <- cov.matrix <- create.cov.matrix(r0, R)
  cond.mean <- rep(0, K * J)
  if (completed.stages > 1) {
    cond.cov.matrix <- create.cond.cov.matrix(
      cov.matrix, K,
      completed.stages - 1
    )
    cond.mean <- create.cond.mean(cov.matrix, K, completed.stages - 1,
      zscores = zscores[[completed.stages - 1]]
    )
  }

  # adjust upper boundaries in light of observed sample sizes:
  for (i in 1:(2^K - 1)) {
    treatments <- intersect(selected.trts[[completed.stages]], get.hyp(i))
    if ((length(treatments > 0)) && (alpha.star[[i]][J] > 0) &&
      (alpha.star[[i]][J] < 1)) {
      for (j in completed.stages:J) {
        try(
          new.u <- uniroot(excess.alpha, c(0, 10),
            alpha.star = alpha.star[[i]][completed.stages:J],
            treatments = treatments, cov.matrix = cond.cov.matrix,
            lower.boundary = l[[i]][completed.stages:J],
            upper.boundary = u[[i]][completed.stages:J],
            selection.method = selection.method, K = K,
            stage = j - completed.stages + 1, z.means = cond.mean
          )$root,
          silent = TRUE
        )
        if (is.null(new.u)) {
          stop("upper boundary not between 0 and 10")
        }
        u[[i]][j] <- round(new.u, 2)
      }
      l[[i]][J] <- u[[i]][J]
    }
  }
  if (J > completed.stages) {
    cond.cov.matrix <- create.cond.cov.matrix(cov.matrix, K, completed.stages)
    cond.mean <- create.cond.mean(
      cov.matrix, K, completed.stages,
      zscores[[completed.stages]]
    )
  }
  for (i in 1:(2^K - 1)) { # get conditional errors
    treatments <- intersect(selected.trts[[completed.stages]], get.hyp(i))
    if ((length(treatments > 0)) && (alpha.star[[i]][J] > 0) &&
      (alpha.star[[i]][J] < 1)) {
      max.z <- max(zscores[[completed.stages]][treatments])
      best.treatment <-
        treatments[which.max(zscores[[completed.stages]][treatments])]
      if (max.z <= u[[i]][completed.stages]) {
        alpha.star[[i]][completed.stages] <- 0
      }
      if (max.z > u[[i]][completed.stages]) {
        alpha.star[[i]][completed.stages:J] <- 1
        if (J > completed.stages) {
          l[[i]][(completed.stages + 1):J] <-
            u[[i]][(completed.stages + 1):J] <- -Inf
        }
      } else if (max.z <= l[[i]][completed.stages]) {
        alpha.star[[i]][completed.stages:J] <- 0
        if (J > completed.stages) {
          l[[i]][(completed.stages + 1):J] <-
            u[[i]][(completed.stages + 1):J] <- Inf
        }
      } else if (selection.method == "select.best") {
        for (j in (completed.stages + 1):J) {
          alpha.star[[i]][j] <- excess.alpha(
            cut.off = u[[i]][j],
            alpha.star = rep(0, J - completed.stages),
            treatments = best.treatment, cov.matrix = cond.cov.matrix,
            lower.boundary = l[[i]][(completed.stages + 1):J],
            upper.boundary = u[[i]][(completed.stages + 1):J],
            selection.method = selection.method, K = K,
            stage = j - completed.stages, z.means = cond.mean
          ) +
            alpha.star[[i]][j - 1]
        }
      } else {
        for (j in (completed.stages + 1):J) {
          alpha.star[[i]][j] <- excess.alpha(
            cut.off = u[[i]][j],
            alpha.star = rep(0, J - completed.stages), treatments = treatments,
            cov.matrix = cond.cov.matrix, lower.boundary =
              l[[i]][(completed.stages + 1):J],
            upper.boundary = u[[i]][(completed.stages + 1):J],
            selection.method = selection.method, K = K,
            stage = j - completed.stages, z.means = cond.mean
          )
        }
      }
    }
  }
  if (is.matrix(nfuture)) {
    sample.sizes[(completed.stages + 1):J, ] <- nfuture
    if (!all(diff(sample.sizes) >= 0)) {
      stop("total sample size per arm cannot decrease between stages.")
    }
    R <- sample.sizes[, -1] / sample.sizes[1, 1]
    r0 <- sample.sizes[, 1] / sample.sizes[1, 1]
    cov.matrix <- create.cov.matrix(r0, R)
    cond.cov.matrix <- create.cond.cov.matrix(cov.matrix, K, completed.stages)
    cond.mean <- create.cond.mean(cov.matrix, K, completed.stages,
      zscores = zscores[[completed.stages]]
    )
  }
  if (J > completed.stages) {
    for (i in 1:(2^K - 1)) {
      treatments <- intersect(selected.trts[[completed.stages + 1]], get.hyp(i))
      if ((length(treatments > 0)) && (alpha.star[[i]][J] > 0) &&
        (alpha.star[[i]][J] < 1)) {
        for (j in (completed.stages + 1):J) {
          try(
            new.u <- uniroot(excess.alpha, c(0, 10),
              alpha.star = alpha.star[[i]][(completed.stages + 1):J],
              treatments = treatments, cov.matrix = cond.cov.matrix,
              lower.boundary = l[[i]][(completed.stages + 1):J],
              upper.boundary = u[[i]][(completed.stages + 1):J],
              selection.method = selection.method, K = K,
              stage = j - completed.stages, z.means = cond.mean
            )$root,
            silent = TRUE
          )
          if (is.null(new.u)) {
            stop("upper boundary not between 0 and 10")
          }
          u[[i]][j] <- round(new.u, 2)
        }
        l[[i]][J] <- u[[i]][J]
      }
    }
  }

  res <- NULL
  res$l <- l
  res$u <- u
  res$sample.sizes <- sample.sizes
  res$K <- K
  res$J <- J
  res$alpha.star <- alpha.star
  res$selection <- selection.method
  res$zscores <- zscores
  res$selected.trts <- selected.trts
  class(res) <- "MAMS.stepdown"

  return(res)
}

###############################################################################
###################### tite.mams function #####################################
###############################################################################
#' @title Function to design multi-arm multi-stage studies with time-to-event
#' endpoints
#' @description The function determines (approximately) the boundaries of a
#' multi-arm multi-stage study with time-to-event endpoints for a given boundary
#' shape and finds the required number of events.
#' @param hr Interesting treatment effect on the scale of hazard ratios
#' (default=2).
#' @param hr0 Uninteresting treatment effect on the scale of hazard ratios
#' (default=1.2).
#' @param K Number of experimental treatments (default=4).
#' @param J Number of stages (default=2).
#' @param alpha One-sided familywise error rate (default=0.05).
#' @param power Desired power (default=0.9).
#' @param r Vector of allocation ratios (default=1:2).
#' @param r0 Vector ratio on control (default=1:2).
#' @param ushape Shape of upper boundary. Either a function specifying the shape
#'  or one of "pocock", "obf" (the default), "triangular" and "fixed".
#' @param lshape Shape of lower boundary. Either a function specifying the shape
#'  or one of "pocock", "obf", "triangular" and "fixed" (the default).
#' @param ufix Fixed upper boundary (default=NULL). Only used if shape="fixed".
#' @param lfix Fixed lower boundary (default=0). Only used if shape="fixed".
#' @param nstart Starting point for finding the sample size (default=1).
#' @param nstop Stopping point for finding the sample size (default=NULL).
#' @param sample.size Logical if sample size should be found as well
#' (default=TRUE).
#' @param Q Number of quadrature points per dimension in the outer integral
#' (default=20).
#' @param parallel if TRUE (default), allows parallelisation of the computation
#' via a user-defined strategy specified by means of the function
#' future::plan(). If not set differently, the default strategy is sequential,
#' which corresponds to a computation without parallelisation.
#' @param print if TRUE (default), indicate at which stage the computation is.
#' @details This function finds the (approximate) boundaries and sample size of
#' a multi-arm multi-stage study with time-to-event endpoints with K active
#' treatments plus control in which all promising treatments are continued at
#' interim analyses as described in Magirr et al (2012). It is a wrapper around
#' the basic mams function to facilitate its use with time-to-event endpoints,
#' following ideas of Jaki & Magirr (2013). Note that the sample size is
#' calculated as the required number of events, from which the total sample
#' size can be estimated (e.g., Whitehead 2001). See ?mams for further details
#' on the basic methodology.
#' @returns An object of the class MAMS containing the following components:
#'\item{l}{Lower boundary.}
#'\item{u}{Upper boundary.}
#'\item{n}{Sample size on control in stage 1.}
#'\item{N}{Maximum total sample size.}
#'\item{K}{Number of experimental treatments.}
#'\item{J}{Number of stages in the trial.}
#'\item{alpha}{Familywise error rate.}
#'\item{alpha.star}{Cumulative familywise error rate spent by each analysis.}
#'\item{power}{Power under least favorable configuration.}
#'\item{rMat}{Matrix of allocation ratios. First row corresponds to control
#' while subsequent rows are for the experimental treatments.}
#' @author Philip Pallmann, Dominic Magirr
#' @export
#' @examples
#' \dontrun{
#' ## An example 2-stage design with triangular efficacy and futility boundaries
#' tite.mams(hr=2, hr0=1.5, K=3, J=2, alpha=0.05, power=0.9,
#'           r=1:2, r0=1:2, ushape="triangular", lshape="triangular")
#' }
tite.mams <- function(hr = 1.5, hr0 = 1.1, K = 4, J = 2, alpha = 0.05,
                          power = 0.9, r = 1:2, r0 = 1:2, ushape = "obf",
                          lshape = "fixed", ufix = NULL, lfix = 0, nstart = 1,
                          nstop = NULL, sample.size = TRUE, Q = 20,
                          parallel = TRUE, print = TRUE) {
  if (hr0 < 1) {
    stop("The uninteresting effect must be 1 or larger.")
  }

  if (hr <= hr0) {
    stop("The interesting effect must be larger than the uninteresting one.")
  }

  p <- pnorm(log(hr) / sqrt(2))
  p0 <- pnorm(log(hr0) / sqrt(2))

  mams(
    K = K, J = J, alpha = alpha, power = power, r = r, r0 = r0, p = p, p0 = p0,
    delta = NULL, delta0 = NULL, sd = NULL, ushape = ushape, lshape = lshape,
    ufix = ufix, lfix = lfix, nstart = nstart, nstop = nstop,
    sample.size = sample.size, Q = Q, type = "tite", parallel = parallel,
    print = print
  )
}

###############################################################################
###################### print MAMS.stepdown function ###########################
###############################################################################
#' Generic print function for class MAMS.stepdown.
#' @details
#' print produces a brief summary of an object from class MAMS.stepdown
#' including boundaries and requires sample size if initially requested.
#' @param x An output object of class MAMS
#' @param digits Number of significant digits to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @return Text output.
#' @export
print.MAMS.stepdown <- function(x, digits=max(3, getOption("digits") - 4),
                                                                        ...) {
# find the nth intersection hypothesis (positions of 1s in binary n)
    get.hyp <- function(n) {
        indlength = ceiling(log(n)/log(2)+.0000001)
        ind = rep(0,indlength)
        newn=n

        for (h in seq(1,indlength)){
            ind[h] = (newn / (2^(h-1))) %% 2
            newn = newn - ind[h]*2^(h-1)
        }
        seq(1,indlength)[ind==1]
    }

    cat(paste("Design parameters for a ", x$J, " stage trial with ", x$K,
                                                " treatments\n\n",sep=""))
    res <- t(x$sample.sizes)
    colnames(res)<-paste("Stage",1:x$J)
    rownames(res) <- c("Cumulative sample size  (control):",
            paste("Cumulative sample size per stage (treatment ", 1:x$K, "):"))

    print(res)
    cat(paste("\nMaximum total sample size: ",
              sum(x$sample.sizes[x$J,]),"\n\n"))

    for (i in 1:length(x$l)){

        cat(paste("\nIntersection hypothesis H_{",
                              paste(get.hyp(i), collapse = " "), "}:","\n\n"))

        res <- matrix(NA,nrow=3,ncol=x$J)
        colnames(res)<-paste("Stage",1:x$J)
        rownames(res) <- c("Conditional error", "Upper boundary",
                                                "Lower boundary")
        res[1,] <- x$alpha.star[[i]]
        res[2,] <- x$u[[i]]
        res[3,] <- x$l[[i]]

        print(res)

    }
}

###############################################################################
###################### summary MAMS.stepdown function #########################
###############################################################################
#' Generic summary function for class MAMS.stepdown.
#' @details
#' print produces a brief summary of an object from class MAMS.stepdown
#' including boundaries and requires sample size if initially requested.
#' @param object An output object of class MAMS
#' @param digits Number of significant digits to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @return Text output.
#' @export
summary.MAMS.stepdown <-function(object,
                          digits=max(3, getOption("digits") - 4), ...) {

  print(object)

}

###############################################################################
###################### plot MAMS.stepdown function ############################
###############################################################################
#' Plot method for MAMS.stepdown objects
#' @description produces as plot of the boundaries.
#' @param x An output object of class MAMS.stepdown
#' @param col A specification for the default plotting color (default=`NULL`).
#' See `par` for more details.
#' @param pch Either an integer specifying a symbol or a single character to be
#' used as the default in plotting points (default=`NULL`). See `par`
#' for more details.
#' @param lty A specification for the default line type to be used between
#' analyses (default=`NULL`). Setting to zero suppresses plotting of the
#' lines. See `par` for more details.
#' @param main An overall title for the plot (default=`NULL`).
#' @param xlab A title for the x axis (default=`"Analysis"`).
#' @param ylab A title for the y axis (default=`"Test statistic"``).
#' @param ylim A title for the y axis (default=`"Test statistic"`).
#' @param type Type of plot to be used (default=`NULL`). See `plot`
#' for more details.
#' @param las A specification of the axis labeling style. The default `1`
#' ensures the labels are always horizontal. See `?par` for details.
#' @param bty Should a box be drawn around the legend? The default \code{"n"}
#' does not draw a box, the alternative option \code{"o"} does.
#' @param ... Further arguments passed to or from other methods.
#' @return Graphic output.
#' @author Thomas Jaki, Dominique-Laurent Couturier
#' @export
#' @examples \dontrun{
#' # 2-stage design with triangular boundaries
#' res <- mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2,
#'  p=0.65, p0=0.55, ushape="triangular", lshape="triangular", nstart=30)
#'
#' plot(res)
#' }
plot.MAMS.stepdown <- function(x, col=NULL, pch=NULL, lty=NULL, main=NULL,
xlab="Analysis", ylab="Test statistic", ylim=NULL, type=NULL, bty="n", las=1,
                                                                          ...) {
# find the nth intersection hypothesis (positions of 1s in binary n)
    get.hyp <- function(n) {
        indlength = ceiling(log(n)/log(2)+.0000001)
        ind = rep(0,indlength)
        newn=n

        for (h in seq(1,indlength)){
            ind[h] = (newn / (2^(h-1))) %% 2
            newn = newn - ind[h]*2^(h-1)
        }
        seq(1,indlength)[ind==1]
    }

    if (is.null(type)) type<-"p"
    if (is.null(bty)) bty<-"n"
    if (is.null(pch)) pch<-1
    if (is.null(las)) las<-1
    if (is.null(col)) col<-1:length(x$l)
    if (length(col) != length(x$l)) {
    stop("There must be as many colours as hypotheses.")
    }
    if (is.null(lty)) lty<-2
    if (is.null(ylim)) {

        lmin <- min(unlist(lapply(x$l,
                                    function(a) min(a[(a!=Inf) & (a!=-Inf)]))))
        if (!is.null(x$zscores)) lmin <- min(lmin,
                              min(unlist(x$zscores)[unlist(x$zscores) != -Inf]))
        umax <- max(unlist(lapply(x$u,
                                    function(a) max(a[(a!=Inf) & (a!=-Inf)]))))
        r <- umax - lmin
        ylim <- c(lmin - r/6, umax + r/6)

    }

    matplot(1:x$J, cbind(x$l[[1]], x$u[[1]]), type=type, pch=pch, main=main,
              col=0, ylab=ylab, xlab=xlab, ylim=ylim, axes=FALSE, las=las, ...)
    mtext(1:x$J,side=1,at=1:x$J)
    #  axis(side=2)
    axis(side=2,at=seq(-10,10,1), las=las)
    lines(x$u[[1]],lty=lty)
    lines(x$l[[1]][1:(x$J)],lty=lty)

    completed.stages <- length(x$zscores)
    if (completed.stages > 0) {
        for (i in 1:completed.stages){
            for (k in 1:x$K){
                points(i, x$zscores[[i]][k], col = 2 ^ (k - 1), pch = 3)
            }
        }
    }



    legend.text <- NULL

    #if (length(col) < length(x$l)) col <- rep(col, length(x$l))

    for (i in 1:length(x$l)){
        legend.text <- c(legend.text, paste("H_{",
                                      paste(get.hyp(i), collapse = " "), "}"))
        legend.col <- c(col, i)
        if ((x$alpha.star[[i]][x$J] > 0) && (x$alpha.star[[i]][x$J] < 1)) {

            matpoints(1:x$J, cbind(x$l[[i]], x$u[[i]]), type=type, pch=pch,
                  col=col[i], ylab=ylab, xlab=xlab, ylim=ylim, axes=FALSE, ...)

            lines(x$u[[i]],lty=lty, col=col[i])
            lines(x$l[[i]][1:(x$J)],lty=lty, col=col[i])
        }



    }

    legend("bottomright", legend=legend.text, bty=bty, lty=lty, col=col)
}
