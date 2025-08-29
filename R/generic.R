#' Shows changes and news
#'
#' Functions showing changes since previous versions.
#'
#' Displays the changes and news given in the NEWS file of the package.
#' @usage MAMSNews()
#' @return Screen output.
#' @author Thomas Jaki
#' @export
#' @examples \dontrun{MAMSNews()}
MAMSNews <- function() {
  file.show(system.file("NEWS", package="MAMS"))
}

#' Function to design multi-arm multi-stage studies with normal endpoints
#' @name mams
#' @aliases mams MAMS
#'
#'The function determines the boundaries of a multi-arm multi-stage study
#'for a given boundary shape and finds the required number of subjects.
#' @usage mams(obj = NULL, K=NULL, J=NULL, alpha=NULL, power=NULL, r=NULL, 
#'                      r0=NULL, p=NULL, p0=NULL, delta=NULL, delta0=NULL, 
#'                      sd=NULL,ushape=NULL, lshape=NULL, ufix=NULL, lfix=NULL, 
#'                      nstart=NULL, nstop=NULL, sample.size=NULL, Q=NULL,
#'                      type=NULL, parallel=NULL, print=NULL, nsim=NULL, 
#'                      H0=NULL, method=NULL)
#' @param obj object of class `MAMS`
#' @param K Number of experimental treatments (default=4).
#' @param J Number of stages (default=2).
#' @param alpha One-sided familywise error rate (default=0.05).
#' @param power Desired power (default=0.9).
#' @param r Vector of allocation ratios (default=1:2).
#' @param r0 Vector ratio on control (default=1:2).
#' @param p Interesting treatment effect on the probability scale.
#' See Details (default=0.75).
#' @param p0 Uninteresting treatment effect on the probability scale.
#' See Details (default=0.5). Note that this parameter affects the sample size 
#' estimation to ensure that Treatment 1 is best only when selecting the 
#' simultaneous stopping rule (i.e., 'method = "simultaneous"') and, for all 
#' methods, as values for treatment arms 2 to K in the simulations under H1.
#' @param delta Interesting treatment effect on the traditional scale.
#' See Details (default=NULL).
#' @param delta0 Uninteresting treatment effect on the traditional scale.
#' See Details (default=NULL). Note that this parameter affects the sample size 
#' estimation to ensure that Treatment 1 is best only when selecting the 
#' simultaneous stopping rule (i.e., 'method = "simultaneous"') and, for all 
#' methods, as values for treatment arms 2 to K in the simulations under H1.
#' @param sd Standard deviation, assumed to be known. 
#' See Details (default=NULL).
#' @param ushape Shape of upper boundary. Either a function specifying the shape
#' or one of "pocock", "obf" (the default), "triangular" and "fixed".
#' See details.
#' @param lshape Shape of lower boundary. Either a function specifying the shape
#' or one of "pocock", "obf", "triangular" and "fixed" (the default). 
#' See details.
#' @param ufix Fixed upper boundary (default=NULL). Only used if shape="fixed".
#' @param lfix Fixed lower boundary (default=0). Only used if shape="fixed".
#' @param nstart Starting point for finding the sample size (default=1).
#' @param nstop Stopping point for finding the sample size (default=NULL).
#' @param sample.size Logical if sample size should be found as well 
#' (default=TRUE).
#' @param Q Number of quadrature points per dimension in the outer integral
#' (default=20).
#' @param type Will be changed automatically by the wrappers tite.mams()
#' (to "tite") and ordinal.mams() (to "ordinal") to customise the output.
#' @param parallel if TRUE (default), allows parallelization of the computation
#' via a user-defined strategy specified by means of the function 
#' future::plan().
#'  If not set differently, the default strategy is sequential, which 
#' corresponds to a computation without parallelization.
#' @param print if TRUE (default), indicate at which stage the computation is.
#' @param nsim a scalar indicating the number of simulations (default = 50'000, 
#' minimum = 1000) 
#' @param H0 if TRUE (default), the simulation also considers the case with all
#' effect sizes set to 0.
#' @param method Type of the desired design: `simultaneous`(default) for 
#' simultaneous stopping rules, `sep` for separate stopping, and `dtl`
#' for drop-the-losers design.
#' @returns An object of the class MAMS containing the following components:
#' \item{l}{Lower boundary.}
#' \item{u}{Upper boundary.}
#' \item{n}{Sample size on control in stage 1.}
#' \item{N}{Maximum total sample size.}
#' \item{K}{Number of experimental treatments.}
#' \item{J}{Number of stages in the trial.}
#' \item{alpha}{Familywise error rate.}
#' \item{alpha.star}{Cumulative familywise error rate spent by each analysis.}
#' \item{power}{Power under least favorable configuration.}
#' \item{rMat}{Matrix of allocation ratios. First row corresponds to control
#' while subsequent rows are for the experimental treatments.}
#' \item{sim}{a list indicating, for each hypothesis of interest (null and/or 
#' alternative), the expected sample size and standard deviation per group
#'  (ess), the cumulated probability of efficacy (efficacy) and futility 
#' (futility) per treatment arm and look}
#' \item{input}{the list of all input parameters except K and J}
#' @details
#' This function finds the boundaries and sample size of a multi-arm multi-stage
#'  study with K active treatments plus control in which all promising
#'  treatments are continued at interim analyses as described in
#'  Magirr et al (2012). At each interim analysis the test statistics are
#'  compared to the lower (futility) bound and any treatment whose corresponding
#'  test statistic falls below that bound is discontinued. Similarly if any test
#'  statistic exceeds the upper (efficacy) bound the null hypothesis
#'  corresponding to that treatment can be rejected and superiority of that
#'  treatment over control claimed. At the same time the study is stopped.
#'  If at least one test statistic exceeds the lower bound and none exceeds the
#'  upper bound the study is continued and further patients are recruited to all
#'  remaining experimental treatments plus control.
#'
#'  The design is found under the least favorable configuration, which requires
#'  an interesting treatment effect `p` that if present we would like to find
#'  with
#'  high probability and an uninteresting effect `p0`. Both `p` and `p0` are
#'  parameterized as \eqn{P(X_k > X_0 ) = p}{P(X_k > X_0 ) = p}, that is the
#'  probability of a randomly selected person on treatment k observing a
#'  better outcome than a random person on control. For `p=0.5` the experimental
#'  treatment and control perform equally well.
#'  The advantage of this parameterization is that no knowledge about the
#'  variance is required. To convert traditional effect sizes, 
#'  \eqn{\delta}{delta}
#'  to this format use 
#'  \eqn{p=\Phi(\frac{\delta}{\sqrt{2}\sigma})}{Phi(delta/(2^0.5*sigma))}.
#'  Alternatively, the interesting and uninteresting effect size can also be
#'  specified directly on the traditional scale of `delta` and `delta` with an
#'  additional specification of the standard deviation `sd` assumed to be known.
#'
#'  The shape of the boundaries (`ushape`, `lshape`) are either using the 
#'   predefined shapes following Pocock (1977), O'Brien & Fleming (1979) or the 
#'   triangular Test (Whitehead, 1997) using options "`pocock"`, `"obf"`or 
#'   `"triangular"` respectively, are constant (option `"fixed"`) or supplied 
#'   in as a function.
#'   If a function is passed it should require exactly one argument specifying
#'   the number of stages and return a vector of the same length. The lower
#'   boundary shape is required to be non-decreasing while the upper boundary
#'   shape needs to be non-increasing. If a fixed lower boundary is used, `lfix`
#'   must be smaller than \eqn{\Phi^{-1}(1-\alpha)/2}{Phi(1-alpha)/2}
#'   to ensure that it is smaller than the upper boundary.
#'
#'   The default starting point for finding the sample size is `nstart=1`, and
#'   the default point where the search is stopped (when `nstop=NULL`) is 3 
#'   times the sample size of the corresponding fixed single-stage design.
#'
#'   Computation of designs with more than four stages are very time consuming
#'   and not advised. The parameter `sample.size` controls whether the required
#'   sample size is computed as well. Setting to `FALSE` approximately halves 
#'   the computation time.
#'
#'   For designs with more than 2 stages, parallelization of the computation by
#'   means of the packages \pkg{future} and \pkg{future.apply} lead to decreased
#'   computation times when choosing a parallelization strategy like, for
#'   example, `multicore` (using separate forked R processes, available to
#'   unix/osx users) or `multisession` (using separate R sessions, available
#'   to all users) (refer to \link[future:plan]{future::plan()} for detail).
#'
#' @author Thomas Jaki, Dominic Magirr, Dominique-Laurent Couturier
#' and Nikita Mozgunov
#' @references
#' Jaki T., Pallmann P., and Magirr D. (2019), \emph{The R Package MAMS for
#' Designing Multi-Arm Multi-Stage Clinical Trials}, \bold{Journal of
#' Statistical Software}, 88(4), 1-25. Link:
#' \doi{10.18637/jss.v088.i04}
#'
#' Magirr D., Jaki T., and Whitehead J. (2012), \emph{A generalized Dunnett test
#' for multi-arm multi-stage clinical studies with treatment selection},
#' \bold{Biometrika}, 99(2), 494-501. Link:
#' \doi{10.1093/biomet/ass002}
#'
#' Pocock S.J. (1977), \emph{Group sequential methods in the design and analysis
#' of clinical trials}, \bold{Biometrika}, 64(2), 191-199.
#'
#' O'Brien P.C., and Fleming T.R. (1979), \emph{A multiple testing procedure for
#' clinical trials}, \bold{Biometrics}, 35(3), 549-556.
#'
#' Whitehead J. (1997), \emph{The Design and Analysis of Sequential
#' Clinical Trials}, \bold{Wiley}: Chichester, UK.
#' Wason J, Stallard N, Bowden J, Jennison C. A multi-stage drop-the-losers 
#' design for multi-arm clinical trials. Statistical Methods in Medical 
#' Research. 2017;26(1):508-524. doi:10.1177/0962280214550759
#' @seealso [new.bounds], [ordinal.mams], [tite.mams], [MAMS].
#' @keywords design
#' @export
#' @examples
#' \dontrun{
#' ## A fixed sample (single stage) design specified on the p scale
#' mams(K=4, J=1, alpha=0.05, power=0.9, r=1, r0=1, p=0.65, p0=0.55)
#'
#' ## The same design specified on the delta scale
#'
#' mams(K=4, J=1, alpha=0.05, power=0.9, r=1, r0=1, p=NULL, p0=NULL,
#'      delta=0.545, delta0=0.178, sd=1)
#'
#' ## An example in Table 1 of Magirr et al (2012)
#' # 2-stage design with O'Brien & Fleming efficacy and zero futility boundary
#'
#' mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2, p=0.65, p0=0.55,
#'      ushape="obf", lshape="fixed", lfix=0, nstart=40)
#' 
#' ## An example of separate stopping rules
#' # 2-stage design with O'Brien & Fleming efficacy and zero futility boundary
#'
#' mams(method = "sep",K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2, 
#'        p=0.65, p0=0.55, ushape="obf", lshape="fixed", lfix=0, nstart=40)
#' 
#' # An example of running drop-the-losers design 
#' # `K` should be defined as vector length of J defining allocation arms per 
#' # stages with final element equal to 1.
#' mams(method = "dtl", K=c(4,1), J=2, alpha=0.05,  
#'          power=0.9, r=1:2, r0=1:2, p=0.65, p0=0.55, ushape="obf", 
#'          lshape="fixed", lfix=0, nstart=40)
#' 
#' 
#'
#' # Note that these examples may take a few minutes to run
#'
#' ## 3-stage design with Triangular efficacy and futility boundary
#' mams(K=4, J=3, alpha=0.05, power=0.9, r=1:3, r0=1:3, p=0.65, p0=0.55,
#'   ushape="triangular", lshape="triangular", nstart=30)
#'
#' ## Different allocation ratios between control and experimental treatments.
#' ## Twice as many patients are randomized to control at each stage.
#'   mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=c(2, 4), p=0.65, 
#'   p0=0.55, ushape="obf", lshape="fixed", lfix=0, nstart=30)
#'
#'
#'   ##
#'   ## example considering different parallelization strategies
#'   ##
#'
#'
#'   # parallel = FALSE (future framework not used)
#'   set.seed(1)
#'   system.time(
#'   print(mams(K=4, J=3, alpha=0.05, power=0.9, r=1:3, r0=1:3,
#'           p=0.65, p0=0.55, ushape="triangular", lshape="triangular",
#'           nstart=30, parallel = FALSE))
#'           )
#'  # parallel = TRUE (default) with default strategy (sequential computation)
#'  plan(sequential)
#'  set.seed(1)
#'  system.time(
#'  print(mams(K=4, J=3, alpha=0.05, power=0.9, r=1:3, r0=1:3,
#'           p=0.65, p0=0.55, ushape="triangular", lshape="triangular", 
#'           nstart=30))
#'           )
#'  # parallel = TRUE(default) with multisession strategy (parallel computation)
#'  plan(multisession)
#'  set.seed(1)
#'  system.time(
#'  print(mams(K=4, J=3, alpha=0.05, power=0.9, r=1:3, r0=1:3, 
#'           p=0.65, p0=0.55, ushape="triangular", lshape="triangular",
#'            nstart=30))
#'           )
#'           plan("default")
#'           }
mams <- function(
  obj = NULL, K = NULL, J = NULL, alpha = NULL, power = NULL, r = NULL, 
  r0 = NULL, p = NULL, p0 = NULL, delta = NULL, delta0 = NULL, sd = NULL,
  ushape = NULL, lshape = NULL, ufix = NULL, lfix = NULL, 
  nstart = NULL, nstop = NULL, sample.size = NULL, Q = NULL,
  type = NULL, parallel = NULL, print = NULL, nsim = NULL, 
  H0 = NULL, method = NULL) {


  if (!is.null(obj)) {
  obj = unpack_object(obj)
  if (!inherits(obj, "MAMS")) {
      stop("Only works on MAMS objects generated with MAMS-3.0.0 and above")
    }
  }


  if (!is.null(attributes(obj)$altered)) {
  stop("This object was altered by sim function.")
  }

  # Store the match call for reference
  mc = match.call()
  # Default parameters in a list
  defaults <- list(
    K=4, J=2, alpha=0.05, power=0.9, r=NULL, r0=NULL, p=0.75, p0=0.5,
    delta=NULL, delta0=NULL, sd=NULL, ushape="obf", lshape="fixed",
    ufix=NULL, lfix=0, nstart=1, nstop=NULL, sample.size=TRUE, Q=20,
    type="normal", parallel=TRUE, print=TRUE, nsim=50000, H0=TRUE, 
    method="simultaneous"
  )
 
  # Convert the provided arguments into a list
  user_defined <- list(
    K=K, J=J, alpha=alpha, power=power, r=r, r0=r0, p=p, p0=p0, 
    delta=delta, delta0=delta0, sd=sd, ushape=ushape, lshape=lshape, 
    ufix=ufix, lfix=lfix, nstart=nstart, nstop=nstop, sample.size=sample.size, 
    Q=Q, type=type, parallel=parallel, print=print, nsim=nsim, H0=H0, 
    method=method
  )
    user_defined <- Filter(Negate(is.null), user_defined)

    if (any(!is.null(user_defined[["p"]]), !is.null(user_defined[["p0"]]))) {
    defaults$delta <- NULL
    defaults$delta0 <- NULL
    defaults$sd <- NULL
  } else if (any(!is.null(user_defined$delta), !is.null(user_defined$delta0), 
              !is.null(user_defined$sd))) {
    defaults$p <- NULL
    defaults$p0 <- NULL
  }


  # If obj is provided, use its values
  if (!is.null(obj)) {
    if (!inherits(obj, "MAMS")) {
      stop("Only works on MAMS objects generated with MAMS-3.0.0 and above")
    }
    obj_params <- obj
    obj_params$method <- attr(obj, "method") 

  if (!is.null(user_defined[["method"]])) {
    if (obj_params[["method"]] != user_defined[["method"]] & 
                                  obj_params[["method"]] == "dtl") {
      obj_params[["K"]] = obj_params[["K"]][1]
      warning(paste("K is set to", obj_params[["K"]][1]))
    }
  }
    if (any(!is.null(user_defined[["p"]]), !is.null(user_defined[["p0"]]))) {
    obj_params$delta <- NULL
    obj_params$delta0 <- NULL
    obj_params$sd <- NULL
  } else if (any(!is.null(user_defined$delta), !is.null(user_defined$delta0), 
              !is.null(user_defined$sd))) {
    obj_params$p <- NULL
    obj_params$p0 <- NULL
  }

    final_params <- modifyList(obj_params, user_defined)
    
  } else {
    final_params <- modifyList(defaults, user_defined)
  }

if (!is.null(final_params$delta)) {
    if (is.null(final_params$delta0)) {
    final_params$delta0  <- 0
      warning("delta0 set to 0")
    }
    if (is.null(final_params$sd)) {
      final_params$sd  <- 1
      warning("sd set to 1")
    }
}

# Check if 'J' was provided by the user
  user_J <- "J" %in% names(user_defined)
  
  if (user_J) {
    if (!"r" %in% names(user_defined)) {
      final_params$r <- 1:final_params$J
    }
    if (!"r0" %in% names(user_defined)) {
      final_params$r0 <- 1:final_params$J
    }
  } else {
    if (is.null(final_params$r)) {
      final_params$r <- 1:final_params$J
    }
    if (is.null(final_params$r0)) {
      final_params$r0 <- 1:final_params$J
    }
  }
  obj <- structure(final_params, class = "MAMS", mc = mc)
  attr(obj, "method") <- final_params$method

  obj <- switch(obj$method,
                "simultaneous" = mams.check.simultaneous(obj),
                "dtl"          = mams.check.dtl(obj),
                "sep"          = mams.check.sep(obj),
                cat("\nThis is an object with an unknown method\n")
  )
  obj <- switch(obj$method,
                "simultaneous" = mams.fit.simultaneous(obj),
                "dtl"          = mams.fit.dtl(obj),
                "sep"          = mams.fit.sep(obj),
                cat("\nThis is an object with an unknown method\n")
  )
  return(obj)
}


#' Simulating multi-arm multi-stage designs
#'
#' @param obj An object of class MAMS
#' @description
#' The function simulates multi-arm multi-stage designs and estimates power and
#' expected sample size.
#' @usage mams.sim(obj=NULL,nsim=NULL, nMat=NULL,
#' u=NULL, l=NULL, pv=NULL, deltav=NULL, sd=NULL, ptest=NULL,
#' parallel=NULL, H0=NULL, K = NULL)
#' @param obj an object of class `MAMS`. The parameters/design of the
#' considered in the output of the `mams()` function are considered as
#' reference for the simulation. If other parameters are given, their values
#' *override* the parameters of the `MAMS` object
#' @param nsim Number of simulations (default=`50000``).
#' @param nMat Jx(K+1) dimensional matrix of observed/expected sample sizes.
#' Rows correspond to stages and columns to arms. First column is control
#' (default: `NULL`).
#' @param u Vector of previously used upper boundaries (default=`NULL`).
#' @param l Vector of previously used upper boundaries (default=`NULL`).
#' @param pv Vector of size K of true treatment effects on the probability
#' scale. See Details (default=`NULL`).
#' @param deltav Vector of size K of true treatment effects on the traditional
#' scale. See Details (default=`NULL`).
#' @param sd Standard deviation. See Details (default=`NULL`).
#' @param ptest Vector of treatment numbers for determining power.
#' For example, c(1, 2) will count rejections of one or both hypotheses for
#' testing treatments 1 and 2 against control (default=`1`).
#' @param parallel if `TRUE` (default), allows parallelization of the
#' computation via a user-defined strategy specified by means of the function
#' \code{\link[future]{plan}}. If not set differently, the default
#' strategy is `sequential`, which corresponds to a computation without
#' parallelization.
#' @param H0 if `TRUE` (default), the simulation also considers the case with
#' all effect sizes set to 0.
#' @param K Allocation for treatment arms (used only with method = "dtl")
#' @returns An object containing the following components:
#' \item{l}{Lower boundary.}
#' \item{u}{Upper boundary.}
#' \item{n}{Sample size on control in stage 1.}
#' \item{N}{Maximum total sample size.}
#' \item{K}{Number of experimental treatments.}
#' \item{J}{Number of stages in the trial.}
#' \item{rMat}{Matrix of allocation ratios. First row corresponds to control 
#' and second row to experimental treatments.}
#' \item{nsim}{Number of simulation runs.}
#' \item{typeI}{The proportion any hypothesis is rejected.}
#' \item{power}{The proportion the first hypothesis is rejected and the 
#'  corresponding test statistic is largest.}
#' \item{ptest}{The vector `ptest`.}
#' \item{prop.rej}{The proportion of times at least one of the hypothesis 
#' specified by `ptest` is rejected.}
#' \item{exss}{The expected sample size.}
#' @export
#' @author Thomas Jaki, Dominic Magirr and Dominique-Laurent Couturier
#' @details
#' This function simulates multi-arm multi-stage studies for a given matrix of
#' sample sizes and boundaries given by the vectors `u` and  `l`.
#' The effect difference between each experimental treatment and control is
#' given by `pv` and is parameterized as
#'  \eqn{P(X_k > X_0 ) = p}{P(X_k > X_0 ) = p}.
#' That is the probability of a randomly selected person on treatment k
#' observing a better outcome than a random person on control. 
#' For `pv=rep(0.5,4)` the experimental treatments and control perform equally 
#' well (i.e. the global null hypothesis is true). 
#' The advantage of this parameterization is that no knowledge about the 
#' variance is required. To convert traditional effect sizes, 
#' \eqn{\delta}{delta} to this format use
#' \eqn{p=\Phi(\frac{\delta}{\sqrt{2}\sigma})}{Phi(delta/(2^0.5*sigma))}.
#' Alternatively, the effect size can also be specified directly on the
#' traditional scale of `deltav` with an additional specification of
#' the standard deviation `sd`.
#'
#' he function returns the probability of rejecting any hypothesis (`typeI`),
#' the power to reject the first hypothesis when the first treatment has the
#' largest estimated effect, the proportion of rejections of the hypothesis
#' specified by `ptest` (`prop.rej`) as well as the expected
#' sample size.
#' @references
#' Jaki T., Pallmann P., and Magirr D. (2019), \emph{The R Package MAMS for
#' Designing Multi-Arm Multi-Stage Clinical Trials}, \bold{Journal of
#' Statistical Software}, 88(4), 1-25. Link:
#' \doi{10.18637/jss.v088.i04}
#'
#' Magirr D., Jaki T., and Whitehead J. (2012), \emph{A generalized Dunnett test
#' for multi-arm multi-stage clinical studies with treatment selection},
#' \bold{Biometrika}, 99(2), 494-501. Link:
#' \doi{10.1093/biomet/ass002}
#' @seealso [mams], [MAMS].
#' @keywords design
#' @examples
#' \dontrun{
#'# Note that some of these examples may take a few minutes to run
#'
#'# 2-stage design with O'Brien & Fleming efficacy and zero futility boundary
#'# with equal sample size per arm and stage. Design can be found using
#'# mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2, ushape="obf", 
#'      # lshape="fixed",
#'      # lfix=0, p=0.65, p0=0.55)
#'
#'# under global null hypothesis (using the pv scale)
#'mams.sim(nsim=10000, nMat=matrix(c(44, 88), nrow=2, ncol=5), 
#'              u=c(3.068, 2.169),
#'              l=c(0.000, 2.169), pv=rep(0.5, 4), ptest=1)
#'
#'# under global null hypothesis (using the deltav scale)
#'mams.sim(nsim=10000, nMat=matrix(c(44, 88), nrow=2, ncol=5), 
#'         u=c(3.068, 2.169),
#'         l=c(0.000, 2.169), pv=NULL, deltav=rep(0, 4), sd=1, ptest=1)
#'
#'# under LFC
#'mams.sim(nsim=10000, nMat=matrix(c(44, 88), nrow=2, ncol=5), 
#'         u=c(3.068, 2.169),
#'         l=c(0.000, 2.169), pv=c(0.65, 0.55, 0.55, 0.55), ptest=1:2)
#'
#'# when all treatments doing similarly well
#'mams.sim(nsim=10000, nMat=matrix(c(44, 88), nrow=2, ncol=5),
#'         u=c(3.068, 2.169),
#'         l=c(0.000, 2.169), pv=c(0.63, 0.62, 0.60, 0.61), ptest=4)
#'
#'##
#'## example considering different parallelisation strategies
#'##
#'
#'# parallel = FALSE (future framework not used)
#'set.seed(1)
#'system.time(
#'  print(mams.sim(nsim=25000, nMat=matrix(c(44, 88), nrow=2, ncol=5), 
#'                 u=c(3.068, 2.169),
#'                 l=c(0.000, 2.169), pv=c(0.65, 0.55, 0.55, 0.55),
#'                 ptest=1:2, parallel=FALSE))
#')
#'# parallel = TRUE (default) with default strategy (sequential computation)
#'plan(sequential)
#'set.seed(1)
#'system.time(
#'  print(mams.sim(nsim=25000, nMat=matrix(c(44, 88), nrow=2, ncol=5), 
#'                 u=c(3.068, 2.169),
#'                 l=c(0.000, 2.169), pv=c(0.65, 0.55, 0.55, 0.55), ptest=1:2))
#')
#'# parallel = TRUE (default) with multisession strategy (parallel computation)
#'plan(multisession)
#'set.seed(1)
#'system.time(
#'  print(mams.sim(nsim=25000, nMat=matrix(c(44, 88), nrow=2, ncol=5),
#'                 u=c(3.068, 2.169),
#'                 l=c(0.000, 2.169), pv=c(0.65, 0.55, 0.55, 0.55), ptest=1:2))
#')
#'plan("default")
#' }
mams.sim <- function(obj=NULL,nsim=NULL, nMat=NULL,
                    u=NULL, l=NULL, pv=NULL, deltav=NULL, sd=NULL, ptest=NULL,
                    parallel=NULL, H0=NULL, K = NULL) {
  if (!is.null(obj)) {

  switch(attr(obj, "method"),
          "simultaneous" = mams.sim.simultaneous(obj=obj, nsim=nsim, 
                                      nMat=nMat, u=u, l=l, pv=pv, 
                                      deltav=deltav, sd=sd, ptest=ptest, 
                                      parallel = parallel, H0=H0),
          "dtl"          = mams.sim.dtl(obj = obj, nsim = nsim, nMat =nMat, 
                                            u = u, l = l, pv = pv, 
                                            deltav = deltav, sd = sd, 
                                            ptest = ptest, H0 = H0, K = K),
          "sep"          = mams.sim.sep(obj=obj, nsim=nsim, nMat=nMat, u=u, 
                                          l=l, pv=pv, deltav=deltav, sd=sd, 
                                          ptest=ptest, H0=H0),
          cat("\nThis is an object with an unknown method\n")
  )
} else { 
mams.sim.simultaneous(obj=NULL, nsim, nMat, u, l, pv, 
                          deltav, sd, ptest, parallel, H0)
}
}

#' Generic print function for class MAMS.
#' @details
#' print produces a brief summary of an object from class MAMS 
#' including boundaries and requires sample size if initially requested.
#' @param x An output object of class MAMS
#' @param digits Number of significant digits to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @return Text output.
#' @author Thomas Jaki, Dominic Magirr, Philip Pallmann
#' @references
#' Magirr D, Jaki T, Whitehead J (2012) A generalized Dunnett test for multi-arm
#' multi-stage clinical studies with treatment selection. 
#' Biometrika, 99(2), 494-501. Stallard N, Todd S (2003) Sequential designs for 
#' phase III clinical trials
#' incorporating treatment selection. Statistics in Medicine, 22(5), 689-703.
#' Magirr D, Stallard N, Jaki T (2014) Flexible sequential designs for multi-arm
#' clinical trials. Statistics in Medicine, 33(19), 3269-3279.
#' @seealso \link{mams}, \link{stepdown.mams}.
#' @keywords classes
#' @export
#' @examples \dontrun{
#' # 2-stage design with triangular boundaries
#' res <- mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2, 
#'              p=0.65,p0=0.55,
#'              ushape="triangular", lshape="triangular", nstart=30)
#'
#' print(res)
#' }
print.MAMS <- function(x, digits=max(3, getOption("digits") - 4), ...) {
  switch(attr(x, "method"),
                "simultaneous" = mams.print.simultaneous(x, digits, ...),
                "dtl"          = mams.print.dtl(x,digits,...),
                "sep"          = mams.print.sep(x,digits,...),
                cat("\nThis is an object with an unknown method\n")
  )
}

#' Generic summary function for class MAMS.
#' @description Produces a detailed summary of an object from class MAMS
#' @param object An output object of class MAMS
#' @param digits Number of significant digits to be printed.
#' @param extended TRUE or FALSE
#' @param ... Further arguments passed to or from other methods.
#' @return Text output.
#' @author Dominique-Laurent Couturier
#' @export
#' @examples \dontrun{
#' # 2-stage design with triangular boundaries
#' res <- mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2,
#'              p=0.65, p0=0.55,
#'              ushape="triangular", lshape="triangular", nstart=30)
#'
#' summary(res)
#' }
summary.MAMS <- function(object, digits=max(3, getOption("digits") - 4),
                            extended=FALSE,...) {

  switch(attr(object, "method"),
                "simultaneous" = mams.summary.simultaneous(object,
                                                  digits, extended=FALSE, ...),
                "dtl"          = mams.summary.dtl(object,digits,...),
                "sep"          = mams.summary.sep(object,digits,...),
                cat("\nThis is an object with an unknown method\n")
  )
}

#' Plot method for MAMS objects
#' @description produces as plot of the boundaries.
#' @param x An output object of class MAMS
#' @param ask A logical indicating if R should wait for the next plot to be 
#' displayed.
#' @param which A vector indicating which plot(s) to define. `1` displays
#' the efficacy and futility limits per stage, `2` displays the efficacy
#' nd futility probabilities per stage, `1:2` (default) displays both.
#' @param new A logical indicating if the new plot of the futility and efficacy
#' limits should be displayed (default=`TRUE``).
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
plot.MAMS <- function(x, ask=TRUE, which=1:2, new = TRUE, col=NULL, 
                          pch=NULL, lty=NULL, main=NULL, xlab="Analysis", 
                          ylab="Test statistic", ylim=NULL, type=NULL, 
                          las=1, ...) {

  switch(attr(x, "method"),
                "simultaneous" = mams.plot.simultaneous(x,...),
                "dtl"          = cat("There is no plot for this method"),
                "sep"          = mams.plot.sep(x,...),
                cat("\nThis is an object with an unknown method\n")
  )
}
##################### Service functions #######################################
# Function to pack the object
pack_object <- function(obj) {
  # Define the top-level keys that should stay at the top level
  top_level_keys <- c("l", "u", "n", "N", "K", "J", "alpha", "alpha.star", 
                      "power", "rMat", "sim")
  
  input <- list()
  packed_obj <- list()
  
  # Copy top-level elements directly to the packed object
  for (key in names(obj)) {
    if (key %in% top_level_keys) {
      packed_obj[[key]] <- obj[[key]]
    } else {
      # Move everything else to input
      input[[key]] <- obj[[key]]
    }
  }
  
  if ("par" %in% names(obj) && is.list(obj$par)) {
    input$par <- obj$par
  }
  
  packed_obj$input <- input
  
  obj$attributes  <- attributes(obj)
  attrs = obj$attributes[names(obj$attributes) != "names"]

  for (i in names(attrs)) {
  attr(packed_obj, names(attrs[i]))  <- attrs[[i]]
  }
  
    return(packed_obj)
}
# Function to `unpack` the object
unpack_object <- function(packed_obj) {
  unpacked_obj <- list()
  
  for (key in names(packed_obj)) {
    if (key != "input") {
      unpacked_obj[[key]] <- packed_obj[[key]]
    }
  }
  
  if ("input" %in% names(packed_obj) && is.list(packed_obj$input)) {
    for (key in names(packed_obj$input)) {
      unpacked_obj[[key]] <- packed_obj$input[[key]]
    }
  }
  
  attrs =  attributes(packed_obj)[2:length(attributes(packed_obj))]

  for (i in names(attrs)) {
  attr(unpacked_obj, names(attrs[i])) <- attrs[[i]]
  }
attributes(unpacked_obj)
  return(unpacked_obj)
}
