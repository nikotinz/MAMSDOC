#' @name MAMS
#' @aliases MAMS-package MAMS
#' @title Designing Multi-Arm Multi-Stage Studies
#' @description
#' This package allows to design multi-arm multi-stage (MAMS)
#' studies with asymptotically normal endpoints and known variance.
#' It considers normal, binary, ordinal and time-to-event endpoints 
#' in which either the single best treatment or all promising treatments 
#' are continued at the interim analyses. 
#'
#' @details
#' Currently implemented functions are:
#' \itemize{
#'   \item \code{\link[MAMS:mams]{mams()}}: a function allowing to design multi-arm multi-stage studies with normal endpoints,
#'   \item \code{\link[MAMS:new.bounds]{new.bounds()}}: a function allowing to update the lower and upper boundaries of a multi-arm multi-stage study, typically initially defined by \code{\link[MAMS:mams]{mams()}}, based on observed sample sizes,
#'   \item \code{\link[MAMS:mams.sim]{mams.sim()}}: a function allowing to simulate multi-arm multi-stage studies given chosen boundaries and sample size, and estimates power and expected sample size,
#'   \item \code{\link[MAMS:stepdown.mams]{stepdown.mams()}}: a function allowing to find stopping boundaries for a 2- or 3-stage (step-down) multiple-comparisons-with-control test,
#'   \item \code{\link[MAMS:stepdown.update]{stepdown.update()}}: a function allowing to update the stopping boundaries of a multi-arm multi-stage study, typically initially defined by \code{\link[MAMS:stepdown.mams]{stepdown.mams()}}, at an interim analysis as well as allowing for unplanned treatment selection and/or sample-size reassessment,
#'   \item \code{\link[MAMS:ordinal.mams]{ordinal.mams()}}: a function allowing to design multi-arm multi-stage studies with ordinal or binary endpoints,
#'   \item \code{\link[MAMS:tite.mams]{tite.mams()}}: a function allowing to design multi-arm multi-stage studies with time-to-event endpoints.
#' }
#'
#' We refer to Jaki et al (2019) for an overview of the package as well as to Magirr et al (2012) and Magirr et al (2014) for theoretical details.
#'
#' \bold{Parallelisation}
#'
#' Since version 2.0.0, \pkg{MAMS} relies on the package \pkg{future.apply} for parallel computation. The package \pkg{future.apply} is part of the \pkg{future} parallelisation framework that requires users to define their parallelisation strategy by means of the function \code{\link[future:plan]{future::plan()}}. This function takes several options like, for example, \code{sequential} (default strategy corresponding to a computation without parallelisation), \code{multicore} (using separate forked \pkg{R} processes, available to unix/osx users) and \code{multisession} (using separate \pkg{R} sessions, available to all users). We refer to Bengtsson H. (2022) for an overview of the \pkg{future} framework.
#'
#' Note that, for the functions of \pkg{MAMS} to be available to workers defined by \code{\link[future:plan]{future::plan()}}, \pkg{MAMS} has to be installed at a location available under \code{\link{.libPaths}} (by default, \pkg{R} installs packages in the directory corresponding to the first element of \code{\link{.libPaths}}).  
#'
#' \bold{Reproducibility}
#'
#' Results of the \pkg{MAMS} package for studies involving more than 2 stages are seed-dependent (as the Gaussian quadrature integration of the multivariate normal distribution relies on probabilities estimated by means of the randomised Quasi-Monte-Carlo procedure of Genz and Bretz in \code{\link[mvtnorm:pmvnorm]{mvtnorm::pmvnorm()}}).
#'
#' Results are reproducible if a seed is set before the evaluation of a function of the \pkg{MAMS} package (typically by means of the function \code{\link{set.seed}}).
#'
#' When \code{parallel=TRUE}, the \pkg{future} package assigns independent streams of L'Ecuyer pseudo-random numbers to each parallelised task, allowing results to be reproducible when a seed is set, even when using a different parallelisation strategy and/or a different number of workers. When \code{parallel=FALSE}, the random number generation is handled by base \pkg{R} directly instead of by the \pkg{future} package, so that, if the number of stages is larger than 2, evaluations using the same seed will not lead to the same exact results with \code{parallel=FALSE} and \code{parallel=TRUE}. 
#'
#' @author
#' Thomas Jaki, Dominique-Laurent Couturier, Dominic Magirr, Nikita Mozgunov, Philip Pallmann
#'
#' Maintainer: Thomas Jaki \email{thomas.jaki@pm.me}.
#'
#' @references
#' Jaki T., Pallmann P. and Magirr D. (2019), \emph{The R Package MAMS for Designing Multi-Arm Multi-Stage Clinical Trials}, \bold{Journal of Statistical Software}, 88(4), 1-25. Link: \doi{10.18637/jss.v088.i04}
#'
#' Magirr D., Jaki T. and Whitehead J. (2012), \emph{A generalized Dunnett test for multi-arm multi-stage clinical studies with treatment selection}, \bold{Biometrika}, 99(2), 494-501. Link: \doi{10.1093/biomet/ass002}
#'
#' Magirr D., Stallard N. and Jaki T. (2014), \emph{Flexible sequential designs for multi-arm clinical trials}, \bold{Statistics in Medicine}, 33(19), 3269-3279. Link: \doi{10.1002/sim.6183}
#'
#' Bengtsson H. (2022), \emph{A Unifying Framework for Parallel and Distributed Processing in R using Futures}, to appear in \bold{The R Journal}. Link: \href{https://journal.r-project.org/archive/2021/RJ-2021-048/index.html}{accepted version}
#'
#' @keywords package MAMS
#' @useDynLib MAMS, .registration=TRUE
#' @import mvtnorm
#' @import utils
#' @import future
#' @import future.apply
#' @import cli
#' @importFrom graphics axis legend lines matplot matpoints mtext points abline
#' @importFrom graphics layout par text title
#' @importFrom utils globalVariables
#' @importFrom stats dnorm pnorm qnorm rnorm uniroot qbinom quantile var 
#' @importFrom stats setNames
#' @importFrom grDevices dev.new devAskNewPage gray
#' @importFrom methods is
"_PACKAGE"
