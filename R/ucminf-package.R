################################################################################
#                                                                              #
# ucminf package doc                                                           #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# ucminf package overview                                                      #
# Algorithms:  quasi-Newton type +                                             #
#              BFGS updating of the inverse Hessian +                          #
#              soft line search with a trust region                            #
#------------------------------------------------------------------------------#

#' ucminf: General-Purpose Unconstrained Non-Linear Optimization
#'
#' The \pkg{ucminf} package provides an algorithm for general-purpose
#' unconstrained non-linear optimization.
#'
#' @name ucminf-package
#'
#' @aliases ucminf-package
#'
# @docType package
#' _UCMINF
#'
#' @useDynLib ucminf, mfopt, .registration = TRUE, .fixes = "C_"
#'
#' @section Bugreport: Any bug or suggestion can be reported using the
#' \code{ucminf} tracker facilities at:
#' \url{https://github.com/hdakpo/ucminf/issues}
#'
#' @author K Hervé Dakpo, Hans Bruun Nielsen, and Stig Bousgaard Mortensen
#'
# @importFrom numDeriv hessian
# @importFrom calculus jacobian
NULL
