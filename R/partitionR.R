#' Partitioning community CV
#' 
#' `PartitionR()` is a function used to partition the temporal coefficient
#' of variation of a community into the mean coefficient of variation of species
#' and three stabilizing components: the dominance, asynchrony and averaging effects
#' (see Details).
#' 
#' @usage partitionR(z, n = 1)
#' 
#' @param z A `matrix` containing repeated measurements of species abundances. 
#' The `matrix` must contain numerical values only, with years in rows and species in
#' columns. Remove any extra column.
#' @param n Only species appearing more than `n` years (`integer`, defaults to 1) are used in the calculations.
#' 
#' @return Returns an object of class `'comstab'`. It contains the mean coefficient of variation of species (\eqn{\overline{CV}}),
#'  the dominance (\eqn{ \Delta}), asynchrony (\eqn{ \Psi}) and averaging (\eqn{ \omega}) effects.
#'  
#' @details The analytic framework is described in details in Segrestin *et al.* (in prep).
#' In short, the partitioning relies on the following equation: \deqn{CV=\overline{CV} \Delta \Psi \omega} 
#' where \eqn{CV} is the community coefficient of variation (reciprocal of community stability), 
#' \eqn{\overline{CV}} is the mean species coefficient of variation, \eqn{ \Delta} is the dominance effect,
#' \eqn{ \Psi} is the asynchrony effect, and \eqn{ \omega} is the averaging effect.
#'
#' @references Segrestin *et al.* (in prep) A unified framework for partitioning the drivers of stability of ecological communities
#' 
#' @examples
#' require(stats)
#' 
#' # Simulates a custom community time series using 'comTS()':
#' z <- comTS(nsp = 10, ny = 30, even = 0.6, mvs = 1.5, sync = "0")
#' 
#' # Runs the partitioning of the community coefficient of variation:
#' partitionR(z)
#' 
#' @importFrom stats var
#' 
#' @author Jules Segrestin, \email{jsegrestin@@gmail.com}
#' 
#' @export

partitionR <- function(z, n = 1){
  
  if(!is.matrix(z)) stop("Error: z is not a matrix")
  if(!is.numeric(z)) stop("Error: non-numerical values in z")
  if(dim(z)[1] == 1) stop("Error: single-row matrix")
  if(dim(z)[2] == 1) warning("This analysis is not relevant for single-species communities")
  if(!is.numeric(n)) stop("n must be numeric")
  
  # Remove absent species
  z <- z[, colSums(z) > 0, drop = F]
  
  # Number of years each species is recorded
  ni <- apply(X = z, MARGIN = 2, FUN = function(x) sum(x > 0))
  z <- z[, ni > n] #remove transient species appearing only n years
  
  # Total CV
  varsum <- stats::var(rowSums(z))
  meansum <- mean(rowSums(z))
  CV <- sqrt(varsum)/meansum
  
  # Mean CV
  vari <- apply(X = z, MARGIN = 2, FUN = stats::var)
  meani <- colMeans(z)
  CVi <- sqrt(vari) / meani
  CVbar <- mean(CVi)
  
  # Delta
  sumvar <- sum(vari)
  sumsd <- sum(sqrt(vari))
  CVtilde <- sumsd / meansum
  Delta <- CVtilde / CVbar
  
  # Asynchrony
  rootPhi <- sqrt(varsum) / sumsd
  beta <- log(1/2) / (log(sumvar/(sumsd^2)))
  Psi <- rootPhi^beta
  CVpsi <- CVtilde * Psi
  omega <- rootPhi / Psi
  
  res <- c(CVbar, Delta, Psi, omega)
  names(res) <- c("CVbar", "Delta", "Psi", "omega")
  class(res) <- "comstab"
  return(res)
}

#' @export
print.comstab <- function(x, ...){
  cat("\nPartitionning of the community temporal variability (CV)")
  cat("\nSee Segrestin et al. (2023)")
  cat("\n")
  cat(paste0("\nCV = ", round(prod(x), 2),
             "\nMean CV = ", round(x[1], 2),
             "\nDominance effect = ", round(x[2], 2),
             "\nAsynchrony effect = ", round(x[3], 2),
             "\nAveraging effect = ", round(x[4], 2)))
  cat("\n")
  cat("\nRelatives contributions:")
  rel <- log(x[2:4]) / log(prod(x[2:4]))
  cat(paste0("\n% Dominance = ", round(rel[1], 2)),
      paste0("\n% Asynchrony = ", round(rel[2], 2)),
      paste0("\n% Averaging = ", round(rel[3], 2)))
}