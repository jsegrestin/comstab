#' Partitioning of the temporal CV of ecological communities
#' 
#' `PartitionR()` is a function used to partition the temporal coefficient of variation of a community
#'  into the variability of the average species and three stabilizing effects: the dominance, asynchrony and averaging effects
#' (see Details).
#' 
#' @usage partitionR(z, ny = 1)
#' 
#' @param z A `matrix` containing repeated measurements of species abundances. 
#' The `matrix` must contain numerical values only, with years in rows and species in
#' columns. Remove any extra column.
#' @param ny Only species appearing more than `ny` years (`integer`, defaults to 1) are used in the calculations.
#' 
#' @return Returns an object of class `'comstab'`.
#' @return An object of class `'comstab'` is a list containing the following components:
#'  * `'CVs'` a named vector of calculated coefficient of variations. `CVe` is the CV of an average species,
#'  `CVtilde` is the mean of species CVs weighted by their relative abundances, `CVa` is the expected community CV if 
#'   the community was stabilized by species asynchrony only, and `CVc` is the observed community CV.
#'  * `'Stabilization'` a named vector of the stabilizing effects. `tau` is the total stabilization, `Delta` is
#'  the dominance effect, `Psi` is the asynchrony effect, and `omega` is the averaging effect.
#'  * `'Relative'` a named vector of the relative contributions of each stabilizing effect to the total stabilization.
#'  `Delta_cont`, `Psi_cont`, and `omega_cont` are the relative contribution of respectively, the dominance, asynchrony, and averaging effects to the total stabilization.
#'  Returns a vector of NAs if any Stabilizing effect is higher than 1.
#'  
#' @details The analytic framework is described in details in Segrestin *et al.* (2024).
#' In short, the partitioning relies on the following equation: \deqn{CV_{com} = CV_e \Delta \Psi \omega} 
#' where \eqn{CV_{com}} is the community coefficient of variation (reciprocal of community stability), 
#' \eqn{CV_e} is the expected community CV when controlling for the dominance structure and species temporal synchrony,
#' \eqn{ \Delta} is the dominance effect, \eqn{ \Psi} is the asynchrony effect, and \eqn{ \omega} is the averaging effect.
#'
#' @references Segrestin *et al.* (2024) A unified framework for partitioning the drivers of stability of ecological communities. Global Ecology and Biogeography, https://doi.org/10.1111/geb.13828
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
#' @importFrom stats var coef lm cor.test setNames
#' 
#' @author Jules Segrestin, \email{jsegrestin@@gmail.com}
#' 
#' @export

partitionR <- function(z, ny = 1){
  
  if(!is.matrix(z)) stop("Error: z is not a matrix")
  if(!is.numeric(z)) stop("Error: non-numerical values in z")
  if(dim(z)[1] == 1) stop("Error: single-row matrix")
  if(!is.numeric(ny)) stop("ny must be numeric")
  
  # Remove absent species
  z <- z[, colSums(z) > 0, drop = F]
  
  # Number of years each species is recorded
  nyi <- apply(X = z, MARGIN = 2, FUN = function(x) sum(x > 0))
  z <- z[, nyi > ny] #remove transient species appearing only ny years
  n <- ncol(z)
  
  # Total CV
  varsum <- stats::var(rowSums(z))
  meansum <- mean(rowSums(z))
  CV <- sqrt(varsum)/meansum
  if(CV == 0) stop("The community CV is zero. This analysis does not apply to perfectly stable communities.")
  
  if(dim(z)[2] == 1) {
    
    warning("This analysis is not relevant for single-species communities. All stabilizing effects were fixed to 1.")
    
    # outputs
    CVs <- stats::setNames(object = c(CV, CV, CV, CV),
                           nm = c("CVe", "CVtilde", "CVa", "CVc"))
    Stabilization <- stats::setNames(object = c(1, 1, 1, 1),
                                     nm = c("tau", "Delta", "Psi", "omega"))
    Relative <- stats::setNames(object = rep(NA, 3),
                                nm = c("Delta_cont", "Psi_cont", "omega_cont"))
    res <- list(CVs = CVs, Stabilization = Stabilization, Relative = Relative)
    class(res) <- "comstab"
    return(res)
  }
  
  if(dim(z)[2] > 1) {
    
    # Expected community CV if all species had even abundances
    vari <- apply(X = z, MARGIN = 2, FUN = stats::var)
    meani <- colMeans(z)
    CVi <- sqrt(vari) / meani
    TPL <- stats::coef(stats::lm(log10(CVi) ~ log10(meani)))
    CVe <- 10^TPL[1] * (meansum / n) ^ TPL[2]
    
    if(dim(z)[2] > 2){
      testcor <- stats::cor.test(log10(CVi), log10(meani))$p.value > 0.05
      if (testcor) warning("No significant power law between species CVs and abundances.")
    }
    
    # Dominance effect
    sumsd <- sum(sqrt(vari))
    CVtilde <- sumsd / meansum
    Delta <- CVtilde / CVe
    if(Delta > 1) warning("Destabilizing effect of dominants. Relative effects cannot be computed.")
    
    # Compensatory dynamics
    sdsum <- sqrt(varsum)
    rootPhi <- sdsum / sumsd
    
    # Asynchrony effect
    sumvar <- sum(vari)
    beta <- log10(1/2) / (log10(sumvar / (sumsd^2)))
    Psi <- rootPhi^beta
    
    # averaging effect
    omega <- rootPhi / Psi
    if(omega > 1) warning("Community diversity is lower than the null diversity. Relative effects cannot be computed.")
    
    # total stabilization
    tau <- Delta * Psi * omega
    
    # outputs
    CVs <- stats::setNames(object = c(CVe, CVtilde, CVtilde * Psi, CV),
                           nm = c("CVe", "CVtilde", "CVa", "CVc"))
    Stabilization <- stats::setNames(object = c(tau, Delta, Psi, omega),
                                     nm = c("tau", "Delta", "Psi", "omega"))
    if(any(Stabilization > 1)){
      Relative <- stats::setNames(object = rep(NA, 3),
                                  nm = c("Delta_cont", "Psi_cont", "omega_cont"))
    } else {
      Relative <- stats::setNames(object = c(log10(Delta)/log10(tau), log10(Psi)/log10(tau), log10(omega)/log10(tau)),
                                  nm = c("Delta_cont", "Psi_cont", "omega_cont"))
    }
    
    res <- list(CVs = CVs, Stabilization = Stabilization, Relative = Relative)
    class(res) <- "comstab"
    return(res)
  }
}

#' @export
print.comstab <- function(x, ...){
  cat("\nPartitionning of the community temporal variability (CV)")
  cat("\nSee Segrestin et al. (2024)")
  cat("\n")
  cat(paste0("Community CV = ", round(x$CVs["CVc"], 2),
             "\nTotal stabilization = ", round(x$Stabilization["tau"], 2),
             "\nDominance effect = ", round(x$Stabilization["Delta"], 2),
             "\nAsynchrony effect = ", round(x$Stabilization["Psi"], 2),
             "\nAveraging effect = ", round(x$Stabilization["omega"], 2)))
  cat("\n")
  cat("\nRelatives contributions:")
  cat(paste0("\n% Dominance = ", round(x$Relative["Delta_cont"], 2)),
      paste0("\n% Asynchrony = ", round(x$Relative["Psi_cont"], 2)),
      paste0("\n% Averaging = ", round(x$Relative["omega_cont"], 2)))
}