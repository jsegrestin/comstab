#' Computing Taylor's power law
#' 
#' `TPL()` is a function used to compute the power function between
#' species temporal means and variances.
#' 
#' @usage TPL(z, ny = 1, method = c("OLS", "SMA"))
#' 
#' @param z A `matrix` containing repeated measurements of species abundances. 
#' The `matrix` must contain numerical values only, with years in rows and species in
#' columns. Remove any extra column.
#' @param ny Only species appearing more than `ny` years (`integer`, defaults to 1) are used in the calculations.
#' @param method Regression method used to estimate the scaling coefficient `(factor)`: 
#' * `OLS` (defaults) : Ordinary least square
#' * `SMA`: Standardized major axis, sometimes called reduced major axis
#' 
#' @return An object of class `'TPL'` is a list containing the following components:
#' * `'data'` a `data.frame` containing species temporal means and variance
#' * `'test'` the result of the Pearson's correlation test between log10-transformed species means and variances (see `cor.test` in the `stats` package)
#' * `a` the intercept of the power law
#' * `b` the scaling coefficient of the power law
#' 
#' @examples
#' require(smatr)
#' 
#' # Simulates a custom community time series using 'comTS()':
#' z <- comTS(nsp = 10, ny = 30, even = 0.6, mvs = 1.5, sync = "0")
#' 
#' # Computes Taylor's power law:
#' TPL(z)
#' 
#' @author Jules Segrestin, \email{jsegrestin@@gmail.com}
#' @importFrom stats var coef lm cor.test
#' @importFrom smatr sma
#' 
#' @export

TPL <- function(z, ny = 1, method = c("OLS", "SMA")){
  
  if(!is.matrix(z)) stop("Error: z is not a matrix")
  if(!is.numeric(z)) stop("Error: non-numerical values in z")
  if(any(z < 0)) stop("Error: negative values in z")
  if(dim(z)[1] == 1) stop("Error: single-row matrix")
  if(!is.numeric(ny)) stop("ny must be numeric")
  method <- match.arg(method)
  
  # Replace Nas with zeros
  z[is.na(z)] <- 0
  # Remove absent species
  z <- z[, colSums(z) > 0, drop = FALSE]
  # Remove transient species appearing only ny years
  nyi <- apply(X = z, MARGIN = 2, FUN = function(x) sum(x > 0))
  z <- z[, nyi > ny, drop = FALSE] 
  
  if(dim(z)[2] == 1) stop("Only one species found. Taylor's power law cannot be computed")
  
  meani <- colMeans(z)
  vari <- apply(X = z, MARGIN = 2, FUN = stats::var)
  
  test <- stats::cor.test(log10(meani), log10(vari))
  
  if(method == "OLS") TPL <- stats::coef(stats::lm(log10(vari) ~ log10(meani)))
  if(method == "SMA") TPL <- smatr::sma(vari ~ meani, log = "xy")$coef[[1]][, 1]
  
  res <- list(data = data.frame(mean = meani, var = vari), test = test, 
              method = method, a = TPL[1], b = TPL[2])
  class(res) <- "TPL"
  return(res)
}

#' @export
print.TPL <- function(x, ...){
  cat("\nTaylor's power law")
  cat("\n")
  print(x$test)
  cat(paste0(x$method, "-based scaling coefficient b = ", round(x$b, 2)))
}