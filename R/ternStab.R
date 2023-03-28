#' Plotting the relative contribution of stabilizing effects
#' 
#' `ternStab()` is a graph function used to represent the relative contributions
#' of the three stabilizing effects ("Dominance", "Asynchrony" and "Averaging") on a
#' ternary plot.
#' 
#' @usage ternStab(x, ..., point = TRUE, add = FALSE)
#' 
#' @param x object of class `'comstab'`.
#' @param ... other parameters to be passed through to plotting functions.
#' @param point plot the community on the ternary plot (`logical`, defaults to TRUE)
#' @param add add the community on the current plot window (`logical`, defaults to FALSE)
#' 
#' @examples
#' require(Ternary)
#' 
#' # Simulates a custom community time series using 'comTS()':
#' z <- comTS(nsp = 10, ny = 30, even = 0.6, mvs = 1.5, sync = "0")
#' 
#' # Runs the partitioning of the community coefficient of variation:
#' x <- partitionR(z)
#' 
#' # Plots the relative contributions
#' par(mar = c(0, 0, 0, 0))
#' ternStab(x)
#' 
#' # Adds a second community on the ternary plot
#' z2 <- comTS(nsp = 15, ny = 30, even = .7, mvs = 1.1, sync = "1")
#' x2 <- partitionR(z2)
#' ternStab(x2, add = TRUE, col = "red")
#' 
#' @author Jules Segrestin, \email{jsegrestin@@gmail.com}
#' @import Ternary
#' @export

ternStab <- function(x, ..., point = TRUE, add = FALSE){
  if(!inherits(x, 'comstab')) stop("x must be an object of class 'comstab'.")
  
  # calculate the relative contributions of each stabilizing components
  rel <- log(x[2:4]) / log(prod(x[2:4]))
  
  # plot the result
  if(!add){
    Ternary::TernaryPlot(alab = "Dominance contribution \U2192", 
                         blab = "Asynchrony contribution \U2192",
                         clab = "\U2190 Averaging contribution",
                         point = "up", grid.minor.lines = 0, 
                         grid.lty = "solid", grid.lwd = .5,
                         axis.rotate = FALSE,
                         padding = .1, ...)
  }
  if(point){
    Ternary::TernaryPoints(rel, ...)
  }
}