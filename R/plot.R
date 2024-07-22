#' Plotting a `comstab` object
#' 
#' Plotting method for object inheriting from class "comstab".
#' 
#' @usage 
#' \method{plot}{comstab}(x, \dots, xlab = "", ylab = "Log scale", cex.comp = 1)
#' 
#' @param x object of class `'comstab'`.
#' @param ... other parameters to be passed through to plotting functions.
#' @param xlab a label for the x axis, removed by defaults.
#' @param ylab a label for the y axis, defaults to 'Log scale'.
#' @param cex.comp A numerical value giving the label size of stabilizing components.
#' This is an absolute measure, not scaled by par("cex").
#' 
#' @return No return value, graphical function.
#' 
#' @examples
#' require(graphics)
#' 
#' # Simulates a custom community time series using 'comTS()':
#' z <- comTS(nsp = 10, ny = 30, even = 0.6, mvs = 1.5, sync = "0")
#' 
#' # Runs the partitioning of the community coefficient of variation:
#' x <- partitionR(z)
#' 
#' # Plots the result
#' plot(x)
#' 
#' @author Jules Segrestin, \email{jsegrestin@@gmail.com}
#' @import graphics
#' @export

plot.comstab <- function(x, ..., xlab = "", ylab = "Log scale", cex.comp = 1){
  dots <- list(...)
  y <- x$CVs
  graphics::plot(x = 1:4, y = y, t = "b", log = "y", xaxt = "n",
                 xlab = xlab, ylab = ylab, ...)
  graphics::axis(side = 1, at = 1:4, labels = c(expression(CV[e]), expression(widetilde(CV)), expression(CV[a]), expression(CV[com])), ...)
  graphics::mtext(text = c("Dominance", "Asynchrony", "Averaging"), side = 1, at = 1.5:3.5, line = -1, cex = cex.comp)
  
  if("cex.axis" %in% names(dots)){
    graphics::mtext(text = c(expression(Delta), expression(psi), expression(omega)), side = 1, at = 1.5:3.5, line = -1-cex.comp, cex = dots$cex.axis)
  }else{
    graphics::mtext(text = c(expression(Delta), expression(psi), expression(omega)), side = 1, at = 1.5:3.5, line = -1-cex.comp, cex = par()$cex.axis)
  }
}

#' Plotting a `TPL` object
#' 
#' Plotting method for object inheriting from class "TPL".
#' 
#' @usage 
#' \method{plot}{TPL}(x, \dots, xlab = "Mean", ylab = "Variance", log = "xy")
#' 
#' @param x object of class `'comstab'`.
#' @param ... other parameters to be passed through to plotting functions.
#' @param xlab a label for the x axis, defaults to 'Mean".
#' @param ylab a label for the y axis, defaults to 'Variance'.
#' @param log Defaults to 'xy' to show the log-log relationship.
#' 
#' @return No return value, graphical function.
#' 
#' @examples
#' require(graphics)
#' 
#' # Simulates a custom community time series using 'comTS()':
#' z <- comTS(nsp = 10, ny = 30, even = 0.6, mvs = 1.5, sync = "0")
#' 
#' # Computes Taylor's power law:
#' x <- TPL(z)
#' 
#' # Plots the result
#' plot(x)
#' 
#' @author Jules Segrestin, \email{jsegrestin@@gmail.com}
#' @import graphics
#' @export
#' 

plot.TPL <- function(x, ..., xlab = "Mean", ylab = "Variance", log = "xy"){
  xx <- x$data[, 1]
  yy <- x$data[, 2]
  graphics::plot(x = xx, y = yy, xlab = xlab, ylab = ylab, log = log, ...)
  graphics::mtext(text = paste0("cor = ", round(x$test$estimate, 2),
                                "; Pval ", ifelse(x$test$p.value < 0.001, "< 0.001", 
                                                  paste0("= ", round(x$test$p.value, 3)))), line = -1)
  graphics::mtext(text = paste0("b = ", round(x$b, 2)), line = -2)
  
  if(x$test$p.value < 0.05) graphics::abline(a = x$a, b = x$b, ...)
}