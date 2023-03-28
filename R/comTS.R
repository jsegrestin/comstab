#' Simulating time series for ecological communities
#' 
#' `comTS()` is a function used to simulate `com`munity `T`ime `S`eries 
#' based on custom parameters.
#' 
#' @details
#' 
#' The simulation produces temporal series of abundances of `'nsp'` species for `'ny'` years.
#' The mean abundance of each species is determined by a geometric series of `'nsp'` values using a constant
#' `'even'` ratio between successive terms. Hence, a high `'even'` value means a community with even species
#' abundances while a low `'even'` means a strongly dominated community. Species temporal variances are
#' calculated following Taylor’s power law, using with a `'mvs'` scaling coefficient.
#' Finally, for each species, the `'ny'` abundance values are sampled from a normal distribution
#' with the corresponding species parameters, and bounded to positive values. Since each species is simulated independently of others,
#' the resulting simulation can be used to describe a community where the `'nsp'` species fluctuate independently (`'sync'` = 0).
#' To simulate an overall positive synchrony (`'sync'` = 1), the temporal sequences 
#' of each species are sorted to maximize the number of years with all species 
#' having values above or below their respective means (one random selection among many possible combinations).
#' A stabilizing negative synchrony (`'sync'` = -1) is simulated by sorting the temporal sequences
#' of each species to maximize the number of years where successively abundant species
#' have values above and below their respective means (one random selection among many possible combinations).
#' 
#' The simulation use a simplistic approach where species fluctuations are not related to any underlying environmental factor
#' nor demographic parameters. Consequently, the temporal order of the simulated abundances for each species cannot be considered realistic.
#' Nevertheless, this simplification has little influence on the analyses performed in this `R` package since none of the indices 
#' calculated depend on the temporal order of individual series, but rather describe the overall variability and temporal coordination of species.
#' 
#' @usage comTS(nsp, ny, even, mvs, sync = c("-1", "0", "1"))
#' 
#' @param nsp Number of species in the community `(integer)`.
#' @param ny Length of the time series in years `(integer)`.
#' @param even Parameter of the geometric rank-abundance curve 
#' ranging between 0 and 1 `(numeric)` .
#' @param mvs Scaling coefficient of the mean-variance 
#' relationship ranging between 1 and 2 `(numeric)`.
#' @param sync Level of synchrony between species `(factor)`:
#'  * `"0"`: independant fluctuations.
#'  * `"+1"`: positive synchrony.
#'  * `"-1"`: anti-synchronous fluctuations.
#' 
#' @return A `matrix` of `'ny'` rows and `'nsp'` columns, containing numerical values of species abundances.
#' 
#' @examples 
#' require(stats)
#' 
#' comTS(nsp = 10, ny = 30, even = 0.6, mvs = 1.5, sync = "0")
#' 
#' @author Jules Segrestin, \email{jsegrestin@@gmail.com}
#' @importFrom stats rnorm
#' @export

comTS <- function(nsp, ny, even, mvs, sync = c("-1", "0", "1")){
  
  if(!is.numeric(nsp)) stop("argument 'nsp' must be numeric")
  if(!is.numeric(ny)) stop("argument 'ny' must be numeric")
  if(!is.numeric(even)) stop("argument 'even' must be numeric")
  if(even <= 0 | even >= 1) stop("argument 'even' must range between 0 and 1")
  if(!is.numeric(mvs)) stop("argument 'mvs' must be numeric")
  if(mvs  <= 1 | mvs >= 2) stop("argument 'mvs' must range between 1 and 2")
  match.arg(sync)
  
  # geometric rank abundance curve
  x <- 100
  for (i in 1:(nsp-1)) x <- c(x, even * x[i])
  x <- x / sum(x) * 100
  
  # Mean-Variance scaling intercept
  K <- .1
  
  # sampling in normal distribution for each species
  sim <- sapply(x, function(m) stats::rnorm(ny, m, sqrt(K * m^mvs)))
  
  # sorting species values according to their synchrony level
  if(sync == "+1") {
    sim <- apply(sim, MARGIN = 2, function(x) c(x[x >= mean(x)], x[x < mean(x)]))
    sim <- sim[sample(1:ny, ny, replace = F),]
  }
  
  if(sync == "-1"){
    pos <- apply(sim[, seq(1, nsp, 2)], MARGIN = 2, function(x) c(x[x >= mean(x)], x[x < mean(x)]))
    neg <- apply(sim[, seq(2, nsp, 2)], MARGIN = 2, function(x) c(x[x < mean(x)], x[x >= mean(x)]))
    sim <- cbind(pos, neg)
    sim <- sim[sample(1:ny, ny, replace = F), order(colMeans(sim), decreasing = T)]
  }
  
  #removing negative values
  sim[sim < 0] <- 0
  return(sim)
}