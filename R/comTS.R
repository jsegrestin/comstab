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
#' the resulting simulation can be used to describe a community where the `'nsp'` species fluctuate independently (`'sync'` = "0").
#' To simulate an overall positive synchrony (`'sync'` = "1"), the temporal sequences 
#' of each species are sorted to maximize the number of years with all species 
#' having values above or below their respective median (one random selection among many possible combinations).
#' A stabilizing negative synchrony (`'sync'` = "-1") is simulated by sorting the temporal sequences
#' of each species to maximize the number of years where successively abundant species
#' have values above and below their respective median (one random selection among many possible combinations).
#' High positive synchrony (`'sync'` = "2") and high negative synchrony (`'sync'` = "-2") are generated using a similar approach but
#' sorting values based on four quartiles instead of using the median only.
#' 
#' The simulation uses a simplistic approach where species fluctuations are not related to any underlying environmental factor
#' nor demographic parameters. Consequently, the temporal order of the simulated abundances for each species cannot be considered realistic.
#' Nevertheless, this simplification has little influence on the analyses performed in this `R` package since none of the indices 
#' calculated depend on the temporal order of individual series, but rather describe the overall variability and temporal coordination of species.
#' 
#' @usage comTS(nsp, ny, even, mvs, sync = c("-2", "-1", "0", "1", "2"))
#' 
#' @param nsp Number of species in the community `(integer)`.
#' @param ny Length of the time series in years `(integer > 5)`.
#' @param even Parameter of the geometric rank-abundance curve 
#' ranging between 0 and 1 `(numeric)` .
#' @param mvs Scaling coefficient of the mean-variance 
#' relationship ranging between 1 and 2 `(numeric)`.
#' @param sync Level of synchrony between species `(factor)`:
#'  * `"0"`: independant fluctuations.
#'  * `"1"`: positive synchrony.
#'  * `"2"`: high positive synchrony.
#'  * `"-1"`: anti-synchronous fluctuations.
#'  * `"-2"`: high anti-synchronous fluctuations.
#' 
#' @return 
#' A `matrix` of `'ny'` rows and `'nsp'` columns, containing numerical values of species abundances.
#' The parameters used to compute species values (even, mvs, and sync) are stored as attributes of the matrix.
#' 
#' @examples 
#' require(stats)
#' 
#' comTS(nsp = 10, ny = 30, even = 0.6, mvs = 1.5, sync = "0")
#' 
#' @author Jules Segrestin, \email{jsegrestin@@gmail.com}
#' @importFrom stats runif rnorm var lm anova quantile
#' @export

comTS <- function(nsp, ny, even, mvs, sync = c("-2", "-1", "0", "1", "2")){
  
  if(!is.numeric(nsp)) stop("argument 'nsp' must be numeric")
  if(!is.numeric(ny)) stop("argument 'ny' must be numeric")
  if(ny < 5) stop("'ny' must be at least 5")
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
  K <- round(stats::runif(n = 1, min = 0.001, max = 1), 3)
  
  # sampling in normal distribution for each species
  sim <- sapply(X = x, function(m) stats::rnorm(n = ny, mean = m, sd = sqrt(K * m^mvs)))
  
  #removing negative values
  sim <- pmax(sim, 0)
  
  #check TPL
  check_TPL <- function(sim, K, mvs){
    nsp <- ncol(sim)
    x <- colMeans(sim)
    var_pred <- K * x^mvs
    var_obs <- apply(X = sim, MARGIN = 2, FUN = stats::var)
    TPL_m <- rep(log10(x), 2)
    TPL_v <- log10(c(var_pred, var_obs))
    TPL_type <- rep(c("pred", "obs"), each = nsp)
    TPL <- stats::lm(formula = TPL_v ~ TPL_m * TPL_type)
    TPL <- stats::anova(object = TPL)
    return(TPL$`Pr(>F)`[3] < .001)
  }
  check <- check_TPL(sim = sim, K = K, mvs = mvs)
  
  loop <- 1
  while(check & loop < 100){
    K <- round(stats::runif(n = 1, min = 0.001, max = 1), 3)
    sim <- sapply(X = x, function(m) stats::rnorm(n = ny, mean = m, sd = sqrt(K * m^mvs)))
    sim <- pmax(sim, 0)
    check <- check_TPL(sim = sim, K = K, mvs = mvs)
    loop <- loop + 1
  }
  if (loop == 100) stop("fail to simulate a community matching the input parameters: decrease species richness or increase community evenness")
  
  # sorting values according to the species synchrony level
  split_quantile <- function(x, probs) {
    q <- stats::quantile(x = x, probs = probs, names = F)
    if(any(duplicated(q))){
      q <- stats::quantile(x = x, probs = c(0, .5, 1), names = F)
      if(any(duplicated(q))) {
        return(sort(x, decreasing = T))
      } else {
        unlist(split(x = x, f = cut(x = x, breaks = q, include.lowest = TRUE))) 
      }
    } else {
      unlist(split(x = x, f = cut(x = x, breaks = q, include.lowest = TRUE)))
    }
  }
  
  if(sync == "1") {
    sim <- apply(X = sim, MARGIN = 2, FUN = function(x) split_quantile(x = x, probs = c(0, .5, 1)))
    sim <- sim[sample(1:ny, ny, replace = F), ]
  }
  
  if(sync == "2") {
    sim <- apply(X = sim, MARGIN = 2, FUN = function(x) split_quantile(x = x, probs = c(0, .25, .5, .75, 1)))
    sim <- sim[sample(x = 1:ny, size = ny, replace = F), ]
  }
  
  if(sync == "-1"){
    pos <- apply(X = sim[, seq(1, nsp, 2), drop = FALSE], MARGIN = 2, FUN = function(x) split_quantile(x = x, probs = c(0, .5, 1)))
    neg <- apply(X = sim[, seq(2, nsp, 2), drop = FALSE], MARGIN = 2, FUN = function(x) rev(split_quantile(x = x, probs = c(0, .5, 1))))
    sim <- cbind(pos, neg)
    sim <- sim[sample(x = 1:ny, size = ny, replace = F), order(colMeans(sim), decreasing = T)]
  }
  
  if(sync == "-2"){
    pos <- apply(X = sim[, seq(1, nsp, 2), drop = FALSE], MARGIN = 2, FUN = function(x) split_quantile(x = x, probs = c(0, .25, .5, .75, 1)))
    neg <- apply(X = sim[, seq(2, nsp, 2), drop = FALSE], MARGIN = 2, FUN = function(x) rev(split_quantile(x = x, probs = c(0, .25, .5, .75, 1))))
    sim <- cbind(pos, neg)
    sim <- sim[sample(x = 1:ny, size = ny, replace = F), order(colMeans(sim), decreasing = T)]
  }
  
  attr(sim, "even") <- c(even = even)
  attr(sim, "mvs") <- c(intercept = K, slope = mvs)
  attr(sim, "sync") <- c(sync = sync)
  return(sim)
}