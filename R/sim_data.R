#' Generate Simulated Data
#'
#' It returns simulated data. This is corresponding to high dimensional linear model simulation settings in the manuscript.
#'
#' @param seed Random seed.
#'
#' @param n Number of subjects.
#'
#' @param p Number of variables.
#'
#' @param beta Coefficients.
#'
#' @return A list object containing the simulated data.
#'
#' @author Chong Wu and Wei Pan
#'
#' @references
#' Wu, C., Xu, G., Shen, X., & Pan, W. (2018+). An adaptive test for high-dimensional generalized linear models with application to detect gene-environment interactions, Submitted.
#'
#' @examples
#'
#' # Generate the data (codes for the simulations in the manuscript)
#' n = 30
#' signal.r = 0
#' nInformative = 3
#' p = 40
#' seed = 1
#' s = 0.01
#' non.zero = floor((p/2) * s)
#' alpha = c(rep(0,p/2 - non.zero), runif(non.zero,-signal.r,signal.r))
#' beta = c(rep(2,nInformative), rep(0,(p/2- 3)), alpha)
#' dat = sim_data(seed, n = n, p = p, beta = beta)
#'
#' X = dat$X
#' Y = dat$Y
#'
sim_data <- function(seed,n,p, beta) {
    
    set.seed(seed)
    
    Sigma <- autocorr.mat(p)
    
    X <- matrix(rnorm(n*p), n, p)
    X <- t(t(chol(Sigma))%*%t(X))
    
    X <- X - rep(1, nrow(X)) %*% t(colMeans(X)) #center it

    error <- rnorm(n,0,1)
    error <- error -mean(error)
    Y <- X%*%beta + error
    
    out = list(X = X, Y = Y)
    out
}

