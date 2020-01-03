#' Adaptive Interaction Sum of Powered Score (aiSPU) Tests
#'
#' It returns p-values of the iSPU tests and aiSPU test.
#'
#' @param Y Response or phenotype data. It can be either a binary trait or a quantitative trait. A vector with length n (number of subjects).
#'
#' @param X Variables of interest; each row for a subject, and each column
#'     for predictor).A matrix with dimension n by p.
#'
#' @param cov Ancillary covariates without penalization, such as age, gender, and etc.
#'
#' @param cov2 High-dimensional covariates with penalization.
#'
#' @param pow Power used in iSPU test. A vector of the powers.
#'
#' @param model Use "gaussian" for a quantitative trait, and use "binomial" for a binary trait.
#'
#' @param n.perm The number of boostrap replications
#'
#' @param penalty Penalty for cov2. The default is truncated lasso (tlp). We recommend using tlp in real data applications.
#'
#' @param tau Tunning parameters for tlp.
#'
#' @param standardize Logical flag for x variable standardization, prior to fitting the model sequence. The coefficients are always returned on the original scale. Default is standardize=TRUE. If variables are in the same units already, you might not wish to standardize.
#'
#' @param dfmax Limit the maximum number of variables in the model. Useful when the dimension of cov2 is high.
#'
#' @param pmax Limit the maximum number of variables ever to be nonzero.
#'
#' @param resample Methods for calculating p-values. The default is the asymptotics-based method ("asy").
#'
#' @return A list object, Ts : test statistics for the iSPU tests (in the order of the specified pow) and finally for the aSPU test.
#'         pvs : p-values for the iSPU and aiSPU tests.
#'
#' @author Chong Wu and Wei Pan
#'
#' @references
#' Wu, C., Xu, G., Shen, X., & Pan, W. (2018+). A Regularization-Based Adaptive Test for High-Dimensional Generalized Linear Models, Submitted.
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
#' cov = NULL
#' X.tmp = X
#' cov2 = X.tmp[,1:(p/2)]
#' X = X.tmp[,(p/2 + 1):p]
#' aispu(Y, X,cov = NULL, cov2, pow = c(1:6, Inf), model= "gaussian",penalty = "tlp", n.perm = 10,resample = "boot")
#' aispu(Y, X,cov = NULL, cov2, pow = c(1:6, Inf), model= "gaussian",penalty = "tlp", n.perm = 10,resample = "asy")
aispu <- function(Y, X,cov = NULL, cov2, pow = c(1:6, Inf), model= c("gaussian", "binomial"),n.perm = 1000, penalty = c("tlp","lasso","ridge","net","mcp","SCAD"), tau = 0.1, standardize = FALSE, dfmax = 1000, pmax = 1000, resample = c("asy","asy-boot","boot"),bandwidth = 3) {
    
    X <- as.matrix(X)
    model <- match.arg(model)
    resample <- match.arg(resample)
    
    if(resample == "boot") {
        res = aispu_boot(Y = Y, X = X, cov = cov, SNP = cov2, pow = pow,model = model, n.perm = n.perm, penalty = penalty, tau = tau, standardize = standardize, dfmax = dfmax, pmax = pmax)
        return(res)
    } else if(resample == "asy-boot"){
        res = aispu_asy(Y = Y, X = X, cov = cov, SNP = cov2, pow = pow,model = model, n.perm = n.perm, penalty = penalty, tau = tau, standardize = standardize, dfmax = dfmax, pmax = pmax)
        return(res)
        
    } else {
        res = apvalband_aSPU(Y = Y, X = X, cov = cov, SNP = cov2, pow = pow,bandwidth = bandwidth, model = model,  penalty = penalty, tau = tau, standardize = standardize, dfmax = dfmax, pmax = pmax)
        return(res)
        
    }
}

