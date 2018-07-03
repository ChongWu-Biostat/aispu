#' Three Step Procedure
#'
#' It returns p-values of the NST and ST. Can only apply to linear model. Not applicable for binary outcome.
#'
#' @param Y.f Response or phenotype data. It can be either a binary trait or a quantitative trait. A vector with length n (number of subjects).
#'
#' @param X.f Variables of interest; each row for a subject, and each column
#'     for predictor).A matrix with dimension n by p.
#'
#' @param sub.size The size for the first pruning step. The authors recommend set sub.size = 0.3 * n.
#'
#' @param test.set The indicator for the variable of interest.
#'
#' @param M the number of bootstrap replications
#'
#' @return The non-studentized and studentized statistics and corresponding p-values.
#'
#' @author Zhang, X., & Cheng, G.
#'
#' @references
#' Zhang, X., & Cheng, G. (2017). Simultaneous inference for high-dimensional linear models. Journal of the American Statistical Association, 112(518), 757-768.
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
#' sub.size <- n*0.3
#' test.set = (p/2 + 1):p
#'
#' X = dat$X
#' Y = dat$Y
#'
#' cov = NULL
#" nperm = 100
#' #ST(X, Y, sub.size, test.set,M = nperm) #not run, time-consuming
ST <- function(X.f, Y.f, sub.size, test.set, M=500)
{
n <- dim(X.f)[1]
p <- dim(X.f)[2]

n1 <- sub.size
n0 <- n-floor(n1)
S1 <- sample(1:n, floor(n1), replace=FALSE)
X.sub <- X.f[S1,]
Y.sub <- Y.f[S1]
cvfit <- cv.glmnet(X.sub, Y.sub, intercept=FALSE)
cf <- as.numeric(coef(cvfit, s="lambda.min"))[-1]
set1 <- (1:p)[abs(cf)>0]
resi <- Y.sub-X.sub%*%cf
beta.m <- t(scale(X.sub[,-set1])) %*% resi #standardize
screen.set <- sort(order(abs(beta.m),decreasing=TRUE)[1:(n0-1-length(set1))])
a <- (1:p)[-set1]
screen.set <- union(a[screen.set],set1)
X <- X.f[-S1,screen.set]
Y <- Y.f[-S1]

node <- score.nodewiselasso(X, wantTheta=TRUE, verbose=FALSE, lambdaseq="quantile",
parallel=FALSE, ncores=2, oldschool = FALSE, lambdatuningfactor = 1)
Theta <- node$out
Gram<-t(X)%*%X/n0

sreg <- scalreg(X,Y)
beta.hat <- sreg$coefficients
sigma.sq <- sum((Y-X%*%beta.hat)^2)/(n0-sum(abs(beta.hat)>0))
test.set.i <- intersect(screen.set,test.set)
index <- screen.set%in%test.set.i

Omega <- diag(Theta%*%Gram%*%t(Theta))*sigma.sq
beta.db <- beta.hat+Theta%*%t(X)%*%(Y-X%*%beta.hat)/n0
margin.st <- sqrt(n0)*abs(beta.db[index])/sqrt(Omega[index])
margin.nst <- sqrt(n0)*abs(beta.db[index])
stat.st <- max(margin.st)
stat.nst <- max(margin.nst)

stat.boot.st <- stat.boot.nst <- rep(NA,M)
for(i in 1:M)
 {
 e <- rnorm(n0)
 xi.boot <- Theta[index,]%*%t(X)%*%e*sqrt(sigma.sq)/sqrt(n0)
 stat.boot.nst[i] <- max(abs(xi.boot))
 stat.boot.st[i] <- max(abs(xi.boot/sqrt(Omega[index])))
 }

p.nst = sum(stat.boot.nst >= stat.nst) / M
p.st = sum(stat.boot.st >= stat.st) / M

result <- list(stat.nst, p.nst, stat.st, p.st)
names(result) <- c("non-studentized test","p-non-studentized test","studentized test","p-studentized test")

return(result)
}

# An example



