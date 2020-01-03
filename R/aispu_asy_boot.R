#require(mvtnorm)

aispu_asy <- function(Y, X,cov = NULL, SNP, pow = c(1:6, Inf), model= c("gaussian", "binomial"),n.perm = 1000, penalty = c("tlp","lasso","ridge","net","mcp","SCAD"), tau = 0.1, standardize = FALSE, dfmax = 1000, pmax = 1000){
    
    model = match.arg(model)
    
    n <- dim(X)[1]
    p <- dim(X)[2]
    
    if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)
    k <- ncol(X)
    
    if(is.null(cov)) {
        penalty.factor = rep(1,dim(SNP)[2])
        cov = SNP
    } else {
        penalty.factor = c(rep(0,dim(cov)[2]),rep(1,dim(SNP)[2]))
        cov = cbind(cov, SNP)
    }
    
    if(is.null(cov)) {
        r <- Y - mean(Y)
        U <- t(X) %*% r
        U <- U / n
        X.new = sweep(X,MARGIN=1,r,`*`)
    } else {
        if(penalty =="tlp") {
            tau1 = tau[1]
            cv = cv.glmTLP(cov,Y,family = model, tau = tau1, penalty.factor = penalty.factor,  standardize = standardize, dfmax = dfmax,pmax = pmax)
            cv = cbind(cv$cvm,cv$lambda, cv$tau)
            colnames(cv) = c("cv","lambda","tau")
            
            if(length(tau) > 1) {
                for(i in 2:length(tau)) {
                    tau1 = tau[i]
                    cv.tmp = cv.glmTLP(cov,Y,family = model, tau = tau1, penalty.factor = penalty.factor, standardize = standardize, dfmax = dfmax,pmax = pmax)
                    cv.tmp = cbind(cv.tmp$cvm,cv.tmp$lambda, cv.tmp$tau)
                    colnames(cv.tmp) = c("cv","lambda","tau")
                    cv = rbind(cv,cv.tmp)
                }
            }
            
            cv = cv[cv[,1] ==min(cv[,1]),]
            if(class(cv) == "matrix") {
                cv = cv[cv[,3] ==min(cv[,3]),]
            }
            if(class(cv) == "matrix")  {
                cv = cv[ceiling(dim(cv)[1]/2),]
            }
            
            lambda.tmp = cv[2]
            tau.tmp = cv[3]
            
            fit1 = glmTLP(cov,Y,family = model,lambda = lambda.tmp, penalty.factor = penalty.factor,  tau = tau.tmp, standardize = standardize, dfmax = dfmax,pmax = pmax)
            
            yfits <- predict(fit1,cov,type = "response")
            yresids <- Y - yfits
            
        } else if (penalty == "lasso") {
            cvfit <- cv.glmnet(cov, Y,family = model, penalty.factor = penalty.factor,  standardize = standardize, dfmax = dfmax,pmax = pmax,alpha = 1)
            lambda.tmp = cvfit$lambda.min
            yfits <- predict(cvfit, newx = cov, s = "lambda.min",type = "response")
            yresids <- Y - yfits
        } else if (penalty == "ridge") {
            cvfit <- cv.glmnet(cov, Y,family = model, penalty.factor = penalty.factor,  standardize = standardize, dfmax = dfmax,pmax = pmax,alpha = 0)
            lambda.tmp = cvfit$lambda.min
            yfits <- predict(cvfit, newx = cov, s = "lambda.min",type = "response")
            yresids <- Y - yfits
        } else if (penalty == "net") {
            cvfit <- cv.glmnet(cov, Y,family = model,penalty.factor = penalty.factor,  standardize = standardize, dfmax = dfmax,pmax = pmax,alpha = 0.5)
            lambda.tmp = cvfit$lambda.min
            yfits <- predict(cvfit, newx = cov, s = "lambda.min",type = "response")
            yresids <- Y - yfits
        } else if (penalty == "mcp") {
            cvfit = cv.ncvreg(cov, Y, family=model,penalty.factor = penalty.factor,  penalty="MCP")
            lambda.tmp =cvfit$lambda.min
            fit = ncvreg(cov,Y, family = model,penalty.factor = penalty.factor,  penalty = "MCP",lambda = lambda.tmp)
            yfits = predict(fit,cov,type = "response")
            yresids <- Y - yfits
        } else if (penalty == "scad") {
            cvfit = cv.ncvreg(cov, Y, family=model,penalty.factor = penalty.factor, penalty="SCAD")
            lambda.tmp =cvfit$lambda.min
            fit = ncvreg(cov,Y, family = model,penalty.factor = penalty.factor, penalty = "SCAD",lambda = lambda.tmp)
            yfits = predict(fit,cov,type = "response")
            yresids <- Y - yfits
        }
        
        X.new <- sweep(X,MARGIN=1,yresids,`*`)
        U <- t(X) %*% yresids
        U <- U / n
        
        sam.cov = cov(X.new)
        
        diag.sam.cov <- diag(sam.cov)
        diag.sam.cov[diag.sam.cov <= 10^(-10)] <- 10^(-10)
        
        # test stat's:
        Ts <- max(U^2/diag.sam.cov)
    }
    
    ## calculate the expecatation, variance, and covariance of L(gamma)
    parametric.boot <- matrix(NA, n.perm, dim(U)[1])
    T0s = matrix(0, nrow=n.perm, ncol=1)
    
    Y0 = Y
    
    for (b in 1:n.perm){
        set.seed(b)
        if (is.null(cov)) {
            Y0 <- sample(Y, length(Y))
            ##  Null score vector:
            parametric.boot[b,]<- t(X) %*% (Y0-mean(Y0))
        } else {
            
            ## with nuisance parameters:
            if ( model == "gaussian") {
                Y0 <- yfits + sample(yresids, n, replace = F )
            } else {
                for(i in 1:n) Y0[i] <- sample(c(1,0), 1, prob=c(yfits[i], 1-yfits[i]) )
            }
            
            if (penalty == "tlp") {
                fit1 = glmTLP(cov,Y0,family = model,lambda = lambda.tmp, tau = tau.tmp, penalty.factor = penalty.factor, standardize = standardize, dfmax = dfmax,pmax = pmax)
                yfits0 <- predict(fit1,cov,type = "response")
            } else if (penalty == "lasso") {
                fit <- glmnet(cov, Y0,family = model,lambda = lambda.tmp,penalty.factor = penalty.factor, standardize = standardize, dfmax = dfmax,pmax = pmax,alpha = 1)
                yfits0 <- predict(fit, newx = cov,type = "response")
            } else if (penalty == "ridge") {
                fit <- glmnet(cov, Y0,family = model,lambda = lambda.tmp,penalty.factor = penalty.factor, standardize = standardize, dfmax = dfmax,pmax = pmax,alpha = 0)
                yfits0 <- predict(fit, newx = cov,type = "response")
            } else if (penalty == "net") {
                fit <- glmnet(cov, Y0,family = model,lambda = lambda.tmp, penalty.factor = penalty.factor, standardize = standardize, dfmax = dfmax,pmax = pmax,alpha = 0.5)
                yfits0 <- predict(fit, newx = cov,type = "response")
            } else if (penalty == "mcp") {
                fit = ncvreg(cov,Y0, family = model,penalty.factor = penalty.factor,  penalty = "MCP",lambda = lambda.tmp)
                yfits0 = predict(fit,cov,type = "response")
            } else if (penalty == "scad") {
                fit = ncvreg(cov,Y0, family = model,penalty.factor = penalty.factor, penalty = "SCAD",lambda = lambda.tmp)
                yfits0 = predict(fit,cov,type = "response")
            }
            U0<-t(X) %*% (Y0 - yfits0)
            U0 <- U0/n
            T0s[b, 1] = max(U0^2/diag.sam.cov)
            
            parametric.boot[b,]<- U0
            
        }
    }
    
    ##observed statistics
    pval <- numeric(length(pow) + 1)
    L <- numeric(length(pow))
    if(sum(pow==Inf)==0) {
        L.e <- numeric(length(pow))
        L.var <- numeric(length(pow))
        
        stan.L <- numeric(length(pow))
        
        boot.test.stat <- matrix(NA,n.perm,length(pow))
        
        for (b in 1:length(pow)) {
            boot.test.stat[,b] <-rowSums(parametric.boot^{pow[b]})
        }
    } else {
        L.e <- numeric(length(pow) - 1)
        L.var <- numeric(length(pow) - 1)
        
        stan.L <- numeric(length(pow) - 1)
        
        boot.test.stat <- matrix(NA,n.perm,length(pow)-1)
        
        for (b in 1:(length(pow)-1)) {
            boot.test.stat[,b] <-rowSums(parametric.boot^{pow[b]})
        }
    }
    
    L.e <- colMeans(boot.test.stat)
    boot.var <- var(boot.test.stat)
    
    L.var <- diag(boot.var)
    boot.cor <- cor(boot.test.stat)
    
    ########################################
    
    for(i in 1:length(pow)){
        if(pow[i] != Inf){
            L[i] <- sum(U^(pow[i]))
            stan.L[i] <- (L[i] - L.e[i]) / sqrt(L.var[i])
            if(pow[i] %% 2 == 1) pval[i] <- 2 * (1 - pnorm(abs(stan.L[i])))
            if(pow[i] %% 2 == 0) pval[i] <- 1 - pnorm(stan.L[i])
        } else {
            pval.inf <- round(sum(abs(Ts)<=abs(T0s)) / n.perm, digits = 8)
            pval[i] <- pval.inf
        }
    }
    
    f.pow <- pow[pow != Inf]
    odd.ga <- f.pow[f.pow %% 2 == 1]
    odd.ga.id <- which(f.pow %% 2 == 1)
    even.ga <- f.pow[f.pow %% 2 == 0]
    even.ga.id <- which(f.pow %% 2 == 0)
    n.odd.ga <- length(odd.ga)
    n.even.ga <- length(even.ga)
    R_O <- matrix(NA, n.odd.ga, n.odd.ga)
    R_E <- matrix(NA, n.even.ga, n.even.ga)
    
    R_O <- boot.cor[odd.ga.id,odd.ga.id]
    R_E <- boot.cor[even.ga.id,even.ga.id]
    
    TO <- max(abs(stan.L[odd.ga.id]))
    TE <- max(stan.L[even.ga.id])
    pval_O <- 1 - pmvnorm(lower = -rep(TO, n.odd.ga), upper = rep(TO, n.odd.ga), mean = rep(0, n.odd.ga), sigma = R_O)
    pval_E <- 1 - pmvnorm(lower = rep(-Inf, n.even.ga), upper = rep(TE, n.even.ga), mean = rep(0, n.even.ga), sigma = R_E)
    
    if(sum(pow==Inf)==0) {
        pval.min <- min(c(pval_O, pval_E))
        pval[length(pow) + 1] <- 1 - (1 - pval.min)^2
        names(pval) <- c(paste("iSPU(", pow,")", sep = ""), "aiSPU")
    } else {
        
        pval.min <- min(c(pval_O, pval_E, pval.inf))
        pval[length(pow) + 1] <- 1 - (1 - pval.min)^3
        names(pval) <- c(paste("iSPU(", pow, ")", sep = ""), "aiSPU")
    }
    out = list(pvs = pval, Ts = L)
    return(out)
}
