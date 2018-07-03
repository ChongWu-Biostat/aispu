aispu_boot <- function(Y, X, cov=NULL, SNP, pow=c(1:6, Inf), model=c("gaussian", "binomial"), n.perm=1000, penalty = c("tlp","lasso","ridge","net","mcp","SCAD"),  tau = 0.1, standardize = FALSE, dfmax = 1000, pmax = 1000){
    
    model = match.arg(model)
    
    n <- length(Y)
    if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)
    k <- ncol(X)
    
    if(is.null(cov)) {
        penalty.factor = rep(1,dim(SNP)[2])
        cov = SNP
    } else {
        penalty.factor = c(rep(0,dim(cov)[2]),rep(1,dim(SNP)[2]))
        cov = cbind(cov, SNP)
    }
    
    if (is.null(cov)){
        ## NO nuisance parameters:
        Xg <- XUs <- X
        U <- t(Xg) %*% (Y-mean(Y))
        yresids <- Y-mean(Y)
        yfits <- rep(mean(Y), n)
    } else {
        ## with nuisance parameters:
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
            
            
            coef.save <- c(fit1$a0,as.matrix(fit1$beta))

            yfits <- predict(fit1,cov,type = "response")
            yresids <- Y - yfits
            
        } else if (penalty == "lasso") {
            cvfit <- cv.glmnet(cov, Y,family = model, penalty.factor = penalty.factor,  standardize = standardize, dfmax = dfmax,pmax = pmax,alpha = 1)
            lambda.tmp = cvfit$lambda.min
            yfits <- predict(cvfit, newx = cov, s = "lambda.min",type = "response")
            yresids <- Y - yfits
            
            coef.tmp <- predict(cvfit,s = "lambda.min",type = "coefficients")
            coef.save <- as.matrix(coef.tmp)
            
        } else if (penalty == "ridge") {
            cvfit <- cv.glmnet(cov, Y,family = model, penalty.factor = penalty.factor,  standardize = standardize, dfmax = dfmax,pmax = pmax,alpha = 0)
            lambda.tmp = cvfit$lambda.min
            yfits <- predict(cvfit, newx = cov, s = "lambda.min",type = "response")
            yresids <- Y - yfits
            
            coef.tmp <- predict(cvfit,s = "lambda.min",type = "coefficients")
            coef.save <- as.matrix(coef.tmp)
            
        } else if (penalty == "net") {
            cvfit <- cv.glmnet(cov, Y,family = model,penalty.factor = penalty.factor,  standardize = standardize, dfmax = dfmax,pmax = pmax,alpha = 0.5)
            lambda.tmp = cvfit$lambda.min
            yfits <- predict(cvfit, newx = cov, s = "lambda.min",type = "response")
            yresids <- Y - yfits
            
            coef.tmp <- predict(cvfit,s = "lambda.min",type = "coefficients")
            coef.save <- as.matrix(coef.tmp)
            
        } else if (penalty == "mcp") {
            cvfit = cv.ncvreg(cov, Y, family=model,penalty.factor = penalty.factor,  penalty="MCP")
            lambda.tmp =cvfit$lambda.min
            fit = ncvreg(cov,Y, family = model,penalty.factor = penalty.factor,  penalty = "MCP",lambda = lambda.tmp)
            yfits = predict(fit,cov,type = "response")
            yresids <- Y - yfits
            
            coef.save <- as.matrix(fit$beta)

        } else if (penalty == "scad") {
            cvfit = cv.ncvreg(cov, Y, family=model,penalty.factor = penalty.factor, penalty="SCAD")
            lambda.tmp =cvfit$lambda.min
            fit = ncvreg(cov,Y, family = model,penalty.factor = penalty.factor, penalty = "SCAD",lambda = lambda.tmp)
            yfits = predict(fit,cov,type = "response")
            yresids <- Y - yfits
            
            coef.save <- as.matrix(fit$beta)
        }
        
        Us <- XUs <- X
        U <- t(XUs) %*% (Y - yfits)
        
    }
    
    X.new = sweep(X,MARGIN=1,yresids,`*`)
    sam.cov = cov(X.new)
    
    diag.sam.cov <- diag(sam.cov)
    diag.sam.cov[diag.sam.cov <= 10^(-10)] <- 10^(-10)


    # test stat's:
    Ts <- rep(0, length(pow))
    for(j in 1:length(pow)){
        if (pow[j] < Inf)
        Ts[j] = sum(U^pow[j]) else Ts[j] = max(U^2/diag.sam.cov)
    }
    
    # bootstrap:
    T0s = matrix(0, nrow=n.perm, ncol=length(pow))
    Y0 = Y
    for(b in 1:n.perm){
        if (is.null(cov)) {
            Y0 <- sample(Y, length(Y))
            #########Null score vector:
            U0 <- t(Xg) %*% (Y0-mean(Y0))
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
            U0<-t(XUs) %*% (Y0 - yfits0)
        }
        
        # test stat's:
        for(j in 1:length(pow))
        if (pow[j] < Inf)
        T0s[b, j] = sum(U0^pow[j]) else T0s[b, j] = max(U0^2/diag.sam.cov)
        
    }
    
    # bootstrap-based p-values:
    #pPerm0 <- apply( matrix( rep(abs(Ts),n.perm), nrow = n.perm, byrow = T) < abs(T0s), 2, mean)
    
    pPerm0 = rep(NA,length(pow))
    for ( j in 1:length(pow))
    {
        pPerm0[j] = round( sum(abs(Ts[j])<=abs(T0s[,j])) / n.perm, digits = 8)
        P0s = ( ( n.perm - rank( abs(T0s[,j]) ) ) + 1 ) / (n.perm)
        if (j == 1 ) minp0  = P0s else minp0[which(minp0>P0s)] = P0s[which(minp0>P0s)]
    }
    
    Paspu <- (sum(minp0 <= min(pPerm0)) + 1) / (n.perm+1)
    pvs <- c(pPerm0, Paspu)
    
    Ts <- c(Ts, min(pPerm0))
    names(Ts) <- c(paste("iSPU(", pow,")", sep=""), "aiSPU")
    names(pvs) = names(Ts)
    
    list(Ts = Ts, pvs = pvs,coef = coef.save)
}



