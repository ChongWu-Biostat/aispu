require(mvtnorm)

best.band <- function(sam, bandwidth, cv.fold= 5, norm= "F") {
    p <- dim(sam)[2]
    n <- dim(sam)[1]
    fold.size <- round(n/cv.fold)
    sam.idx <- sample(1:n, size = n, replace = FALSE)
    n.bandwidth <- length(bandwidth)
    diff.norm <- matrix(0, cv.fold, n.bandwidth)
    for(i in 1:cv.fold){
        if(i == cv.fold){
            temp.idx <- sam.idx[((i - 1)*fold.size + 1):n]
        }else{
            temp.idx <- sam.idx[((i - 1)*fold.size + 1):(i*fold.size)]
        }
        sam.train <- sam[-temp.idx,]
        sam.test <- sam[temp.idx,]
        sam.train.cov <- cov(sam.train)
        sam.test.cov <- cov(sam.test)
        
        for(j in 1:n.bandwidth){
            sam.train.cov.band <- sam.train.cov
            sam.train.cov.band[abs(row(sam.train.cov.band) - col(sam.train.cov.band)) > bandwidth[j]] <- 0
            diff.norm[i, j] <- norm(sam.train.cov.band - sam.test.cov, type = norm)
        }
    }
    diff.norm <- colMeans(diff.norm)
    best <- which(diff.norm == min(diff.norm))
    best <- best[1]
    return(bandwidth[best])
}


best.band2 <- function(sam, bandwidth, cv.fold= 5, norm= "F") {
    p <- dim(sam)[2]
    n <- dim(sam)[1]
    fold.size <- round(n/cv.fold)
    sam.idx <- sample(1:n, size = n, replace = FALSE)
    n.bandwidth <- length(bandwidth)
    diff.norm <- matrix(0, cv.fold, n.bandwidth)
    for(i in 1:cv.fold){
        if(i == cv.fold){
            temp.idx <- sam.idx[((i - 1)*fold.size + 1):n]
        }else{
            temp.idx <- sam.idx[((i - 1)*fold.size + 1):(i*fold.size)]
        }
        sam.train <- sam[-temp.idx,]
        sam.test <- sam[temp.idx,]
        sam.train.cov <- t(sam.train) %*% sam.train / (dim(sam.train)[1]-1)
        
        sam.test.cov <- t(sam.test) %*% sam.test / (dim(sam.test)[1]-1)
        
        for(j in 1:n.bandwidth){
            sam.train.cov.band <- sam.train.cov
            sam.train.cov.band[abs(row(sam.train.cov.band) - col(sam.train.cov.band)) > bandwidth[j]] <- 0
            diff.norm[i, j] <- norm(sam.train.cov.band - sam.test.cov, type = norm)
        }
    }
    diff.norm <- colMeans(diff.norm)
    best <- which(diff.norm == min(diff.norm))
    best <- best[1]
    return(bandwidth[best])
}


best.band3 <- function(sam, cov.new, bandwidth, cv.fold= 5, norm= "F"){
    p <- dim(sam)[2]
    n <- dim(sam)[1]
    fold.size <- round(n/cv.fold)
    sam.idx <- sample(1:n, size = n, replace = FALSE)
    n.bandwidth <- length(bandwidth)
    diff.norm <- matrix(0, cv.fold, n.bandwidth)
    for(i in 1:cv.fold){
        if(i == cv.fold){
            temp.idx <- sam.idx[((i - 1)*fold.size + 1):n]
        }else{
            temp.idx <- sam.idx[((i - 1)*fold.size + 1):(i*fold.size)]
        }
        sam.train <- sam[-temp.idx,]
        sam.test <- sam[temp.idx,]
        cov.train <- cov.new[-temp.idx,]
        cov.test <- cov.new[temp.idx,]
        
        V11 <- t(cov.train) %*% cov.train / (dim(cov.train)[1]-1)
        V12 <- t(cov.train) %*% sam.train / (dim(cov.train)[1]-1)
        V22 <- t(sam.train)%*% sam.train / (dim(sam.train)[1]-1)
        sam.train.cov <- V22 - t(V12) %*% solve(V11) %*% V12
        
        
        V11 <- t(cov.test) %*% cov.test / (dim(cov.test)[1]-1)
        V12 <- t(cov.test) %*% sam.test / (dim(cov.test)[1]-1)
        V22 <- t(sam.test)%*% sam.test / (dim(sam.test)[1]-1)
        sam.test.cov <- V22 - t(V12) %*% solve(V11) %*% V12
        
        for(j in 1:n.bandwidth){
            sam.train.cov.band <- sam.train.cov
            sam.train.cov.band[abs(row(sam.train.cov.band) - col(sam.train.cov.band)) > bandwidth[j]] <- 0
            diff.norm[i, j] <- norm(sam.train.cov.band - sam.test.cov, type = norm)
        }
    }
    diff.norm <- colMeans(diff.norm)
    best <- which(diff.norm == min(diff.norm))
    best <- best[1]
    return(bandwidth[best])
}


SPU_E <- function(gamma, n, cov) {
    if(gamma %% 2 == 1){
        e <- 0
    }
    if(gamma %% 2 == 0){
        ga.half <- gamma/2
        e <- (factorial(gamma)*sum(diag(cov)^ga.half))/(2^ga.half * factorial(ga.half) * n^{ga.half})
    }
    return(e)
}


SPU_E2 <- function(gamma, n, cov) {
    out = 0
    if(gamma %% 2 == 0) {
        ga.half <- gamma/2
        e <- (factorial(gamma)*diag(cov)^ga.half)/(2^ga.half * factorial(ga.half) * n^{ga.half})
        out = e%*%e
    }
    return(out)
}

SPU_E3 <- function(s,t, n, cov) {
    out = 0
    p <- dim(cov)[1]
    gamma = s
    ga.half <- gamma/2
    e1 <- (factorial(gamma)*diag(cov)^ga.half)/(2^ga.half*factorial(ga.half) * n^{ga.half})
    gamma = t
    ga.half <- gamma/2
    e2 <- (factorial(gamma)*diag(cov)^ga.half)/(2^ga.half*factorial(ga.half) * n^{ga.half})
    out = e1 %*% e2
    
    return(out)
}


c_space <- function(s, t) {
    c1.vec <- c2.vec <- c3.vec <- NULL
    for(c3 in 0:min(c(t, s))){
        for(c1 in 0:floor((t - c3)/2)) {
            for(c2 in 0:floor((s - c3)/2)) {
                if( t - 2*c1 - c3 == 0 & s - 2*c2 - c3  == 0){
                    c1.vec <- c(c1.vec, c1)
                    c2.vec <- c(c2.vec, c2)
                    c3.vec <- c(c3.vec, c3)
                }
            }
        }
    }
    return(data.frame(c1 = c1.vec, c2 = c2.vec, c3 = c3.vec))
}



SPU_Var <- function(gamma,n, cov) {
    p <- dim(cov)[1]
    if(gamma == 1){
        var <- (1/n)*(t(rep(1, p)) %*% cov %*% rep(1, p))
        var <- as.numeric(var)
    } else if(gamma %%2 ==0) {
        P1 <- SPU_E(2*gamma, n, cov)
        P2 <- -SPU_E2(gamma, n, cov)
        C <- c_space(gamma, gamma)
        C = C[C[,3] >0,]
        n.case <- dim(C)[1]
        P3 <- 0
        
        diags <- diag(cov)
        p <- length(diags)
        mat1 <- matrix(rep(diags, p), p, p, byrow = FALSE)
        mat2 <- matrix(rep(diags, p), p, p, byrow = TRUE)
        for(i in 1:n.case){
            c1 <- C$c1[i]
            c2 <- C$c2[i]
            c3 <- C$c3[i]
            
            mat <- mat1^(c1)*mat2^(c2)*cov^(c3)
            diag(mat) <- 0
            N <- (factorial(gamma))^2*sum(mat)
            D <- (n^gamma * factorial(c1)*factorial(c2)*factorial(c3)*2^(c1 + c2))
            P3 <- P3 + N/D
        }
        var <- P1 + P2 + P3
    } else {
        P1 = SPU_E(2* gamma, n, cov)
        C <- c_space(gamma, gamma)
        n.case <- dim(C)[1]
        P3 <- 0
        
        diags <- diag(cov)
        p <- length(diags)
        mat1 <- matrix(rep(diags, p), p, p, byrow = FALSE)
        mat2 <- matrix(rep(diags, p), p, p, byrow = TRUE)
        for(i in 1:n.case){
            c1 <- C$c1[i]
            c2 <- C$c2[i]
            c3 <- C$c3[i]
            
            mat <- mat1^(c1)*mat2^(c2)*cov^(c3)
            diag(mat) <- 0
            N <- (factorial(gamma))^2*sum(mat)
            D <- (n^gamma * factorial(c1)*factorial(c2)*factorial(c3)*2^(c1 + c2))
            P3 <- P3 + N/D
        }
        var <- P1 + P3
    }
    return(var)
}

SPU_Cov <- function(s,t,n, cov) {
    p <- dim(cov)[1]
    var = 0
    if(s %% 2 ==0) {
        P1 <- SPU_E(s+t, n, cov)
        P2 <- -SPU_E3(s,t, n, cov)
        C <- c_space(s, t)
        C = C[C[,3] >0,]
        n.case <- dim(C)[1]
        P3 <- 0
        
        diags <- diag(cov)
        p <- length(diags)
        mat1 <- matrix(rep(diags, p), p, p, byrow = FALSE)
        mat2 <- matrix(rep(diags, p), p, p, byrow = TRUE)
        for(i in 1:n.case){
            c1 <- C$c1[i]
            c2 <- C$c2[i]
            c3 <- C$c3[i]
            
            mat <- mat1^(c1)*mat2^(c2)*cov^(c3)
            diag(mat) <- 0
            N <- factorial(s) * factorial(t) * sum(mat)
            D <- (n^((s+t)/2) * factorial(c1)*factorial(c2)*factorial(c3)*2^(c1 + c2))
            P3 <- P3 + N/D
        }
        var <- P1 + P2 + P3
    } else {
        P1 <- SPU_E(s+t, n, cov)
        C <- c_space(s, t)
        n.case <- dim(C)[1]
        P3 <- 0
        
        diags <- diag(cov)
        p <- length(diags)
        mat1 <- matrix(rep(diags, p), p, p, byrow = FALSE)
        mat2 <- matrix(rep(diags, p), p, p, byrow = TRUE)
        for(i in 1:n.case){
            c1 <- C$c1[i]
            c2 <- C$c2[i]
            c3 <- C$c3[i]
            
            mat <- mat1^(c1)*mat2^(c2)*cov^(c3)
            diag(mat) <- 0
            N <- factorial(s) * factorial(t) * sum(mat)
            D <- (n^((s+t)/2) * factorial(c1)*factorial(c2)*factorial(c3)*2^(c1 + c2))
            P3 <- P3 + N/D
        }
        var <- P1 + P3
    }
    return(var)
}

apvalband_aSPU <- function(Y, X,cov = NULL,SNP, pow = c(1:6, Inf),bandwidth, cv.fold = 5, model= "gaussian", cov.est.type = "Method1" ,norm = "F", penalty = c("tlp","lasso","ridge","net","mcp","SCAD"), tau = 0.05, standardize = FALSE, dfmax = 1000, pmax = 1000, nlambda = 10 ){
    
    tuning.time = proc.time()[3]
    
    n = dim(X)[1]
    p = dim(X)[2]
    
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
        r = Y - mean(Y)
        U = t(X) %*% r
        U = U / n
        X.new = sweep(X,MARGIN=1,r,`*`)
    } else {
        if(penalty =="tlp") {
            tau1 = tau[1]
            cv = cv.glmTLP(cov,Y,family = model, tau = tau1, penalty.factor = penalty.factor,  standardize = standardize, dfmax = dfmax,pmax = pmax, nlambda = nlambda)
            cv = cbind(cv$cvm,cv$lambda, cv$tau)
            colnames(cv) = c("cv","lambda","tau")
            
            if(length(tau) > 1) {
                for(i in 2:length(tau)) {
                    tau1 = tau[i]
                    cv.tmp = cv.glmTLP(cov,Y,family = model, tau = tau1, penalty.factor = penalty.factor, standardize = standardize, dfmax = dfmax,pmax = pmax,nlambda = nlambda)
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
            
        }
        
        r = yresids
        
        X.new = sweep(X,MARGIN=1,r,`*`)
        U = t(X) %*% r
        U = U / n
    }
    
    tuning.time = proc.time()[3] - tuning.time
    other.time = proc.time()[3]
    
    if(missing(bandwidth)) bandwidth <- seq(from = 0, to = p, by = floor(p/50))
    if(any(bandwidth < 0)){
        cat("Negative values specified in bandwidth are removed.\n")
        bandwidth <- bandwidth[bandwidth < 0]
    }
    if(any(bandwidth != floor(bandwidth))){
        cat("Non-integers specified in bandwidth are converted to their integer parts.")
        bandwidth <- floor(bandwidth)
    }
    
    if(cov.est.type == "Method1") {
        
        sam.cov = cov(X.new)
        output.opt.bw <- TRUE
        if(length(bandwidth) > 1){
            optim.bandwidth <- best.band(X.new,bandwidth, cv.fold, norm)
        }
        if(length(bandwidth) == 1){
            optim.bandwidth <- bandwidth
        }
        if(optim.bandwidth > 0){
            cov.est <- sam.cov
            cov.est[abs(row(cov.est) - col(cov.est)) > optim.bandwidth] <- 0
        }
        if(optim.bandwidth == 0){
            cov.est <- diag(diag(sam.cov))
        }
        
    } else if(cov.est.type == "Method2") {
        
        sam.cov  = t(X.new)%*% X.new /(dim(X.new)[1]-1)
        output.opt.bw <- TRUE
        if(length(bandwidth) > 1){
            optim.bandwidth <- best.band2(X.new,bandwidth, cv.fold, norm)
        }
        if(length(bandwidth) == 1){
            optim.bandwidth <- bandwidth
        }
        if(optim.bandwidth > 0){
            cov.est <- sam.cov
            cov.est[abs(row(cov.est) - col(cov.est)) > optim.bandwidth] <- 0
        }
        if(optim.bandwidth == 0){
            cov.est <- diag(diag(sam.cov))
        }
    } else {
        
        cov.new = sweep(cov,MARGIN=1,r,`*`)
        V11 = t(cov.new) %*% cov.new / (dim(cov.new)[1]-1)
        V12 = t(cov.new) %*% X.new / (dim(cov.new)[1]-1)
        V22 = t(X.new)%*% X.new / (dim(X.new)[1]-1)
        sam.cov  = V22 - t(V12) %*% solve(V11) %*% V12
        
        output.opt.bw <- TRUE
        
        if(length(bandwidth) > 1){
            optim.bandwidth <- best.band3(X.new, cov.new, bandwidth, cv.fold, norm)
        }
        if(length(bandwidth) == 1){
            optim.bandwidth <- bandwidth
        }
        if(optim.bandwidth > 0){
            cov.est <- sam.cov
            cov.est[abs(row(cov.est) - col(cov.est)) > optim.bandwidth] <- 0
        }
        if(optim.bandwidth == 0) {
            cov.est <- diag(diag(sam.cov))
        }
    }
    
    ##observed statistics
    pval <- numeric(length(pow) + 1)
    L <- numeric(length(pow))
    L.e <- numeric(length(pow) - 1)
    L.var <- numeric(length(pow) - 1)
    stan.L <- numeric(length(pow) - 1)
    
    for(i in 1:length(pow)){
        if(pow[i] != Inf){
            L[i] <- sum(U^(pow[i]))
            L.e[i] <- SPU_E(pow[i], n, cov.est)
            L.var[i] <- SPU_Var(pow[i], n, cov.est)
            stan.L[i] <- (L[i] - L.e[i]) / sqrt(L.var[i])
            if(pow[i] %% 2 == 1) pval[i] <- 2 * (1 - pnorm(abs(stan.L[i])))
            if(pow[i] %% 2 == 0) pval[i] <- 1 - pnorm(stan.L[i])
        }else{
            diag.sam.cov <- diag(sam.cov)
            diag.sam.cov[diag.sam.cov <= 10^(-10)] <- 10^(-10)
            L[i] <- n * max(U^2/diag.sam.cov)
            stan.L.inf <- L[i] - (2*log(p) - log(log(p)))
            pval[i] <- pval.inf <- 1 - exp(-exp(-stan.L.inf/2)/sqrt(pi))
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
    diag(R_O) <- 1
    diag(R_E) <- 1
    for(s in odd.ga){
        for(t in odd.ga){
            if(s < t){
                L.cov <- SPU_Cov(s, t, n, cov.est)
                io <- which(odd.ga == s)
                jo <- which(odd.ga == t)
                i <- which(pow == s)
                j <- which(pow == t)
                R_O[io, jo] <- L.cov/sqrt(L.var[i]*L.var[j])
                R_O[jo, io] <- R_O[io, jo]
            }
        }
    }
    
    for(s in even.ga){
        for(t in even.ga){
            if(s < t){
                L.cov <- SPU_Cov(s, t, n, cov.est)
                ie <- which(even.ga == s)
                je <- which(even.ga == t)
                i <- which(pow == s)
                j <- which(pow == t)
                R_E[ie, je] <- L.cov/sqrt(L.var[i]*L.var[j])
                R_E[je, ie] <- R_E[ie, je]
            }
        }
    }
    TO <- max(abs(stan.L[odd.ga.id]))
    TE <- max(stan.L[even.ga.id])
    pval_O <- 1 - pmvnorm(lower = -rep(TO, n.odd.ga), upper = rep(TO, n.odd.ga), mean = rep(0, n.odd.ga), sigma = R_O)
    pval_E <- 1 - pmvnorm(lower = rep(-Inf, n.even.ga), upper = rep(TE, n.even.ga), mean = rep(0, n.even.ga), sigma = R_E)
    pval.min <- min(c(pval_O, pval_E, pval.inf))
    pval[length(pow) + 1] <- 1 - (1 - pval.min)^3
    names(pval) <- c(paste("SPU_", pow, sep = ""), "aSPU")
    
    other.time = proc.time()[3] - other.time
    
    time = c(tuning.time,other.time)
    
    out = list(pvs = pval, Ts = L,optim.bandwidth = optim.bandwidth,runtime = time)
    return(out)
}

