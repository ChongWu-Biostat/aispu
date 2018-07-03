autocorr.mat <- function(p = 1000, rho = 0.3) {
    mat <- diag(p)
    res = rho^abs(row(mat)-col(mat))
    res[abs(row(mat)-col(mat)) >=3] = 0
    return(res)
    
}
