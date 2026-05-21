#' README
#' ------
#' - author: Liang-Cheng Chen
#' - created_at: 2026-04-03
#' - desc:
#'   - implement logdet and logmarglik functions to calculate log-likelihood

# setwd(".")

# Problem 1 -----------------------------------------------------------------------------------

logdet <- function(R) {
    #' `logdet` function calculates the log of the determinant of the matrix.
    #'
    #'  Args:
    #'  ----
    #'  - R (matrix): the target matrix
    #'
    #'  Returns:
    #'  --------
    #'  - output (num): the log of the det of R

    eigs <- eigen(R)$values

    if (any(eigs <= 0)) {
        # handle non-positive eigenvalues case
        stop("matrix R is not positive definite")
        return(NA_real_)
    }
    res <- sum(log(eigs))

    return(res)
}


# Problem 2 -----------------------------------------------------------------------------------

logmarglik <- function(data, A) {
    #' `logmarglik` calculates the marginal log likelihood given the data and features.
    #'
    #' Args:
    #' -----
    #' - data (matrix | data.frame-like): (n x p) matrix-like data
    #' - A (vector): vector selecting target columns of data
    #'
    #' Return:
    #' -------
    #' - logmarglik (num): the log likelihood value

    m <- as.matrix(data) # NOTE: to accomodate both data.table and data.frame
    n <- nrow(m)
    k <- length(A)


    D1 <- m[, 1]
    DA <- m[, A] # NOTE: data.table will fail with this syntax

    MA <- diag(k) + t(DA) %*% DA

    temp <- (1 + t(D1) %*% D1 - t(D1) %*% DA %*% solve(MA) %*% t(DA) %*% D1)
    loglik <- lgamma((n+k+2)/2) - lgamma((k + 2)/2) - 0.5*logdet(MA) - ((n+k+2)/2) * log(temp)

    return(as.numeric(loglik))
}




# Testing Area --------------------------------------------------------------------------------

# # load data
# df <- data.table::fread("hw1/erdata.txt")
# 
# # log-likelihood calculation
# logmarglik(df, A = c(2, 5, 10))
# logmarglik(as.data.frame(df), A = c(2, 5, 10))
# logmarglik(as.matrix(df), A = c(2, 5, 10))
