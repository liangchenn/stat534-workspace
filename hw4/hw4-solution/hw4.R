#' README
#' ------
#' - author: Liang-Cheng Chen
#' - created_at: 2026-05-11


# Packages ------------------------------------------------------------------------------------
library(data.table)
library(MASS) # for multivariate normal
library(snow) # parallel


# Problem 1 -----------------------------------------------------------------------------------

#' Goal: function computing the Lapalace approximation of the marginal likelihood
#' Setup: univariate logistic regression model


getLaplaceApprox <- function(response, explanatory, data, betaMode) {
    
    # data
    xy <- get_xy(response, explanatory, data)
    y <- xy$y
    x <- xy$x
    
    # hessian
    hess_mode <- get_post_hessian(x=x, y=y, betas=betaMode)
    neg_hess <- -hess_mode
    
    # log deteriminant
    log_det <- as.numeric(
        determinant(neg_hess, logarithm = TRUE)$modulus
    )
    
    log_marginal_lik <- log(2 * pi) + post_loglik(x=x, y=y, betas=betaMode) - 0.5 * log_det
    
    return(log_marginal_lik)
}



# Problem 2 -----------------------------------------------------------------------------------

#' Goal: estimate coefficients betas with Metropolis-Hasting Sampling

getPosteriorMeans <- function(response, explanatory, data, betaMode, niter) {
    
    # data
    xy <- get_xy(response, explanatory, data)
    y <- xy$y
    x <- xy$x
    
    # setups
    hess_mode <- get_post_hessian(x, y, betas=betaMode)
    proposal_cov <- -solve(hess_mode)
    
    beta_chain <- matrix(NA_real_, nrow = niter, ncol = length(betaMode))
    colnames(beta_chain) <- paste0("beta", 1:length(betaMode))
    
    # initialization
    beta_curr <- betaMode
    logpost_curr <- post_loglik(x=x, y=y, betas=beta_curr)
    
    # main sampling
    for (iter in 1:niter) {
        # sample proposal dist.
        beta_prop <- as.numeric(MASS::mvrnorm(
            n = 1,
            mu = beta_curr,
            Sigma = proposal_cov
        ))
        
        logpost_prop <- post_loglik(x, y, betas=beta_prop)
        log_accept_prob <- (logpost_prop - logpost_curr)
        
        if (log(runif(1)) <= log_accept_prob) {
            # case: accept w/ prob.
            beta_curr <- beta_prop
            logpost_curr <- logpost_prop
        }
        # case: reject
        beta_chain[iter, ] <- beta_curr
    }
    
    return(beta_chain)
}



# Problem 3 -----------------------------------------------------------------------------------

#' Goal: bayes logistic model with 1 variable, 1 outcome

bayesLogistic <- function(apredictor, response, data, NumberOfIterations) {
    
    # get beta mode
    betaMode <- get_betas_mode_NR(
        response = response,
        explanatory = apredictor,
        data = data
    )
    
    # calc. loglik
    log_marginal_lik <- getLaplaceApprox(
        response = response,
        explanatory = apredictor,
        data = data,
        betaMode = betaMode
    )
    
    # chain
    beta_chain <- getPosteriorMeans(
        response = response,
        explanatory = apredictor,
        data = data,
        betaMode = betaMode,
        niter = NumberOfIterations
    )
    beta_bayes <- colMeans(beta_chain)
    
    # mle
    dat_mat <- as.matrix(data)
    y <- dat_mat[, response]
    x <- dat_mat[, apredictor]
    mle_fit <- glm(y ~ x, family = binomial())
    mle_coef <- coef(mle_fit)
    
    # results
    res <- list(
        apredictor=apredictor,
        logmarglik=log_marginal_lik,
        beta0bayes=unname(beta_bayes[1]),
        beta1bayes=unname(beta_bayes[2]),
        beta0mle=unname(mle_coef[1]),
        beta1mle=unname(mle_coef[2])
    )
    
    return(res)   
}


# Algo: Newton-Ralphson -----------------------------------------------------------------------

get_betas_mode_NR <- function(
        response, explanatory, data, 
        tol=1e-4, max_iter=100, verbose=FALSE
    ) {
    
    # get data
    xy <- get_xy(response, explanatory, data)
    x <- xy$x
    y <- xy$y
    
    # initialization
    betas <- c(0, 0)
    
    # main
    for (iter in 1:max_iter) {
        grad <- get_post_gradient(x=x, y=y, betas=betas)
        hess <- get_post_hessian(x=x, y=y, betas=betas)
        
        betas_new <- betas - solve(hess, grad) # NOTE: better than solve(hess) %*% grad
        
        # early stop condition
        if (max(abs(betas_new - betas)) < tol) {
            if (verbose) {cat(sprintf("Converge at iter=%d", iter))}
            return(betas_new)
        }
        
        # update parameters
        betas <- betas_new
    }
    cat("NR Algo did not converge")
    return(betas)
}


# Main Parallel Function

main <- function(datafile, NumberOfIterations, clusterSize) {
    
    # load data
    data <- read.table(datafile, header = FALSE)
    
    # extract variables
    response <- ncol(data)
    lastPredictor <- ncol(data) - 1
    
    # cluster setups
    cluster <- makeCluster(clusterSize, type="SOCK")
    on.exit(stopCluster(cluster), add = TRUE)
    
    clusterEvalQ(
        cluster,
        {library(MASS)}
    ) # run expression on each cluster
    
    clusterExport(cluster, c(
        "get_xy",
        "log1pexp",
        "post_loglik",
        "get_post_gradient",
        "get_post_hessian",
        "get_betas_mode_NR",
        "getLaplaceApprox",
        "getPosteriorMeans",
        "bayesLogistic"
    ), envir = environment()) # export helpers to each cluster
    
    
    # Main estimation
    results <- clusterApply(
        cluster,
        1:lastPredictor,
        bayesLogistic,
        response,
        data,
        NumberOfIterations
    )
    
    # NOTE: really bad design
    # print out the results
    # for(i in 1:lastPredictor)
    # {
    #     cat('Regression of Y on explanatory variable ',results[[i]]$apredictor,
    #         ' has log marginal likelihood ',results[[i]]$logmarglik,
    #         ' with beta0 = ',results[[i]]$beta0bayes,' (',results[[i]]$beta0mle,')',
    #         ' and beta1 = ',results[[i]]$beta1bayes,' (',results[[i]]$beta1mle,')',
    #         '\n')    
    # }
    
    return(results)
}






# Helpers -------------------------------------------------------------------------------------

# 1. get proper data (X, y)
get_xy <- function(response, explanatory, data) {
    m <- as.matrix(data)
    return(list(
        x=m[, explanatory],
        y=m[, response]
    ))
}


# 2. log-likelihood functions

log1pexp <- function(x) {
    # if (x > 0) {
    #     return(x + log1p(exp(-x)))
    # }
    # return(log1p(exp(x)))
    #NOTE: need to vectorize func, above is not
    ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}


post_loglik <- function(x, y, betas) {
    
    beta0 <- betas[1]
    beta1 <- betas[2]
    
    temp <- beta0 + beta1 * x
    
    loglik <- sum(y * temp - log1pexp(temp))
    log_prior <- -log(2 * pi) - 0.5 * sum(betas**2)
    log_post <- log_prior + loglik
    
    return(log_post)
}

get_post_gradient <- function(x, y, betas) {
    
    beta0 <- betas[1]
    beta1 <- betas[2]
    
    temp <- beta0 + beta1 * x
    pi_ <- plogis(temp)
    
    grad0 <- sum(y - pi_) - beta0
    grad1 <- sum((y - pi_) * x) - beta1
    
    return(c(grad0, grad1))
}

get_post_hessian <- function(x, y, betas) {
    
    beta0 <- betas[1]
    beta1 <- betas[2]
    
    temp <- beta0 + beta1 * x
    pi_ <- plogis(temp)
    w <- pi_ * (1 - pi_)
    
    h00 <- -sum(w) - 1
    h01 <- -sum(w * x)
    h11 <- -sum(w * x^2) - 1
    
    hess <- matrix(c(h00, h01,
             h01, h11), nrow = 2, byrow = TRUE)
    
    return(hess)
}






# Testing -------------------------------------------------------------------------------------
# 
# df <- fread("hw2/534binarydatasmall.txt")
# 
# # NR algo
# get_betas_mode_NR(ncol(df), 1, df, verbose = T, tol = 1e-6)
# 
# # MH algo
# beta_chain <- getPosteriorMeans(
#     response = ncol(df),
#     explanatory = c(1),
#     data = df,
#     betaMode = get_betas_mode_NR(ncol(df), 1, df, verbose = T, tol = 1e-6),
#     niter = 10000
# )
# 
# # bayes logistic
# res <- bayesLogistic(1, ncol(df), df, 10000)
# 
# # full data
# df <- fread("hw2/534binarydata.txt")
# s <- Sys.time()
# res <- bayesLogistic(1, ncol(df), df, 10000)
# e <- Sys.time()
# 
# res_df <- rbindlist(res)
# 
# # viz.
# library(ggplot2)
# ggplot(res_df)+
#     aes(x=beta0mle, y=beta0bayes)+
#     # geom_line()+
#     geom_point()+
#     geom_abline(slope = 1, intercept = 0, color='salmon')
# 
# ggplot(res_df)+
#     aes(x=beta1mle, y=beta1bayes)+
#     # geom_line()+
#     geom_point()+
#     geom_abline(slope = 1, intercept = 0, color='salmon')
