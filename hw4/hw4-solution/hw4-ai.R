library(MASS)
library(snow)



# Problem 1 -----------------------------------------------------------------------------------

AI_getLaplaceApprox <- function(response, explanatory, data, betaMode) {
    # Extract variables
    y <- data[, response]
    x <- data[, explanatory]
    n <- length(y)
    
    b0 <- betaMode[1]
    b1 <- betaMode[2]
    
    # Calculate pi_i (inverse logit) [cite: 15, 21]
    # Using pi_ to avoid conflict with R's built-in pi [user request]
    linear_pred <- b0 + b1 * x
    pi_ <- exp(linear_pred) / (1 + exp(linear_pred))
    
    # Calculate log-likelihood l(beta|D) [cite: 25]
    # We use a small epsilon to avoid log(0)
    eps <- 1e-16
    log_lik <- sum(y * log(pi_ + eps) + (1 - y) * log(1 - pi_ + eps))
    
    # Calculate l*(beta) 
    # l*(beta) = -log(2*pi) - 0.5*(beta0^2 + beta1^2) + l(beta|D)
    l_star <- -log(2 * base::pi) - 0.5 * (b0^2 + b1^2) + log_lik
    
    # Calculate the Hessian matrix D^2 l*(beta) [cite: 38, 39]
    # From log-likelihood derivatives[cite: 26], the second derivatives are:
    # d2l/db0^2 = -sum(pi * (1-pi))
    # d2l/db0db1 = -sum(x * pi * (1-pi))
    # d2l/db1^2 = -sum(x^2 * pi * (1-pi))
    # The prior adds -1 to the diagonal elements of the Hessian
    
    weight <- pi_ * (1 - pi_)
    h00 <- -sum(weight) - 1
    h01 <- -sum(x * weight)
    h11 <- -sum(x^2 * weight) - 1
    
    hessian <- matrix(c(h00, h01, h01, h11), nrow = 2)
    
    # Calculate log marginal likelihood [cite: 52, 55]
    # log(2*pi) + l_star - 0.5 * log(det(-Hessian))
    det_neg_hessian <- det(-hessian)
    log_marg_lik <- log(2 * base::pi) + l_star - 0.5 * log(det_neg_hessian)
    
    return(log_marg_lik)
}



# Problem 2 -----------------------------------------------------------------------------------

AI_getPosteriorMeans <- function(response, explanatory, data, betaMode, niter) {
   
    # Extract data [cite: 7]
    y <- data[, response]
    x <- data[, explanatory]
    
    # 1. Setup the Proposal Covariance Matrix [cite: 65]
    # We need the Hessian at the mode to define the proposal distribution
    b0_mode <- betaMode[1]
    b1_mode <- betaMode[2]
    
    # Calculate pi_i at the mode [cite: 15, 21]
    lin_pred_mode <- b0_mode + b1_mode * x
    pi_mode <- exp(lin_pred_mode) / (1 + exp(lin_pred_mode))
    
    # Compute Hessian D^2 l*(beta) at the mode [cite: 39]
    weight_mode <- pi_mode * (1 - pi_mode)
    h00 <- -sum(weight_mode) - 1
    h01 <- -sum(x * weight_mode)
    h11 <- -sum(x^2 * weight_mode) - 1
    hessian_mode <- matrix(c(h00, h01, h01, h11), nrow = 2)
    
    # Proposal covariance is the negative inverse of the Hessian [cite: 65]
    prop_sigma <- solve(-hessian_mode)
    
    # 2. Define the l*(beta) function for log-posterior density [cite: 32]
    get_l_star <- function(beta_vec) {
        b0 <- beta_vec[1]
        b1 <- beta_vec[2]
        
        lp <- b0 + b1 * x
        # Using pi_ to avoid built-in pi constant [user instruction]
        pi_ <- exp(lp) / (1 + exp(lp))
        
        # Log-likelihood [cite: 25]
        # Add small epsilon to avoid log(0)
        log_lik <- sum(y * log(pi_ + 1e-16) + (1 - y) * log(1 - pi_ + 1e-16))
        
        # Return l*(beta) [cite: 32]
        return(-log(2 * base::pi) - 0.5 * (b0^2 + b1^2) + log_lik)
    }
    
    # 3. Metropolis-Hastings Algorithm [cite: 57, 75]
    # Initialize storage and starting values [cite: 58, 94]
    samples <- matrix(0, nrow = niter, ncol = 2)
    samples[1, ] <- betaMode
    current_l_star <- get_l_star(betaMode)
    
    for (k in 2:niter) {
        # Generate candidate from bivariate normal [cite: 62, 63]
        candidate <- MASS::mvrnorm(n = 1, mu = samples[k - 1, ], Sigma = prop_sigma)
        candidate_l_star <- get_l_star(candidate)
        
        # Calculate log acceptance probability [cite: 70, 82]
        log_alpha <- candidate_l_star - current_l_star
        
        # Acceptance step [cite: 78, 81, 84, 86]
        if (log(runif(1)) <= log_alpha) {
            samples[k, ] <- candidate
            current_l_star <- candidate_l_star
        } else {
            samples[k, ] <- samples[k - 1, ]
        }
    }
    
    # 4. Calculate posterior means 
    beta0_mean <- mean(samples[, 1])
    beta1_mean <- mean(samples[, 2])
    
    return(c(beta0 = beta0_mean, beta1 = beta1_mean))
}



# Problem 3 -----------------------------------------------------------------------------------

# Helpers
# Find the mode of the posterior using Newton-Raphson
getMode <- function(response, explanatory, data) {
    y <- data[, response]
    x <- data[, explanatory]
    beta <- c(0, 0) # Initial values (0,0) [cite: 44]
    eps <- 0.0001
    diff <- 1
    
    while(diff > eps) {
        beta_old <- beta
        lp <- beta[1] + beta[2] * x
        pi_ <- exp(lp) / (1 + exp(lp)) # Inverse logit [cite: 21]
        
        # Gradient of l* [cite: 26, 37]
        grad <- c(sum(y - pi_) - beta[1], sum(y * x - pi_ * x) - beta[2])
        
        # Hessian of l* [cite: 38, 39]
        weight <- pi_ * (1 - pi_)
        h00 <- -sum(weight) - 1
        h01 <- -sum(x * weight)
        h11 <- -sum(x^2 * weight) - 1
        hessian <- matrix(c(h00, h01, h01, h11), nrow = 2)
        
        # Update step [cite: 46]
        beta <- beta - solve(hessian) %*% grad
        diff <- max(abs(beta - beta_old)) # Convergence check [cite: 47]
    }
    return(as.vector(beta))
}

# baye logistic
AI_bayesLogistic <- function(apredictor, response, data, NumberOfIterations) {
    # 1. Find the mode via Newton-Raphson
    beta_mode <- getMode(response, apredictor, data)
    
    # 2. Compute Laplace Approximation (Log Marginal Likelihood) [cite: 89]
    log_ml <- AI_getLaplaceApprox(response, apredictor, data, beta_mode)
    
    # 3. Compute Posterior Means via Metropolis-Hastings [cite: 90, 91]
    post_means <- AI_getPosteriorMeans(response, apredictor, data, beta_mode, NumberOfIterations)
    
    # 4. Optional: Calculate MLEs using standard logistic regression
    fit_mle <- glm(data[, response] ~ data[, apredictor], family = "binomial")
    mles <- coef(fit_mle)
    
    # Return the results as a list [cite: 135, 136]
    return(list(
        apredictor = apredictor,
        logmarglik = log_ml,
        beta0bayes = post_means[1],
        beta1bayes = post_means[2],
        beta0mle = mles[1],
        beta1mle = mles[2]
    ))
}


# main
library(snow) # Required for parallel computing [cite: 140]

AI_main <- function(datafile, NumberOfIterations, clusterSize) {
    data <- read.table(datafile, header = FALSE)
    
    response_idx <- ncol(data)
    lastPredictor <- ncol(data) - 1
    
    cluster <- makeCluster(clusterSize, type = "SOCK")
    
    clusterExport(
        cluster,
        c("getMode", "AI_getLaplaceApprox", "AI_getPosteriorMeans", "AI_bayesLogistic")
    )
    
    results <- clusterApply(
        cluster,
        1:lastPredictor,
        AI_bayesLogistic,
        response_idx,
        data,
        NumberOfIterations
    )
    
    for (i in 1:lastPredictor) {
        cat(
            "Regression of Y on explanatory variable", results[[i]]$apredictor,
            "has log marginal likelihood", results[[i]]$logmarglik,
            "with beta0 =", results[[i]]$beta0bayes, "(", results[[i]]$beta0mle, ")",
            "and beta1 =", results[[i]]$beta1bayes, "(", results[[i]]$beta1mle, ")",
            "\n"
        )
    }
    
    stopCluster(cluster)
    return(results)
}


# Example execution [cite: 142]
# res <- AI_main('534binarydata.txt', 10000, 10)