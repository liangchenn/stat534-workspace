# This script is for AI helper functions
# Helper to calculate AIC for a specific set of explanatory variables
# isValidLogisticRCDD <- function(response, explanatory, data) {return(TRUE)}

AI_get_AIC <- function(response, explanatory, data) {
    if (length(explanatory) == 0) {
        model <- glm(data[, response] ~ 1, family = binomial(link = "logit"))
    } else {
        model <- glm(data[, response] ~ as.matrix(data[, explanatory]), 
                     family = binomial(link = "logit"))
    }
    return(AIC(model))
}

# Helper to find all valid neighboring models
AI_get_valid_neighbors <- function(current_vars, all_vars, response, data) {
    valid_neighbors <- list()
    
    # Potential neighbors: Addition
    to_add <- setdiff(all_vars, current_vars)
    for (v in to_add) {
        neighbor <- sort(c(current_vars, v))
        if (isValidLogisticRCDD(response, neighbor, data)) {
            valid_neighbors[[length(valid_neighbors) + 1]] <- neighbor
        }
    }
    
    # Potential neighbors: Deletion
    if (length(current_vars) > 0) {
        for (i in 1:length(current_vars)) {
            neighbor <- current_vars[-i]
            if (isValidLogisticRCDD(response, neighbor, data)) {
                valid_neighbors[[length(valid_neighbors) + 1]] <- neighbor
            }
        }
    }
    
    return(valid_neighbors)
}

# main function design
AI_M3search <- function(response, data, n_iter) {
    # 1. Initialization
    all_vars <- setdiff(1:ncol(data), response)
    
    # Randomly generate a valid starting model
    # We'll try random subsets until one is valid
    valid_start <- FALSE
    while (!valid_start) {
        start_size <- sample(0:length(all_vars), 1)
        current_model <- sort(sample(all_vars, start_size))
        if (isValidLogisticRCDD(response, current_model, data)) {
            valid_start <- TRUE
        }
    }
    
    best_model <- current_model
    best_aic <- AI_get_AIC(response, best_model, data)
    
    # 2. Iteration Loop
    for (r in 1:n_iter) {
        # Step 1 & 2: Get valid neighbors of current model
        nbd_curr <- AI_get_valid_neighbors(current_model, all_vars, response, data)
        
        # Step 3: Uniformly sample a candidate model M_prime
        idx_prime <- sample(seq_along(nbd_curr), 1)
        m_prime <- nbd_curr[[idx_prime]]
        
        # Step 4: Get valid neighbors of the candidate model
        nbd_prime <- AI_get_valid_neighbors(m_prime, all_vars, response, data)
        
        # Step 5 & 6: Calculate acceptance probabilities
        # p = -AIC - log(#nbd)
        aic_prime <- AI_get_AIC(response, m_prime, data)
        aic_curr <- AI_get_AIC(response, current_model, data)
        
        p_prime <- -aic_prime - log(length(nbd_prime))
        p_curr <- -aic_curr - log(length(nbd_curr))
        
        # Step 7 & 8: Metropolis-Hastings Acceptance Step
        accepted <- FALSE
        if (p_prime > p_curr) {
            accepted <- TRUE
        } else {
            u <- runif(1)
            if (log(u) < (p_prime - p_curr)) {
                accepted <- TRUE
            }
        }
        
        if (accepted) {
            current_model <- m_prime
            # Update global best if the new current model is better than any seen before
            if (aic_prime < best_aic) {
                best_model <- m_prime
                best_aic <- aic_prime
            }
        }
        # If not accepted, current_model and best_model remain unchanged
    }
    
    # 3. Return results
    return(list(
        "best model" = best_model,
        "best AIC" = best_aic
    ))
}