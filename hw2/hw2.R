#' README
#' ------
#' - author: Liang-Cheng Chen
#' - created_at: 2026-04-20
#' - desc:
#'   - implement kic for logistic regression given model

# setwd(".")

# Problem 1 -----------------------------------------------------------------------------------

getLogisticKIC <- function(response, explanatory, data, type = "A") {
    #' 
    #' Args:
    #' -----
    #' - response (int): the response variable index.
    #' - explanatory (vector[int]): the vector of explanatory variables indexes.
    #' - data (data.frame like): the target data in data.frame-like object.
    #' - type (str): {A, B}, denoting kic, BIC, respectively. Default to use AIC.
    #' 
    #' Return
    #' ------
    #' - KIC (num): KIC value.
    
    # arg check
    if (!type %in% c("A", "B")) {stop(sprintf("Type %s is not allowed. Should be 'A' or 'B'", type))}
    
    # accomodate different data inputs
    data <- as.matrix(data)
    
    # logistic regression model
    if(!length(explanatory)){
        # empty variable case
        model <- glm(data[, response] ~ 1, family = binomial(link='logit'))
    } else {
        # regular case
        model <- glm(data[,response] ~ data[, as.numeric(explanatory)],
                       family = binomial(link='logit'))
    }
    
    # KIC calc.
    multiplier <- ifelse(type == "A", 2, log(nrow(data)))
    deviance_ <- model$deviance
    kic <- deviance_ + multiplier * (1 + length(explanatory))
    
    return(kic)
}

getLogisticAIC <- function(response, explanatory, data) {
    #TODO: there's warnings about the glm estimation, currently suppress, worth investigating.
    return(getLogisticKIC(response, explanatory, data, type = "A") |> suppressWarnings())
}

getLogisticBIC <- function(response, explanatory, data) {
    #TODO: there's warnings about the glm estimation, currently suppress, worth investigating.
    return(getLogisticKIC(response, explanatory, data, type = "B") |> suppressWarnings())
}



# Problem 2 -----------------------------------------------------------------------------------

forwardSearch <- function(response, data, KICfunc) {
    #'
    #' Args
    #' ----
    #' - response (int): response var index.
    #' - data (df | dt | matrix): target data in data.frame like object.
    #' - KICfunc (function): function to compute AIC/BIC value.
    #'
    #' Returns
    #' -------
    #' - A list with:
    #'   - model: selected explanatory variable indices.
    #'   - kic: criterion value of the final model.
    #'   - path: order of variables added during forward search.
    
    # initial setups
    all_vars <- setdiff(1:ncol(data), response)
    n_all_vars <- length(all_vars)
    
    explanatory <- c()
    
    
    # initial kic
    kic_curr <- KICfunc(
        response = response, 
        explanatory = explanatory,
        data = data
    )
    
    while((length(explanatory) < n_all_vars)) {
        
        # add new variable from remaining set
        available_vars_idx <- setdiff(all_vars, explanatory)
        kics <- sapply(available_vars_idx, function(e) {KICfunc(response, union(explanatory, e), data)})
        
        if (kic_curr < min(kics)) {
            # case: stop when there no kic improvement
            break
        }
        
        # regular case: add variable to model that has min kic
        explanatory <- c(explanatory, available_vars_idx[which.min(kics)])
        kic_curr <- min(kics)
    }
    
    return(list(model=sort(explanatory), kic=kic_curr, path=explanatory))
}


forwardSearchAIC <- function(response, data) {
    return(forwardSearch(response, data, getLogisticAIC))
}

# forwardSearchBIC <- function(response, data) {
#     return(forwardSearch(response, data, getLogisticBIC))
# }


# Problem 3 -----------------------------------------------------------------------------------

backwardSearch <- function(response, data, KICfunc) {
    #'
    #' Args
    #' ----
    #' - response (int): response var index.
    #' - data (df | dt | matrix): target data in data.frame like object.
    #' - KICfunc (function): function to compute AIC/BIC value.
    #'
    #' Returns
    #' -------
    #' - A list with:
    #'   - model: selected explanatory variable indices.
    #'   - kic: criterion value of the final model.
    #'   - path: order of variables added during backward search.
    
    # initial setups
    all_vars <- setdiff(1:ncol(data), response)
    n_all_vars <- length(all_vars)
    explanatory <- all_vars
    path <- c()
    
    # initial kic
    kic_curr <- KICfunc(
        response = response, 
        explanatory = explanatory,
        data = data
    )
    
    while((length(explanatory) > 0)) {
        
        # remove variable from the explanatory set
        available_vars_idx <- explanatory
        kics <- sapply(available_vars_idx, function(e) {KICfunc(response, setdiff(explanatory, e), data)})
        
        if (kic_curr < min(kics)) {
            # case: stop when there no kic improvement
            break
        }
        
        # regular case: remove variable from model that has min kic
        explanatory <- setdiff(explanatory, available_vars_idx[which.min(kics)])
        path <- c(path, available_vars_idx[which.min(kics)])
        kic_curr <- min(kics)
    }
    
    return(list(model=sort(explanatory), kic=kic_curr, path=path))
}


backwardSearchAIC <- function(response, data) {
    return(backwardSearch(response, data, getLogisticAIC))
}

# backwardSearchBIC <- function(response, data) {
#     return(backwardSearch(response, data, getLogisticBIC))
# }


# Problem 4 -----------------------------------------------------------------------------------

# foward search w/ BIC
forwardSearchBIC <- function(response, data) {
    return(forwardSearch(response, data, getLogisticBIC))
}

# backward search w/ BIC
backwardSearchBIC <- function(response, data) {
    return(backwardSearch(response, data, getLogisticBIC))
}



# Testing Area --------------------------------------------------------------------------------

# test small data
test_df <- fread("hw2/534binarydatasmall.txt")

forwardSearchAIC(response = 11, data = test_df)
backwardSearchAIC(response = 11, data = test_df)

forwardSearchBIC(response = 11, data = test_df)
backwardSearchBIC(response = 11, data = test_df)
# 
# 
# # full data
# df <- fread("hw2/534binarydata.txt")
# 
# forwardSearchAIC(response = 61, data = df)
# backwardSearchAIC(response = 61, data = df)
# 
# forwardSearchBIC(response = 61, data = df)
# backwardSearchBIC(response = 61, data = df)


# forwardSearchAIC(61, df)
# backwardSearchAIC(61, df)



