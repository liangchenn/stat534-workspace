#' README
#' ------
#' - author: Liang-Cheng Chen
#' - created_at: 2026-04-27
#' - desc:
#'   - implement MC3 algorithm



# Packages ------------------------------------------------------------------------------------

library(rcdd)
library(data.table)


# Problem 1 -----------------------------------------------------------------------------------

# main function

MC3search <- function(response, data, n_iter) {
    
    s <- Sys.time()
    # available variables set
    available_vars <- setdiff(1:ncol(data), response)
    
    # generate subset randomly
    while(T) {
        subset_vars <- .generate_subset(available_vars)
        k <- length(subset_vars)
        # check if logistic valid, if not keep sampling
        if (isValidLogisticRCDD(response = response, explanatory = subset_vars, data = data)) {
            # if valid, proceed
            break
        }
    }
    
    # initial model
    curr_model <- subset_vars
    # initial best model
    curr_best <- curr_model
    
    curr_aic <- getLogisticAIC(response, curr_model, data)
    best_aic <- curr_aic
    
    # MC iterations
    for (i in 1:n_iter) {
        
        # print(sprintf("...Running for iteration %d / %d", i, n_iter))
        
        
        # step 1 & 2: find valid neighbor variables
        valid_neighbors <- .find_valid_neighbor(response = response, explanatory = curr_model, data = data)

        # step 3: uniformly sample a model from valid neighbors
        model_candidate <- valid_neighbors[[sample(1:length(valid_neighbors), 1)]] # it's a list
        
        # step 4: find neighbor of curr model candidate
        neighbors_of_candidate <- .find_valid_neighbor(response = response, explanatory = model_candidate, data = data)
        
        # step 5: calc. p for MA'
        candidate_aic <- getLogisticAIC(response, model_candidate, data)
        p_candidate <- -candidate_aic - log(length(neighbors_of_candidate))
        
        # step 6: calc. p for MAr
        p_curr <- -curr_aic - log(length(valid_neighbors))
        
        # step 7 & 8: current model picking w/ accept-rejecting sampling
        # NOTE: more readable to have accept case organized
        accept <- FALSE
        
        if (p_candidate > p_curr) {
            accept <- TRUE
        } else {
            u <- runif(1)
            if (log(u) < p_candidate - p_curr) {
                accept <- TRUE
            }
        }
        
        if (accept) {
            # update current model if accepting new candidate
            curr_model <- model_candidate
            curr_aic <- candidate_aic
        }
        
        if (curr_aic < best_aic) {
            curr_best <- curr_model
            best_aic <- curr_aic
        }

        
    } # end of iteration loop
    e <- Sys.time()
    # print(sprintf("time elapsed: %.4f s", e-s))
    return(list(
        bestAIC = best_aic,
        bestAICvars = sort(curr_best)
    ))
}


# helper functions

.generate_subset <- function(variables) {
    # step 1. decide a size k of variables
    k <- sample(0:length(variables), 1) # NOTE: should include empty case
    # step 2. sample k variable from variable set
    subset_vars <- sample(variables, size = k)
    
    return(subset_vars)
}


isValidLogisticRCDD <- function(response, explanatory, data) {
    if(0==length(explanatory))
    {
        #we assume that the empty logistic regresion is valid
        return(TRUE);
    }
    logisticreg = suppressWarnings(glm(data[,response] ~ as.matrix(data[,as.numeric(explanatory)]),family=binomial(link=logit),x=TRUE));
    tanv = logisticreg$x;
    tanv[data[,response] == 1, ] <- (-tanv[data[,response] == 1, ]);
    vrep = cbind(0, 0, tanv);
    #with exact arithmetic; takes a long time
    #lout = linearity(d2q(vrep), rep = "V");

    lout = linearity(vrep, rep = "V");
    return(length(lout)==nrow(data));
}


.find_valid_neighbor <- function(response, explanatory, data) {
    # find current available vars set
    available_vars <- setdiff(1:ncol(data), response)
    remaining_set <- setdiff(available_vars, explanatory)
    # add one var
    plus_one <- lapply(remaining_set, function(var) {union(explanatory, var)})
    # minus one var
    minus_one <- lapply(explanatory, function(var) {setdiff(explanatory, var)})
    # NOTE: here if the minus one is empty, return numeric(0) with len = 0
    
    # find valid neighbors
    valid_neighbors <- Filter(
        function(nbd) {
            isValidLogisticRCDD(response = response, explanatory = nbd, data = data)
        },
        c(plus_one, minus_one)
    )

    return(valid_neighbors)
}



.getLogisticKIC <- function(response, explanatory, data, type = "A") {
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
    return(.getLogisticKIC(response, explanatory, data, type = "A") |> suppressWarnings())
}



# Testing -------------------------------------------------------------------------------------

# df <- fread("hw2/534binarydatasmall.txt")
# df <- fread("hw2/534binarydata.txt")
# 
# setDF(df)
# 
# 
# MC3search(ncol(df), df, n_iter = 25)
# 
# s <- proc.time()
# MC3search(ncol(df), df, n_iter = 10)
# e <- proc.time()
# 
# library(parallel)
# set.seed(9527)
# results <- mclapply(
#     X = 1:2,
#     FUN = function(i) {
#         MC3search(ncol(df), df, n_iter = 25)
#     },
#     mc.cores = detectCores() - 1
# )