# STAT/BIOSTAT 534 Homework 4 solution
# Bayesian univariate logistic regression with Laplace approximation and
# Metropolis-Hastings sampling.

# Numerically stable log(1 + exp(x)).
log1pexp <- function(x) {
  ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

# Numerically stable inverse-logit exp(x) / (1 + exp(x)).
inv_logit <- function(x) {
  ifelse(x >= 0, 1 / (1 + exp(-x)), exp(x) / (1 + exp(x)))
}

# Build the two-column design matrix for one univariate logistic regression.
# Column 1 is the intercept; column 2 is the selected explanatory variable.
design_matrix <- function(x) {
  cbind(Intercept = 1, X = as.numeric(x))
}

# Log-likelihood for logistic regression:
# sum_i y_i * eta_i - log(1 + exp(eta_i)), where eta_i = beta0 + beta1*x_i.
log_likelihood <- function(beta, x, y) {
  X <- design_matrix(x)
  eta <- drop(X %*% beta)
  sum(y * eta - log1pexp(eta))
}

# l*(beta0, beta1) from the homework PDF. This is the log-likelihood plus
# the log N(0,1) prior densities for beta0 and beta1.
log_posterior_kernel <- function(beta, x, y) {
  -log(2 * pi) - 0.5 * sum(beta^2) + log_likelihood(beta, x, y)
}

# Gradient of l*. The likelihood contribution is X' * (y - pi); the prior
# contribution is -beta.
gradient_log_posterior <- function(beta, x, y) {
  X <- design_matrix(x)
  eta <- drop(X %*% beta)
  pi_hat <- inv_logit(eta)
  drop(t(X) %*% (y - pi_hat)) - beta
}

# Hessian matrix of l*. For logistic regression, the likelihood Hessian is
# -X' W X, where W has diagonal entries pi_i * (1 - pi_i). The prior adds -I.
hessian_log_posterior <- function(beta, x, y) {
  X <- design_matrix(x)
  eta <- drop(X %*% beta)
  pi_hat <- inv_logit(eta)
  w <- pi_hat * (1 - pi_hat)
  -crossprod(X, X * w) - diag(2)
}

# Find the posterior mode using Newton-Raphson, starting from (0, 0).
newton_raphson_mode <- function(x, y, epsilon = 1e-4, max_iter = 100) {
  beta <- c(0, 0)
  current_value <- log_posterior_kernel(beta, x, y)

  for (iter in seq_len(max_iter)) {
    grad <- gradient_log_posterior(beta, x, y)
    hess <- hessian_log_posterior(beta, x, y)

    # The PDF update is beta_new = beta - H^{-1} grad.
    step <- tryCatch(solve(hess, grad), error = function(e) MASS::ginv(hess) %*% grad)
    candidate <- as.numeric(beta - step)

    # Small step-halving safeguard: Newton steps should improve near the mode,
    # but this prevents numerical overshoot for difficult predictors.
    candidate_value <- log_posterior_kernel(candidate, x, y)
    shrink <- 1
    while (!is.finite(candidate_value) && shrink > 1e-8) {
      shrink <- shrink / 2
      candidate <- as.numeric(beta - shrink * step)
      candidate_value <- log_posterior_kernel(candidate, x, y)
    }

    while (is.finite(candidate_value) &&
           candidate_value < current_value &&
           shrink > 1e-8) {
      shrink <- shrink / 2
      candidate <- as.numeric(beta - shrink * step)
      candidate_value <- log_posterior_kernel(candidate, x, y)
    }

    if (all(abs(candidate - beta) < epsilon)) {
      beta <- candidate
      break
    }

    beta <- candidate
    current_value <- candidate_value
  }

  names(beta) <- c("beta0", "beta1")
  beta
}

# Problem 1: log of the Laplace approximation to the marginal likelihood.
# Formula: log(2*pi) + l*(mode) - 0.5*log(det(-H(mode))).
laplace_log_marginal_likelihood <- function(x, y, mode = NULL) {
  if (is.null(mode)) {
    mode <- newton_raphson_mode(x, y)
  }

  hess <- hessian_log_posterior(mode, x, y)
  precision <- -hess
  log_det_precision <- as.numeric(determinant(precision, logarithm = TRUE)$modulus)

  log(2 * pi) + log_posterior_kernel(mode, x, y) - 0.5 * log_det_precision
}

# Draw one sample from a bivariate normal distribution. MASS is part of the
# recommended R distribution and is available on this machine.
rmvnorm2 <- function(mean, sigma) {
  as.numeric(MASS::mvrnorm(n = 1, mu = mean, Sigma = sigma))
}

# Problem 2: Metropolis-Hastings sampler and posterior mean estimates.
metropolis_hastings <- function(x, y, NumberOfIterations, mode = NULL) {
  if (is.null(mode)) {
    mode <- newton_raphson_mode(x, y)
  }

  hess_at_mode <- hessian_log_posterior(mode, x, y)
  proposal_cov <- solve(-hess_at_mode)

  samples <- matrix(NA_real_, nrow = NumberOfIterations, ncol = 2)
  colnames(samples) <- c("beta0", "beta1")

  current <- as.numeric(mode)
  current_log_post <- log_posterior_kernel(current, x, y)
  accepted <- 0L

  for (k in seq_len(NumberOfIterations)) {
    proposed <- rmvnorm2(current, proposal_cov)
    proposed_log_post <- log_posterior_kernel(proposed, x, y)

    # Work on the log scale: accept if log(U) <= l*(proposal) - l*(current).
    log_acceptance_ratio <- proposed_log_post - current_log_post
    if (log(runif(1)) <= log_acceptance_ratio) {
      current <- proposed
      current_log_post <- proposed_log_post
      accepted <- accepted + 1L
    }

    samples[k, ] <- current
  }

  list(
    samples = samples,
    beta_mean = colMeans(samples),
    acceptance_rate = accepted / NumberOfIterations
  )
}

# Optional check against frequentist MLEs from glm().
compute_mle <- function(x, y) {
  fit <- glm(y ~ x, family = binomial())
  coef_values <- coef(fit)
  names(coef_values) <- c("beta0", "beta1")
  coef_values
}

# This is the function requested in the starter code. It runs one univariate
# Bayesian logistic regression: Y against predictor column apredictor.
bayesLogistic <- function(apredictor, response, data, NumberOfIterations) {
  x <- as.numeric(data[[apredictor]])
  y <- as.numeric(data[[response]])

  mode <- newton_raphson_mode(x, y)
  log_marg_lik <- laplace_log_marginal_likelihood(x, y, mode)
  mh <- metropolis_hastings(x, y, NumberOfIterations, mode)
  mle <- compute_mle(x, y)

  list(
    apredictor = apredictor,
    logmarglik = log_marg_lik,
    beta0bayes = mh$beta_mean[1],
    beta1bayes = mh$beta_mean[2],
    beta0mle = mle[1],
    beta1mle = mle[2],
    beta0mode = mode[1],
    beta1mode = mode[2],
    acceptance_rate = mh$acceptance_rate
  )
}

# Choose the local dataset name. The PDF says 534binarydata.txt, but this
# directory contains 534binarydata-1.txt.
default_datafile <- function() {
  if (file.exists("534binarydata.txt")) {
    "534binarydata.txt"
  } else if (file.exists("534binarydata-1.txt")) {
    "534binarydata-1.txt"
  } else {
    stop("Could not find 534binarydata.txt or 534binarydata-1.txt.")
  }
}

# Parallel driver for Problem 3. If the snow package is installed, it uses
# snow. Otherwise it falls back to base R's parallel package with the same
# clusterApply-style workflow so the script remains runnable.
main <- function(datafile = default_datafile(), NumberOfIterations = 10000,
                 clusterSize = 10) {
  data <- read.table(datafile, header = FALSE)

  response <- ncol(data)
  lastPredictor <- ncol(data) - 1
  predictors <- seq_len(lastPredictor)

  if (clusterSize <= 1) {
    results <- lapply(
      predictors, bayesLogistic,
      response = response, data = data,
      NumberOfIterations = NumberOfIterations
    )
  } else if (requireNamespace("snow", quietly = TRUE)) {
    cluster <- snow::makeCluster(clusterSize, type = "SOCK")
    on.exit(snow::stopCluster(cluster), add = TRUE)
    snow::clusterEvalQ(cluster, library(MASS))
    snow::clusterExport(
      cluster,
      c("log1pexp", "inv_logit", "design_matrix", "log_likelihood",
        "log_posterior_kernel", "gradient_log_posterior",
        "hessian_log_posterior", "newton_raphson_mode",
        "laplace_log_marginal_likelihood", "rmvnorm2",
        "metropolis_hastings", "compute_mle", "bayesLogistic"),
      envir = environment()
    )
    results <- snow::clusterApply(
      cluster, predictors, bayesLogistic,
      response, data, NumberOfIterations
    )
  } else {
    cluster <- parallel::makeCluster(clusterSize, type = "PSOCK")
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterEvalQ(cluster, library(MASS))
    parallel::clusterExport(
      cluster,
      c("log1pexp", "inv_logit", "design_matrix", "log_likelihood",
        "log_posterior_kernel", "gradient_log_posterior",
        "hessian_log_posterior", "newton_raphson_mode",
        "laplace_log_marginal_likelihood", "rmvnorm2",
        "metropolis_hastings", "compute_mle", "bayesLogistic"),
      envir = environment()
    )
    results <- parallel::clusterApply(
      cluster, predictors, bayesLogistic,
      response, data, NumberOfIterations
    )
  }

  for (i in predictors) {
    cat(
      "Regression of Y on explanatory variable ", results[[i]]$apredictor,
      " has log marginal likelihood ", results[[i]]$logmarglik,
      " with beta0 = ", results[[i]]$beta0bayes, " (", results[[i]]$beta0mle, ")",
      " and beta1 = ", results[[i]]$beta1bayes, " (", results[[i]]$beta1mle, ")",
      " acceptance rate = ", results[[i]]$acceptance_rate,
      "\n",
      sep = ""
    )
  }

  invisible(results)
}

# When this file is executed with Rscript, run the complete homework job.
# When it is sourced from another R session, only define the functions.
if (sys.nframe() == 0) {
  set.seed(534)
  main(default_datafile(), NumberOfIterations = 10000, clusterSize = 10)
}
