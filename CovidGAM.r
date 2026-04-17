#######################################
### PROJECT 3: SMOOTH DECONVOLUTION ###
#######################################

### Contributions of each member ###

# Isaac's Contributions: Worked on optimization of the smoothing parameter, lambda, by means 
# of the minimize_BIC, best_fit and BIC functions. Also, I helped refactor code for improved 
# readability and parameter flexibility, so it could potentially be configured for different 
# datasets.

# Eoghan's Contributions: Worked on the following functions: init_model_matrices,
# pnll, gll, fd, np_boot, the finite-difference comparison and plotting the bootstrap
# CI.

# Seán's Contributions: Worked on the following functions: init_model_matrices, pnll, gll, 
# and plotted model with lambda = 5e-5. Helped with final plot also.  


# Github repo link: https://github.com/YohoCodes/StatProj4_ISE.git


### Overview & high level outline ###

# In this project we first fit a deconvolution model to Covid-19 death data for 
# 150 days in 2020 given an arbitrary lambda. Lambda is a smoothing parameter
# for our penalised negative log likelihood for the Poisson distribution that models
# deaths. The smoothing parameter allows us to prevent overfitting our data. 
# To balance model fit and complexity, we find a value of lambda that minimises the 
# model BIC. The fitting is implemented using the "BFGS" method built in to optim 
# to minimise the penalised negative log-likelihood, using exact derivatives.

# Essentially, we are trying to recover the trajectory of the infection curve of 
# the cases that lead to the observed deaths in this time period. We are able to
# do this as we have good data to suggest that if d is the interval from infection
# to death, then log(d) ~ N(3.151, 0.469^2). There will be some uncertainty 
# associated with this estimated infection curve, so we perform non-parametric
# bootstrapping to estimate a 95% confidence interval for daily new infection rate.

###############################################
### Please Set your Working Directory Below ###
###############################################

#setwd('') # <<<<<<<< User's Working Directory

# Importing library for spline design
library(splines)

##########################
### Defining functions ###
##########################


init_model_matrices <- function(K, t, spline_range, min_window, max_window, edur, sdur) {
  
  # The purpose of this function is to calculate 3 matrices used in our deconvolution
  # model discussed below
  
  # The inputs for this function are as follows:
  
    # K: Number of splines used to cover the infection period considered
    # t: The days in the dataset (days on which deaths are observed)
    # spline_range: The days on which infections can occur leading to observed deaths
                  # These are the days we model using cubic splines
    # min_window: The minimum number of days before an individual's death that we consider
                  # as a possible day of infection leading to that death
    # max_window: The maximium number of days before an individual's death that we consider
                # as a possible day of infection leading to that death
    # edur: Expected log interval from infection to death
    # sdur: Standard deviation of log days from infection to death
    
  # This function outputs 3 matrices, X, X_tilde and S, that represent the model 
  # matrix for deaths, the model matrix for infections and the penalty matrix respectively
  
  
  # Creating a total of K+4 knots to create K cubic splines
  # (K-2) of these knots cover the interval over which infections are evaluated
  inner_knots <- seq(min(spline_range), max(spline_range), length.out = K-2)
  # Calculating the interval between knots
  int_length <- inner_knots[2] - inner_knots[1]
  # Adding 3 knots to either side of the interval "spline_range", giving us (K+4) total knots.
  outer_left_knots <- seq(min(spline_range) - (3*int_length), min(spline_range) - int_length, by = int_length)
  outer_right_knots <- seq(max(spline_range) + int_length, max(spline_range) + (3*int_length), by = int_length)
  knot_vals <- c(outer_left_knots, inner_knots, outer_right_knots)
  # The columns of X_tilde representing 80 unscaled splines
  X_tilde <- splineDesign(knots = knot_vals, x = spline_range)
  
  # n (amount of rows in X) is the number of days on which deaths were observed
  n <- length(t) 
  d <- 1:max_window
  
  # pd is the probability vector for days from infection to death, where pd_j is the probability of
  # a death on day i being caused by infection on day (i-j)
  pd <- dlnorm(d,edur,sdur); pd <- pd/sum(pd) 
  
  # Initialising an empty model matrix
  X <- matrix(0, nrow = n, ncol = K)
  
  # For each i, a different window of X_tilde is evaluated to create a row of model matrix X
  for (i in 1:n) {
    # Window width increases until i = max_window - min_window + 1
    window_width <- min(min_window - 1 + i, max_window) 
    # Vector of all j indices from 1 to the window width
    j <- 1:window_width
    # Rows of X_tilde to be summed after multiplying each of them by pd[j]
    rows <- min_window + i - j
    # Each row of X_tilde within the current window is multiplied by pd_j 
    # Columns are then summed to create row i of our model matrix X
    X[i, ] <- colSums(X_tilde[rows, , drop = FALSE] * pd[j])
  }
  
  # A square matrix S that constructs the penalty term in the final log likelihood function
  S <- crossprod(diff(diag(K),diff=2))
  
  return(list(X_tilde = X_tilde,X = X,S = S))
}


pnll <- function(gamma, X, y, S, lambda = 5e-5, wb=1) {
  
  # The purpose of this function is to calculate the penalised negative log likelihood
  # of the Poisson distribution that models deaths
  
  # The inputs for this function are as follows:
  
    # gamma: Reparameterisation of our spline coefficient vector beta. gamma_k is exp(beta_k),
           # used to ensure that the estimated daily new infection rate f(t) is positive
    # X: Our model matrix for deaths, from "init_model_matrices"
    # y: Our observed deaths (variable "nhs")
    # S: Penalty matrix, from "init_model_matrices"
    # lambda: A smoothing parameter for our penalised negative log likelihood (PNLL), default
            # value of 0.00005
    # wb: A weight parameter used later when "pnll" is called within the non-parametric 
        # bootstrapping function "np_boot", default = 1
  
  # This function outputs our penalised negative log-likelihood, a scalar, evaluated at gamma, a vector
  
  # Transforming our gammas to betas for calculations
  beta <- exp(gamma)
  # Computing expected number of deaths per day
  mu <- as.vector(X %*% beta)
  # Computing negative log-likelihood for Poisson-distributed deaths 
  # "wb" has no effect when taking its default value
  nll <- -sum(wb * (y * log(mu) - mu)) 
  # Defining a smoothing penalty which will be added to the log-likelihood:
  P <- 0.5 * lambda * (beta %*% S %*% beta)
  # Outputting the penalized negative log-likelihood as a scalar
  as.numeric(nll + P)
}


gll <- function(gamma, X, y, S, lambda=5e-5, wb=1) {
  
  # The purpose of this function is to calculate the derivative vector of the PNLL w.r.t. gamma
  
  # The inputs for this function are as follows:
  
    # gamma: Reparameterisation of our spline coefficient vector beta. gamma_k is exp(beta_k),
           # used to ensure that the estimated daily new infection rate f(t) is positive
    # X: Our model matrix for deaths, from "init_model_matrices"
    # y: Our observed deaths (variable "nhs")
    # S: A square matrix S that constructs the penalty term in the final log likelihood function
    # lambda: A smoothing parameter for our penalised negative log likelihood (PNLL), default
            # value of 0.00005
    # wb: A weight parameter used later when "pnll" is called within the non-parametric 
        # bootstrapping function "np_boot", default = 1
  
  # This function outputs the gradient vector of our PNLL
  
  # Converting gamma to beta
  beta <- exp(gamma)
  # Expected number of deaths per day
  mu <- as.vector(X %*% beta)
  # Computing the first term, call it v, in the F matrix computation
  v <- wb * ((y / mu) - 1)
  # Computing the gradient of the PNLL w.r.t. gamma which is the derivative of log-likelihood
  # w.r.t. gamma + derivative of the penalty term w.r.t gamma 
  # We have ordered the computation to make it as efficient as possible, using element-wise
  # multiplication and regular matrix multiplication where appropriate
  grad <- - beta * (t(X) %*% v) + lambda * (beta * (S %*% beta))
  as.vector(grad)
}


fd <- function(gamma, X, y, S, lambda = 5e-5, eps = 1e-7) {
  
  # This function produces a finite differencing estimate of the derivative vector,
  # which we can then compare to the output of our "gll" function.

  # The inputs for this function are as follows:
  
    # gamma: Reparameterisation of our spline coefficient vector beta. gamma_k is exp(beta_k),
           # used to ensure that the estimated daily new infection rate f(t) is positive
    # X: Our model matrix for deaths, from "init_model_matrices"
    # y: Our observed deaths (vatiable "nhs")
    # S: A square matrix S that constructs the penalty term in the final log likelihood function
    # lambda: A smoothing parameter for our penalised negative log likelihood (PNLL), default
            # value of 0.00005
    # wb: A weight parameter used later when "pnll" is called within the non-parametric 
        # bootstrapping function "np_boot", default = 1
    # eps: the finite difference interval that determines the precision of the gradient
  
  # This function outputs the finite differencing estimate of the PNLL gradient vector
  
  # Initialising a zero vector for the finite difference gradient, equal in length to "gamma"
  fd_vec <- numeric(length(gamma))
  # Finite difference uses F(y0+eEi)−F(y0)
  # Calculation for F(y0)
  pnll0 <- pnll(gamma, X, y, S, lambda) # pnll for inputted gamma
  
  # Looping to compute partial derivative for each parameter
  for (i in 1:length(gamma)) {
    # Making a copy so we don't overwrite gamma
    gamma1 <- gamma
    # Increase gamma1[i] by eps
    gamma1[i] <- gamma1[i] + eps
    # Compute resulting PNLL at perturbed Y
    pnll1 <- pnll(gamma1, X, y, S, lambda)
    # Apply finite difference formula: Creating (F(y0+eEi)−F(y0))/ e
    fd_vec[i] <- (pnll1 - pnll0)/eps
  }
  
  fd_vec
}


BIC <- function(lambda, beta_hat, y, X, S) {
  
  # This function evaluates the Bayesian information criterion. It's purpose is to provide
  # a reference for the goodness of fit for a model by penalizing model complexity

  # The inputs for this function are as follows:
  
    # lambda: A smoothing parameter for our penalised negative log likelihood (PNLL), default
            # value of 0.00005
    # beta_hat: The parameter estimated that minimizes the penalized negative log likelihood
    # y: Our observed deaths (variable "nhs")
    # X: Our model matrix for deaths, from "init_model_matrices"
    # S: A square matrix S that constructs the penalty term in the final log likelihood function

  # This function outputs a scalar value that quantifies the goodness of fit for the model.
  # Lower BIC values correspond with a better goodness of fit
  
  # n is the number of data points
  n <- length(y)
  # First, mu must be calculated based on the gamma_hat estimate
  mu <- as.vector(X %*% beta_hat)
  # Compute W as vector for more efficient computation
  W_v <- y/mu^2
  # Hessian without and with smoothing
  H0 <- t(X * W_v) %*% X
  H_lam <- H0 + lambda*S
  # Using Cholesky decomposition to compute inverse efficiently
  R <- chol(H_lam)
  H_lam_inv <- backsolve(R, forwardsolve(t(R),diag(ncol(H_lam))))
  # Computing trace of H_lam_inv %*% H0 efficiently
  EDF <- sum(H_lam_inv * t(H0))
  # Negative log likelihood for Poisson
  nll <- -sum(y * log(mu) - mu) # Removed lgamma(y+1)
  # BIC Equation
  2*nll + log(n) * EDF
}


minimize_BIC <- function(range=seq(-13,-7,length=50), BIC, par, fn, gr, y, X, S) {
  
  # This function takes the BIC, evaluated over a specified range of lambda inputs, and finds the 
  # minimum output value. The function returns the corresponding lambda that yielded that BIC value
  
  # The inputs for this function are as follows:

    # range: The range of log(lambda) values to be tested in the BIC function
    # BIC: Function to be minimized where the minimum output corresponds with the optimum lambda
    # par: The initial vector that optim uses as a starting point for optimization
    # fn: The function to be minimized by optim which helps this function find the corresponding fit 
        # for a given lambda value being tested - not to be confused with the BIC function being 
        # optimized in this function
    # gr: The gradient of the function being optimized by optim - not the gradient of the BIC
    # y: The data points corresponding to deaths on each julian day
    # X: Our model matrix for deaths, from "init_model_matrices"
    # S: A square matrix S that constructs the penalty term in the final log likelihood function

  # This function outputs the optimized lambda value that minimizes the BIC criterion. It also
  # returns a vector of all the BIC values returned, so the BIC curve can be examined
  
  # Exponentiating the range to get the lambda values
  lambda <- exp(range) # Lambda up to this point has been a scalar, but this is a vector
  BICs <- numeric(length(lambda)) # Creating a vector to store the BIC values

  for (i in seq_along(lambda)) { # Looping through the indices of each lambda value
    # Fitting the model with the current lambda value
    fit_i <- optim(par=par, fn=fn, gr=gr, y=y, X=X, S=S, lambda=lambda[i], method='BFGS')
    # Retrieving the appropriate estimate for beta vector by exponentiating the vector of gammas
    beta_hat <- exp(fit_i$par)
    # Calculating the BIC value for the current lambda value
    BICs[i] <- BIC(lambda[i], beta_hat=beta_hat, y=y, X=X, S=S)
  }
  # Returning the lambda value that minimizes the BIC on the chosen range
  return(list(lambda=lambda[which.min(BICs)], BICs=cbind(lambda, BICs)))
}



best_fit <- function(par, fn, gr, y, X, S, lambda) {
  
  # The purpose of this function is to find the optimal gamma vector using the
  # 'BFGS' method of optim. This will allow us to estimate "mu_hat" and "f_hat",
  # our predicted death and infection curves respectively.
  
  # The inputs for this function are as follows:

      # par: The initial gamma values to be optimised over
      # fn: The objective function we wish to minimise (pnll function)
      # gr: A function to return the gradient (gll function)
      # y: The data points corresponding to deaths on each Julian day
      # X: Our model matrix for deaths, from "init_model_matrices"
      # S: A square matrix S that constructs the penalty term in the final log likelihood function
  
  # This function returns a list containing the optimized spline coefficients (beta_hat), 
  # the predicted deaths for each day (mu_hat), and the predicted infections for each day (f_hat)
  
  # Using optim with the 'BFGS' method to optimize the provided function
  fit <- optim(par=par, fn=fn, gr=gr, y=y, X=X, S=S, lambda=lambda, method='BFGS')
  gamma_hat <- fit$par # Retrieving optimised gamma
  beta_hat <- exp(gamma_hat) # Converting gamma estimate to beta estimate
  mu_hat <- as.vector(X_config %*% beta_hat) # Modelling the expected death counts
  f_hat <- as.vector(X_tilde_config %*% beta_hat) # Modelling the estimate for infections
  return(list(beta_hat=beta_hat, mu_hat=mu_hat, f_hat=f_hat))
}


np_boot <- function(nb=200, par, fn, gr, y, X, X_tilde, S, lambda) {

  # The purpose of this function is to perform non-parametric bootstrapping of our
  # "n" day-death pairs in order to produce a 95% confidence interval to reflect the
  # uncertainty in our estimate f_hat.
  
  # The inputs for this function are as follows:
    # nb: The number of bootstrap replicates to be performed, default = 200
    # par: Initial gamma vector supplied to optim
    # fn: Objective function minimized by optim
    # gr: Gradient of the objective function supplied to optim
    # y: Observed daily death counts, as used throughout the model
    # X: Death model matrix returned by init_model_matrices
    # X_tilde: Infection model matrix returned by init_model_matrices
    # S: Penalty matrix returned by init_model_matrices
    # lambda: Smoothing parameter for the penalised likelihood
  
  # The outputs of this function is a list containing 2 vectors: "f_ci_lower" and 
  # "f_ci_upper". These vectors are of equal length and their corresponding entries
  # represent the 95% confidence interval for the daily new infection rate for each
  # day.
  
  # n is the number of data points in the dataset
  n <- length(y)
  # The number of rows in boot_f is the same as the number of predicted values for the infection curve
  nrows_boot_f <- nrow(X_tilde)
  # initialized empty matrix to store bootstraps
  boot_f <- matrix(0, nrow = nrows_boot_f, ncol = nb)

  for (i in 1:nb) { # bootstrap loop
    # non-parametric bootstrap sample, tabulating to count the occurrences of each day-death pair
    # Bootstrap resampling is equivalent to re-weighting the terms in the log-likelihood. These are our weights
    wb <- tabulate(sample(n,replace=TRUE),n)
    # Re-optimising gamma, but this time using the newly calculated weights
    opt_boot <- optim(par=gamma_config, fn=fn, gr=gr, y=y, wb=wb, X=X, S=S, lambda=lambda, method="BFGS")
    # Converting optimal gamma values to betas
    beta_hat_boot <- exp(opt_boot$par)

    # Computing the predicted values for all days using the new betas from this bootstrap
    # The matrix boot_f has one column per bootstrap replicate, with each row representing an infection day
    boot_f[, i] <- X_tilde_config %*% beta_hat_boot
  }

  # Taking the 2.5th percentile and 97.5th percentile of each of those 200 bootstrap resamples 
  # to get the empirical 95% confidence interval. MARGIN = 1 means we apply the function row by row.
  f_ci_lower <- apply(boot_f, 1, function(x) quantile(x, 0.025))
  f_ci_upper <- apply(boot_f, 1, function(x) quantile(x, 0.975))
  
  return(list(f_ci_lower=f_ci_lower, f_ci_upper=f_ci_upper))
}


##########################################################
### Importing Data and Configuring Function Parameters ###
##########################################################

########################## Dataset Config #########################
# Project is designed to be configurable, making it easier to 
# adjust for different datsets with minimal code augmentation

# Insert name of dataset as string
dataset_name <- "engcov.txt"
# Insert concatenated vector in form (predictor, response)
cols <- c("julian","nhs")

# Saved dataset to global environment
data <- read.table(dataset_name, header = TRUE)
filtered_data <- data[,cols]

t_config <- as.vector(data[[cols[1]]]) # Predictor Vector
y_config <- as.vector(data[[cols[2]]]) # Response Vector

########################## Parameter Config ########################

# Number of splines to use in the model
K_config <- 80
# Minimum window of days for evaluation
min_window_config <- 30
# Maximum windows of days for evaluation
max_window_config <- 80
# Expected log of days from infected to death
edur_config <- 3.151
# Standard deviation of log of days from infected to death
sdur_config <- .469

# Automatically sets spline range so it's min_window_config days longer than the dataset
spline_range_config <- seq(
  from = min(t_config) - min_window_config,
  to = max(t_config)
)

# Initial vector values to test gradient function and optimize model
gamma_config <- rep(0, K_config)


#######################################################################
### Initializing Model Matrices and Saving Them as Global Variables ###
#######################################################################

# Model matrices are initialized based on the parameter config above
model_matrices <- init_model_matrices(K=K_config, 
                                      t=t_config, 
                                      spline_range=spline_range_config, 
                                      min_window=min_window_config,
                                      max_window=max_window_config,
                                      edur=edur_config, 
                                      sdur=sdur_config)

# Saving model matrices to global environment
X_tilde_config <- model_matrices$X_tilde
X_config <- model_matrices$X
S_config <- model_matrices$S


###########################################################
### Testing Gradient Function Using Finite Differencing ###
###########################################################

# Comparing gradient function to finite differencing algorithm
fd_out <- fd(gamma_config, X=X_config, y=y_config, S=S_config)
grad_out <- gll(gamma_config, X=X_config, y=y_config, S=S_config)
# If the values of fd_out and grad_out are very similar, the computed gradient by "gll" is accurate
head(cbind(fd_out, grad_out))


#############################################################################################
### Finding Best Fit for the Smooth Deconvolution Model and Computing Confidence Interval ###
#############################################################################################

# Initializing the model with an arbitrary lambda value
M_0 <- best_fit(par=gamma_config, fn=pnll, gr=gll, y=y_config, X=X_config, S=S_config, lambda=5e-5)

# Finding optimal lambda value by minimizing BIC
gamma_config <- log(M_0$beta_hat) # Changing gamma_config to our new estimate since convergence will be reached quicker
search_range <- seq(-13,-7,length=50) # Each value in this grid search range is the log of lambda
# This function below minimizes the BIC and returns both the corresponding lambda and the BIC values for each lambda
min_BIC <- minimize_BIC(range=search_range, BIC=BIC, par=gamma_config, fn=pnll, gr=gll, y=y_config, X=X_config, S=S_config)

lambda_hat <- min_BIC$lambda # lambda value that minimized the BIC

# Model with lambda value that minimizes BIC
M_1 <- best_fit(par=gamma_config, fn=pnll, gr=gll, y=y_config, X=X_config, S=S_config, lambda=lambda_hat)

# Computing confidence interval about the entire infection curve (f) using a non-parametric bootstrap
bootstrap <- np_boot(nb=200, par=gamma_config, fn=pnll, gr=gll, y=y_config, X=X_config, 
                     X_tilde=X_tilde_config, S=S_config, lambda=lambda_hat)


##################################################################
### Initial Plotting - Observing Adjust Lambda's Effect on Fit ###
##################################################################

library(ggplot2)

arb_lambda <- ggplot(data = NULL, aes(x = t_config, y = data.frame(M_1$f_hat))) +
  # Plotting the data points for deaths
  geom_point(data = data.frame(t = t_config, 
                               y = y_config), 
             aes(x = t, y = y, color = "Observed Deaths"), size = 1.5) + 
  # Adding the fitted death curve
  geom_line(data = data.frame(t_config = t_config, 
                              mu_hat = M_0$mu_hat), 
            aes(x = t_config, y = mu_hat, color = 'Fitted Death Curve'), linewidth = 1.1) +
  # Adding the fitted infection curve
  geom_line(data = data.frame(spline_range_config = spline_range_config,  
                              f_hat = M_0$f_hat), 
            aes(x = spline_range_config, 
                y = f_hat, color = 'Fitted Infection Curve'), linewidth = 1.1) +
  # Creating legend
  scale_color_manual(
    name = "Legend",
    values = c(
      "Observed Deaths" = "darkgreen",
      "Fitted Death Curve" = "black",
      "Fitted Infection Curve" = "blue")
  ) +
  labs(x = "Day of year",y = "Counts", title = expression("Daily Deaths, Fitted Death Curve and Predicted Infections - using arbitrary " * lambda)) +
  
  theme_minimal()

arb_lambda



#####################################################################################
### Final Plotting - Using Optimized Lambda and Overlaying the Confidence interval ###
#####################################################################################



Op_with_CI <- ggplot(data = NULL, aes(x = t_config, y = data.frame(M_1$f_hat))) + 
  # Plotting the data points for deaths
  geom_point(data = data.frame(t = t_config, 
                               y = y_config), 
             aes(x = t, y = y, color = "Observed Deaths"), size = 1.5) + 
  # Adding the fitted death curve
  geom_line(data = data.frame(t_config = t_config, 
                              mu_hat = M_1$mu_hat), 
            aes(x = t_config, y = mu_hat, color = 'Fitted Death Curve'), linewidth = 1.1) +
  # Adding the fitted infection curve
  geom_line(data = data.frame(spline_range_config = spline_range_config, 
                              f_hat = M_1$f_hat), aes(x = spline_range_config, 
                                                      y = f_hat, color = 'Fitted Infection Curve'), linewidth = 1.1) + 
  # Adding the upper and lower CI
  geom_line(data = data.frame(spline = spline_range_config,
                              lower = unlist(bootstrap[[1]])), 
            aes(x = spline_range_config, y = lower, color = 'CI Bounds'), 
            linetype = "longdash", linewidth = .75, alpha = 0.75) +
  geom_line(data = data.frame(spline = spline_range_config, 
                              upper = unlist(bootstrap[[2]])), 
            aes(x = spline_range_config, y = upper), color = 'red', linetype = "longdash", linewidth = .75, alpha = 0.75) +
  # Creating legend
  scale_color_manual(
    name = "Legend",
    values = c(
      "Observed Deaths" = "darkgreen",
      "Fitted Death Curve" = "black",
      "Fitted Infection Curve" = "blue",
      "CI Bounds" = "red")
  ) +
  labs(x = "Day of year",y = "Counts", 
       title = expression("Daily Deaths, Fitted Death Curve, Predicted Infections (with CI) - using optimised " * lambda)) +
  
  theme_minimal()

Op_with_CI 
