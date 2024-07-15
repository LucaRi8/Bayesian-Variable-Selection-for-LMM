# Metropolis within gibbs sampler for spike and slab priors in linear mixed model

library(Matrix)
library(emulator)
library(fastDummies)
library(mvnfast)
library(tidyverse)

# computation of gamma distribution (given sigma_u, sigma_e, u)
pi_gamma <- function(gamma, sigma_e, theta, Z, data, response, v, u)
{
  # v parameter of the prior variance (in D0)
  # n and p are sample size and number of parameter
  # assume prior beta0 = 0 
  # the intercept is always included
  y <- as.matrix((data[, response]))
  data <- data %>% dplyr::select(-response)
  p <- length(gamma)
  p_u <- ncol(Z)
  n <- nrow(data)
  # model-matrix construction
  pos_gamma <- which(gamma == 1)
  data <- data.frame("Intercept" = 1, data[, pos_gamma])
  nonnumeric_feature <- names(which(!sapply(data, is.numeric)))
  if(ncol(data) == 1){
    X <- as.matrix(data)
  } else if(length(nonnumeric_feature) != 0) {
    X <- as.matrix(dummy_cols(data, select_columns = nonnumeric_feature, 
                              remove_selected_columns = TRUE, 
                              remove_first_dummy = T))
  } else X <- as.matrix(data)
  
  # distribution calculation
  D <- diag(ncol(X)*v)
  D_inv <- diag(ncol(X)*(1/v))
  XtX <- crossprod(X) + D_inv
  XtX_inv <- chol2inv(chol(XtX))
  epsilon_u <- y - Z %*% u
  S <- quad.tform(XtX_inv, as.vector(crossprod(epsilon_u, X)))
  
  m_like <- -0.5*determinant(D, logarithm = T)$modulus + 
    sum(gamma)*log(theta) + (p-sum(gamma))*log(1-theta) +
    0.5*determinant(XtX_inv, logarithm = T)$modulus + 1/(2*sigma_e) * S -
    1/(2*sigma_e) * crossprod(epsilon_u)
  
  return(list(mlike = m_like, XtX = XtX, D = D, D_inv = D_inv, 
              XtX_inv = XtX_inv, epsilon_u = epsilon_u, X = X, gamma = gamma))
}


# gibbs sampler
gibbs_spikeslab <- function(R, beta_init, sigma_init, u_init, gamma_init, data, Z, response, theta,
                            ni_e, ni_u, lambda_e, lambda_u, v)
{
  y <- as.matrix((data[, response]))
  p_b <- length(gamma_init)
  p_u <- ncol(Z)
  p <- length(gamma_init)
  n <- nrow(data)
  par_init <- c(u_init, sigma_init)
  # matrix in which store the simulated quantities
  parout <- matrix(NA, R+1, length(par_init))
  gammaout <- matrix(NA, R+1, length(gamma_init))
  beta_list <- list()
  # put the starting points
  gammaout[1, ] <- gamma_init
  parout[1, ] <- par_init
  beta_list[[1]] <- beta_init
  # gamma distribution at starting points
  m_0 <- pi_gamma(gamma_init, sigma_e = par_init[length(par_init)-1], 
                  theta = theta, Z, data, response, v, par_init[1:p_u])
  gamma0 <- gamma_init 
  
  ZtZ <- crossprod(Z)
  # start the metropolis within gibbs sampler
  for(i in 1:R)
  {
    # propose a new gammma 
    gamma1 <- gamma0
    pos <- sample(1:p, size = 1)
    ifelse(gamma1[pos]==0, gamma1[pos] <- 1, gamma1[pos] <- 0)
  
    m_1 <- pi_gamma(gamma1, sigma_e = parout[i, length(par_init)-1], 
                    theta = theta, Z, data, response, v, parout[i, 1:p_u])
    m_0 <-  pi_gamma(gamma0, sigma_e = parout[i, length(par_init)-1], 
                     theta = theta, Z, data, response, v, parout[i, 1:p_u])
    # metropolis acceptance rate
    alpha <- exp(m_1$mlike - m_0$mlike) 
    if(as.numeric(runif(1) < alpha))
    {
      gamma0 <- gamma1
      m_0 <- m_1
    }
    gammaout[i+1, ] <- gamma0
    
    # sample beta parameters
    sigma_e <- parout[i, length(par_init)-1]
    epsilon_u <- y - Z %*% parout[i, 1:p_u]
    beta_list[[i+1]] <- mvnfast::rmvn(1, mu = m_0$XtX_inv %*% crossprod(m_0$X, epsilon_u),
                                              sigma = m_0$XtX_inv * sigma_e)
    colnames(beta_list[[i+1]]) <- colnames(m_0$X)
    
    # sample u parameters
    sigma_u <- parout[i, length(par_init)]
    B <- Matrix(diag(ncol(ZtZ)) * (sigma_e / sigma_u))
    epsilon_b <- y - m_0$X %*% as.vector(beta_list[[i+1]])
    parout[i+1, 1:p_u] <- mvnfast::rmvn(1, mu = chol2inv(chol(ZtZ + B)) %*% crossprod(Z, epsilon_b),
                                  sigma = sigma_e * chol2inv(chol(ZtZ + B)))
    
    # sample sigma_e parameters 
    epsilon_tot <- epsilon_b - Z %*% parout[i+1, 1:p_u]
    beta_cor <- as.vector(beta_list[[i+1]])
    parout[i+1, length(par_init)-1] <- 1/rgamma(1, (n/2 + ni_e),
                                                drop(crossprod(epsilon_tot)/2 + 
                                                   quad.form(m_0$D_inv, beta_cor)/2 + 
                                                   lambda_e))
    # sample sigma_u parameters
    u_cor <- parout[i+1, 1:p_u]
    parout[i+1, length(par_init)] <- 1/rgamma(1, ni_u + p_u/2, (crossprod(u_cor)/2 + lambda_u))
    
    print(i)
  }
  return(list(beta = beta_list, parameters = parout, gamma = gammaout))
}
