# RJMCMC for bayesian model selection
library(Matrix)
library(emulator)
library(fastDummies)

# compute the marginal likelihood of gamma given y and X
m_like <- function(gamma, y, modelMAT, theta, lambda, ni, v, pos_list, n, p) {
  # v parameter of the prior variance (in D0)
  # n and p are sample size and number of parameter
  # assume prior beta0 = 0 e D0 = diag(p)
  # the intercept is always included
  if(colnames(modelMAT)[1] == "(Intercept)" )
  {
    pos_gamma <- c(1, which(gamma == 1)+1)
  } else pos_gamma <- which(gamma == 1)
  pos_in_data <- do.call(c, lapply(pos_gamma, function(i) pos_list[[i]]))
  X <- Matrix(modelMAT[, pos_in_data])
  D0 <- diag(ncol(X)) * v
  yty <- crossprod(y)
  XtX <- crossprod(X) + diag(rep(1/v, ncol(X)))
  XtX_inv <-chol2inv(chol(XtX))
  S <- yty -  crossprod(y, X) %*% XtX_inv %*% t(crossprod(y, X))
  
  gamma_post <- 0.5* determinant(XtX_inv, logarithm = T)$modulus -
    (n/2 + ni) * log(lambda + S/2) -
    (1/2)*log(v**sum(gamma)) + dbinom(sum(gamma), size = p, prob = theta, log = T)
  
  return(gamma_post)
}


# RJ MCMC function
RJ_MCMC <- function(R, y, data, lambda, ni, theta = 0.5, gamma_init, v = 1, pos_list, verbose) {
  # lambda, ni iperparameter for sigma
  # zeta, phi iperparameter for theta
  # theta prob for gamma distribution
  n <- nrow(data)
  p <- length(gamma_init)
  gamma_sim <- matrix(NA, R+1, length(gamma_init)) # store simulated value
  m_sim <- numeric(R+1) # store m(gamma| X, y)
  # init marginal likelihood
  m_0 <- m_like(gamma_init, y, data, theta, lambda, ni, v, pos_list, n, p)
  # store the inizialisation parameter
  gamma_sim[1, ] <- gamma_init
  m_sim[1] <- m_0
  gamma0 <- gamma_init 
  for(i in 1:R){
    gamma1 <- gamma0
    pos <- sample(1:p, size = 1)
    ifelse(gamma1[pos]==0, gamma1[pos] <- 1, gamma1[pos] <- 0)
    m_1 <- m_like(gamma1, y, data, theta, lambda, ni, v, pos_list, n, p)
    alpha <- exp(m_1 - m_0)
    if(as.numeric(runif(1) < alpha))
    {
      gamma0 <- gamma1
      m_0 <- m_1
    }
    gamma_sim[i+1, ] <- gamma0
    m_sim[i+1] <- m_0
    if(verbose == T) {
      if(i %% 100 == 0) print(i)
    }
  }
  return(list(gamma = gamma_sim, lm = m_sim))
}