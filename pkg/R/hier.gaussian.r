# The object of parameters must have

# Likelihood and gradient for 

##' Example of a function computing the log likelihood.
##' @param y
##' @param priors 
##' @param lp compute the log likelihood
##' @param lgrad compute the log gradient
##' @return log posterior or log gradient

loglikelihood <- function(y,theta,grad=FALSE) {
  llk <- sum(apply(y,1,function(x) dmvnorm(x,theta,diag(2),log=TRUE)))
  if (grad) {
    attr(llk,"grad") <- lgradient(y,theta)
  }
  return(llk)
}

##' Log prior for theta_j (i.e. for a single case j)
lprior <- function(theta,mu,priors,grad=FALSE) {
  a <- priors$sigma$alpha
  b <- priors$sigma$beta
  lprior <- -.5 * log(2*pi) + a*log(b) - log(gamma(a)) + log(gamma(a+.5)) - 
    (a+.5)*log(.5*(theta - mu)^2 + b)
  lprior <- sum(lprior)
  if (grad) {
    attr(lprior,"grad") <- -2*(a+.5)*(theta - mu)/((theta-mu)^2 + 2*b)
  }
  return(lprior)
}

lgradient <- function(y,theta,priors) {

}


lposterior <- function(y,theta,mu,priors=list(theta=list(mu=0,sigma=1),sigma=list(alpha=2,beta=.1)),grad=FALSE){
  llk <- loglikelihood(y,theta,grad)
  lp <- lprior(theta,mu,priors,grad)
  lpost <- llk+lp
  if (grad) {
    attr(lpost,"grad") <- attr(llk,"grad") + attr(lp,"grad")
  }
  return(lpost)
}

lposterior.all <- function(ys,value,priors) {
  J <- length(ys)
  sum(sapply(1:J,function(j) {
    lposterior(ys[[j]],value$theta[j,],value$mu,priors)
  })) +
    sum(dgamma(value$sigma,priors$sigma$alpha,priors$sigma$beta,log=TRUE))
}

##' Example of an mcmc function
##' @param ys list of vectors
##' @param method mcmc method to use
##' @param param 
##' @param beta current
##' @return list with parameters at each iteration and a vector of log posteriors at each iteration

mcmc.hier.gaussian <- function(ys,method,theta=NULL,niter=100,burnin=10,priors=list(mu=list(mu=0,sigma=1),sigma=list(alpha=2,beta=.1)),verbose=TRUE,...) {
  
  value <- list()
  J <- length(ys)
  value$mu <- rnorm(2,priors$mu$mu,priors$mu$sigma)
  value$sigma <- rgamma(2,priors$sigma$alpha,priors$sigma$beta)
  value$theta <- t(sapply(1:J,function(j) {
    c(rnorm(1,value$mu[1],value$sigma[1]),
      rnorm(1,value$mu[2],value$sigma[2]))
  }))
  
  values <- list()
  P <- ncol(ys[[1]])
  lps <- rep(0,niter)
  for (iter in 1:niter) {
    
    lps[iter] <- lposterior.all(ys,value,priors)
    
    # Gibbs sample upper parameters
    value <- gibbs.sigma.hier.gaussian(value,priors)
    if (verbose) print(head(value$mu))
    
    value <- gibbs.mu.hier.gaussian(value,priors)
    if (verbose) print(head(value$sigma))
    
    # MCMC step on lower parameters
    lower <- lapply(1:J,function(j) {
      lp <- function(val) {
        value$theta[j,] <- val
        lposterior(ys[[j]],value$theta[j,],value$mu,priors)
      }
      current <- value$theta[j,]
      method(current, lp)
    })
    value$theta <- do.call(rbind,lower)
    
    # Save progress
    if (iter > burnin) {
      values[[iter]] <- value
    }
    if (verbose) {
      cat(lps[iter],"\n")
    }
  }
  return(list(value=value,lps=lps,values=values))
}


rinvchisq <- function(n=1,df,scale=1) {
  df * scale / rchisq(n,df)
}

rinvgamma <- function(n,alpha,beta) {  # using beta=10
  1/rgamma(n,alpha,beta)
}

gibbs.mu.hier.gaussian <- function(value,priors) {
  P <- ncol(value$theta)
  J <- nrow(value$theta)
  value$mu <- rnorm(P,colMeans(value$theta),value$sigma/sqrt(J))
  return(value) 
}

gibbs.sigma.hier.gaussian <- function(value,priors) {
  P <- ncol(value$theta)
  J <- nrow(value$theta)
  a <- priors$sigma$alpha
  b <- priors$sigma$beta
  mu.p <- matrix(value$mu,nr=J,nc=P,byrow=TRUE)
  s2  <- .5*colSums((value$theta - mu.p)^2)
  value$sigma <- sqrt(rinvgamma(P,a + P/2, b + s2))
  return(value)
}

