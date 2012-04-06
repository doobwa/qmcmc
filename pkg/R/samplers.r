##' General purpose Metropolis-Hastings sampler
##' @param current
##' @param likelihood
##' @param olp log density of current
##' @param sd  standard deviation for Gaussian proposal density
##' @return new state and its log density as an attribute

mh <- function(current,lposterior,olp=NULL,sd=.1) {
  if (is.null(olp)) {
    olp <- attr(current,"log.density") <- lposterior(current)
  }
  cand <- current + rnorm(length(current),0,sd)
  clp <- lposterior(cand)
  if (clp - olp > log(runif(1))) {
    current <- cand
    attr(current,"log.density") <- clp
  }
  return(current)
}

##' General purpose slice sampler (Radford Neal's code)
##' @param current parameter vector
##' @param lposterior function returning the log posterior
##' @param olp log density of current
##' @param sd  standard deviation for Gaussian proposal density
##' @return new state and its log density as an attribute
slice <- function(current,lposterior,olp=NULL,m=20) {
  if (is.null(olp)) {
    olp <- attr(current,"log.density") <- lposterior(current)
  }
  lpost <- function(val) {
    current[p] <- val
    lposterior(current)
  }
  for (p in 1:length(current)) {
    val <- uni.slice.alt(current[p],lpost,gx0=olp,m=m)
    current[p] <- val
    olp <- attr(val,"log.density")
  }
  attr(current,"log.density") <- olp
  return(current)
}

##' General purpose HMC (Radford Neal's code)
##' @param current parameter vector
##' @param lposterior function returning the log posterior
##' @param epsilon step size
##' @param L number of steps
##' 
## Changed U and grad_U to not be negative llk, etc.
## Added log.density attributes
hmc <- function (current_q, lposterior, epsilon, L) 
{
  U      <- function(x) -lposterior(x)
  grad_U <- function(x) -lposterior(x,lp=FALSE)
    q = current_q
    p = rnorm(length(q), 0, 1)
    current_p = p
    p = p - epsilon * grad_U(q)/2
    for (i in 1:L) {
        q = q + epsilon * p
        if (i != L) 
            p = p - epsilon * grad_U(q)
    }
    p = p - epsilon * grad_U(q)/2
    p = -p
    current_U = U(current_q)
    current_K = sum(current_p^2)/2
    proposed_U = U(q)
    proposed_K = sum(p^2)/2
    attr(current_q,"log.density") <- -current_U
    attr(q,"log.density") <- -proposed_U
    if (runif(1) < exp(current_U - proposed_U + current_K - proposed_K)) {
        return(q)
    }
    else {
        return(current_q)
    }  
}
