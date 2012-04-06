library(mvtnorm)

# Set up example of hierarchical normal data
# y_ij ~ N_2(theta_j,1)
# theta_jp ~ N(mu_p,sigma_p^2)
# sigma_p ~ InvGamma(alpha,beta)

set.seed(1)
mu <- c(0,0)
sigma <- c(3,3)
J <- 5
theta <- t(sapply(1:J,function(j) {
  rmvnorm(1,mu,diag(sigma))
}))
priors <- list(mu=list(mu=0,sigma=3),sigma=list(alpha=2,beta=.1))

n <- rep(4,J)                
ys <- lapply(1:J,function(j) {
  rmvnorm(n[j],theta[1,],diag(2))
})

value <- list(theta=theta,mu=mu,sigma=sigma)

test_that("can gibbs sample mu and sigma",{
  gibbs.mu.hier.gaussian(value,priors)
  gibbs.sigma.hier.gaussian(value,priors)
})

loglikelihood(ys[[1]],theta[1,],grad=FALSE)
lposterior(ys[[1]],theta[1,],mu)
lprior(theta[1,],mu,priors,grad=TRUE)

method <- slice
lower <- lapply(1:J,function(j) {
  lp <- function(val) {
    value$theta[j,] <- val
    lposterior(ys[[j]],value$theta[j,],value$mu,priors)
  }
  current <- theta[j,]
  method(current, lp)
})
lower <- do.call(rbind,lower)#t(sapply(lower,function(x) x$final$pv))

lposterior.all(ys,value,priors)

test_that("can sample lower level",{

  value$pv <- theta[1,]
  lprior(value,grad=TRUE)

})

fit <- mcmc.hier.gaussian(ys,slice)
lposterior.all(ys,fit$value,priors)
