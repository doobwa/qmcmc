library(mvtnorm)
library(MCMCpack)

# Set up example of hierarchical normal data
# y_ij ~ N_2(theta_j,1)
# theta_jp ~ N(mu_p,sigma_p^2)
# sigma_p ~ InvGamma(alpha,beta)


# Explore prior on sigma^2
xs <- seq(.01,10,by=.1)
plot(xs,dinvgamma(xs,2,2),type="l")

set.seed(1)
mu <- c(0,0)
sigma <- c(3,3)
J <- 10
theta <- t(sapply(1:J,function(j) {
  rmvnorm(1,mu,diag(sigma))
}))
priors <- list(mu=list(mu=0,sigma=3),sigma=list(alpha=2,beta=2))

# Generate data
n <- rep(10,J)                
ys <- lapply(1:J,function(j) {
  rmvnorm(n[j],theta[j,],diag(2))
})

value <- truth <- list(theta=theta,mu=mu,sigma=sigma)

# Compare true values to empirical estimates
colMeans(ys[[5]])
truth$theta[5,]

test_that("can gibbs sample mu and sigma",{
  gibbs.mu.hier.gaussian(value,priors)
  gibbs.sigma.hier.gaussian(value,priors)
})

mse <- function(a,b) mean((a-b)^2)
niter <- 100

test_that("gibbs sampling mu works",{
  value <- truth
  values <- list()
  for (iter in 1:niter) {
    value <- gibbs.mu.hier.gaussian(value,priors)
    values[[iter]] <- value
  }
  mu.samples <- do.call(rbind,lapply(values,function(v) v$mu))
                                        # should be closer to 0
  err <- mse(colMeans(mu.samples),colMeans(truth$theta))
  expect_that(err < .01, is_true())
})

test_that("gibbs sampling sigma works",{
  value <- truth
  values <- list()
  for (iter in 1:niter) {
    value <- gibbs.sigma.hier.gaussian(value,priors)
    values[[iter]] <- value
  }
  sigma.samples <- do.call(rbind,lapply(values,function(v) v$sigma))
  colMeans(sigma.samples)
  apply(truth$theta,2,sd)^2
})

# Make sure basic functions run
loglikelihood.hier.gaussian(ys[[1]],theta[1,],grad=FALSE)
lposterior.hier.gaussian(ys[[1]],theta[1,],mu)
lprior.hier.gaussian(theta[1,],mu,sigma,priors,grad=FALSE)

# Compare prior on theta with sigma and with sigma integrated out
par(mfrow=c(1,2))
thetas <- seq(-3,3,by=.01)
lps <- sapply(thetas,function(x) {
  lprior.hier.gaussian(c(x,theta[1,2]),mu,sigma,priors)
})
plot(thetas,lps,type="l")

lps <- sapply(thetas,function(x) {
  lprior.nosigma.hier.gaussian(c(x,theta[1,2]),mu,priors)
})
plot(thetas,lps,type="l",col="red")

# Sample just a single lower level theta
j <- 1
cs <- TRUE
value <- truth
values <- list()
method <- slice
for (iter in 1:50) {
  lp <- function(val) {
    value$theta[j,] <- val
    lposterior.hier.gaussian(ys[[j]],value$theta[j,],value$mu,value$sigma,priors,collapse.sigma=cs)
  }
  current <- value$theta[j,]
  value$theta[j,] <- method(current, lp)
  values[[iter]] <- value$theta[j,,drop=FALSE]
}
values <- do.call(rbind,values)
plot(values)
points(truth$theta[j,,drop=FALSE],col="red")

# Plot sampled values to true values
par(mfrow=c(1,2))
plot(values[,1],type="l")
abline(h=truth$theta[j,1],col="red")
plot(values[,2],type="l")
abline(h=truth$theta[j,2],col="red")

# Sample all lower level
value <- truth
values <- list()
method <- slice
for (iter in 1:50) {
  lower <- lapply(1:J,function(j) {
    lp <- function(val) {
      value$theta[j,] <- val
      lposterior.hier.gaussian(ys[[j]],value$theta[j,],value$mu,value$sigma,priors)
    }
    current <- value$theta[j,]
    method(current, lp)
  })
  value$theta <- do.call(rbind,lower)#t(sapply(lower,function(x) x$final$pv))
  values[[iter]] <- value$theta
}
values <- melt(values)

th <- melt(truth$theta)
qplot(L1,value,data=values,geom="line") + geom_hline(data=th,aes(yintercept=value),colour="red") + facet_grid(X1~X2)

# Fit all at once
fit <- mcmc.hier.gaussian(ys,slice,priors=priors)
sigmas <- do.call(rbind,lapply(fit$values,function(f) f$sigma))
plot(sigmas)

lposterior.all(ys,fit$value,priors) # TODO: broken
lposterior.all(ys,truth,priors)

