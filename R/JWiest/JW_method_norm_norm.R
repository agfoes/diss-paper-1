library(rjags)
library(coda)

jw_model <- "
  data {
    for (i in 1:n) {
      d[i] <- 1
    }
  }
  
  model {
    for (i in 1:n) {
      y[i] ~ dbin(py[i], m)
      d[i] ~ dinterval(z[i], cuts[i, 1:2])
      z[i] ~ dgamma(shape, shape/ez[i])
      logit(py[i]) <- inprod(x[i, 1:p], beta[1:p]) + delt * z[i]
      log(ez[i]) <- inprod(x[i, 1:p], alph[1:p])
    }
    
    shape ~ dgamma(0.001, 0.001)
    delt ~ dnorm(0, 0.001)
    for (j in 1:p) {
      beta[j] ~ dnorm(0, 0.001)
      alph[j] ~ dnorm(0, 0.001)
    }
  }
"

set.seed(123)

n <- 100
m <- 50

# true parameters
beta0 <- 0
beta1 <- 1
delt <- 1

gamma0 <- 0
gamma1 <- 1
shape <- 0.25

# fully observed covariate
x <- runif(n, -1, 1)

# coarsened covariate
ez <- exp(gamma0 + gamma1 * x)
z_true <- rgamma(n, shape = shape, rate = shape / ez)

# outcome
eta <- beta0 + beta1 * x + delt * z_true
py <- plogis(eta)

y <- rbinom(n, size = m, prob = py)

eps <- 0.01
upper_inf <- 1e4

## 50% 
observed <- rep(TRUE, n)
idx_cen <- sample(
  seq_len(n),
  size = round(0.5*n),
  replace = FALSE
)

observed[idx_cen] <- FALSE
q <- median(z_true)

cuts <- matrix(NA, nrow = n, ncol = 2)
z_init <- numeric(n)
for (i in seq_len(n)) {
  if (observed[i]) {
    cuts[i, ] <- c(max(0, z_true[i] - eps), z_true[i] + eps)
  } else if (z_true[i] <= q) {
    cuts[i, ] <- c(0, q)
  } else {
    cuts[i, ] <- c(q, upper_inf)
  }
  
  lower <- cuts[i, 1]
  upper <- cuts[i, 2]
  if (upper == upper_inf) {
    z_init[i] <- lower + 1
  } else {
    z_init[i] <- (lower + upper) / 2
  }
  
}

# model matrix
X <- cbind(1, x)
mydata <- list(y = y,
               x = X,
               cuts = cuts,
               m = m,
               n = n,
               p = ncol(X))
myinit <- list(z = z_init)
mymodel <- jags.model(textConnection(jw_model),
                      data = mydata,
                      inits = myinit,
                      n.chains = 1)

update(mymodel, n.iter = 5000)
mysample <- coda.samples(mymodel,
                         variable.names = c("beta", "delt", "alph", "shape"),
                         n.iter = 50000)
summary(mysample)
effectiveSize(mysample)
autocorr.diag(mysample)
plot(mysample)





