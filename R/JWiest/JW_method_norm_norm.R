library(rjags)
library(coda)

jw_model <- "{
  data {
    dimx <- dim(x)
    n <- dim1[1]
    p <- dimx[2]
    
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
}"


