profile.ci.for.icc.sem <- function(mod, level = 0.95){
  mle.dev <- devfun2(mod)

  mle.pars <- getME(mod, 'theta')
  mle.pars <- c(mle.pars*sigma(mod), sigma(mod))^2

  delta.mle <- 2*sum(mle.pars[2:3])

  delta <- delta.mle
  mle.dev <- devfun2(mod)

  mle.dev.SqSc <- function(par){
    mle <- mle.dev(sqrt(par))

    return(mle)
  }

  profile.sem <- function(delta){
    heq <- function(par){
      return(par[2] + par[3] - delta/2)
    }

    heq.jac <- function(par){
      return(matrix(c(0, 1, 1), nrow = 1))
    }

    hin <- function(par){
      return(c(par[1], par[2], par[3]))
    }

    hin.jac <- function(par){
      return(matrix(c(1, 0, 0,
                      0, 1, 0,
                      0, 0, 1), byrow = TRUE))
    }

    constrOptim.nl(par = mle.pars, fn = mle.dev.SqSc, heq = heq, heq.jac = heq.jac, hin = hin, hin.jac = hin.jac, control.outer = list(trace = 0))$value
  }

  profile.sem <- Vectorize(profile.sem)

  deltas <- c(delta.mle)
  zs <- c(0)

  z.max <- qnorm(0.99)
  D.z <- z.max / 8

  deltas <- c(deltas, 1.01*deltas)
  zs <- c(zs, sqrt(profile.sem(tail(deltas, 1)) - deviance(mod)))

  i <- 2
  while (tail(zs, 1) < z.max){
    DzDpsi <- (zs[i] - zs[i - 1])/(deltas[i] - deltas[i-1])

    dpsi <- D.z / DzDpsi

    deltas <- c(deltas, tail(deltas, 1) + dpsi)
    zs <- c(zs, sqrt(profile.sem(tail(deltas, 1)) - deviance(mod)))

    i <- i + 1
  }

  i <- 2
  while (zs[1] > -z.max){
    DzDpsi <- (zs[i-1] - zs[i])/(deltas[i-1] - deltas[i])

    dpsi <- D.z / DzDpsi

    deltas <- c(deltas[1] - dpsi, deltas)
    zs <- c(-sqrt(profile.sem(deltas[1]) - deviance(mod)), zs)

    i <- i + 1
  }

  cc <- pchisq(zs^2, 1)
  cc.spline <- splinefun(sqrt(deltas), cc)

  cd <- pnorm(zs)
  cd.spline <- splinefun(sqrt(deltas), cd)

  cq.spline <- splinefun(cd, sqrt(deltas))

  alpha = 1 - level

  return(list(ci = cq.spline(c(alpha/2, 1-alpha/2)), delta.mle = sqrt(delta.mle), delta = sqrt(deltas), cc = cc.spline, cd = cd.spline, cq = cq.spline))
}
