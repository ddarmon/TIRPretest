compute.reliability.measures <- function(dat, level = 0.95, B = 500){

  X <- matrix(dat$measure, nrow = length(unique(dat$id)))

  ICC.out <- psych::ICC(X)

  icc.11 <- as.numeric(c(ICC.out$results[1, 2], ICC.out$results[1, 7:8]))
  icc.21 <- as.numeric(c(ICC.out$results[2, 2], ICC.out$results[2, 7:8]))

  mod.icc11 <- lmer(measure ~ 1 + (1 | id), data = dat)

  confint.icc11.prof <- confint(mod.icc11, level = level, method = 'profile', quiet = TRUE)
  confint.icc11.boot <- confint(mod.icc11, level = level, method = 'boot', nsim = B, quiet = TRUE)

  ci.icc11.sem.prof <- sqrt(2)*confint.icc11.prof[2, ]
  ci.icc11.mdd.prof <- qnorm(0.975)*sqrt(2)*confint.icc11.prof[2, ]

  ci.icc11.sem.boot <- sqrt(2)*confint.icc11.boot[2, ]
  ci.icc11.mdd.boot <- qnorm(0.975)*sqrt(2)*confint.icc11.boot[2, ]

  sig.icc11 <- sigma(mod.icc11)
  psi <- unname(sig.icc11 * getME(mod.icc11, "theta"))

  icc11.sem <- sqrt(2)*sig.icc11
  icc11.mdd <- icc11.sem*qnorm(0.975)

  results <- matrix(NA, nrow = 5, ncol = 3)

  results[1, ] <- icc.11

  results[2, 1] <- icc11.sem
  results[2, 2:3] <- ci.icc11.sem.prof

  results[3, 1] <- icc11.sem
  results[3, 2:3] <- ci.icc11.sem.boot

  results[4, 1] <- icc11.mdd
  results[4, 2:3] <- ci.icc11.mdd.prof

  results[5, 1] <- icc11.mdd
  results[5, 2:3] <- ci.icc11.mdd.boot

  colnames(results) <- c('Estimate', sprintf('%g%% LB', 100*level), sprintf('%g%% UB', 100*level))
  rownames(results) <- c('ICC(1, 1)', 'SEM (Prof.)', 'SEM (Boot)', 'MDD (Prof.)', 'MDD (Boot)')

  ret <- list(icc11 = results)

  mySumm <- function(.) {
    s <- sigma(.)
    sd.components <- unname(s * getME(., "theta"))

    c(beta =getME(., "beta"), sigma = s, sig0 = sd.components, diff.var = 2*(s^2 + sd.components[2]))
  }

  mod.icc21 <- lmer(measure ~ 1 + (1 | id) + (1 | visit), data = dat)

  params <- mySumm(mod.icc21)

  boot.out <-  bootMer(mod.icc21, mySumm, nsim = B)#, .progress = 'txt')

  diff.boot.ci <- boot::boot.ci(boot.out, index=5, type=c("norm", "basic", "perc"))

  ci.icc21.sem <- sqrt(diff.boot.ci$percent[4:5])

  ci.icc21.mdd <- qnorm(0.975)*sqrt(diff.boot.ci$percent[4:5])

  icc21.sem <- sqrt(params[5])
  icc21.mdd <- qnorm(0.975)*icc21.sem

  results <- matrix(nrow = 3, ncol = 3)

  results[1, ] <- icc.21

  results[2, 1] <- icc21.sem
  results[2, 2:3] <- ci.icc21.sem
  results[3, 1] <- icc21.mdd
  results[3, 2:3] <- ci.icc21.mdd

  colnames(results) <- c('Estimate', sprintf('%g%% LB', 100*level), sprintf('%g%% UB', 100*level))
  rownames(results) <- c('ICC(2, 1)', 'SEM (Boot)', 'MDD (Boot)')

  ret$icc21 <- results

  return(ret)
}
