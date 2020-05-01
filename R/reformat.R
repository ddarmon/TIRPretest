# For a generic data matrix
format.for.icc <- function(X){
  nv <- ncol(X) # number of visits per subject
  ns <- nrow(X) # number of subjects

  x <- X

  dim(x) <- NULL # flattens X in column-major order

  x.df <- data.frame(id = rep(1:ns, times = nv), visit = rep(1:nv, each = ns), measure = x)

  return(x.df)
}

# For the specific Excel format used by NKI
restack.df <- function(dat, measure){
  ids <- rep(dat$Subject_ID, times = 3)
  visits <- rep(1:3, each = nrow(dat))
  measurements <- c()

  for (i in 1:3){
    col.name <- sprintf('%s%g', measure, i)
    measurements <- c(measurements, dat[[col.name]])
  }

  dat.stacked <- data.frame(id = ids, visit = visits, measure = measurements)

  return(dat.stacked)
}
