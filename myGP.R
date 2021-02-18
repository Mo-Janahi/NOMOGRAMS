myGP <- function (X, Z, d, g, dK = FALSE) 
{
  n <- nrow(X)
  m <- ncol(X)
  if (is.null(n)) 
    stop("X must be a matrix")
  if (length(Z) != n) 
    stop("must have nrow(X) = length(Z)")
  if (length(d) == 1) 
    d <- rep(d, m)
  else if (length(d) != m) 
    stop("must have length(d) = ncol(X)")
  out <- .C("newGPsep_R", m = as.integer(m), n = as.integer(n), 
            X = as.double(t(X)), Z = as.double(Z), d = as.double(d), 
            g = as.double(g), dK = as.integer(dK), gpsepi = integer(1), 
            PACKAGE = "laGP")
  return(out)
}
