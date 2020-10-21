# NOT RUN {
## examining a particular laGP call from the larger analysis provided
## in the aGP documentation

## A simple 2-d test function used in Gramacy & Apley (2014);
## thanks to Lee, Gramacy, Taddy, and others who have used it before
f2d <- function(x, y=NULL)
{
  if(is.null(y)) {
    if(!is.matrix(x) && !is.data.frame(x)) x <- matrix(x, ncol=2)
    y <- x[,2]; x <- x[,1]
  }
  g <- function(z)
    return(exp(-(z-1)^2) + exp(-0.8*(z+1)^2) - 0.05*sin(8*(z+0.1)))
  z <- -g(x)*g(y)
}

## build up a design with N=~40K locations
x <- seq(-2, 2, by=0.02)
X <- as.matrix(expand.grid(x, x))
Z <- f2d(X)

## optional first pass of nearest neighbor
Xref <- matrix(c(-1.725, 1.725), nrow=TRUE)
out <- laGP(Xref, 6, 50, X, Z, method="nn")

## second pass via ALC, ALCOPT, MSPE, and ALC-ray respectively,
## conditioned on MLE d-values found above
out2 <- laGP(Xref, 6, 50, X, Z, d=out$mle$d)
out2.alcopt <- laGP(Xref, 6, 50, X, Z, d=out2$mle$d, method="alcopt")
out2.mspe <- laGP(Xref, 6, 50, X, Z, d=out2$mle$d, method="mspe")
out2.alcray <- laGP(Xref, 6, 50, X, Z, d=out2$mle$d, method="alcray")

## look at the different designs
plot(rbind(X[out2$Xi,], X[out2.mspe$Xi,]), type="n",
     xlab="x1", ylab="x2", main="comparing local designs")
points(Xref[1], Xref[2], col=2, cex=0.5)
text(X[out2$Xi,], labels=1:50, cex=0.6)
text(X[out2.alcopt$Xi,], labels=1:50, cex=0.6, col="forestgreen")
text(X[out2.mspe$Xi,], labels=1:50, cex=0.6, col="blue")
text(X[out2.alcray$Xi,], labels=1:50, cex=0.6, col="red")
legend("right", c("ALC", "ALCOPT", "MSPE", "ALCRAY"),
       text.col=c("black", "forestgreen", "blue", "red"), bty="n")

## compare computational time
c(nn=out$time, alc=out2$time, alcopt=out2.alcopt$time, 
  mspe=out2.mspe$time, alcray=out2.alcray$time)

# }
# NOT RUN {
## Joint path sampling: a comparison between ALC-ex, ALC-opt and NN

## defining a predictive path
wx <- seq(-0.85, 0.45, length=100)
W <- cbind(wx-0.75, wx^3+0.51)

## three comparators from Sun, et al. (2017)
## larger-than-default "close" argument to capture locations nearby path
p.alc <- laGPsep(W, 6, 100, X, Z, close=10000, lite=FALSE)
p.alcopt <- laGPsep(W, 6, 100, X, Z, method="alcopt", lite=FALSE)
## note that close=10*(1000+end) would be the default for method = "alcopt"
p.nn <- laGPsep(W, 6, 100, X, Z, method="nn", close=10000, lite=FALSE)

## time comparison
c(alc=p.alc$time, alcopt=p.alcopt$time, nn=p.nn$time)

## visualization
plot(W, type="l", xlab="x1", ylab="x2", xlim=c(-2.25,0), ylim=c(-0.75,1.25), lwd=2)
points(X[p.alc$Xi,], col=2, cex=0.6)
lines(W[,1]+0.25, W[,2]-0.25, lwd=2)
points(X[p.nn$Xi,1]+0.25, X[p.nn$Xi,2]-0.25, pch=22, col=3, cex=0.6)
lines(W[,1]-0.25, W[,2]+0.25, lwd=2) 
points(X[p.alcopt$Xi,1]-0.25, X[p.alcopt$Xi,2]+0.25, pch=23, col=4, cex=0.6)
legend("bottomright", c("ALC-opt", "ALC-ex", "NN"), pch=c(22, 21, 23), col=c(4,2,3))
# }
