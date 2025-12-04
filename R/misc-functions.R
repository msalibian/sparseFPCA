#' Operator squared norm distance
#'
#' This function computes the squared "operator
#' norm distance" between its arguments. This is
#' closely related to the squared Frobenius norm
#' of the difference of the matrices.
#'
#' @param gammax the matrix representation of the first operator
#' @param gamma0 the matrix representation of the second operator
#'
#' @return the squared Frobenius norm of the difference of
#' the matrices, divided by the number of
#' elements.
#'
#' @export
norma <- function(gammax, gamma0){
  norm <- NA
  kx <- dim(gammax)[2]
  k0 <- dim(gamma0)[2]
  if(kx!=k0) { print("cuidado!!") }
  if(kx==k0) {
    gamaresta<-gammax-gamma0
    #  gamarestavec<-c(lowerTriangle(gamaresta,diag=TRUE),upperTriangle(gamaresta,diag=FALSE))
    # norm<-sum(gamarestavec*gamarestavec) /(kx*kx)
    norm<- sum(gamaresta^2)/(kx*kx)
  }
  return( norm )
}


#' All possible subsets
#'
#' This function returns all possible subsets of a given set.
#'
#' This function returns all possible subsets of a given set.
#'
#' @param n the size of the set of which subsets will be computed
#' @param k the size of the subsets
#' @param set an optional vector to be taken as the set
#'
#' @return A matrix with subsets listed in its rows
#'
#' @export
subsets <- function(n, k, set = 1:n) {
  if(k <= 0) NULL else if(k >= n) set
  else rbind(cbind(set[1], Recall(n - 1, k - 1, set[-1])), Recall(n - 1, k, set[-1]))
}

#' @export
clean <- function(X, h=.02) {
  N <- length(X$x)
  for (i in 1:N){
    t <- X$pp[[i]]
    ind <- c(1, which(diff(t) > h) + 1)
    X$pp[[i]] <- X$pp[[i]][ind]
    X$x[[i]] <- X$x[[i]][ind]
  }
  return(X)
}

#' @export
mu.f <- function(p) {
  tmp <- 5 + sin(p*pi*4)*10*exp(-2*p) + 5 * sin(p*pi/3) + 2 * cos(p*pi/2)
  return(tmp)
}



#' Generate irregular data
#'
#' Generate sparse irregular longitudinal samples.
#'
#' Generate sparse irregular longitudinal samples. The number of observations
#' per curve is \code{rbinom(nb.n, nb.p)}. It uses a Karhunen-Loeve
#' representations uing the mean function \code{mu.f} and (2)
#' eigenfunctions (in \code{phi}). The scores are \code{rnorm(0, 25)}
#' and \code{rnorm(0, 1)}.
#'
#' @param n the number of samples
#' @param mu.c the mean
#' @param nb.n the maximum number of observations per curve
#' @param nb.p the probability of success when computing the number of observations
#' per curve
#' @param mu.f the mean function
#' @param phi a list of length 2 containing two "eigenfunctions"
#'
#' @return A list with two components:
#' \item{x}{A list with the \code{n} vectors of observations (one per curve)}
#' \item{pp}{A corresponding list with the vectors of times at which observations
#' were taken}
#'
#' @export
data.coor <- function(n, mu.c=10, nb.n=10, nb.p=.25) {
  # generate irregular data
  mi <- 0
  ma <- 1
  pp <- x <- vector('list', n)
  np <- rnbinom(n, nb.n, nb.p)
  phi <- vector('list', 2)
  phi[[1]] <- function(t) sqrt(2)*cos(2*pi*t)
  phi[[2]] <- function(t) sqrt(2)*sin(2*pi*t)
  for(i in 1:n) {
    npi <- np[i]
    pp[[i]] <- sort( runif(npi, min=mi, max=ma) )
    xis <- rnorm(2, sd=.5) * c(5, 1) #
    p1 <- phi[[1]](pp[[i]]) * xis[1]
    p2 <- phi[[2]](pp[[i]]) * xis[2]
    x[[i]] <- mu.c + mu.f(pp[[i]]) + p1 + p2 # + rnorm(npi, sd=.5)
  }
  return(list(x=x, pp=pp))
}

#' Epanechnikov kernel
#'
#' This function computes the Epanechnikov kernel.
#'
#' This function computes the Epanechnikov kernel.
#'
#' @param x a real number
#'
#' @return 3/4(1-x^2) if x <= 1, 0 otherwise
#'
#' @export
k.epan <- function(x) {
  a <- 0.75*(1-x^2)
  k.epan <- a*(abs(x)<=1)
  return(k.epan)
}

# gg <- function(aa) sqrt( 25/4 * phi[[1]](aa)^2 + phi[[2]](aa)^2/4)
# ggs <- function(aa,bb) 25/4 * phi[[1]](aa) * phi[[1]](bb) + phi[[2]](aa) * phi[[2]](bb)/4

#' @export
psi <- function(r, k=1.345)
  pmin(k, pmax(-k, r))

#' @export
rho <- function(a, cc=3) {
  tmp <- 1-(1-(a/cc)^2)^3
  tmp[abs(a/cc)>1] <- 1
  return(tmp)
}

#' @export
rhoprime <- function(a, cc=3) {
  tmp <- 6*(a/cc)*(1-(a/cc)^2)^2
  tmp[abs(a/cc)>1] <- 0
  return(tmp)
}

#' @export
tukey.weight <- function(a) {
  tmp <- 6*(1-a^2)^2
  tmp[ abs(a) > 1 ] <- 0
  return(tmp)
}

# tukey.weight <- function(a) {
#   tmp <- 1/abs(a)
#   tmp[ abs(a) < 1 ] <- 1
#   return(tmp)
# }




#' @export
funwchiBT <- function(u,cBT){
  c <- cBT
  U <- abs(u)
  Ugtc <- (U > c)
  w <- (3-3*(u/c)^2+(u/c)^4)/(c*c)
  w[Ugtc] <- 1/(U[Ugtc]*U[Ugtc])
  w
}



#' Generate regular data
#'
#' Generate regular longitudinal samples.
#'
#' Generate regular longitudinal samples. The observations
#' per curve are \code{seq(0, 1, length=np)}. It uses a Karhunen-Loeve
#' representations uing the mean function \code{mu.f} and (2)
#' eigenfunctions (in \code{phi}). The scores are \code{rnorm(0, 25)}
#' and \code{rnorm(0, 1)}.
#'
#' @param n the number of samples
#' @param mu.c the mean
#' @param np the number of observations per curve
#' @param mu.f the mean function
#' @param phi a list of length 2 containing two "eigenfunctions"
#'
#' @return A list with two components:
#' \item{x}{A list with the \code{n} vectors of observations (one per curve)}
#' \item{pp}{A corresponding list with the vectors of times at which observations
#' were taken}
#'
#' @export
data.coor2 <- function(n, mu.c=10, np=100, mu.f, phi) {
  # generate irregular data
  mi <- 0
  ma <- 1
  pp <- x <- vector('list', n)
  # np <- rnbinom(n, nb.n, nb.p)
  for(i in 1:n) {
    pp[[i]] <- seq(mi, ma, length=np)
    xis <- rnorm(2, sd=.5) * c(5, 1) #
    p1 <- phi[[1]](pp[[i]]) * xis[1]
    p2 <- phi[[2]](pp[[i]]) * xis[2]
    x[[i]] <- mu.c + mu.f(pp[[i]]) + p1 + p2 # + rnorm(npi, sd=.5)
  }
  return(list(x=x, pp=pp))
}

#' @export
matrixx <- function(X, mh) {
  # build all possible pairs Y_{ij}, Y_{il}, j != l
  # and the corresponding times t_{ij}, t_{i,l}, j != l
  # They are returned in "m" and "mt" below
  #
  n <- length(X$x)
  M <- MT <- NULL
  for (i in 1:n){
    comb <- subsets(length(X$x[[i]]),2)
    if(class(comb)[1] != 'matrix') comb <- matrix(comb, byrow=TRUE, ncol=2)
    Maux <- matrix(X$x[[i]][comb]-mh[[i]][comb],ncol=2)
    MTaux <- matrix(X$pp[[i]][comb],ncol=2)
    M <- rbind(M, Maux, cbind(Maux[,2], Maux[,1]))
    MT <- rbind(MT, MTaux, cbind(MTaux[,2], MTaux[,1]))
    #
    # I'm not sure whether we need both pairs for each j != l,
    # i.e. (Y_{ij}, Y_{il}) and (Y_{il}, Y_{ij}) check the two lines of code above
    #
#     tmp <- cbind(1:length(X$x[[i]]), 1:length(X$x[[i]]) )
#     M <- rbind(M, matrix(X$x[[i]][tmp]-mh[[i]][tmp],ncol=2))
#     MT <- rbind(MT, matrix(X$pp[[i]][tmp],ncol=2))
  }
  return(list(m=M,mt=MT))
}


matrixx2 <- function(X, mh){
  # this version keeps the diagonal
  n <- length(X$x)
  M <- MT <- NULL
  for (i in 1:n){
    comb <- subsets(length(X$x[[i]]),2)
    if(class(comb)[1] != 'matrix') comb <- matrix(comb, byrow=TRUE, ncol=2)
    Maux <- matrix(X$x[[i]][comb]-mh[[i]][comb],ncol=2)
    MTaux <- matrix(X$pp[[i]][comb],ncol=2)
    M <- rbind(M, Maux, cbind(Maux[,2], Maux[,1]))
    MT <- rbind(MT, MTaux, cbind(MTaux[,2], MTaux[,1]))
        tmp <- cbind(1:length(X$x[[i]]), 1:length(X$x[[i]]) )
        M <- rbind(M, matrix(X$x[[i]][tmp]-mh[[i]][tmp],ncol=2))
        MT <- rbind(MT, matrix(X$pp[[i]][tmp],ncol=2))
  }
  return(list(m=M,mt=MT))
}



kernel1 <- function(s,t)
{
  return((1/2)*(1/2)^(0.9 * abs(s-t)))
}

kernel2 <- function(s,t)
{
  return(10*min(s,t))
}

kernel3 <- function(s,t)
{
  return((1/2)*exp(-0.9*log(2)*abs(s-t)^2)  )
}


# phin <- function(t,n){sqrt(2)*sin(0.5*(2*n-1)*pi*t)}

#' @export
data.inf.C4 <- function(n, p, ep1, Sigma, mu.c=30,sigma.c=0.01,k=5) {

  ### Functional data
  mi <- 0
  ma <- 1
  np <- p
  pp <- seq(mi, ma, length=np)
  phin <- function(t,n){sqrt(2)*sin(0.5*(2*n-1)*pi*t)}

  direction.cont <- phin(pp,k)
  x <- my.mvrnorm(n, rep(0,p), Sigma)
  Z <-rnorm(n,mu.c,sigma.c)

  outl <- rbinom(n, size=1, prob=ep1)

  # add outliers
  if(ep1 > 0) {
  ncu <- sum(outl)
  for (i in 1:n)
    x[i,] <- x[i,] + outl[i]*Z[i]*direccion.cont
  }
  return(list(x=x, outl=outl, pp=pp))
}


# mu.f2 <- function(p) {
#   tmp <- 5 + cos(p*pi*4)*10*exp(-2*p)
#   return(tmp)
# }

#' @export
mu.f2 <- function(p) {
  tmp <- p + exp(-2*p)
  return(tmp*5)
}

#' @export
data.coor3 <- function(n, mu.c=10, nb.n=10, nb.p=.25) {
  # generate irregular data
  mi <- 0
  ma <- 1
  pp <- x <- vector('list', n)
  np <- rnbinom(n, nb.n, nb.p)
  phin <- function(t,n){sqrt(2)*sin(0.5*(2*n-1)*pi*t)}
  for(i in 1:n) {
    npi <- np[i]
    pp[[i]] <- sort( runif(npi, min=mi, max=ma) )
    xis <- rnorm(3, sd=.5) * c(5, 1, .5) #
    p1 <- phin(pp[[i]], 1) * xis[1]
    p2 <- phin(pp[[i]], 2) * xis[2]
    p3 <- phin(pp[[i]], 3) * xis[3]
    x[[i]] <- mu.f2(pp[[i]]) + p1 + p2 + p3 # + rnorm(npi, sd=.5)
  }
  return(list(x=x, pp=pp))
}

# ggs2 <- function(aa,bb) 25/4 * phin(aa,1) * phin(bb,1) + phin(aa,2) * phin(bb,2)/4

#' @export
data.coor4 <- function(n, mu.c=10, nb.n=10, nb.p=.25, k = 5, eps=0.1) {
  # generate irregular data
  mi <- 0
  ma <- 1
  pp <- x <- vector('list', n)
  np <- pmax( rnbinom(n, nb.n, nb.p), 2)
  outs <- rbinom(n, size=1, prob=eps)
  # lambdas <- c(6, 1, .5, .2, .1) / 2
  lambdas <- c(3, 1, 0.25, 0.1, 0.05)
  ps <- matrix(0, max(np), k)
  xxis <- matrix(0, n, k)
  #if(k > length(lambdas)) lambdas <- c(lambdas, rep(.1, k - length(lambdas)))
#   if(k > length(lambdas)) lambdas <- c(lambdas, .1^((1:(k - length(lambdas)))+1))
#   else lambdas <- lambdas[1:k]
  if(k > length(lambdas)) { lambdas <- c(lambdas, .1^((1:(k - length(lambdas)))+1))
  } else {
    lambdas <- lambdas[1:k]
  }
  for(i in 1:n) {
    npi <- np[i]
    pp[[i]] <- sort( runif(npi, min=mi, max=ma) )
    xis <- rnorm(k) * lambdas[1:k]
    if (outs[i] == 1) {
      #xis <- c(rnorm(1,sd=1/16), rt(n=1,df=2), rt(n=1,df=2))
      xis <- rnorm(3, mean=mu.c, sd=1)
      if(k > 3) xis <- c(xis, rep(0, k - 3))
      }
    xis <- xis[1:k]
    for(j in 1:k) ps[1:npi, j] <- phin(pp[[i]], j) * xis[j]
    x[[i]] <- mu.f2(pp[[i]]) + rowSums(ps[1:npi, ]) # + rnorm(npi, sd=.5)
    xxis[i, ] <- xis
  }
  return(list(x=x, pp=pp, outs=outs, xis=xxis, k=k, lambdas=lambdas^2))
}

#' @export
data.coor5 <- function(n, mu.c=10, nb.n=10, nb.p=.25, eps=0.1) {
  # generate irregular data
  mi <- 0
  ma <- 1
  pp <- x <- vector('list', n)
  np <- rnbinom(n, nb.n, nb.p)
  ncontam=rbinom(1,n,eps)
  phin <- function(t,n){sqrt(2)*sin(0.5*(2*n-1)*pi*t)}
  for(i in 1:n){
    npi <- np[i]
    pp[[i]] <- sort( runif(npi, min=mi, max=ma) )
    xis <- rnorm(3, sd=.5) * c(5, 1, .5)
    #
    p1 <- phin(pp[[i]], 1) * xis[1]
    p2 <- phin(pp[[i]], 2) * xis[2]
    p3 <- phin(pp[[i]], 3) * xis[3]
    if (rbinom(1,1,eps) == 1) x[[i]] <- mu.f2cont(pp[[i]]) + p1 + p2 + p3

    else x[[i]] <- mu.f2(pp[[i]]) + p1 + p2 + p3 # + rnorm(npi, sd=.5)
  }
  return(list(x=x, pp=pp, ncontam=ncontam))
}

#' @export
mu.f2cont <- function(p) {
  tmp <- 25*sqrt(p)
  return(tmp)
}

# gtshat.new <-function(X,t0,s0,h,mh,matx,cc=1.56,eps=1e-6){
#
#   # now add an intercept...
#
#   n=length(X$x)
#   B=rnorm(1);
#   err=1+eps
#   j=0
#
#   # sigma.hat <- sqrt(ggs(t0,t0) - ggs(t0, s0)^2 / ggs(s0, s0))
#
#   # print(sigma.hat)
#   if(missing(mh)) mh <- mu.hat3(X=X,h=h,cc=cc,ep=eps)
#
#   # XX=mu.hat3(X,mh,h=h,cc=cc,ep=eps)
#   if(missing(matx)) matx <- matrixx(X=X, mh=mh)
#
#   M <- matx$m
#   MT <- matx$mt
#
#
#   w2 <- k.epan((MT[,1]-t0)/h)
#   w3 <- k.epan((MT[,2]-s0)/h)
#
#   we <- w2*w3
#
#   M <- M[ we > 0, ]
#   we <- we[ we > 0]
#   theta0 <- repmed.w(x=M[,1], y=M[,2], w=we)
#   sigma.hat <- mad.w( M[,2] - theta0[1] - theta0[2] * M[,1], w=we)
#
#   ## Need to include the intercept in the iterations
#   ## below...
#   z <- cbind(rep(1, length(M[,1])), M[,1])
#   alpha <- theta0[1]
#   beta <- theta0[2]
#   while ( ( (j <- j+1) < 1000 ) && (err > eps) ){
#
#     re <- (M[,2] - alpha - beta*M[,1] )/ sigma.hat
#     # w1 <- psi((M[,1]-B*M[,2])/sigma.hat,cc)/(M[,1]-B*M[,2])
#     w1 <- psi(re,cc)/re
#     w1[ is.na(w1) ] <- 1 # missing values are "psi(0)/0", that lim is 1
#     w <- w1 * we
#     tmp <- solve( t(z) %*% (z * w), t(z) %*% (M[,2]*w))
#     # B2=sum(w*M[,1]*M[,2])/sum(w*M[,2]^2)
#     B2 <- tmp[2]
#     # print(c(j, mean(rho((M[,1]-B*M[,2])/sigma.hat,cc)*k.epan((MT[,1]-t0)/h)*k.epan((MT[,2]-s0)/h))))
#     err <- abs(B2/B - 1)
#     #print(c(B, B2))
#     beta <- B <- B2
#     alpha <- tmp[1]
#   }
#   return(beta)
# }


#' @export
gtshat.med <- function(X,t0,s0,h,mh,matx,cc=1.56,eps=1e-6){
  # median of slopes
  n <- length(X$x)
  err <- 1+eps
  j <- 0
  if(missing(mh)) mh <- mu.hat3(X=X,h=h,cc=cc,ep=eps)
  if(missing(matx)) matx <- matrixx(X=X, mh=mh)
  M <- matx$m
  MT <- matx$mt
  w2 <- k.epan((MT[,1]-t0)/h)
  w3 <- k.epan((MT[,2]-s0)/h)
  we <- w2*w3
  M <- M[ we > 0, ]
  we <- we[ we > 0]
  # get rid of points on the diagonal
#   xs <- M[,1]
#   ys <- M[,2]
#   tmp <- (xs != ys)
#   B <- median( M[tmp,2] / M[tmp,1 ])
#   sigma.hat <- mad( M[tmp,2] - B * M[tmp,1])
  B <- median( M[,2] / M[,1 ])
#   sigma.hat <- mad( M[,2] - B * M[,1])
# #   B <- med.w( x=M[,2] / M[,1 ], w=we)
# #   sigma.hat <- mad.w( x=M[,2] - B * M[,1], w=we)
#   while ( ( (j <- j+1) < 1000 ) && (err > eps) ){
#     w1 <- psi((M[,2]-B*M[,1])/sigma.hat,cc)/(M[,2]-B*M[,1])
#     w1[ is.na(w1) ] <- 1/sigma.hat # missing values are "psi(0)/0", that lim is 1/sigma.hat
#     w <- w1 * we
#     B2 <- sum(w*M[,1]*M[,2])/sum(w*(M[,1]^2))
#     err <- abs(B2/B - 1)
#     B <- B2
#   }
  return(B)
}

#' @export
gtshat <- function(X, t0, s0, h, mh, matx, cc=3.443689, eps=1e-6){ # 3.443689 4.685065 1.345
  n <- length(X$x)
  err <- 1+eps
  j <- 0
  if(missing(mh)) mh <- mu.hat3(X=X,h=h,cc=cc,ep=eps)
  if(missing(matx)) matx <- matrixx(X=X, mh=mh)
  M <- matx$m
  MT <- matx$mt
  w2 <- k.epan((MT[,1]-t0)/h)
  w3 <- k.epan((MT[,2]-s0)/h)
  we <- w2*w3
  M <- M[ we > 0, ]
  we <- we[ we > 0]
  if( length(we)==0 ) return(NA)
  # get rid of points on the diagonal
  #   xs <- M[,1]
  #   ys <- M[,2]
  #   tmp <- (xs != ys)
  #   B <- median( M[tmp,2] / M[tmp,1 ])
  #   sigma.hat <- mad( M[tmp,2] - B * M[tmp,1])
  B <- median( M[,2] / M[,1 ])
  sigma.hat <- mad( M[,2] - B * M[,1])
    # B <- med.w( x=M[,2] / M[,1 ], w=we)
    # sigma.hat <- mad.w( x=M[,2] - B * M[,1], w=we)
  while ( ( (j <- j+1) < 1000 ) && (err > eps) ){
    # w1 <- rhoprime( (M[,2]-B*M[,1]) / sigma.hat, cc) / (M[,2]-B*M[,1])
    # w1[ is.na(w1) ] <- 1/sigma.hat # missing values are "psi(0)/0", that lim is 1/sigma.hat
    w1 <- tukey.weight( ( (M[,2]-B*M[,1]) / sigma.hat ) / cc )
    w <- w1 * we
    B2 <- sum(w*M[,1]*M[,2])/sum(w*(M[,1]^2))
    err <- abs(B2/B - 1)
    if( is.na(err) || is.nan(err) ) return(NA)
    B <- B2
  }
  return(B)
}

# gtshat.lin <-function(X,t0,s0,h,mh,matx,cc=1.56,eps=1e-6){
#   n <- length(X$x)
#   err <- 1+eps
#   j <- 0
#   if(missing(mh)) mh <- mu.hat3(X=X,h=h,cc=cc,ep=eps)
#   if(missing(matx)) matx <- matrixx(X=X, mh=mh)
#   M <- matx$m
#   MT <- matx$mt
#   w2 <- k.epan((MT[,1]-t0)/h)
#   w3 <- k.epan((MT[,2]-s0)/h)
#   we <- w2*w3
#   M <- M[ we > 0, ]
#   MT <- MT[ we > 0, ]
#   we <- we[ we > 0]
#   B <- median( M[,2] / M[,1 ])
#   sigma.hat <- mad( M[,2] - B * M[,1])
#   #   B <- med.w( x=M[,2] / M[,1 ], w=we)
#   #   sigma.hat <- mad.w( x=M[,2] - B * M[,1], w=we)
#   t <- cbind(MT[,1] - t0, MT[,2] - s0)
#   tt <- cbind(rep(1, nrow(t)), t)*M[,1]
#   beta <- rlm(x=tt, y=M[,2])$coef
#
#   while ( ( (j <- j+1) < 1000 ) && (err > eps) ){
#       re <- as.vector(M[,2] - tt %*% beta) / sigma.hat
#       w <- ( psi(re,cc)/re )
#       w[ is.nan(w) ] <- 1
#       w <- w * we
#       beta.n <- solve( t( tt * w ) %*% tt , t(tt * w) %*% M[,2] )
#       err <- sum( (beta.n - beta)^2 )
#       beta <- beta.n
#     }
#   return(beta[1])
# }


# gtshat2 <-function(X,t0,s0,h,mh,matx,cc=1.56,eps=1e-6){
#   # this version uses the diagonal
#   n <- length(X$x)
#   err <- 1+eps
#   j <- 0
#   if(missing(mh)) mh <- mu.hat3(X=X,h=h,cc=cc,ep=eps)
#   if(missing(matx)) matx <- matrixx(X=X, mh=mh)
#   M <- matx$m
#   MT <- matx$mt
#   w2 <- k.epan((MT[,1]-t0)/h)
#   w3 <- k.epan((MT[,2]-s0)/h)
#   we <- w2*w3
#   M <- M[ we > 0, ]
#   we <- we[ we > 0]
#   # get rid of points on the diagonal
#     xs <- M[,1]
#     ys <- M[,2]
#     tmp <- (xs != ys)
#     B <- median( M[tmp,2] / M[tmp,1 ])
#     sigma.hat <- mad( M[tmp,2] - B * M[tmp,1])
# #   B <- median( M[,2] / M[,1 ])
# #   sigma.hat <- mad( M[,2] - B * M[,1])
#   #   B <- med.w( x=M[,2] / M[,1 ], w=we)
#   #   sigma.hat <- mad.w( x=M[,2] - B * M[,1], w=we)
#   while ( ( (j <- j+1) < 1000 ) && (err > eps) ){
#     w1 <- psi((M[,2]-B*M[,1])/sigma.hat,cc)/(M[,2]-B*M[,1])
#     w1[ is.na(w1) ] <- 1/sigma.hat # missing values are "psi(0)/0", that lim is 1/sigma.hat
#     w <- w1 * we
#     B2 <- sum(w*M[,1]*M[,2])/sum(w*(M[,1]^2))
#     err <- abs(B2/B - 1)
#     B <- B2
#   }
#   return(B)
# }

#' @export
gtshat.mscale <-function(X,t0,s0,h,mh,matx,cc=1.56,eps=1e-6){
  n <- length(X$x)
  err <- 1+eps
  j <- 0
  if(missing(mh)) mh <- mu.hat3(X=X,h=h,cc=cc,ep=eps)
  if(missing(matx)) matx <- matrixx(X=X, mh=mh)
  M <- matx$m
  MT <- matx$mt
  w2 <- k.epan((MT[,1]-t0)/h)
  w3 <- k.epan((MT[,2]-s0)/h)
  we <- w2*w3
  M <- M[ we > 0, ]
  we <- we[ we > 0]
  B <- median( M[,2] / M[,1 ])
  sigma.hat <- mscale.w( M[,2] - B * M[,1], w = we)
  while ( ( (j <- j+1) < 1000 ) && (err > eps) ){
    w1 <- psi((M[,2]-B*M[,1])/sigma.hat,cc)/(M[,2]-B*M[,1])
    w1[ is.na(w1) ] <- 1/sigma.hat # missing values are "psi(0)/0", that lim is 1/sigma.hat
    w <- w1 * we
    B2 <- sum(w*M[,1]*M[,2])/sum(w*(M[,1]^2))
    err <- abs(B2/B - 1)
    B <- B2
  }
  return(B)
}

#' @export
gtthat <- function(X, t0, h, muhat, b=.5, cc=1.54764, initial.sc, max.it=300, eps=1e-10) {
  # find the scale, full iterations
  # muhat should be a list as returned by mu.hat3
  i <- 0
  err <- 1 + eps
  t <- unlist(X$pp)
  x <- unlist(X$x)
  if(missing(muhat))
    muhat <- mu.hat3(X=X,h=h,cc=cc,ep=eps)
  muhat <- unlist(muhat)
  if(missing(initial.sc)) sc <- mad(x - muhat) else sc <- initial.sc
  while( ( (i <- i+1) < max.it ) && (err > eps) ) {
    kerns <- k.epan((t-t0)/h)
    sc2 <- sqrt(sc^2 * sum(kerns * rho((x-muhat)/sc,cc)) / (b * sum(kerns)))
    err <- abs(sc2/sc - 1)
    if( is.na(err) || is.nan(err) ) return(NA)
    #   print(summary(kerns))
    #   print(sc2)
    #   print(sc)
    # }
    sc <- sc2
  }
  return(list(gtt=sc, muhat=muhat))
}

#' @export
matrixx.rot <- function(X) {
  tmp <- relist(rep(0, length(unlist(X$x))), X$x)
  return( matrixx(X=X, mh=tmp) )
}



#' @export
uhat <- function(X, t0, h=0.1, cc=1.56, ep=1e-6, max.it=100){
  x <- unlist(X$x)
  t <- unlist(X$pp)
  s <- localMAD(X,t0,h)
  oldu <- (u <- median(x)) + 100*ep
  it <- 0
  kw <- k.epan((t-t0)/h)
  while( ((it <- it + 1) < max.it ) && ( (abs(oldu) - u) > ep ) ){
    w <- psi((x-u)/s,cc)/(x-u)
    w[ is.nan(w) ] <- 1
    w <- w*kw
    oldu <- u
    u <- sum(w*x)/sum(w)
  }
  return(u)
}

#' @export
uhat.lin <- function(X, t0, h=0.1, cc=1.345, ep=1e-6, max.it=100){
  x <- unlist(X$x)
  t <- unlist(X$pp)
  s <- localMAD(X,t0,h)
  # oldu <- (u <- median(x)) + 100*ep
  it <- 0
  tt <- cbind(rep(1, length(t)), t-t0)
  wk <- k.epan((t-t0)/h)
  beta <- MASS::rlm(x=tt[wk>0,], y=x[wk>0])$coef
  while( ((it <- it + 1) < max.it ) ){
    re <- as.vector(x - tt %*% beta)/s
    w <- ( psi(re,cc)/re )
    w[ is.nan(w) ] <- 1
    w <- w * wk
    beta.n <- solve( t( tt * w ) %*% tt, t(tt * w) %*% x )
    if( any( is.na(beta.n) | is.nan(beta.n) ) ) return(NA)
    if( sum( (beta.n - beta)^2 ) < ep ) it = max.it
    beta <- beta.n
  }
  return(beta.n[1])
}



# mu.hat <- function(X,h=0.1,cc=1.56,ep=1e-6) {
#   n=length(X$x)
#   us <- vector('list', n)
#
#   for(i in 1:n){
#     tmp=rep(0,length(X$x[[i]]))
#     tt=X$pp[[i]]
#
#     for (j in 1:length(tt)){
#       tmp[j]=uhat(X=X, t0=tt[j], h=h, cc=cc, ep=ep)
#       }
#     us[[i]]<- tmp
#   }
#
#   return(list(x=X$x,pp=X$pp,us=us))
#
# }
#

#' @export
mu.hat2 <- function(X,h=0.1,cc=1.56,ep=1e-6) {
  tt <- unlist(X$pp)
  nt <- length(tt)
  us <- rep(0, nt)
  for(i in 1:nt)
    us[i] <- uhat(X=X, t0=tt[i], h=h, cc=cc, ep=ep)
  return(us)
}

#' @export
mu.hat2.lin <- function(X,h=0.1,cc=1.56,ep=1e-6) {
  tt <- unlist(X$pp)
  nt <- length(tt)
  us <- rep(0, nt)
  for(i in 1:nt)
    us[i] <- uhat.lin(X=X, t0=tt[i], h=h, cc=cc, ep=ep)
  return(us)

}

#' @export
localMAD <- function(X,t0,h) {
  x <- unlist(X$x)
  t <- unlist(X$pp)
  window <- which(t < t0+h & t > t0-h)
  return( mad(x[window]) )
}

#' @export
mu.hat3 <- function(X,h=0.1,cc=1.56,ep=1e-6) {
  us=relist(mu.hat2(X=X,h=h,cc=cc,ep=ep),X$x)
  return(us)
  #return(list(x=X$x,pp=X$pp,u=us))
}

#' @export
mu.hat3.lin <- function(X,h=0.1,cc=1.56,ep=1e-6) {
  us=relist(mu.hat2.lin(X=X,h=h,cc=cc,ep=ep),X$x)
  return(us)
  #return(list(x=X$x,pp=X$pp,u=us))
}



# true.mu.hat3 <- function(X) {
#   us=relist(true.mu.hat2(X),X$x)
#   return(us)
#   #return(list(x=X$x,pp=X$pp,u=us))
# }
#
#
# true.mu.hat2 <- function(X) {
#
#   tt <- unlist(X$pp)
#   nt <- length(tt)
#   us <- rep(0, nt)
#   for(i in 1:nt)
#     us[i] <- mu.f2(tt[i])
#   return(us)
#
# }



#' @export
mu.hat4 <- function(X, cc=1.56, ep=1e-6) {
  # all curves observed at the same time points
  pp <- X$pp[[1]]
  np <- length(pp)
  muh <- rep(0, np)
  for(i in 1:np) muh[i] <- m.location(sapply(X$x, function(a,i) a[i], i=i), cc=cc, ep=ep)
  return(muh)
}

#' @export
m.location <- function(x,cc=1.56,ep=1e-6, max.it=500){

  s <- mad(x)
  u <- median(x)

  for (i in 1:max.it){

    w <- psi((x-u)/s,cc)/(x-u)
    w[abs(x-u)<1e-8] <- 1
    oldu <- u
    u <- sum(w*x)/sum(w)
    if(abs(u-oldu)/abs(oldu)<ep){
      break
    }
  }
  return(u)
}



kern1 <- function(s,t)
{
  return((1/2)*(1/2)^(0.9 * abs(s-t)))
}

kern2 <- function(s,t)
{
  return(10*pmin(s,t))
}

kern3 <- function(s,t)
{
  return((1/2)*exp(-0.9*log(2)*abs(s-t)^2)  )
}

# phin <- function(t,n){sqrt(2)*sin(0.5*(2*n-1)*pi*t)}
# phin <- function(t,n){sqrt(2)*sin(0.5*(2*n-1)*pi*t)}

#' @export
my.mvrnorm <- function(n, mu, si) {
  p <- length(mu)
  ul <- svd(si)
  u <- ul$u
  l <- diag(sqrt(ul$d))
  x <- scale(matrix(rnorm(n*p), n, p) %*% u %*% l %*% t(u), center = -mu, scale=FALSE)
  attr(x, 'scaled:center') <- NULL
  if(n == 1) return(drop(x)) else return(x)

}

#' @export
data.coor6 <- function(n, mu.c=10, nb.n=10, nb.p=.25, eps=0.1, cont.k=5, cont.mu = 30) {
  # generate irregular data
  # infinite dimensional process
  mi <- 0
  ma <- 1
  pp <- x <- vector('list', n)
  np <- pmax(rnbinom(n, nb.n, nb.p), 2) # no curve with a single obs
  ncontam <- rbinom(1,n,eps)
  phin <- function(t,n){sqrt(2)*sin(0.5*(2*n-1)*pi*t)}
  for(i in 1:n) {
    npi <- np[i]
    pp[[i]] <- sort( runif(npi, min=mi, max=ma) )
    Sigma1 <- outer(pp[[i]], pp[[i]], kern2)
    x[[i]] <- my.mvrnorm(1, rep(0,npi), Sigma1)
#     if (i <= ncontam) xis <- c(rnorm(1,sd=1/8), rt(n=1,df=2), rnorm(1,sd=5))
#     else xis <- rnorm(3, sd=.5) * c(5, 1, .5)
#     p1 <- phin(pp[[i]], 1) * xis[1]
#     p2 <- phin(pp[[i]], 2) * xis[2]
#     p3 <- phin(pp[[i]], 3) * xis[3]
    if( i <= ncontam) x[[i]] <- x[[i]] + rnorm(1, cont.mu, .5) * phin(pp[[i]], cont.k)
    x[[i]] <- x[[i]] + mu.f2(pp[[i]])
  }
  return(list(x=x, pp=pp, ncontam=ncontam))
}

#' @export
data.muller <- function(n, eps, mu.c=12, ns=2:4, max.tries=100) {
  # generate irregular data
  # ns = vector of possible number of obs per trajectory
  mu <- function(t) t + sin(t)
  phi <- vector('list', 2)
  phi[[1]] <- function(t) -cos(pi*t/10)/sqrt(5)
  phi[[2]] <- function(t) sin(pi*t/10)/sqrt(5)
  lambdas <- c(4, 1)  # these are the variances of the scores
  # set.seed(seed)
  tries <- 1
  repeated <- TRUE
  while( (tries < max.tries) & (repeated) ) {
    tts <- seq(0, 10, length=51) + rnorm(51, sd=sqrt(.1))
    tts <- sort(pmin( pmax(tts, 0), 10)) # we'll sample from tts[2:50]
    repeated <- ( length(tts[2:50]) != length(unique(tts[2:50])) )
    tries <- tries + 1
  }
  if( (tries == max.tries) & repeated ) stop('Repeated times')
  # ns <- 2:5
  pp <- x <- vector('list', n)
  np <- sample(ns, n, repl=TRUE)
  xxis <- matrix(0, n, 2)
  outs <- rbinom(n, size=1, prob=eps)
  for(i in 1:n) {
    npi <- np[i]
    pp[[i]] <- sort( tts[ sample(2:50, npi) ] )
    xis <- rnorm(2) * sqrt( lambdas )
    if(outs[i]==1) xis[2] <- rnorm(n=1, mean=mu.c, sd=1)
    x[[i]] <- mu(pp[[i]]) + phi[[1]](pp[[i]])*xis[1] + phi[[2]](pp[[i]])*xis[2]
    xxis[i, ] <- xis
  }
  return(list(x=x, pp=pp, xis=xxis, mu=mu, phis=phi, lambda=lambdas, outs=outs))
}

export <- function(X, file='data-efpac.csv') {
  N=length(X$x)
  tmp <- sapply(X$x, length)
  mpt=max(tmp)
  XT=matrix(NaN,nrow=2*N,ncol=mpt)
  for(i in 1:N) XT[i,1:(tmp[i])] <- X$x[[i]]
  for(i in (N+1):(2*N)) XT[i,1:(tmp[i-N])] <- X$pp[[i-N]]
  write.csv(XT, file=file, row.names=F)
}

#' @export
repmed <- function(x, y) {
  n <- length(x)
  tmp <- rep(0,n)
  for(i in 1:n) {
    hh <- (x != x[i])
    tmp[i] <- median( (y[hh]-y[i])/(x[hh]-x[i]) )
  }
  b <- median(tmp)
  a <- median( y - b * x )
  return(c(a,b))
}

#' @export
med.w <- function(x, w=rep(1,n)) {
  oo <- order(x)
  n <- length(x)
  w <- w[oo] / sum(w)
  tmp <- min( (1:n)[ cumsum(w) >= .5] )
  return( (x[oo])[tmp] )
}

#' @export
repmed.w <- function(x, y, w) {
  n <- length(x)
  tmp <- rep(0,n)
  for(i in 1:n) {
    hh <- (x != x[i])
    tmp[i] <- med.w( x=(y[hh]-y[i])/(x[hh]-x[i]), w=w[hh] )
  }
  b <- med.w(x=tmp, w)
  a <- med.w(x=y - b * x, w=w)
  return(c(a,b))
}

#' @export
mad.w <- function(x, w)
{
  mu <- med.w(x=x, w=w)
  return(med.w(x=abs(x-mu), w=w)*1.4826)
}
#
# x <- rnorm(100)
# y <- (x-1)*(x+2)/5 + rnorm(100, sd=.5*abs(x)/3)
# plot(x,y)
# x0 <- 0.5
# w <- k.epan((x-x0)/.3)
# th <- repmed.w(x=x, y=y, w=w)
# re <- y - th[1] - th[2]*x
# mu.w <- med.w(x=re, w=w)
# (mad.w <- med.w(x=abs(re-mu.w), w=w))
#

# cov.fun.hat <- function(X, h, mh, ma, ncov=50, trace=FALSE) {
#   if(trace) print("Computing cov function")
#   if(missing(mh)) mh <- mu.hat3(X=X, h=h)
#   if(missing(ma)) ma <- matrixx(X=X, mh=mh)
#   mii <- min( ti <- unlist(X$pp) )
#   maa <- max( ti )
#   tt <- seq(mii, maa, length=ncov)
#   ss <- seq(mii, maa, length=ncov)
#   pps <- as.matrix(expand.grid(tt, ss))
#   np <- nrow(pps)
#   ghat <- rep(0, np)
#   for(j in 1:np) {
#     t0 <- pps[j, 1]
#     s0 <- pps[j, 2]
#     if (t0==s0) { ghat[j] <- (gtthat(X=X, t0=t0, h=h, muhat=mh)$gtt)^2
#      #print( gtshat(X=X,t0=t0,s0=s0,h=h,mh=mh,matx=ma,cc=1.56,eps=1e-6))
#     }
#     else {
#       sigma.hat1 <- gtthat(X=X, t0=t0, h=h, muhat=mh)$gtt
#       ghat[j] <- gtshat(X=X,t0=t0,s0=s0,h=h,mh=mh,matx=ma,cc=1.56,eps=1e-6) * sigma.hat1^2
#     }
#   }
#   G <- matrix(ghat,ncov,ncov)
#   G <- ( G + t(G) ) / 2
#   if(trace) print("Done computing cov function")
#   return(list(G=G, grid=pps))
# }

#' @export
cov.fun.hat.mscale <- function(X, h, mh, ma, ncov=50, trace=FALSE) {
  if(trace) print("Computing cov function")
  if(missing(mh)) mh <- mu.hat3(X=X, h=h)
  if(missing(ma)) ma <- matrixx(X=X, mh=mh)
  mii <- min( ti <- unlist(X$pp) )
  maa <- max( ti )
  tt <- seq(mii, maa, length=ncov)
  ss <- seq(mii, maa, length=ncov)
  pps <- as.matrix(expand.grid(tt, ss))
  np <- nrow(pps)
  ghat <- rep(0, np)
  for(j in 1:np) {
    t0 <- pps[j, 1]
    s0 <- pps[j, 2]
    if (t0==s0) ghat[j] <- (gtthat(X=X, t0=t0, h=h, muhat=mh)$gtt)^2
    else {
      sigma.hat1 <- gtthat(X=X, t0=t0, h=h, muhat=mh)$gtt
      ghat[j] <- gtshat.mscale(X=X,t0=t0,s0=s0,h=h,mh=mh,matx=ma,cc=1.56,eps=1e-6) * sigma.hat1^2
    }
  }
  G <- matrix(ghat,ncov,ncov)
  G <- ( G + t(G) ) / 2
  if(trace) print("Done computing cov function")
  return(list(G=G, grid=pps))
}

# cov.fun.hat.lin <- function(X, h, mh, ma, ncov=50, trace=FALSE) {
#   if(trace) print("Computing cov function")
#   if(missing(mh)) mh <- mu.hat3(X=X, h=h)
#   if(missing(ma)) ma <- matrixx(X=X, mh=mh)
#   mii <- min( ti <- unlist(X$pp) )
#   maa <- max( ti )
#   tt <- seq(mii, maa, length=ncov)
#   ss <- seq(mii, maa, length=ncov)
#   pps <- as.matrix(expand.grid(tt, ss))
#   np <- nrow(pps)
#   ghat <- rep(0, np)
#   for(j in 1:np) {
#     t0 <- pps[j, 1]
#     s0 <- pps[j, 2]
#     if (t0==s0) ghat[j] <- (gtthat(X=X, t0=t0, h=h, muhat=mh)$gtt)^2
#     else {
#       sigma.hat1 <- gtthat(X=X, t0=t0, h=h, muhat=mh)$gtt
#       ghat[j] <- gtshat.lin(X=X,t0=t0,s0=s0,h=h,mh=mh,matx=ma,cc=1.56,eps=1e-6) * sigma.hat1^2
#     }
#   }
#   G <- matrix(ghat,ncov,ncov)
#   G <- ( G + t(G) ) / 2
#   if(trace) print("Done computing cov function")
#   return(list(G=G, grid=pps))
# }

#' @export
cov.fun.hat2 <- function(X, h, mh, ma, ncov=50, trace=FALSE) {
  # this function uses the diagonal
  if(trace) print("Computing cov function")
  if(missing(mh)) mh <- mu.hat3(X=X, h=h)
  if(missing(ma)) ma <- matrixx(X=X, mh=mh)
  mii <- min( ti <- unlist(X$pp) )
  maa <- max( ti )
  tt <- seq(mii, maa, length=ncov)
  ss <- seq(mii, maa, length=ncov)
  pps <- as.matrix(expand.grid(tt, ss))
  np <- nrow(pps)
  betahat <- rep(0, np)
  sigmahat <- rep(0, ncov)
  for(j in 1:ncov) sigmahat[j] <- gtthat(X=X, t0=tt[j], h=h, muhat=mh)$gtt
  for(j in 1:np) {
    t0 <- pps[j, 1]
    s0 <- pps[j, 2]
    betahat[j] <- gtshat(X=X,t0=t0,s0=s0,h=h,mh=mh,matx=ma,eps=1e-6)
    # gamma[s0, t0] / gamma[t0, t0]
  }
  G <- betahat <- matrix(betahat, ncov, ncov)
  for(i in 1:ncov)
    for(j in 1:ncov)
      G[i,j] <- betahat[i,j] * sigmahat[i]^2
  G <- ( G + t(G) ) / 2
  if(trace) print("Done computing cov function")
  return(list(G=G, grid=pps))
}

# cov.fun.hat.lin2 <- function(X, h, mh, ma, ncov=50, trace=FALSE) {
#   # this function uses the diagonal
#   if(trace) print("Computing cov function")
#   if(missing(mh)) mh <- mu.hat3(X=X, h=h)
#   if(missing(ma)) ma <- matrixx(X=X, mh=mh)
#   mii <- min( ti <- unlist(X$pp) )
#   maa <- max( ti )
#   tt <- seq(mii, maa, length=ncov)
#   ss <- seq(mii, maa, length=ncov)
#   pps <- as.matrix(expand.grid(tt, ss))
#   np <- nrow(pps)
#   betahat <- rep(0, np)
#   sigmahat <- rep(0, ncov)
#   for(j in 1:ncov) sigmahat[j] <- gtthat(X=X, t0=tt[j], h=h, muhat=mh)$gtt
#   for(j in 1:np) {
#     t0 <- pps[j, 1]
#     s0 <- pps[j, 2]
#     betahat[j] <- gtshat.lin(X=X,t0=t0,s0=s0,h=h,mh=mh,matx=ma,cc=1.56,eps=1e-6)
#   }
#   G <- betahat <- matrix(betahat, ncov, ncov)
#   for(i in 1:ncov)
#     for(j in 1:ncov)
#       G[i,j] <- betahat[i,j] * sigmahat[i]^2
#   G <- ( G + t(G) ) / 2
#   if(trace) print("Done computing cov function")
#   return(list(G=G, grid=pps))
# }

#' @export
mscale.w <- function(r, w, c.BT=1.56, b.BT = 0.5, max.it = 500, eps=1e-8) {
  iiter <- 1
  sigma <- mad(r)
  suma <- sum(w)
  while (iiter < max.it){
    sigmaINI <- sigma
    resid <- r/sigmaINI
    sigma2 <- sum( w*(r^2)*funwchiBT(resid, c.BT), na.rm=T )/( suma*b.BT )
    sigma <- sqrt(sigma2)
    criterio <- abs(sigmaINI-sigma)/sigmaINI
    if(criterio <= eps){
      sigmaBT <- sigma
      iiter <- max.it } else {
        iiter <- iiter+1
        if( iiter == max.it ){
          sigmaBT <- sigma
        }
      }
  }
  return(max(sigmaBT, 1e-10))
}


#' @export
pred <- function(X, muh, cov.fun, tt, ss, k=2, s=20, rho=0) {
  # k = number of eigenfunctions to use in the KL expansion (number of scores)
  # s = number of eigenfunctions to use to compute Sigma_Y
  # (if s < length(Y), then Sigma_Y will be singular)
  # rho = regularization parameter
  # compute eigenfunctions
  eg <- eigen(cov.fun)
  lam <- eg$values[1:max(k,s)]
  ef <- eg$vectors[,1:max(k,s)]
  # eg <- svd(cov.fun)
  # lam <- eg$d[1:max(k,s)]
  # ef <- eg$u[,1:max(k,s)]
  lam[ lam < 0 ] <- 0 # drop negative eigenvalues
  # keep eigenfunctions with "positive" eigenvalues
  s1 <- max( (1:max(k,s))[ lam > 1e-5 ] )
  # standardize eigenfunctions
  normas <- apply(ef, 2, L2.norma.mesh, mesh=tt) # rep(1, max(k,s))
  efn <- scale(ef, center=FALSE, scale=normas)
  # make them into functions
  ff <- vector('list', max(k,s))
  for(i in 1:max(k,s))
    ff[[i]] <- approxfun(tt, efn[,i], method='linear')
  # compute predicted scores
  n <- length(X$x)
  xis <- matrix(NA, n, k)
  pr <- matrix(NA, n, length(tt))
  for(i in 1:n) {
    ti <- X$pp[[i]]
    xic <- X$x[[i]] - muh[[i]]
    phis <- matrix(0, length(ti), max(k,s))
    for(j in 1:max(k,s))
      phis[,j] <- ff[[j]](ti)
    # Sigma_Y
    siy <- phis[,1:s1, drop=FALSE] %*% diag( lam[1:s1] ) %*% t(phis[,1:s1, drop=FALSE])
    # Sigma_Y^{-1} (X - mu)
    rhs <- as.vector( solve(siy + rho * diag(length(ti)), xic ) )
    # scores  = \phi' \lambdas Sigma_Y^{-1} (X - mu)
    if(k > 1) {
      xis[i,] <- t( phis[,1:k, drop=FALSE] %*% diag( lam[1:k] ) ) %*% rhs
    } else {
      xis[i,] <- t( phis[,1, drop=FALSE] * lam[1] ) %*% rhs
    }
    # \hat{X - mu} = \sum phi_j \xi_j
    pr[i,] <- efn[,1:k, drop=FALSE] %*% as.vector( xis[i,] )
  }
  # compute mean function on the grid used to compute
  # the covariance function / matrix
  tobs <- unlist(X$pp)
  oopp <- order(tobs)
  mu.tobs <- unlist(muh)[oopp]
  tobs <- tobs[oopp]
  mu.tt <- approx(tobs, mu.tobs, tt, method='linear')$y
  # \hat{X} = \hat{X - mu} + mu
  pr <- scale(pr, center=-mu.tt, scale=FALSE)
  return(list(xis=xis, pred=pr))
}


# calcula la integral de efe en (0,1) sobre una grila mesh
# 0<mesh[1]<mesh[2]<...<mesh[m]<1

integral <- function(efe,mesh){
  # we need that the first and last element of mesh
  # are the extreme points of the interval where we integrate
  m <- length(mesh)
  t0 <- min(mesh)
  tm1 <- max(mesh)
  primero<- (mesh[1]-t0)*efe[1]
  ultimo<-(tm1-mesh[m])*efe[m]
  sumo<-primero+ultimo
  menos1<-m-1
  for (i in 1:menos1)
  {
    a1<-(efe[i]+efe[i+1])/2
    amplitud<- mesh[i+1]-mesh[i]
    sumo<-sumo+amplitud*a1
  }
  return(sumo)
}

#################################
# Calcula el producto interno entre dos datos (como listas)  sobre una grilla de puntos MESH en (0,1) ordenados
##########################

#' @export
L2.dot.product.mesh <- function(dato1, dato2, mesh)
{
  return(integral(dato1*dato2,mesh))
}



#################################
# Calcula la norma de un dato funcional sobre una grilla de puntos MESH en (0,1) ordenados
# (OJO ES LA NORMA no el cuadrado de la norma!!!)
##########################

#' @export
L2.norma.mesh <- function(dato, mesh)
{
  return(sqrt(L2.dot.product.mesh(dato,dato,mesh)))
}

#' @export
gtshat.Sest <- function(X,t0,s0,h,mh,matx,cc=1.56,eps=1e-6){
  n <- length(X$x)
  err <- 1+eps
  j <- 0
  if(missing(mh)) mh <- mu.hat3(X=X,h=h,cc=cc,ep=eps)
  if(missing(matx)) matx <- matrixx(X=X, mh=mh)
  M <- matx$m
  MT <- matx$mt
  w2 <- k.epan((MT[,1]-t0)/h)
  w3 <- k.epan((MT[,2]-s0)/h)
  we <- w2*w3
  M <- M[ we > 0, ]
  we <- we[ we > 0]
  Sestim<-fastSloc(M, N=50, k=4, best.r=5, bdp=.5)

  cova<-Sestim$covariance

  rho=cova[1,2]/sqrt(cova[1,1]*cova[2,2])
  return(list(covast=cova[1,2],rhost=rho))
}

#' @export
cov.fun.Sest.rho <- function(X, h, mh, ma, ncov=50, trace=FALSE) {
  if(trace) print("Computing cov function")
  if(missing(mh)) mh <- mu.hat3(X=X, h=h)
  if(missing(ma)) ma <- matrixx(X=X, mh=mh)
  mii <- min( ti <- unlist(X$pp) )
  maa <- max( ti )
  tt <- seq(mii, maa, length=ncov)
  ss <- seq(mii, maa, length=ncov)
  pps <- as.matrix(expand.grid(tt, ss))
  np <- nrow(pps)
  ghat <- rep(0, np)
  for(j in 1:np) {
    t0 <- pps[j, 1]
    s0 <- pps[j, 2]
    sigma.t0 <- gtthat(X=X, t0=t0, h=h, muhat=mh)$gtt
    sigma.s0 <- gtthat(X=X, t0=s0, h=h, muhat=mh)$gtt

    if(t0==s0){ghat[j] <- sigma.t0^2}
    else {
      ghat[j] <- gtshat.Sest(X=X,t0=t0,s0=s0,h=h,mh=mh,matx=ma,cc=1.56,eps=1e-6)$rhost *  sigma.t0 * sigma.s0 }

  }
  G <- matrix(ghat,ncov,ncov)
  G <- ( G + t(G) ) / 2
  if(trace) print("Done computing cov function")
  return(list(G=G, grid=pps))
}

#' @export
cov.fun.Sest.cov <- function(X, h, mh, ma, ncov=50, trace=FALSE) {
  if(trace) print("Computing cov function")
  if(missing(mh)) mh <- mu.hat3(X=X, h=h)
  if(missing(ma)) ma <- matrixx(X=X, mh=mh)
  mii <- min( ti <- unlist(X$pp) )
  maa <- max( ti )
  tt <- seq(mii, maa, length=ncov)
  ss <- seq(mii, maa, length=ncov)
  pps <- as.matrix(expand.grid(tt, ss))
  np <- nrow(pps)
  ghat <- rep(0, np)
  for(j in 1:np) {
    t0 <- pps[j, 1]
    s0 <- pps[j, 2]
    sigma.t0 <- gtthat(X=X, t0=t0, h=h, muhat=mh)$gtt
    sigma.s0 <- gtthat(X=X, t0=s0, h=h, muhat=mh)$gtt

    if(t0==s0){ghat[j] <- sigma.t0^2}
    else {
      ghat[j] <- gtshat.Sest(X=X,t0=t0,s0=s0,h=h,mh=mh,matx=ma,cc=1.56,eps=1e-6)$covast }

  }
  G <- matrix(ghat,ncov,ncov)
  G <- ( G + t(G) ) / 2
  if(trace) print("Done computing cov function")
  return(list(G=G, grid=pps))
}

#' @export
mattern.cov <- function(d, nu=5/2, rho=1, sigma=1) {
  tmp <- sqrt(2*nu)*d/rho
  a <- besselK(x=tmp, nu=nu) #, expon.scaled = FALSE)
  b <- tmp^nu
  return( sigma^2 * 2^(1-nu) * b * a / gamma(nu) )
}

#' @export
mu.f3 <- function(p) {
  tmp <- 5*sin( 2*pi*p ) * exp(-p*3)
  return(tmp*2)
}


#' @export
data.mattern <- function(n, mu.c=10, nobs=2:5, k = 5, eps=0.1, phis, lambdas, mean.f = mu.f3) { #nb.n=10, nb.p=.25, k = 5,
  # generate irregular data
  # outliers are in the eigendecomposition
  # lambdas = sqrt of eigenvalues -- sd of scores
  mi <- 0
  ma <- 1
  pp <- x <- vector('list', n)
  np <- sample(nobs, n, replace=TRUE) # pmax( rnbinom(n, nb.n, nb.p), 2)
  outs <- rbinom(n, size=1, prob=eps)
  ps <- matrix(0, max(np), k)
  xxis <- matrix(0, n, k)
  for(i in 1:n) {
    npi <- np[i]
    pp[[i]] <- sort( runif(npi, min=mi, max=ma) )
    xis <- rnorm(k) * lambdas[1:k]
    if (outs[i] == 1) {
      #xis <- c(rnorm(1,sd=1/16), rt(n=1,df=2), rt(n=1,df=2))
      xis[2:3] <- rnorm(2, mean=mu.c, sd=.25)
      # if(k > 3) xis <- c(xis, rep(0, k - 3))
    }
    xis <- xis[1:k]
    for(j in 1:k) ps[1:npi, j] <- phis[[j]](pp[[i]]) * xis[j]
    x[[i]] <- mean.f(pp[[i]]) + rowSums(ps[1:npi, ]) # + rnorm(npi, sd=.5)
    xxis[i, ] <- xis
  }
  return(list(x=x, pp=pp, outs=outs, xis=xxis, k=k, lambdas=lambdas))
}

#' @export
data.mattern2 <- function(n, mu.c=c(20, 25), nobs=3:5, q = length(phis), k=q, eps=0.0, phis, lambdas, mean.f=mu.f3) { #nb.n=10, nb.p=.25,
  # generate irregular data
  # outliers are not in the eigendecomposition
  # lambdas = eigenvalues -- variance of scores
  # q = number of components used to generate
  # k = number of components to return
  mi <- 0
  ma <- 1
  pp <- x <- vector('list', n)
  np <- sample(nobs, n, replace=TRUE) # pmax( rnbinom(n, nb.n, nb.p), 2)
  outs <- rbinom(n, size=1, prob=eps)
  xxis <- matrix(0, n, q)
  ps <- matrix(0, max(np), q) # auxiliar space
  for(i in 1:n) {
    npi <- np[i]
    pp[[i]] <- sort( runif(npi, min=mi, max=ma) )
    xis <- rnorm(q, sd = sqrt(lambdas[1:q]))
    if (outs[i] == 1) xis[2:3] <- rnorm(2, mean=mu.c, sd=1/4) * sqrt(lambdas[2:3])
    for(j in 1:q) ps[1:npi, j] <- phis[[j]](pp[[i]]) * xis[j]
    x[[i]] <- mean.f(pp[[i]]) + rowSums(ps[1:npi, ]) # + rnorm(npi, sd=.5)
    xxis[i, ] <- xis
  }
  return(list(x=x, pp=pp, outs=outs, xis=xxis[, 1:k], k=k,
              q=q, lambdas=lambdas[1:k]))
}

# data.mattern3 <- function(n, mu.c=10, nobs=2:5, k = 5, eps=0.1, phis, lambdas) {
#   # generate irregular data
#   # no. obs. per trajectory are uniform on the set {2, 3, 4, 5}
#   # outliers are in the eigendecomposition
#   # lambdas = sqrt of eigenvalues -- sd of scores
#   mi <- 0
#   ma <- 1
#   pp <- x <- vector('list', n)
#   np <- sample(nobs, n, replace=TRUE) # pmax( rnbinom(n, nb.n, nb.p), 2)
#   outs <- rbinom(n, size=1, prob=eps)
#   ps <- matrix(0, max(np), k)
#   xxis <- matrix(0, n, k)
#   if(k > length(lambdas)) { lambdas <- c(lambdas, .1^((1:(k - length(lambdas)))+1))
#     } else {
#       lambdas <- lambdas[1:k]
#     }
#   for(i in 1:n) {
#     npi <- np[i]
#     pp[[i]] <- sort( runif(npi, min=mi, max=ma) )
#     xis <- rnorm(k) * lambdas[1:k]
#     if (outs[i] == 1) {
#       #xis <- c(rnorm(1,sd=1/16), rt(n=1,df=2), rt(n=1,df=2))
#       xis[2:3] <- rnorm(2, mean=mu.c, sd=.25)
#       # if(k > 3) xis <- c(xis, rep(0, k - 3))
#     }
#     xis <- xis[1:k]
#     for(j in 1:k) ps[1:npi, j] <- phis[[j]](pp[[i]]) * xis[j]
#     x[[i]] <- mu.f3(pp[[i]]) + rowSums(ps[1:npi, ]) # + rnorm(npi, sd=.5)
#     xxis[i, ] <- xis
#   }
#   return(list(x=x, pp=pp, outs=outs, xis=xxis, k=k, lambdas=lambdas))
# }


#' @export
mu.f4 <- function(p) {
  tmp <- mu.f3(p)
  p2 <- p[ p > .5 ]
  tmp[ p > .5 ] <- 5*sin( 2*pi*p2 ) * 2
  return(tmp)
}

#' Trimmed mean of squares
#'
#' This function returns the mean of the smallest
#' (1-alpha)% squares of its arguments.
#'
#' @param x numeric vector
#' @param alpha real number between 0 and 1
#'
#' @return the mean of the smallest (1-alpha)% squared values in \code{x}
#'
#' @export
tm <- function(x, alpha) {
  n <- length(x)
  return( mean( (sort(x^2, na.last=NA))[1:(n - floor(alpha*n))], na.rm=TRUE ) )
}

#' Trimmed mean
#'
#' This function returns the mean of the smallest
#' (1-alpha)% of its arguments.
#'
#' @param x numeric vector
#' @param alpha real number between 0 and 1
#'
#' @return the mean of the smallest (1-alpha)% values in \code{x}
#'
#' @export
tmns <- function(x, alpha) {
  # tmns = trimmed mean, no square
  n <- length(x)
  return( mean( (sort(x, na.last=NA))[1:(n - floor(alpha*n))], na.rm=TRUE ) )
}


# cv.mu <- function(X, k=5, hs=exp(seq(-4, 0, by=.6)), alpha=.4, seed=123 ) {
#   # CV for the mean function
#   lh <- length(hs)
#   n <- length(X$x)
#   if(exists(".Random.seed", where=.GlobalEnv)) old.seed <- .Random.seed
#   set.seed(seed)
#   fl <- sample( (1:n) %% k + 1)
#   tmses <- rep(NA, lh)
#   Xtmp <- vector('list', 2)
#   names(Xtmp) <- c('x', 'pp')
#   for(j in 1:lh) {
#     Xh <- relist(NA, X$x) # predictions go here
#     for(i in 1:k) {
#       this <- (1:n)[ fl != i ]
#       Xtmp$x <- X$x[ this ] # training obs
#       Xtmp$pp <- X$pp[ this ] # training times
#       ps <- unlist(X$pp[ -this ]) # times where we need to predict
#       xhats <- rep(NA, length(ps))
#       for(l in 1:length(ps)) {
#         tmp2 <- try(uhat.lin(X=Xtmp, t0=ps[l], h=hs[j], cc=1.56, ep=1e-6, max.it=100), silent=TRUE)
#         if( class(tmp2) != 'try-error' )
#           xhats[l] <- tmp2
#       }
#       Xh[-this] <- relist(xhats, X$x[ -this ] ) # fill predictions
#     }
#     tmp <- sapply(mapply('-', X$x, Xh), function(a) a^2) # squared residuals, list-wise
#     tmses[j] <- tm(unlist(tmp), alpha=alpha) # 40% trimming, use 60% smallest resids
#   }
#   if(exists('old.seed')) assign('.Random.seed', old.seed, envir=.GlobalEnv)
#   return(list(tmse=tmses, h=hs))
# }

#' @export
cv.mu.par <- function(X, k.cv=5, k = k.cv, hs=exp(seq(-4, 0, by=.6)), alpha=.4, seed=123) {
  # parallel computing version
  # Cross validation for the mean function
  lh <- length(hs)
  n <- length(X$x)
  if(exists(".Random.seed", where=.GlobalEnv)) old.seed <- .Random.seed
  set.seed(seed)
  fl <- sample( (1:n) %% k.cv + 1)
  tmses <- rep(NA, lh)
  Xtmp <- vector('list', 2)
  names(Xtmp) <- c('x', 'pp')
  tmp.par <- foreach(h=hs, .combine=c, .inorder=FALSE, .packages='MASS',
                .export=c('tm', 'uhat.lin', 'localMAD', 'psi', 'k.epan')) %dopar% {
    Xh <- relist(NA, X$x) # predictions go here
    for(i in 1:k.cv) {
      this <- (1:n)[ fl != i ]
      Xtmp$x <- X$x[ this ] # training obs
      Xtmp$pp <- X$pp[ this ] # training times
      ps <- unlist(X$pp[ -this ]) # times where we need to predict
      xhats <- rep(NA, length(ps))
      for(l in 1:length(ps)) {
        tmp2 <- try(uhat.lin(X=Xtmp, t0=ps[l], h=h, cc=1.56, ep=1e-6, max.it=100), silent=TRUE)
        if( !inherits(tmp2, "try-error") )
          xhats[l] <- tmp2
      }
      Xh[-this] <- relist(xhats, X$x[ -this ] ) # fill predictions
    }
    # tmp <- sapply(mapply('-', X$x, Xh), function(a) a^2) # squared residuals, list-wise
    tmp <- mapply('-', X$x, Xh)  # squared residuals, list-wise
    # # tmses <- tm(unlist(tmp), alpha=alpha) # 40% trimming, use 60% smallest resids
    # tmses <- RobStatTM::mscale(unlist(tmp)) #, delta=.3, tuning.chi=2.560841)
    tmp2 <- unlist(tmp)
    if(all(is.na(tmp2))) {
      tmses <- NA } else {
        tmp2 <- tmp2[ !is.na(tmp2) ]
        # tmp2 <- tmp2[ !is.na(tmp2) ]
        # if(length(tmp2) > 0) {
        tmp3 <- RobStatTM::locScaleM(tmp2, psi='bisquare')
        tmses <- tmp3$disper^2 + tmp3$mu^2 # RobStatTM::mscale(tmp2) #, delta=.3, tuning.chi=2.560841)
      }

  }
  if(exists('old.seed')) assign('.Random.seed', old.seed, envir=.GlobalEnv)
  return(list(tmse=tmp.par, h=hs))
}


#' @export
cov.fun.cv.par <- function(X, muh, ncov=50, k.cv=5, hs=exp(seq(-4, 0, by=.6)),
                       alpha=.4, seed=123, k=2, s=20, reg.rho=1e-5) {
  # parallel computing version
  # CV for the covariance function
  # muh = estimated mean function
  mii <- min( ti <- unlist(X$pp) )
  maa <- max( ti )
  tt <- seq(mii, maa, length=ncov)
  lh <- length(hs)
  n <- length(X$x)
  if(exists(".Random.seed", where=.GlobalEnv)) old.seed <- .Random.seed
  set.seed(seed)
  fl <- sample( (1:n) %% k.cv + 1)
  tmspe <- rep(NA, lh)
  Xtest <- Xtrain <- vector('list', 2)
  names(Xtrain) <- c('x', 'pp')
  names(Xtest) <- c('x', 'pp')
  tmp.par <- foreach(h=hs, .combine=c, .inorder=FALSE, .packages=c('MASS', 'mgcv'),
                     .export=c('tm', 'cov.fun.hat2', 'matrixx', 'subsets',
                               'gtthat', 'gtshat', 'mu.hat3',
                               'k.epan', 'psi', 'rho', 'pred.cv', 'tukey.weight',
                               'L2.norma.mesh', 'L2.dot.product.mesh',
                               'integral')) %dopar% {
    Xhat <- relist(NA, X$x)
    for(i in 1:k.cv) {
      this <- (1:n)[ fl != i ] # training set
      Xtrain$x <- X$x[ this ]
      Xtrain$pp <- X$pp[ this ]
      Xtest$x <- X$x[ -this ]
      Xtest$pp <- X$pp[ -this ]
      ma <- matrixx(Xtrain, muh[ this ])
      cov.fun <- try(cov.fun.hat2(X=Xtrain, h=h, mh=muh[ this ],
                      ma=ma, ncov=50, trace=FALSE)) # $G
      if( !inherits(cov.fun, 'try-error') ) {
        if(!any(is.na(cov.fun$G))) {
        uu <- as.vector(cov.fun$G) #
        ttx <- cov.fun$grid #
        cov.fun <- matrix(fitted(gam(uu ~ s(ttx[,1], ttx[,2]))), length(tt), length(tt)) #
        cov.fun <- ( cov.fun + t(cov.fun) ) / 2
        tmp <- try( pred.cv(X=Xtrain, muh=muh[ this ], X.pred=Xtest,
                            muh.pred=muh[ -this ], cov.fun=cov.fun, tt=tt,
                            k=k, s=s, rho=reg.rho) )
        if( !inherits(tmp, 'try-error')) Xhat[ -this ] <- tmp
        }
      }
    }
    # tmp <- sapply(mapply('-', X$x, Xhat), function(a) a^2) # squared residuals, list-wise
    tmp <- mapply('-', X$x, Xhat) # squared residuals, list-wise
    # tmspe <- tm(unlist(tmp), alpha=alpha) # 40% trimming, use 60% smallest resids
    # tmspe <- RobStatTM::mscale(unlist(tmp)) #, delta=.3, tuning.chi=2.560841)
    tmp2 <- unlist(tmp)
    if(any(is.na(tmp2))) {
      tmspe <- NA } else {
    # tmp2 <- tmp2[ !is.na(tmp2) ]
    # if(length(tmp2) > 0) {
        tmp3 <- RobStatTM::locScaleM(tmp2, psi='bisquare')
        tmspe <- tmp3$disper^2 + tmp3$mu^2  # RobStatTM::mscale(tmp2) #, delta=.3, tuning.chi=2.560841)
      }
    # } else {
    #   tmspe <- NA
    # }
  }
  if(exists('old.seed')) assign('.Random.seed', old.seed, envir=.GlobalEnv)
  return(list(tmspe=tmp.par, h=hs, ncov=ncov, k=k, s=s, rho=reg.rho))
}

#' @export
cov.fun.cv <- function(X, muh, ncov=50, k.cv=5, hs=exp(seq(-4, 0, by=.6)),
                       alpha=.4, seed=123, k=2, s=20, reg.rho=1e-5) {
  # CV for the covariance function
  # muh = estimated mean function
  mii <- min( ti <- unlist(X$pp) )
  maa <- max( ti )
  tt <- seq(mii, maa, length=ncov)
  lh <- length(hs)
  n <- length(X$x)
  if(exists(".Random.seed", where=.GlobalEnv)) old.seed <- .Random.seed
  set.seed(seed)
  fl <- sample( (1:n) %% k.cv + 1)
  tmspe <- rep(NA, lh)
  Xtest <- Xtrain <- vector('list', 2)
  names(Xtrain) <- c('x', 'pp')
  names(Xtest) <- c('x', 'pp')
  for(j in 1:lh) {
    Xhat <- relist(NA, X$x)
    for(i in 1:k.cv) {
      this <- (1:n)[ fl != i ] # training set
      Xtrain$x <- X$x[ this ]
      Xtrain$pp <- X$pp[ this ]
      Xtest$x <- X$x[ -this ]
      Xtest$pp <- X$pp[ -this ]
      ma <- matrixx(Xtrain, muh[ this ])
      # cov.fun <- try(cov.fun.hat2(X=Xtrain, h=hs[j], mh=muh[ this ],
      #                             ma=ma, ncov=50, trace=FALSE)$G)
      cov.fun <- try(cov.fun.hat2(X=Xtrain, h=hs[j], mh=muh[ this ],
                                  ma=ma, ncov=50, trace=FALSE))
      if( !inherits(cov.fun, 'try-error') ) {
        if(!any(is.na(cov.fun$G))) {
        uu <- as.vector(cov.fun$G) #
        ttx <- cov.fun$grid #
        cov.fun <- matrix(fitted(gam(uu ~ s(ttx[,1], ttx[,2]))), length(tt), length(tt)) #
        cov.fun <- ( cov.fun + t(cov.fun) ) / 2
        tmp <- try( pred.cv(X=Xtrain, muh=muh[ this ], X.pred=Xtest,
                            muh.pred=muh[ -this ], cov.fun=cov.fun, tt=tt,
                            k=k, s=s, rho=reg.rho) )
        if( !inherits(tmp, 'try-error') ) Xhat[ -this ] <- tmp
        }
      }
    }
    # tmp <- sapply(mapply('-', X$x, Xhat), function(a) a^2) # squared residuals, list-wise
    tmp <- mapply('-', X$x, Xhat) # squared residuals, list-wise
    # tmspe <- tm(unlist(tmp), alpha=alpha) # 40% trimming, use 60% smallest resids
    tmp2 <- unlist(tmp)
    if(any(is.na(tmp2))) {
      tmpse[j] <- NA } else {
    # tmp2 <- tmp2[ !is.na(tmp2) ]
    # if(length(tmp2) > 0) {
        tmp3 <- RobStatTM::locScaleM(tmp2, psi='bisquare')
        tmspe[j] <- tmp3$disper^2 + tmp3$mu^2 # RobStatTM::mscale(tmp2) #, delta=.3, tuning.chi=2.560841)
    }
  }
  if(exists('old.seed')) assign('.Random.seed', old.seed, envir=.GlobalEnv)
  return(list(tmspe=tmspe, h=hs, ncov=ncov, k=k, s=s, rho=reg.rho))
}




# cov.fun.cv <- function(X, muh, ncov=50, k.cv=5, hs=exp(seq(-4, 0, by=.6)),
#                        alpha=.4, seed=123, k=2, s=20, reg.rho=1e-5) {
#   # CV for the covariance function
#   # muh = estimated mean function
#   mii <- min( ti <- unlist(X$pp) )
#   maa <- max( ti )
#   tt <- seq(mii, maa, length=ncov)
#
#   lh <- length(hs)
#   n <- length(X$x)
#   if(exists(".Random.seed", where=.GlobalEnv)) old.seed <- .Random.seed
#   set.seed(seed)
#   fl <- sample( (1:n) %% k.cv + 1)
#   tmspe <- rep(NA, lh)
#   Xtest <- Xtrain <- vector('list', 2)
#   names(Xtrain) <- c('x', 'pp')
#   names(Xtest) <- c('x', 'pp')
#   for(j in 1:lh) {
#     Xhat <- relist(NA, X$x)
#     for(i in 1:k.cv) {
#       this <- (1:n)[ fl != i ] # training set
#       Xtrain$x <- X$x[ this ]
#       Xtrain$pp <- X$pp[ this ]
#       Xtest$x <- X$x[ -this ]
#       Xtest$pp <- X$pp[ -this ]
#       ma <- matrixx(Xtrain, muh[ this ])
#       cov.fun <- cov.fun.hat2(X=Xtrain, h=hs[j], mh=muh[ this ],
#                               ma=ma, ncov=50, trace=FALSE)$G
#       Xhat[ -this ] <- pred.cv(X=Xtrain, muh=muh[ this ], X.pred=Xtest,
#                                muh.pred=muh[ -this ], cov.fun=cov.fun, tt=tt,
#                                k=k, s=s, rho=reg.rho)
#     }
#     tmp <- sapply(mapply('-', X$x, Xhat), function(a) a^2) # squared residuals, list-wise
#     tmspe[j] <- tm(unlist(tmp), alpha=alpha) # 40% trimming, use 60% smallest resids
#     print(c(hs[j], tmspe[j]))
#   }
#   if(exists('old.seed')) assign('.Random.seed', old.seed, envir=.GlobalEnv)
#   return(list(tmspe=tmspe, h=hs, ncov=ncov, k=k, s=s, rho=rho))
# }


#' @export
pred.cv <- function(X, muh, X.pred, muh.pred, cov.fun, tt, k=2, s=20, rho=0) {
  # prediction based on cov.fun
  # but for trajectories not present in X (those in X.pred)
  # X: training set
  # muh: estimated mean for training set
  # X.pred: test set
  # muh.pred: estimated mean for test set

  eg <- eigen(cov.fun)
  lam <- eg$values[1:max(k,s)]
  ef <- eg$vectors[,1:max(k,s)]
  # eg <- svd(cov.fun)
  # lam <- eg$d[1:max(k,s)]
  # ef <- eg$u[,1:max(k,s)]
  lam[ lam < 0 ] <- 0
  s1 <- max( (1:max(k,s))[ lam > 1e-5 ] )
  normas <- apply(ef, 2, L2.norma.mesh, mesh=tt) # rep(1, max(k,s))
  efn <- scale(ef, center=FALSE, scale=normas)
  ff <- vector('list', max(k,s))
  for(i in 1:max(k,s))
    ff[[i]] <- approxfun(tt, efn[,i], method='linear')
  n <- length(X.pred$x)
  xis <- matrix(NA, n, k)
  Xhat <- relist(NA, X.pred$x)
  for(i in 1:n) {
    ti <- X.pred$pp[[i]]
    xic <- X.pred$x[[i]] - muh.pred[[i]]
    phis <- matrix(NA, length(ti), max(k,s))
    for(j in 1:max(k,s))
      phis[,j] <- ff[[j]](ti)
    siy <- phis[,1:s1, drop=FALSE] %*% diag( lam[1:s1] ) %*% t(phis[,1:s1, drop=FALSE])
    rhs <- as.vector( solve(siy + rho * diag(length(ti)), xic ) )
    if(k > 1) {
      xis[i,] <- t( phis[,1:k, drop=FALSE] %*% diag( lam[1:k] ) ) %*% rhs
    } else {
      xis[i,] <- t( phis[,1, drop=FALSE] * lam[1] ) %*% rhs
    }
    Xhat[[i]] <- as.vector( phis[,1:k, drop=FALSE] %*% as.vector( xis[i,] ) ) + muh.pred[[i]]
  }
  return(Xhat)
}


#' @export
pred.cv.whole <- function(X, muh, X.pred, muh.pred, cov.fun, tt, k=2, s=20, rho=0) {
  # prediction based on cov.fun
  # but for trajectories not present in X (those in X.pred)
  # X: training set
  # muh: estimated mean for training set on tt
  # X.pred: test set
  # muh.pred: estimated mean for test set
  eg <- eigen(cov.fun)
  lam <- eg$values[1:max(k,s)]
  ef <- eg$vectors[,1:max(k,s)]
  # eg <- svd(cov.fun)
  # lam <- eg$d[1:max(k,s)]
  # ef <- eg$u[,1:max(k,s)]
  lam[ lam < 0 ] <- 0
  s1 <- max( (1:max(k,s))[ lam > 1e-5 ] )
  normas <- apply(ef, 2, L2.norma.mesh, mesh=tt) # rep(1, max(k,s))
  efn <- scale(ef, center=FALSE, scale=normas)
  ff <- vector('list', max(k,s))
  for(i in 1:max(k,s))
    ff[[i]] <- approxfun(tt, efn[,i], method='linear')
  n <- length(X.pred$x)
  xis <- matrix(NA, n, k)
  Xhat <- relist(NA, X.pred$x)
  for(i in 1:n) {
    ti <- X.pred$pp[[i]]
    xic <- X.pred$x[[i]] - muh.pred[[i]]
    phis <- matrix(NA, length(ti), max(k,s))
    for(j in 1:max(k,s))
      phis[,j] <- ff[[j]](ti) # efn[, i] #
    siy <- phis[,1:s1, drop=FALSE] %*% diag( lam[1:s1] ) %*% t(phis[,1:s1, drop=FALSE])
    rhs <- as.vector( solve(siy + rho * diag(length(ti)), xic ) )
    if(k > 1) {
      xis[i,] <- t( phis[,1:k, drop=FALSE] %*% diag( lam[1:k] ) ) %*% rhs
    } else {
      xis[i,] <- t( phis[,1, drop=FALSE] * lam[1] ) %*% rhs
    }
    # Xhat[[i]] <- as.vector( phis[,1:k, drop=FALSE] %*% as.vector( xis[i,] ) ) + muh
    Xhat[[i]] <- as.vector( efn[, 1:k, drop=FALSE] %*% as.vector( xis[i,] ) ) + muh
  }
  return(Xhat)
}



#' @export
pred.scores <- function(X, muh, X.pred, muh.pred, cov.fun, tt, k=2, s=20, rho=0) {
  # predicted scores based on cov.fun
  # but for trajectories not present in X (those in X.pred)
  # X: training set
  # muh: estimated mean for training set
  # X.pred: test set
  # muh.pred: estimated mean for test set

  eg <- eigen(cov.fun)
  lam <- eg$values[1:max(k,s)]
  ef <- eg$vectors[,1:max(k,s)]
  # eg <- svd(cov.fun)
  # lam <- eg$d[1:max(k,s)]
  # ef <- eg$u[,1:max(k,s)]
  lam[ lam < 0 ] <- 0
  s1 <- max( (1:max(k,s))[ lam > 1e-5 ] )
  normas <- apply(ef, 2, L2.norma.mesh, mesh=tt) # rep(1, max(k,s))
  efn <- scale(ef, center=FALSE, scale=normas)
  ff <- vector('list', max(k,s))
  for(i in 1:max(k,s))
    ff[[i]] <- approxfun(tt, efn[,i], method='linear')
  n <- length(X.pred$x)
  xis <- matrix(NA, n, k)
  # Xhat <- relist(NA, X.pred$x)
  for(i in 1:n) {
    ti <- X.pred$pp[[i]]
    xic <- X.pred$x[[i]] - muh.pred[[i]]
    phis <- matrix(NA, length(ti), max(k,s))
    for(j in 1:max(k,s))
      phis[,j] <- ff[[j]](ti)
    siy <- phis[,1:s1, drop=FALSE] %*% diag( lam[1:s1] ) %*% t(phis[,1:s1, drop=FALSE])
    rhs <- as.vector( solve(siy + rho * diag(length(ti)), xic ) )
    if(k > 1) {
      xis[i,] <- t( phis[,1:k, drop=FALSE] %*% diag( lam[1:k] ) ) %*% rhs
    } else {
      xis[i,] <- t( phis[,1, drop=FALSE] * lam[1] ) %*% rhs
    }
    # Xhat[[i]] <- as.vector( phis[,1:k, drop=FALSE] %*% as.vector( xis[i,] ) ) + muh.pred[[i]]
  }
  return(xis)
}




#' @export
mysparseWiener <- function(n, pts, K, npc, mean.f=function(a) 0, mu.c=0, eps=0) {
  # outliers have 2nd and 3rd scores contaminated
  # n = number of curves
  # K = number of components to use to generate
  # pts = vector of time points where to evaluate
  # npc = vector of integers with the possible number of
  # observations per curve (will be uniform among these)
  # the curves are observed at a subset of points randomly drawn from pts
  phis <- function(a, j) sin(pi*(2*j-1)*a/2)*sqrt(2)
  las <- 1/(pi*(2*(1:K)-1)) # (las = sds of scores)
  tmp <- vector('list', 4)
  names(tmp) <- c('x', 'pp', 'xis', 'lambdas')
  tmp$x <- tmp$pp <- vector('list', n)
  tmp$xis <- matrix(NA, n, K)
  outs <- rbinom(n, size=1, prob=eps)
  a <- outer(pts, 1:K, FUN=phis)
  for(j in 1:n) {
    ii <- sort(sample(length(pts), sample(npc, 1)))
    tmp$xis[j, ] <- rnorm(K, mean=0, sd=las)
    if( outs[j] == 1) tmp$xis[j, 2:3] <- rnorm(2, mean=mu.c, sd=las)
    tmp$x[[j]] <- as.vector( a[ii,, drop=FALSE] %*% tmp$xis[j, ] ) + mean.f(pts[ii])
    tmp$pp[[j]] <- pts[ii]
  }
  tmp$lambdas <- las
  return(tmp)
}

#' @export
mysparseWiener2 <- function(n, q, k, mi=0, ma=1, npc=2:5, mean.f = function(a) 0, mu.c=0, eps=0) {
  # outliers have only the 3rd score contaminated
  # n = number of curves
  # q = number of components to use to generate the process
  # k = number of scores and eigenvalues to return
  # curves are observed at a (uniformly) randomly chosen number of points
  # from U(mi, ma)
  # npc = vector of integers with the possible number of
  # observations per curve (will be uniform among these)
  phis <- function(a, j) sin(pi*(2*j-1)*a/2)*sqrt(2)
  las <- 2/(pi*(2*(1:q)-1)) # (las = sds of scores)
  tmp <- vector('list', 5)
  names(tmp) <- c('x', 'pp', 'xis', 'lambdas', 'outs')
  tmp$x <- tmp$pp <- vector('list', n)
  tmp$xis <- matrix(NA, n, q)
  tmp$outs <- outs <- rbinom(n, size=1, prob=eps)
  # a <- outer(pts, 1:K, FUN=phis)
  for(j in 1:n) {
    pps <- sort(runif(sample(npc, 1), min=mi, max=ma))
    tmp$xis[j, ] <- rnorm(q, mean=0, sd=las)
    if( outs[j] == 1) tmp$xis[j, 3] <- rnorm(1, mean=mu.c, sd=1) * las[3]
    tmp$x[[j]] <- as.vector( outer(pps, 1:q, FUN=phis) %*% tmp$xis[j, ] ) + mean.f(pps)
    tmp$pp[[j]] <- pps
  }
  tmp$lambdas <- las[1:k]^2
  tmp$phis <- phis
  tmp$xis <- tmp$xis[, 1:k]
  return(tmp)
}

#' @export
mysparseWiener3 <- function(n, q, k, mi=0, ma=1, npc=2:5, mean.f = function(a) 0, mu.c=4, eps=0) {
  # contamination is along the 4th eigenfunction
  # (outliers are multiples of \phi_4)
  # n = number of curves
  # q = number of components to use to generate the process
  # k = number of scores and eigenvalues to return
  # curves are observed at a (uniformly) randomly chosen number of points
  # from U(mi, ma)
  # npc = vector of integers with the possible number of
  # observations per curve (will be uniform among these)
  phis <- function(a, j) sin(pi*(2*j-1)*a/2)*sqrt(2)
  las <- 2/(pi*(2*(1:q)-1)) # (las = sds of scores)
  tmp <- vector('list', 5)
  names(tmp) <- c('x', 'pp', 'xis', 'lambdas', 'outs')
  tmp$x <- tmp$pp <- vector('list', n)
  tmp$xis <- matrix(NA, n, q)
  tmp$outs <- outs <- rbinom(n, size=1, prob=eps)
  # a <- outer(pts, 1:K, FUN=phis)
  for(j in 1:n) {
    pps <- sort(runif(sample(npc, 1), min=mi, max=ma))
    tmp$xis[j, ] <- rnorm(q, mean=0, sd=las)
    tmp$x[[j]] <- as.vector( outer(pps, 1:q, FUN=phis) %*% tmp$xis[j, ] ) + mean.f(pps)
    if( outs[j] == 1) tmp$x[[j]] <- phis(a=pps, j=4) * rnorm(1, mean=mu.c, sd=.1)
    tmp$pp[[j]] <- pps
  }
  tmp$lambdas <- las[1:k]^2
  tmp$phis <- phis
  tmp$xis <- tmp$xis[, 1:k]
  return(tmp)
}





#' @export
covresids.old <- function(X, mh, gam.fit) {
    # compute r_{ijl} = [ X_i( t_{ij} ) - beta(t_{ij}, t_{il}) X_i(t_{il}) ]
    n <- length(X$x)
    # les <- unlist( sapply(X$pp, function(a) nrow(subsets(length(a), 2))*2 ) )
    # les2 <- unlist( sapply(X$x, function(a) nrow(subsets(length(a), 2))*2 ) )
    re <- vector('numeric', 0) # sum(les))
    for (i in 1:n) {
      comb <- subsets(length(X$x[[i]]),2)
      if(class(comb)[1] != 'matrix') comb <- matrix(comb, byrow=TRUE, ncol=2)
      Maux <- matrix(X$x[[i]][comb]-mh[[i]][comb],ncol=2)
      MTaux <- matrix(X$pp[[i]][comb],ncol=2)
      M <- rbind(Maux, cbind(Maux[,2], Maux[,1]))
      # MT <- rbind(MT, MTaux, cbind(MTaux[,2], MTaux[,1]))
      aa <- rbind(MTaux, cbind(MTaux[,2], MTaux[,1]))
      tmp.dat <- list(x1 = aa[,1], x2 = aa[,2])
      tmp <- predict(gam.fit, newdata=tmp.dat)
      tmp.dat2 <- list(x1 = tmp.dat$x2, x2 = tmp.dat$x2)
      tmp2 <- predict(gam.fit, newdata=tmp.dat2)
      re <- c(re, as.vector( M[,1] - M[,2] * tmp / tmp2) )
      # print(c(i, les[i], les2[i], length(tmp)))
    }
    return(re)
  }

#' @export
covresids <- function(X, mh, gam.fit) {
  # compute r_{ijl} = [ X_i( t_{ij} ) - beta(t_{ij}, t_{il}) X_i(t_{il}) ]
  n <- length(X$x)
  M <- MT <- NULL
  for (i in 1:n) {
    comb <- subsets(length(X$x[[i]]),2)
    if(class(comb)[1] != 'matrix') comb <- matrix(comb, byrow=TRUE, ncol=2)
    Maux <- matrix(X$x[[i]][comb]-mh[[i]][comb],ncol=2)
    MTaux <- matrix(X$pp[[i]][comb],ncol=2)
    M <- rbind(M, Maux, cbind(Maux[,2], Maux[,1]))
    MT <- rbind(MT, MTaux, cbind(MTaux[,2], MTaux[,1]))
  }
  tmp.dat <- list(x1 = MT[,1], x2 = MT[,2])
  tmp <- predict(gam.fit, newdata=tmp.dat)
  tmp.dat2 <- list(x1 = tmp.dat$x2, x2 = tmp.dat$x2)
  tmp2 <- predict(gam.fit, newdata=tmp.dat2)
  tmp.dat3 <- list(x1 = tmp.dat$x1, x2 = tmp.dat$x1)
  tmp3 <- predict(gam.fit, newdata=tmp.dat3)
  tmp4 <- sqrt( tmp3 - tmp^2 / tmp2 )
  summary(tmp4)
  re <- as.vector( (M[,1] - M[,2] * tmp / tmp2) / tmp4 )
  return(re)
}



#' @export
cov.fun.cv.res.par <- function(X, muh, ncov, k.cv, hs, seed=123) {
  # parallel computing version
  # CV for the covariance function
  # muh = estimated mean function
  mii <- min( ti <- unlist(X$pp) )
  maa <- max( ti )
  tt <- seq(mii, maa, length=ncov)
  lh <- length(hs)
  n <- length(X$x)
  if(exists(".Random.seed", where=.GlobalEnv)) old.seed <- .Random.seed
  set.seed(seed)
  fl <- sample( (1:n) %% k.cv + 1)
  tmspe <- rep(NA, lh)
  Xtest <- Xtrain <- vector('list', 2)
  names(Xtrain) <- c('x', 'pp')
  names(Xtest) <- c('x', 'pp')
  tmp.par <- foreach(h=hs, .combine=c, .inorder=FALSE, .packages=c('MASS', 'mgcv'), .export=c('tm', 'cov.fun.hat2', 'matrixx', 'subsets',
                                'gtthat', 'gtshat', 'mu.hat3', 'rho',
                               'k.epan', 'psi', 'pred.cv', 'covresids',
                               'L2.norma.mesh', 'L2.dot.product.mesh',
                               'integral')) %dopar% {
                                 re <- vector('numeric', 0)
                                 for(i in 1:k.cv) {
                                   this <- (1:n)[ fl != i ] # training set
                                   Xtrain$x <- X$x[ this ]
                                   Xtrain$pp <- X$pp[ this ]
                                   Xtest$x <- X$x[ -this ]
                                   Xtest$pp <- X$pp[ -this ]
                                   ma <- matrixx(Xtrain, muh[ this ])
                                   cov.fun <- try(cov.fun.hat2(X=Xtrain, h=h, mh=muh[ this ],
                                                               ma=ma, ncov=ncov, trace=FALSE)) # $G
                                   if( !inherits(cov.fun, 'try-error') ) {
                                     if(!any(is.na(cov.fun$G))) {
                                       uu <- as.vector(cov.fun$G) #
                                       ttx <- cov.fun$grid #
                                       x1 <- ttx[,1]
                                       x2 <- ttx[,2]
                                       gam.fit <- gam(uu ~ s(x1, x2), family=gaussian())
                                       re <- c(re, covresids(X=Xtrain, mh=muh[ this ], gam.fit=gam.fit) )
                                     }
                                   }
                                 }
                                 # tmp <- sapply(mapply('-', X$x, Xhat), function(a) a^2) # squared residuals, list-wise
                                 if( any(is.na(re)) | (length(re)==0) ) {
                                   tmspe <- NA } else {
                                     tmp3 <- RobStatTM::locScaleM(tmp2, psi='bisquare')
                                     tmspe <- tmp3$disper^2 + tmp3$mu^2 # RobStatTM::mscale(re) #, delta=.3, tuning.chi=2.560841)
                                   }
                               }
  if(exists('old.seed')) assign('.Random.seed', old.seed, envir=.GlobalEnv)
  return(list(tmspe=tmp.par, h=hs, ncov=ncov))
}

#' @export
covresids.new <- function(X, mh, h, seed=123) {
  # compute r_{ijl} = [ X_i( t_{ij} ) - beta(t_{ij}, t_{il}) X_i(t_{il}) ]
  # gtshat <-function(X,t0,s0,h,mh,matx,cc=1.56,eps=1e-6){
  n <- length(X$x)
  M <- MT <- NULL
  X.cv <- vector('list', 2)
  names(X.cv) <- c('x', 'pp')
  re <- vector('numeric', 0)
  for (i in 1:n) {
    comb <- subsets(length(X$x[[i]]),2)
    if(class(comb)[1] != 'matrix') comb <- matrix(comb, byrow=TRUE, ncol=2)
    Maux <- matrix(X$x[[i]][comb]-mh[[i]][comb],ncol=2)
    MTaux <- matrix(X$pp[[i]][comb],ncol=2)
    sigmas <- betas <- rep(NA, nrow(Maux))
    X.cv$x <- X$x[-i]
    X.cv$pp <- X$pp[-i]
    ma.cv <- matrixx(X.cv, mh[-i])
    MT <- ma.cv$mt
    M <- ma.cv$m
    # set.seed(seed + 13*i)
    # j <- sample(nrow(MTaux), 1)
    for(j in 1:nrow(MTaux)) {
      t0 <- MTaux[j, 1]
      s0 <- MTaux[j, 2]
      betas[j] <- gtshat(X=X.cv, t0=t0, s0=s0, h=h, mh=mh[-i], matx=ma.cv)
      we <- k.epan((MT[,1]-t0)/h) * k.epan((MT[,2]-s0)/h)
      # B <- med.w( x=M[,2] / M[,1 ], w=we)
      sigmas[j] <- mad.w( x=M[,2] - betas[j] * M[,1], w=we)
    }
    re <- c(re, (Maux[,1] - Maux[,2] * betas) / sigmas )
    # re <- c(re, (Maux[j ,1] - Maux[j ,2] * betas[j]) / sigmas[j] )
  }
  return(re)
}


# Elliptical-FPCA main function
#
#

#' @export
efpca <- function(X, ncpus=4, opt.h.mu, opt.h.cov, hs.mu=seq(10, 25, by=1), hs.cov=hs.mu,
                  rho.param=1e-3, alpha=.2, k = 3, s = k, trace=FALSE, seed=123, k.cv=5,
                  ncov=50, max.kappa=1e3) {

  # X is a list of 2 named elements "x", "pp"
  # which contain the lists of observations and times for each item
  # X$x[[i]] and X$pp[[i]] are the data for the i-th individual

  # Start cluster
  if( missing(opt.h.mu) || missing(opt.h.cov) ) {
    cl <- makeCluster(ncpus) #, outfile="kk.txt") # stopCluster(cl)
    registerDoParallel(cl)
  }
  # run CV to find smoothing parameters for
  # mean and covariance function
  if(missing(opt.h.mu)) {
    aa <- cv.mu.par(X, alpha=alpha, hs=hs.mu, seed=seed, k.cv=k.cv)
    opt.h.mu <- aa$h[ which.min(aa$tmse) ]
  }
  mh <- mu.hat3.lin(X=X, h=opt.h.mu)
  if(missing(opt.h.cov)) {
    bb <- cov.fun.cv.par(X=X, muh=mh, ncov=ncov, k.cv=k.cv, hs=hs.cov,
                         alpha=alpha, seed=seed, k=k, s=s, reg.rho=rho.param)[1:2]
    opt.h.cov <- bb$h[ which.min(bb$tmspe) ]
  }
  if( exists('cl', inherits=FALSE) ) {
    stopCluster(cl)
  }
  ma <- matrixx(X, mh)
  # Compute the estimated cov function
  cov.fun2 <- cov.fun.hat2(X=X, h=opt.h.cov, mh=mh, ma=ma, ncov=ncov, trace=FALSE)
  # smooth it
  yy <- as.vector(cov.fun2$G)
  xx <- cov.fun2$grid
  tmp <- fitted(mgcv::gam(yy ~ s(xx[,1], xx[,2]), family='gaussian'))
  cov.fun2$G <- matrix(tmp, length(unique(xx[,1])), length(unique(xx[,1])))
  cov.fun2$G <- ( cov.fun2$G + t(cov.fun2$G) ) / 2
  ours <- list(mh=mh, ma=ma, cov.fun2 = cov.fun2)

  # predicted scores, fitted values
  Xpred.fixed <- Xpred <- pred(X=X, muh=ours$mh, cov.fun=ours$cov.fun2$G,
                tt=unique(ours$cov.fun2$grid[,1]),
                ss=unique(ours$cov.fun2$grid[,1]), k=k, s=s, rho=rho.param)

  # # select rho
  # n <- length(X$x)
  # tt.grid <- unique(ours$cov.fun2$grid[,1])
  # sigma2.1 <- vector('list', n) #rep(NA, n)
  # for(j in 1:n)
  #   sigma2.1[[j]] <- (X$x[[j]] - approx(x=tt.grid, y=Xpred$pred[j, ], xout=X$pp[[j]], method='linear')$y)
  # # sigma2.1 <- RobStatTM::mscale(unlist(sigma2.1))^2 #, delta=.3, tuning.chi=2.560841))
  # sigma2.1 <- mscale.long(sigma2.1)^2
  # Xpred <- pred(X=X, muh=ours$mh, cov.fun=ours$cov.fun2$G,
  #               tt=unique(ours$cov.fun2$grid[,1]),
  #               ss=unique(ours$cov.fun2$grid[,1]), k=k, s=s, rho=sigma2.1)
  # sigma2.2 <- vector('list', n) #rep(NA, n)
  # for(j in 1:n)
  #   sigma2.2[[j]] <- (X$x[[j]] - approx(x=tt.grid, y=Xpred$pred[j, ], xout=X$pp[[j]], method='linear')$y)
  # # sigma2.2 <- RobStatTM::mscale(unlist(sigma2.2))^2 #, delta=.3, tuning.chi=2.560841) )
  # sigma2.2 <- mscale.long(sigma2.2)^2
  # rho.param <- sigma2.2

  # select rho with condition number
  # la1 <- svd(ours$cov.fun2$G)$d[1]
  la1 <- eigen(ours$cov.fun2$G)$values[1]
  # rho.param <- uniroot(function(rho, la1, max.kappa) return( (la1+rho)/rho - max.kappa ), la1=la1,
  #                      max.kappa = max.kappa, interval=c(1e-15, 1e15))$root
  rho.param <- la1/(max.kappa-1)
  Xpred <- pred(X=X, muh=ours$mh, cov.fun=ours$cov.fun2$G,
                tt=unique(ours$cov.fun2$grid[,1]),
                ss=unique(ours$cov.fun2$grid[,1]), k=k, s=s, rho=rho.param)
  return(list(cov.fun = ours$cov.fun$G, muh=ours$mh, tt=unique(ours$cov.fun2$grid[,1]),
              ss=unique(ours$cov.fun2$grid[,2]), ma=ma, xis=Xpred$xis, pred=Xpred$pred,
              opt.h.mu=opt.h.mu, opt.h.cov=opt.h.cov, rho.param=rho.param,
              pred.fixed = Xpred.fixed$pred, xis.fixed = Xpred.fixed$xis))
}


# cov.fun.cv.par <- function(X, muh, ncov=50, k.cv=5, hs=exp(seq(-4, 0, by=.6)),
#                            alpha=.4, seed=123, k=2, s=20, reg.rho=1e-5) {
#   # parallel computing version
#   # CV for the covariance function
#   # muh = estimated mean function
#   mii <- min( ti <- unlist(X$pp) )
#   maa <- max( ti )
#   tt <- seq(mii, maa, length=ncov)
#   lh <- length(hs)
#   n <- length(X$x)
#   if(exists(".Random.seed", where=.GlobalEnv)) old.seed <- .Random.seed
#   set.seed(seed)
#   fl <- sample( (1:n) %% k.cv + 1)
#   tmspe <- rep(NA, lh)
#   Xtest <- Xtrain <- vector('list', 2)
#   names(Xtrain) <- c('x', 'pp')
#   names(Xtest) <- c('x', 'pp')
#   tmp.par <- foreach(h=hs, .combine=c, .inorder=FALSE, .packages=c('MASS', 'mgcv'),
#                      .export=c('tm', 'cov.fun.hat2', 'matrixx', 'subsets',
#                                'gtthat', 'gtshat', 'mu.hat3',
#                                'k.epan', 'psi', 'pred.cv',
#                                'L2.norma.mesh', 'L2.dot.product.mesh', 'rho',
#                                'integral')) %dopar% {
#                                  Xhat <- relist(NA, X$x)
#                                  for(i in 1:k.cv) {
#                                    this <- (1:n)[ fl != i ] # training set
#                                    Xtrain$x <- X$x[ this ]
#                                    Xtrain$pp <- X$pp[ this ]
#                                    Xtest$x <- X$x[ -this ]
#                                    Xtest$pp <- X$pp[ -this ]
#                                    ma <- matrixx(Xtrain, muh[ this ])
#                                    # cov.fun <- try(cov.fun.hat2(X=Xtrain, h=h, mh=muh[ this ],
#                                    #                             ma=ma, ncov=50, trace=FALSE)) # $G
#                                    mh <- muh[ this ]
#                                    mii <- min( ti <- unlist(Xtrain$pp) )
#                                    maa <- max( ti )
#                                    tt <- seq(mii, maa, length=ncov)
#                                    ss <- seq(mii, maa, length=ncov)
#                                    pps <- as.matrix(expand.grid(tt, ss))
#                                    np <- nrow(pps)
#                                    betahat <- rep(0, np)
#                                    sigmahat <- rep(0, ncov)
#                                    for(jj in 1:ncov) sigmahat[jj] <- gtthat(X=Xtrain, t0=tt[jj], h=h, muhat=mh)$gtt
#                                    for(jj in 1:np) {
#                                      t0 <- pps[jj, 1]
#                                      s0 <- pps[jj, 2]
#                                      betahat[jj] <- gtshat(X=X$train,t0=t0,s0=s0,h=h,mh=mh,matx=ma,cc=1.56,eps=1e-6)
#                                      # gamma[s0, t0] / gamma[t0, t0]
#                                    }
#                                    G <- betahat <- matrix(betahat, ncov, ncov)
#                                    for(ii in 1:ncov)
#                                      for(jj in 1:ncov)
#                                        G[ii,jj] <- betahat[ii,jj] * sigmahat[ii]^2
#                                    G <- ( G + t(G) ) / 2
#                                    cov.fun <- list(G=G, grid=pps)
#                                    if( class(cov.fun) != 'try-error') {
#                                      if(!any(is.na(cov.fun$G))) {
#                                        uu <- as.vector(cov.fun$G) #
#                                        ttx <- cov.fun$grid #
#                                        cov.fun <- matrix(fitted(gam(uu ~ s(ttx[,1], ttx[,2]))), length(tt), length(tt)) #
#                                        cov.fun <- ( cov.fun + t(cov.fun) ) / 2
#                                        tmp <- try( pred.cv(X=Xtrain, muh=muh[ this ], X.pred=Xtest,
#                                                            muh.pred=muh[ -this ], cov.fun=cov.fun, tt=tt,
#                                                            k=k, s=s, rho=reg.rho) )
#                                        if( class(tmp) != 'try-error') Xhat[ -this ] <- tmp
#                                      }
#                                    }
#                                  }
#                                  # tmp <- sapply(mapply('-', X$x, Xhat), function(a) a^2) # squared residuals, list-wise
#                                  tmp <- mapply('-', X$x, Xhat) # squared residuals, list-wise
#                                  # tmspe <- tm(unlist(tmp), alpha=alpha) # 40% trimming, use 60% smallest resids
#                                  # tmspe <- RobStatTM::mscale(unlist(tmp)) #, delta=.3, tuning.chi=2.560841)
#                                  tmp2 <- unlist(tmp)
#                                  if(any(is.na(tmp2))) {
#                                    tmspe <- NA } else {
#                                      # tmp2 <- tmp2[ !is.na(tmp2) ]
#                                      # if(length(tmp2) > 0) {
#                                      tmspe <- RobStatTM::mscale(tmp2) #, delta=.3, tuning.chi=2.560841)
#                                    }
#                                  # } else {
#                                  #   tmspe <- NA
#                                  # }
#                                }
#   if(exists('old.seed')) assign('.Random.seed', old.seed, envir=.GlobalEnv)
#   return(list(tmspe=tmp.par, h=hs, ncov=ncov, k=k, s=s, rho=reg.rho))
# }


mscale.long <- function(pp, delta = 0.5, tuning.chi = 1.547645, max.it = 100, tol = 1e-06)
{
  # function to compute an M-scale like the longitudinal version of the sample sd
  # solves 1/n \sum_{i=1}^n [ 1/n_i \sum_{j=1}^{n_i} \rho( pp_ij / s1 ) ] = delta
  s0 <- median(abs(unlist(pp)))/0.6745
  err <- tol + 1
  it <- 0
  while ((err > tol) && (it < max.it)) {
    it <- it + 1
    tmp <- mean( sapply(pp, function(a, s0, cc) mean(rho(a/s0, cc = cc)), s0=s0, cc=tuning.chi) )
    s1 <- sqrt( (s0^2 * tmp)/delta )
    err <- abs(s1 - s0)/s0
    s0 <- s1
  }
  return(s0)
}

#' #' @export
#' sparseExp <- function(n, nte, mi=0, ma=1, npc=2:5,
#'                       mean.f = function(a) 0, mu.c=4, eps=0,
#'                       scale = 1, theta = .2*(ma - mi)) {
#'   # contamination ?
#'   # n = number of curves
#'   # nte = size of the U(mi, ma) grid on which to evaluate the process
#'   # (these points will then be sampled for sparsity)
#'   # npc = vector of integers with the possible number of
#'   # observations per curve (will be uniform among these)
#'   # scale, theta = parameters for the covariance
#'   # function which is "scale*exp(-abs(s-t)/theta)"
#'   te <- sort( runif(nte, min=mi, max=ma) )
#'   my.par <- list(scale=scale, theta=theta)
#'   full.dat  <- fda.usc::rproc2fdata(n=n, t=te, sigma="vexponential",
#'                       par.list=my.par)$data
#'   dimnames(full.dat) <- NULL
#'   tmp <- vector('list', 5)
#'   names(tmp) <- c('x', 'pp', 'xis', 'lambdas', 'outs')
#'   tmp$x <- tmp$pp <- vector('list', n)
#'   # tmp$xis <- matrix(NA, n, q)
#'   # tmp$outs <- outs <- rbinom(n, size=1, prob=eps)
#'   for(j in 1:n) {
#'     pps <- sort( sample.int(n=nte, size=sample(npc, 1)) )
#'     tmp$x[[j]] <- as.vector( full.dat[j, pps] ) + mean.f(te[pps])
#'     tmp$pp[[j]] <- te[pps]
#'   }
#'   return(tmp)
#' }

#' @export sparseGaussKernel sparseSquaredExp
sparseGaussKernel <- sparseSquaredExp <- function(n, nte, mi=0, ma=1, npc=2:5,
                      mean.f = function(a) 0, mu.c=30, eps=0,
                      scale = 1, theta = .2*(ma - mi), phi.cont = function(a) NA) {
  # contamination ?
  # n = number of curves
  # nte = size of the U(mi, ma) grid on which to evaluate the process
  # (these points will then be sampled for sparsity)
  # npc = vector of integers with the possible number of
  # observations per curve (will be uniform among these)
  # scale, theta = parameters for the covariance
  # function which is "scale*exp(-abs(s-t)/theta)"

  # contamination is along the eigenfunction in phi.cont
  # with score N(mu.c, 1)
  # eps = proportion of contamination

  # grid of points to evaluate process
  te <- sort( runif(nte, min=mi, max=ma) )
  # G(a, b) = scale^2 * exp( -(a-b)^2 / theta )
  si <- outer(te, te,
              function(a, b, scale, theta) scale^2 * exp( -(a-b)^2 / theta ),
              scale=scale, theta=theta)
  # a square root of si using SVD
  si.svd <- svd(si)
  a <- si.svd$v %*% diag( sqrt( si.svd$d ) )
  # generate n Gaussian vectors with cov matrix si
  p <- length(te)
  full.dat <- matrix(rnorm(n*p), n, p) %*% t(a)
  tmp <- vector('list', 5)
  names(tmp) <- c('x', 'pp', 'xis', 'lambdas', 'outs')
  tmp$x <- tmp$pp <- vector('list', n)
  # tmp$xis <- matrix(NA, n, q)
  tmp$outs <- outs <- rbinom(n, size=1, prob=eps)
  for(j in 1:n) {
    pps <- sort( sample.int(n=nte, size=sample(npc, 1)) )
    tmp$x[[j]] <- as.vector( full.dat[j, pps] ) + mean.f(te[pps])
    if( outs[j] == 1 ) tmp$x[[j]] <- phi.cont(te[pps]) * rnorm(1, mean=mu.c, sd=1)
    tmp$pp[[j]] <- te[pps]
  }
  return(tmp)
}

#' @export sparseGaussKernel2 sparseSquaredExp2
sparseGaussKernel2 <- sparseSquaredExp2 <- function(n, nte, mi=0, ma=1, npc=2:5,
                                                  mean.f = function(a) 0,
                                                  mu.c.1 = 5, mu.c.2 = 5, sigma, eps=0,
                                                  scale = 1, theta = .2*(ma - mi),
                                                  phis) {
  # contamination ?
  # n = number of curves
  # nte = size of the U(mi, ma) grid on which to evaluate the process
  # (these points will then be sampled for sparsity)
  # npc = vector of integers with the possible number of
  # observations per curve (will be uniform among these)
  # scale, theta = parameters for the covariance function which is "scale*exp(-abs(s-t)/theta)"
  # phis = list of first 2 eigenfunctions
  # contamination is on z1 * phis[[1]] + z2 * phis[[2]]
  # where (z1, z2)^T ~ N((mu.c.1, mu.c.2)^T, \Sigma)
  # eps = proportion of contamination

  # grid of points to evaluate process
  te <- sort( runif(nte, min=mi, max=ma) )
  # G(a, b) = scale^2 * exp( -(a-b)^2 / theta )
  si <- outer(te, te,
              function(a, b, scale, theta) scale^2 * exp( -(a-b)^2 / theta ),
              scale=scale, theta=theta)
  # a square root of si using SVD
  si.svd <- svd(si)
  a <- si.svd$v %*% diag( sqrt( si.svd$d ) )
  # generate n Gaussian vectors with cov matrix si
  # p <- length(te)
  full.dat <- matrix(rnorm(n*nte), n, nte) %*% t(a)
  tmp <- vector('list', 5)
  names(tmp) <- c('x', 'pp', 'xis', 'lambdas', 'outs')
  tmp$x <- tmp$pp <- vector('list', n)
  k <- length(phis) # we compute as many scores as eigenfunctions we have
  tmp$xis <- matrix(NA, n, k)
  # to compute scores we need an equispaced grid te...
  # or maybe we can use L2 inner products
  # (but to check this numerically we need the L2 eigenvalues with
  # a non-equispaced grid, which seems difficult)
  tmp$outs <- outs <- rbinom(n, size=1, prob=eps)
  for(j in 1:n) {
    for(e in 1:k) {
      tmp$xis[j, e] <- L2.dot.product.mesh(full.dat[j, ], phis[[e]](te), te)
    }
    pps <- sort( sample.int(n=nte, size=sample(npc, 1)) )
    tmp$x[[j]] <- as.vector( full.dat[j, pps] ) + mean.f(te[pps])
    if( outs[j] == 1 ) {
      z <- t( chol(sigma) ) %*% rnorm(2) + c(mu.c.1, mu.c.2)
      tmp$x[[j]] <- z[1] * phis[[1]](te[pps]) + z[2] * phis[[2]](te[pps])
    }
    tmp$pp[[j]] <- te[pps]
  }
  return(tmp)
}
