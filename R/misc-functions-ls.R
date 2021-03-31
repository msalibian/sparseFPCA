
#' @export
gtshat.ls <-function(X, t0, s0, h, mh, matx, eps=1e-6){
  n <- length(X$x)
  err <- 1+eps
  j <- 0
  M <- matx$m
  MT <- matx$mt
  w2 <- k.epan((MT[,1]-t0)/h)
  w3 <- k.epan((MT[,2]-s0)/h)
  we <- w2*w3
  M <- M[ we > 0, ]
  # MM <- matrix(NA, nrow(M), 3)
  # MM[, 1] <- M[, 1]
  # MM[, 2] <- M[, 1] * (MT[ we > 0, 1]-t0)
  # MM[, 3] <- M[, 1] * (MT[ we > 0, 2]-s0)
  we <- we[ we > 0]
  if( length(we)==0 ) return(NA)
  # B <- as.numeric( coef( lm( M[, 2] ~ MM - 1, weights = we) )[1] )
  B <- as.numeric( coef( lm(M[ ,2] ~ M[, 1] - 1, weights = we) ) )
  # # get rid of points on the diagonal
  # #   xs <- M[,1]
  # #   ys <- M[,2]
  # #   tmp <- (xs != ys)
  # #   B <- median( M[tmp,2] / M[tmp,1 ])
  # #   sigma.hat <- mad( M[tmp,2] - B * M[tmp,1])
  # B <- median( M[,2] / M[,1 ])
  # sigma.hat <- mad( M[,2] - B * M[,1])
  # # B <- med.w( x=M[,2] / M[,1 ], w=we)
  # # sigma.hat <- mad.w( x=M[,2] - B * M[,1], w=we)
  # while ( ( (j <- j+1) < 1000 ) && (err > eps) ){
  #   w1 <- psi((M[,2]-B*M[,1])/sigma.hat,cc)/(M[,2]-B*M[,1])
  #   w1[ is.na(w1) ] <- 1/sigma.hat # missing values are "psi(0)/0", that lim is 1/sigma.hat
  #   w <- w1 * we
  #   B2 <- sum(w*M[,1]*M[,2])/sum(w*(M[,1]^2))
  #   err <- abs(B2/B - 1)
  #   if( is.na(err) || is.nan(err) ) return(NA)
  #   #   print('gtshat')
  #   #   print(B)
  #   #   print(summary(w1))
  #   #   print(summary(w))
  #   #   print(sigma.hat)
  #   # }
  #   B <- B2
  # }
  return(B)
}


#' @export
gtthat.ls <- function(X, t0, h, muhat, max.it=300, eps=1e-10) {
  # find the scale, full iterations
  # muhat should be a list as returned by mu.hat3
  i <- 0
  err <- 1 + eps
  t <- unlist(X$pp)
  x <- unlist(X$x)
  if(missing(muhat))
    muhat <- mu.hat3.lin.ls(X=X, h=h, ep=eps)
  muhat <- unlist(muhat)
  kerns <- k.epan((t-t0)/h)
  sc <- sqrt( sum(kerns * (x - muhat)^2 ) / (sum(kerns) - 1) )
  # if(missing(initial.sc)) sc <- mad(x - muhat) else sc <- initial.sc
  # while( ( (i <- i+1) < max.it ) && (err > eps) ) {
  #   kerns <- k.epan((t-t0)/h)
  #   sc2 <- sqrt(sc^2 * sum(kerns * rho((x-muhat)/sc,cc)) / (b * sum(kerns)))
  #   err <- abs(sc2/sc - 1)
  #   if( is.na(err) || is.nan(err) ) return(NA)
  #   #   print(summary(kerns))
  #   #   print(sc2)
  #   #   print(sc)
  #   # }
  #   sc <- sc2
  # }
  return(list(gtt=sc, muhat=muhat))
}

#' @export
uhat.ls <- function(X, t0, h=0.1, ep=1e-6, max.it=100){
  x <- unlist(X$x)
  t <- unlist(X$pp)
  # s <- localMAD(X,t0,h)
  # oldu <- (u <- median(x)) + 100*ep
  # it <- 0
  kw <- k.epan((t-t0)/h)
  return( sum( kw*x ) / sum(kw) )
  # while( ((it <- it + 1) < max.it ) && ( (abs(oldu) - u) > ep ) ){
  #   w <- psi((x-u)/s,cc)/(x-u)
  #   w[ is.nan(w) ] <- 1
  #   w <- w*kw
  #   oldu <- u
  #   u <- sum(w*x)/sum(w)
  # }
  # return(u)
}

#' @export
uhat.lin.ls <- function(X, t0, h=0.1, ep=1e-6, max.it=100){
  x <- unlist(X$x)
  t <- unlist(X$pp)
  # s <- localMAD(X,t0,h)
  # oldu <- (u <- median(x)) + 100*ep
  it <- 0
  tt <- cbind(rep(1, length(t)), t-t0)
  wk <- k.epan((t-t0)/h)
  return( coef( lm(x ~ I(t - t0), weights=wk) )[1] )
  #
  #   beta <- rlm(x=tt[wk>0,], y=x[wk>0])$coef
  # while( ((it <- it + 1) < max.it ) ){
  #   re <- as.vector(x - tt %*% beta)/s
  #   w <- ( psi(re,cc)/re )
  #   w[ is.nan(w) ] <- 1
  #   w <- w * wk
  #   beta.n <- solve( t( tt * w ) %*% tt, t(tt * w) %*% x )
  #   if( any( is.na(beta.n) || is.nan(beta.n) ) ) return(NA)
  #   if( sum( (beta.n - beta)^2 ) < ep ) it = max.it
  #   beta <- beta.n
  # }
  # return(beta.n[1])
}

#' @export
mu.hat2.ls <- function(X, h=0.1, ep=1e-6) {
  tt <- unlist(X$pp)
  nt <- length(tt)
  us <- rep(0, nt)
  for(i in 1:nt)
    us[i] <- uhat.ls(X=X, t0=tt[i], h=h, ep=ep)
  return(us)
}

#' @export
mu.hat2.lin.ls <- function(X, h=0.1, ep=1e-6) {
  tt <- unlist(X$pp)
  nt <- length(tt)
  us <- rep(0, nt)
  for(i in 1:nt)
    us[i] <- uhat.lin.ls(X=X, t0=tt[i], h=h, ep=ep)
  return(us)

}

#' @export
mu.hat3.ls <- function(X, h=0.1, ep=1e-6) {
  us=relist(mu.hat2.ls(X=X, h=h, ep=ep),X$x)
  return(us)
  #return(list(x=X$x,pp=X$pp,u=us))
}

#' @export
mu.hat3.lin.ls <- function(X, h=0.1, ep=1e-6) {
  us <- relist(mu.hat2.lin.ls(X=X, h=h, ep=ep), X$x)
  return(us)
  #return(list(x=X$x,pp=X$pp,u=us))
}

#' @export
cov.fun.hat2.ls <- function(X, h, mh, ma, ncov=50, trace=FALSE) {
  # this function uses the diagonal
  if(trace) print("Computing cov function")
  if(missing(ma)) ma <- matrixx(X=X, mh=mh)
  mii <- min( ti <- unlist(X$pp) )
  maa <- max( ti )
  tt <- seq(mii, maa, length=ncov)
  ss <- seq(mii, maa, length=ncov)
  pps <- as.matrix(expand.grid(tt, ss))
  np <- nrow(pps)
  betahat <- rep(0, np)
  sigmahat <- rep(0, ncov)
  for(j in 1:ncov) sigmahat[j] <- gtthat.ls(X=X, t0=tt[j], h=h, muhat=mh)$gtt
  for(j in 1:np) {
    t0 <- pps[j, 1]
    s0 <- pps[j, 2]
    betahat[j] <- gtshat.ls(X=X, t0=t0, s0=s0, h=h, mh=mh, matx=ma, eps=1e-6)
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

#' @export
cv.mu.par.ls <- function(X, k.cv=5, k = k.cv, hs=exp(seq(-4, 0, by=.6)), seed=123) {
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
                     .export=c('uhat.lin.ls', 'k.epan')) %dopar% {
                       Xh <- relist(NA, X$x) # predictions go here
                       for(i in 1:k.cv) {
                         this <- (1:n)[ fl != i ]
                         Xtmp$x <- X$x[ this ] # training obs
                         Xtmp$pp <- X$pp[ this ] # training times
                         ps <- unlist(X$pp[ -this ]) # times where we need to predict
                         xhats <- rep(NA, length(ps))
                         for(l in 1:length(ps)) {
                           tmp2 <- try(uhat.lin.ls(X=Xtmp, t0=ps[l], h=h, ep=1e-6, max.it=100), silent=TRUE)
                           if( class(tmp2) != 'try-error' )
                             xhats[l] <- tmp2
                         }
                         Xh[-this] <- relist(xhats, X$x[ -this ] ) # fill predictions
                       }
                       # tmp <- sapply(mapply('-', X$x, Xh), function(a) a^2) # squared residuals, list-wise
                       tmp <- mapply('-', X$x, Xh)  # squared residuals, list-wise
                       # # tmses <- tm(unlist(tmp), alpha=alpha) # 40% trimming, use 60% smallest resids
                       # tmses <- RobStatTM::mscale(unlist(tmp)) #, delta=.3, tuning.chi=2.560841)
                       tmp2 <- unlist(tmp)
                       if(any(is.na(tmp2))) {
                         tmses <- NA } else {
                           # tmp2 <- tmp2[ !is.na(tmp2) ]
                           # if(length(tmp2) > 0) {
                           tmses <- sqrt(mean(tmp2^2)) #RobStatTM::mscale(tmp2) #, delta=.3, tuning.chi=2.560841)
                         }

                     }
  if(exists('old.seed')) assign('.Random.seed', old.seed, envir=.GlobalEnv)
  return(list(tmse=tmp.par, h=hs))
}


#' @export
cov.fun.cv.par.ls <- function(X, muh, ncov=50, k.cv=5, hs=exp(seq(-4, 0, by=.6)),
                           seed=123, k=2, s=20, reg.rho=1e-5) {
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
                     .export=c('cov.fun.hat2.ls', 'matrixx', 'subsets',
                               'gtthat.ls', 'gtshat.ls', 'mu.hat3.lin.ls',
                               'k.epan', 'pred.cv', 'L2.norma.mesh', 'L2.dot.product.mesh',
                               'integral')) %dopar% {
                                 Xhat <- relist(NA, X$x)
                                 for(i in 1:k.cv) {
                                   this <- (1:n)[ fl != i ] # training set
                                   Xtrain$x <- X$x[ this ]
                                   Xtrain$pp <- X$pp[ this ]
                                   Xtest$x <- X$x[ -this ]
                                   Xtest$pp <- X$pp[ -this ]
                                   ma <- matrixx(Xtrain, muh[ this ])
                                   cov.fun <- try(cov.fun.hat2.ls(X=Xtrain, h=h, mh=muh[ this ],
                                                               ma=ma, ncov=50, trace=FALSE)) # $G
                                   if( class(cov.fun) != 'try-error') {
                                     if(!any(is.na(cov.fun$G))) {
                                       uu <- as.vector(cov.fun$G) #
                                       ttx <- cov.fun$grid #
                                       cov.fun <- matrix(fitted(gam(uu ~ s(ttx[,1], ttx[,2]))), length(tt), length(tt)) #
                                       cov.fun <- ( cov.fun + t(cov.fun) ) / 2
                                       tmp <- try( pred.cv(X=Xtrain, muh=muh[ this ], X.pred=Xtest,
                                                           muh.pred=muh[ -this ], cov.fun=cov.fun, tt=tt,
                                                           k=k, s=s, rho=reg.rho) )
                                       if( class(tmp) != 'try-error') Xhat[ -this ] <- tmp
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
                                     tmspe <- sqrt(mean(tmp2^2)) #RobStatTM::mscale(tmp2) #, delta=.3, tuning.chi=2.560841)
                                   }
                                 # } else {
                                 #   tmspe <- NA
                                 # }
                               }
  if(exists('old.seed')) assign('.Random.seed', old.seed, envir=.GlobalEnv)
  return(list(tmspe=tmp.par, h=hs, ncov=ncov, k=k, s=s, rho=reg.rho))
}

# Elliptical-FPCA main function
#
#
#' @export
lsfpca <- function(X, ncpus=4, opt.h.mu, opt.h.cov, hs.mu=seq(10, 25, by=1), hs.cov=hs.mu,
                  rho.param=1e-5, k = 3, s = k, trace=FALSE, seed=123, k.cv=5, ncov=50,
                  max.kappa=1e3) {

  # X is a list of 2 named elements "x", "pp"
  # which contain the lists of observations and times for each item
  # X$x[[i]] and X$pp[[i]] are the data for the i-th individual

  # Start cluster
  if( missing(opt.h.mu) || missing(opt.h.cov) ) {
    cl <- makeCluster(ncpus) # stopCluster(cl)
    registerDoParallel(cl)
  }
  # run CV to find smoothing parameters for
  # mean and covariance function
  if(missing(opt.h.mu)) {
    aa <- cv.mu.par.ls(X, hs=hs.mu, seed=seed, k.cv=k.cv)
    opt.h.mu <- aa$h[ which.min(aa$tmse) ]
  }
  mh <- mu.hat3.lin.ls(X=X, h=opt.h.mu)
  if(missing(opt.h.cov)) {
    bb <- cov.fun.cv.par.ls(X=X, muh=mh, ncov=ncov, k.cv=k.cv, hs=hs.cov, seed=seed,
                            k=k, s=s, reg.rho=rho.param)[1:2]
    opt.h.cov <- bb$h[ which.min(bb$tmspe) ]
  }
  if( exists('cl', inherits=FALSE) ) {
    stopCluster(cl)
  }
  ma <- matrixx(X, mh)
  # Compute the estimated cov function
  cov.fun2 <- cov.fun.hat2.ls(X=X, h=opt.h.cov, mh=mh, ma=ma, ncov=ncov, trace=FALSE)
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
  # # sigma2.1 <- mean(unlist(sigma2.1)^2) # RobStatTM::mscale(unlist(sigma2.1))^2 #, delta=.3, tuning.chi=2.560841))
  # sigma2.1 <- mean( sapply(sigma2.1, function(a) mean(a^2) ) )
  # Xpred <- pred(X=X, muh=ours$mh, cov.fun=ours$cov.fun2$G,
  #               tt=unique(ours$cov.fun2$grid[,1]),
  #               ss=unique(ours$cov.fun2$grid[,1]), k=k, s=s, rho=sigma2.1)
  # sigma2.2 <- vector('list', n) #rep(NA, n)
  # for(j in 1:n)
  #   sigma2.2[[j]] <- (X$x[[j]] - approx(x=tt.grid, y=Xpred$pred[j, ], xout=X$pp[[j]], method='linear')$y)
  # # sigma2.2 <- mean(unlist(sigma2.2)^2) # RobStatTM::mscale(unlist(sigma2.2))^2 #, delta=.3, tuning.chi=2.560841) )
  # sigma2.2 <- mean( sapply(sigma2.2, function(a) mean(a^2) ) )
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
              xis.fixed = Xpred.fixed$xis, pred.fixed = Xpred.fixed$pred,
              opt.h.mu=opt.h.mu, opt.h.cov=opt.h.cov, rho.param=rho.param))
}



#
# cov.fun.cv.new <- function(X, muh, k.cv, hs, hwide, seed=123) {
#   # CV "in the regression problem formulation"
#   # muh = estimated mean function
#   # mii <- min( ti <- unlist(X$pp) )
#   # maa <- max( ti )
#   # tt <- seq(mii, maa, length=ncov)
#   lh <- length(hs)
#   n <- length(X$x)
#   if(exists(".Random.seed", where=.GlobalEnv)) old.seed <- .Random.seed
#   set.seed(seed)
#   fl <- sample( (1:n) %% k.cv + 1)
#   tmspe <- rep(NA, lh)
#   Xtest <- Xtrain <- vector('list', 2)
#   names(Xtrain) <- c('x', 'pp')
#   names(Xtest) <- c('x', 'pp')
#   mspe1 <- mspe <- rep(NA, lh)
#   for(j in 1:lh) {
#     ress1 <- ress <- vector('numeric', 0)
#     for(i in 1:k.cv) {
#       this <- (1:n)[ fl != i ] # training set
#       Xtrain$x <- X$x[ this ]
#       Xtrain$pp <- X$pp[ this ]
#       Xtest$x <- X$x[ -this ]
#       Xtest$pp <- X$pp[ -this ]
#       ma <- matrixx(Xtrain, muh[ this ])
#       ma2 <- matrixx(Xtest, muh[ -this ])
#       beta.cv <- try(betahat.new.ls(X=Xtrain, h=hs[j], mh=muh[ this ],
#                                  ma=ma, ma2=ma2, trace=FALSE))
#       if( class(beta.cv) != 'try-error') {
#         re1 <- re <- rep(NA, length(beta.cv))
#         for(uu in 1:length(beta.cv)) {
#             w2 <- k.epan((ma$mt[,1]-ma2$mt[uu,1])/hwide)
#             w3 <- k.epan((ma$mt[,2]-ma2$mt[uu,1])/hwide)
#             we <- w2*w3
#             M <- ma$m[ we > 0, ]
#             we <- we[ we > 0]
#             B <- NA
#             if( length(we) > 0 ) B <- as.numeric( coef( lm(M[ ,2] ~ M[, 1] - 1, weights = we) ) )
#             supersigma <- sd( as.numeric( M[,2] - B * M[, 1] ) )
#             re[uu] <- (ma2$m[uu ,2] - beta.cv[uu] * ma2$m[uu , 1])/supersigma
#             re1[uu] <- (ma2$m[uu ,2] - beta.cv[uu] * ma2$m[uu , 1])
#         }
#         ress <- c(ress, re)
#         ress1 <- c(ress1, re1)
#       }
#     }
#     mspe1[j] <- mean( ress1^2 )
#     mspe[j] <- mean( ress^2 )
#     print(c(mspe[j], mspe1[j]))
#   }
#   if(exists('old.seed')) assign('.Random.seed', old.seed, envir=.GlobalEnv)
#   return(list(mspe=mspe, mspe1=mspe1))
# }

# betahat.new.ls <- function(X, h, mh, ma, ma2, trace=FALSE) {
#   nn <- dim(ma2$mt)[1]
#   betahat <- rep(NA, nn)
#   # sigmahat <- rep(NA, nn)
#   # for(j in 1:nn) sigmahat[j] <- gtthat(X=X, t0=tt[j], h=h, muhat=mh)$gtt
#   for(j in 1:nn) {
#     t0 <- ma2$mt[j, 1]
#     s0 <- ma2$mt[j, 2]
#     betahat[j] <- gtshat.ls(X=X,t0=t0,s0=s0,h=h,mh=mh,matx=ma,eps=1e-6)
#     # gamma[s0, t0] / gamma[t0, t0]
#   }
#   return(betahat=betahat)
# }
