
data(aids, package='catdata')
X <- vector('list', 2) 
names(X) <- c('x', 'pp')
X$x <- split(aids$cd4, aids$person)
X$pp <- split(aids$time, aids$person)

n <- length(X$x)
shorts <- vector('logical', n)
for(i in 1:n) {
  tmp <- (X$pp[[i]] >= 0)
  X$pp[[i]] <- (X$pp[[i]])[tmp]
  X$x[[i]] <- (X$x[[i]])[tmp]
  if( length(X$pp[[i]]) <= 1 ) shorts[i] <- TRUE
}
X$x <- X$x[!shorts]
X$pp <- X$pp[!shorts]

# length(X$x)
# summary(lens <- sapply(X$x, length))
# table(lens)
# 
# xmi <- min( tmp <- unlist(X$x) )
# xma <- max( tmp )
# ymi <- min( tmp <- unlist(X$pp) )
# yma <- max( tmp ) 
# n <- length(X$x)
# plot(seq(ymi, yma, length=5), seq(xmi, xma,length=5), type='n', xlab='t', ylab='X(t)')
# for(i in 1:n) { lines(X$pp[[i]], X$x[[i]], col='gray', lwd=1, type='b', pch=19, 
#                       cex=1) }
# lens <- sapply(X$x, length)
# set.seed(22)
# ii <- c(sample((1:n)[lens==2], 1), sample((1:n)[lens==5], 1), 
#         sample((1:n)[lens==10], 1))
# for(i in ii) lines(X$pp[[i]], X$x[[i]], col='black', lwd=4, type='b', pch=19, 
#                    cex=1, lty=1)

library(sparseFPCA)
library(doParallel)
library(fdapace)

ncpus <- 4
seed <- 123
rho.param <- 1e-3 
max.kappa <- 1e3
ncov <- 50
k.cv <- 10
k <- 5
s <- k 
hs.mu <- seq(.1, 1.5, by=.1)
hs.cov <- seq(1, 7, length=10)

ours.ls <- lsfpca(X=X, ncpus=ncpus, hs.mu=hs.mu, hs.cov=hs.cov, rho.param=rho.param, 
                  k = k, s = k, trace=FALSE, seed=seed, k.cv=k.cv, ncov=ncov,
                  max.kappa=max.kappa)
ours.r <- efpca(X=X, ncpus=ncpus, hs.mu=hs.mu, hs.cov=hs.cov, rho.param=rho.param,
                alpha=0.2, k = k, s = k, trace=FALSE, seed=seed, k.cv=k.cv, ncov=ncov,
                max.kappa=max.kappa)
myop <- list(error=FALSE, methodXi='CE', dataType='Sparse', 
             userBwCov = 1.5, userBwMu= .3, kernel='epan', verbose=FALSE, nRegGrid=50)
pace <- FPCA(Ly=X$x, Lt=X$pp, optns=myop)

m.ls <- ours.ls$cov.fun
m.r <- ours.r$cov.fun
m.p <- pace$smoothedCov

isSymmetric(m.ls)
isSymmetric(m.r)
isSymmetric(m.p)

round(eigen(m.ls)$values * .109, 3)
round(eigen(m.r)$values * .109, 3)
round(eigen(m.p)$values * .109, 3)

# round(svd(m.ls)$d[1:20] *.109, 3)
# round(svd(m.r)$d[1:20] *.109, 3)
# round(svd(m.p)$d[1:20] *.109, 3)

m.ls.e <- eigen(m.ls)
tmp <- (m.ls.e$values > 0)
tmp2 <- m.ls.e$vectors[, tmp]
m.ls.1 <- tmp2 %*% diag(m.ls.e$values[tmp]) %*% t(tmp2)

m.r.e <- eigen(m.r)
tmp <- (m.r.e$values > 0)
tmp2 <- m.r.e$vectors[, tmp]
m.r.1 <- tmp2 %*% diag(m.r.e$values[tmp]) %*% t(tmp2)

m.p.e <- eigen(m.p)
tmp <- (m.p.e$values > 0)
tmp2 <- m.p.e$vectors[, tmp]
m.p.1 <- tmp2 %*% diag(m.p.e$values[tmp]) %*% t(tmp2)


# ph1 <- svd(m1)$u
# lam <- svd(m1)$d
# ph2 <- svd(m1)$v
# all.equal(m1, ph1 %*% diag(lam) %*% t(ph2))
# all.equal(m1 %*% ph2, ph1 %*% diag(lam) )
# 
# ph1.e <- eigen(m1, symmetric=TRUE)$vectors
# ph1[1:5, 1:5]
# ph1.e[1:5, 1:5]
# 
# all.equal(
#   as.vector(m1 %*% ph1[,3]) / lam[3],
#   as.vector(ph1[,3])
# )
# 
# p <- dim(m1)[1]
# for(j in 1:p) {
#   print( all(
#     abs( ph1[, j] - ph2[, j] ) < (.Machine$double.eps*100)
#     ) || all(
#       abs( ph1[, j] - ph2[, j] ) < (.Machine$double.eps*100)
#     )
#     )
# }


round(eigen(pace$smoothedCov)$values[1:20]*.109, 3)
pace$lambda



ss <- tt <- ours.r$ss
G.r <- ours.r$cov.fun
filled.contour(tt, ss, G.r, main='ROB')
ss <- tt <- ours.ls$ss
G.ls <- ours.ls$cov.fun
filled.contour(tt, ss, G.ls, main='LS')
ss <- tt <- pace$workGrid
G.pace <- pace$smoothedCov
filled.contour(tt, ss, G.pace, main='PACE')

ss <- tt <- ours.ls$ss
G.ls <- ours.ls$cov.fun
filled.contour(tt, ss, G.ls, main='LS')
filled.contour(tt, ss, m.ls.1, main='Pos LS')

ss <- tt <- ours.r$ss
G.r <- ours.r$cov.fun
filled.contour(tt, ss, G.r, main='ROB')
filled.contour(tt, ss, m.r.1, main='Pos ROB')

ss <- tt <- pace$workGrid
G.pace <- pace$smoothedCov
filled.contour(tt, ss, G.pace, main='PACE')
filled.contour(tt, ss, m.p.1, main='Pos PACE')


