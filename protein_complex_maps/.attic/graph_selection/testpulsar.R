library(pulsar)
library(huge)


options(echo=TRUE) # if you want see commands in output file
args <- commandArgs()
print(args)

print(args[1])

set.seed(10010)
p <- 40 ; n <- 1200
dat  <- huge.generator(n, p, "hub", verbose=FALSE, v=.1, u=.3)
lams <- getLamPath(.2, .01, len=40)


hugeargs <- list(lambda=lams, verbose=FALSE)
time1    <- system.time(
                        out.p    <- pulsar(dat$data, fun=huge, fargs=hugeargs, rep.num=20, 
                                                              criterion='stars', seed=10010)
                                    )
fit.p    <- refit(out.p)

out.p
# Mode: serial
# Path length: 40 
# Subsamples:  20 
# Graph dim:   40 
# Criterion:
#   stars... opt: index 15, lambda 0.132
fit.p
# Pulsar-selected refit of huge 
# Path length: 40 
# Graph dim:   40 
# Criterion:
#   stars... sparsity 0.0325

time2 <- system.time(
                     out.b <- pulsar(dat$data, fun=huge, fargs=hugeargs, rep.num=20, criterion='stars',
                                                     lb.stars=TRUE, ub.stars=TRUE, seed=10010))

time2[[3]] < time1[[3]]
# [1] TRUE
opt.index(out.p, 'stars') == opt.index(out.b, 'stars')
# [1] TRUE


library(QUIC)
quicr <- function(data, lambda) {
        S    <- cov(data)
    est  <- QUIC::QUIC(S, rho=1, path=lambda, msg=0, tol=1e-2)
        path <-  lapply(seq(length(lambda)), function(i) {
                                        tmp <- est$X[,,i]; diag(tmp) <- 0
                                                        as(tmp!=0, "lgCMatrix")
                                            })
        est$path <- path
            est
}



quicargs <- list(lambda=lams)
out.q <- pulsar(dat$data, fun=quicr, fargs=quicargs, rep.num=100, criterion='stars',
                                lb.stars=TRUE, ub.stars=TRUE, ncores=2, seed=10010)


library(BatchJobs)
out.batch <- batch.pulsar(dat$data, fun=quicr, fargs=quicargs, rep.num=100,
                                                    criterion='stars', seed=10010
                                                                              #, cleanup=TRUE
                                                                             )

opt.index(out.q, 'stars') == opt.index(out.batch, 'stars')

out.bbatch <- update(out.batch, criterion=c('stars', 'gcd'),
                                          lb.stars=TRUE, ub.stars=TRUE)
