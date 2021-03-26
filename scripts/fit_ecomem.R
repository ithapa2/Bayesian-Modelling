library("splines")
library("rstan")

BTP = read.csv("data/BTP-BI.csv", header=TRUE)
CBS = read.csv("data/CBS-BI.csv", header=TRUE)
TOW = read.csv("data/TOW-BI.csv", header=TRUE)

BTP.sub = BTP[apply(BTP[,c('chron', 'tmin.may', 'pcp.aug')], 1, function(x) all(!is.na(x))),]

CBS.sub = CBS[apply(CBS[,c('chron', 'tmin.may', 'pcp.aug')], 1, function(x) all(!is.na(x))),]

TOW.sub = TOW[apply(BTP[,c('chron', 'tmin.may', 'pcp.aug')], 1, function(x) all(!is.na(x))),]

Y  = as.matrix(data.frame(BTP.sub$chron, CBS.sub$chron, TOW.sub$chron))
X  = as.matrix(data.frame(BTP.sub$pcp.aug, CBS.sub$pcp.aug, TOW.sub$pcp.aug))
d.raw  = as.matrix(data.frame(BTP.sub$tmin.may, CBS.sub$tmin.may, TOW.sub$tmin.may))

test.na = cbind(Y, X, d.raw)
idx.na  = apply(test.na, 1, function(x) any(is.na(x)))

Y = Y[!idx.na,]
X = X[!idx.na,]
d.raw = d.raw[!idx.na,]

d = matrix(NA, nrow=nrow(d.raw), ncol=ncol(d.raw))
d[d.raw<(-6.119)] = 1
d[d.raw>= (-6.119)] = 0


N_years = nrow(Y)
N_sites = ncol(Y)
lag = 6

dat = list(N_years = N_years,
           N_sites = N_sites,
           lag = lag,
           Y = t(Y),
           X = t(X),
           d = t(d))

#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

sm<-stan_model("scripts/ecomem.stan")

fit<-sampling(sm,
              data=dat,
              iter=500,
              chains = 1, 
              cores = 1)#,control = list(adapt_delta=0.95))

plot(fit)

names(fit)
post = extract(fit)

phi = 1/post$tau
phi.mean = mean(phi)

spherical <- function(x, phi) {
  
  # if (x <= 0) {
  #   reject("x negative spherical function undefined; found x = ", x)
  # }
  if (x > phi) { g=0 } else {
    g = (1 - 3.0/2.0 * x/phi + 1.0/2.0 * (x/phi)^3) 
  }
  
  return(g)
}

x = seq(1e-6, lag+2, by=0.001)
mem = unlist(lapply(x, function(x){spherical(x, phi.mean)}))

plot(x, mem, xlab="Year", ylab="Antecedent weight")

gamma_1 = post$gamma_1
quantile(gamma_1, c(0.10, 0.5, 0.90))
