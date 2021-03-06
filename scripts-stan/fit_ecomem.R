library("splines")
library("rstan")

#######################################################################################
## read in the data
#######################################################################################

BTP = read.csv("data/BTP-BI.csv", header=TRUE)
CBS = read.csv("data/CBS-BI.csv", header=TRUE)
TOW = read.csv("data/TOW-BI.csv", header=TRUE)

BTP.sub = BTP[apply(BTP[,c('chron', 'tmin.may', 'pcp.aug', 'pdsi.sep')], 1, function(x) all(!is.na(x))),]

CBS.sub = CBS[apply(CBS[,c('chron', 'tmin.may', 'pcp.aug', 'pdsi.sep')], 1, function(x) all(!is.na(x))),]

TOW.sub = TOW[apply(BTP[,c('chron', 'tmin.may', 'pcp.aug', 'pdsi.sep')], 1, function(x) all(!is.na(x))),]

fire.raw = data.frame(year = seq(1650, 2017, by=1), fire = rep(0))
fire.raw[which(fire.raw$year %in% c(1664, 1804, 1900)), 'fire'] = 1
fire.raw = fire.raw[which(fire.raw$year %in% BTP.sub$year),]

Y  = as.matrix(data.frame(BTP.sub$chron, CBS.sub$chron, TOW.sub$chron))
X1  = as.matrix(data.frame(BTP.sub$pcp.aug, CBS.sub$pcp.aug, TOW.sub$pcp.aug))
X2  = as.matrix(data.frame(BTP.sub$pdsi.sep, CBS.sub$pdsi.sep, TOW.sub$pdsi.sep))
d.raw  = as.matrix(data.frame(BTP.sub$tmin.may, CBS.sub$tmin.may, TOW.sub$tmin.may))
fire.raw = as.matrix(fire.raw$fire)

test.na = cbind(Y, X1, X2, d.raw, fire.raw)
idx.na  = apply(test.na, 1, function(x) any(is.na(x)))

Y = Y[!idx.na,]
X1 = X1[!idx.na,]
X2 = X2[!idx.na,]
d.raw = d.raw[!idx.na,]
fire = fire.raw[!idx.na,]

d = matrix(NA, nrow=nrow(d.raw), ncol=ncol(d.raw))
d[d.raw<(-6.119)] = 1
d[d.raw>= (-6.119)] = 0
# 
# d.raw.hi = max(d.raw)
# d.raw.lo = min(d.raw)
# d = apply(d.raw, 2, function(x){(x-d.raw.lo)/(d.raw.hi -d.raw.lo)})

#d = d.raw

Y = t(Y)

N_years = ncol(Y)
N_sites = nrow(Y)

X1 = t(X1)
X2 = t(X2)

X = array(NA, c(N_sites, N_years, 2))
X[,,1] = X1
X[,,2] = X2

d = t(d)

lag = 6
N_covars = 2

#######################################################################################
## data list
#######################################################################################

dat = list(N_years = N_years,
           N_sites = N_sites,
           lag = lag,
           N_covars = N_covars,
           Y = Y,
           X = X1,
           d = d,
           fire = fire)

#######################################################################################
## sampling
#######################################################################################

#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

sm<-stan_model("scripts/ecomem.stan")

fit<-sampling(sm,
              data=dat,
              iter=1000,
              chains = 1, 
              cores = 1)#,control = list(adapt_delta=0.95))


#######################################################################################
## output
#######################################################################################

spherical <- function(x, phi) {
  
  # if (x <= 0) {
  #   reject("x negative spherical function undefined; found x = ", x)
  # }
  if (x > phi) { g=0 } else {
    g = (1 - 3.0/2.0 * x/phi + 1.0/2.0 * (x/phi)^3) 
  }
  
  return(g)
}


plot(fit)

names(fit)
post = extract(fit)

phi = 1/post$tau
phi.mean = mean(phi)

x = seq(1e-6, lag+2, by=0.001)
mem = unlist(lapply(x, function(x){spherical(x, phi.mean)}))

plot(x, mem, xlab="Year", ylab="Antecedent weight")

phi_fire = 1/post$tau_fire
phi_fire.mean = mean(phi_fire)

x = seq(1e-6, lag+2, by=0.001)
mem.fire = unlist(lapply(x, function(x){spherical(x, phi_fire.mean)}))

plot(x, mem.fire, xlab="Year", ylab="Antecedent weight")


gamma = cbind(post$gamma0, post$gamma1, post$gamma2)
#quantile(gamma1, c(0.10, 0.5, 0.90))
gamma.quants = apply(gamma, 2, quantile, c(0.05, 0.5, 0.95))
gamma.quants

beta = post$beta
beta.quants = apply(beta, 2, quantile, c(0.05, 0.5, 0.95))
beta.quants
