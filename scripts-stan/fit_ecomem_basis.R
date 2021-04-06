library("mgcv")
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

# d = matrix(NA, nrow=nrow(d.raw), ncol=ncol(d.raw))
# d[d.raw<(-6.119)] = 1
# d[d.raw>= (-6.119)] = 0
#
# d.raw.hi = max(d.raw)
# d.raw.lo = min(d.raw)
# d = apply(d.raw, 2, function(x){(x-d.raw.lo)/(d.raw.hi -d.raw.lo)})

d = d.raw

Y = t(Y)

N_years = ncol(Y)
N_sites = nrow(Y)

X1 = t(X1)
X2 = t(X2)

N_covars = 1

X = array(NA, c(N_sites, N_years, N_covars))

X[,,1] = X1
#X[,,2] = X2

d = t(d)

lag = 6

#######################################################################################
## splines
#######################################################################################

# library(splines)
# set.seed(1234)
# num_knots <- lag # true number of knots
# spline_degree <- 3
# num_basis <- num_knots + spline_degree - 1
# knots <- seq(from=0, to=lag, by=1)
# B <- t(bs(knots, knots=knots, df=num_basis, degree=spline_degree, intercept = TRUE))
# num_data = length(X.basis)


t.s = (0:lag)/lag
time = data.frame(t=0:lag,t.s=t.s)
n.knots = lag + 1
foo = mgcv::s(t.s,k=n.knots,bs="cr")

CRbasis = mgcv::smoothCon(foo,
                          data=time,
                          knots=NULL,
                          absorb.cons=TRUE,
                          scale.penalty=TRUE)

RE = diag(ncol(CRbasis[[1]]$S[[1]]))
B = CRbasis[[1]]$X
S = CRbasis[[1]]$S[[1]] +(1E-07)*RE
S_inv = solve(S)

n_basis = ncol(B)
n_knots = nrow(B)


# L = lag
# ### Memory function inputs ###
# # Define basis functions
# bf = list()
# for (j in 1:length(L)){
#   t.s = (0:L[j])/L[j]
#   time = data.frame(t=0:L[j],t.s=t.s)
#   n.knots = L[j] + 1
#   CRbasis = mgcv::smoothCon(mgcv::s(t.s,k=n.knots,bs="cr"),data=time,knots=NULL,absorb.cons=TRUE,
#                             scale.penalty=TRUE)
#   RE = diag(ncol(CRbasis[[1]]$S[[1]]))
#   bf[[j]] = list(S=CRbasis[[1]]$S[[1]]+(1E-07)*RE,
#                  H=CRbasis[[1]]$X,
#                  k=ncol(CRbasis[[1]]$X),
#                  U=chol(CRbasis[[1]]$S[[1]]+(1E-07)*RE))
# }
# 
# B = bf[[1]]$H
# 
# num.basis = 
# num.data = 
# uni.wts = rep(1/(inputs$L[j]+1),inputs$L[j]+1)
# eta = rnorm(inputs$bf[[j]]$k,coef(lm(uni.wts~inputs$bf[[j]]$H-1)),0.1)
#######################################################################################
## data list
#######################################################################################

dat = list(N_years = N_years,
           N_sites = N_sites,
           lag = lag,
           N_covars = N_covars,
           Y = Y,
           X = X,
           d = d,
           fire = fire,
           B = B,
           S_inv = S_inv,
           n_basis = n_basis,
           n_knots = n_knots)

#######################################################################################
## sampling
#######################################################################################

N_iter = 1000

#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

sm<-stan_model("scripts-stan/ecomem_basis.stan")

fit<-sampling(sm,
              data=dat,
              iter=N_iter,
              chains = 1, 
              cores = 1)#,control = list(adapt_delta=0.95))

saveRDS(fit, 'scripts-stan/output/fit_ecomem_basis.RDS')

