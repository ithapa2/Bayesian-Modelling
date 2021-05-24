library("mgcv")
library("rstan")

#######################################################################################
## read in the data
#######################################################################################

BTP = read.csv("data/BTP-BI.csv", header=TRUE)
CBS = read.csv("data/CBS-BI.csv", header=TRUE)
TOW = read.csv("data/TOW-BI.csv", header=TRUE)
# 
# BTP.sub = BTP[apply(BTP[,c('chron', 'tmin.may', 'pcp.aug', 'pdsi.sep')], 1, function(x) all(!is.na(x))),]
# CBS.sub = CBS[apply(CBS[,c('chron', 'tmin.may', 'pcp.aug', 'pdsi.sep')], 1, function(x) all(!is.na(x))),]
# TOW.sub = TOW[apply(BTP[,c('chron', 'tmin.may', 'pcp.aug', 'pdsi.sep')], 1, function(x) all(!is.na(x))),]

fire.raw = data.frame(year = seq(1650, 2017, by=1), fire = rep(0))
fire.raw[which(fire.raw$year %in% c(1664, 1804, 1900)), 'fire'] = 1

## find the start year and end year
## for now require that chronologies and fire record all be the same length with no NA values
## climate data of shorter length NA values will be imputed
year_upper = min(max(BTP$year, na.rm=TRUE), max(CBS$year, na.rm=TRUE), max(TOW$year, na.rm=TRUE))
year_lower = max(min(BTP$year, na.rm=TRUE), min(CBS$year, na.rm=TRUE), min(TOW$year, na.rm=TRUE), min(fire.raw$year[which(fire.raw$fire==1)]))

# subset site and fire data to selected years
BTP.sub = BTP[which(BTP$year %in% seq(year_lower, year_upper)),]
CBS.sub = CBS[which(CBS$year %in% seq(year_lower, year_upper)),]
TOW.sub = TOW[which(TOW$year %in% seq(year_lower, year_upper)),]

fire.raw = fire.raw[which(fire.raw$year %in% BTP.sub$year),]

## define data objects

# chronologies
Y  = as.matrix(data.frame(BTP.sub$chron, CBS.sub$chron, TOW.sub$chron)) 

Y = t(Y)

N_years = ncol(Y)
N_sites = nrow(Y)


# continuous memory var
mem_var = 'tmin.may'
lag = 6

d = t(as.matrix(data.frame(BTP.sub[, mem_var], CBS.sub[, mem_var], TOW.sub[, mem_var])))

# binary memory var
fire.raw = as.matrix(fire.raw$fire)
fire = as.vector(fire.raw)

# covars
covars = c('pcp.aug', 'pdsi.sep')
N_covars = length(covars)

X = array(NA, c(N_sites, N_years, N_covars))
for (i in 1:N_covars){
  X[,,i] = t(as.matrix(data.frame(BTP.sub[,covars[i]], CBS.sub[covars[i]], TOW.sub[covars[i]])))
}

idx.short.na  = which(apply(X, 2, function(x) any(is.na(x))))
# idx.short.na.d  = which(apply(d, 1, function(x) any(is.na(x))))

for (i in 1:N_covars){
  X[,idx.short.na,i] =  matrix(rowMeans(X[,,i], na.rm=TRUE))[,rep(1, length(idx.short.na))]
}

d[idx.short.na,]  = matrix(rowMeans(d, na.rm=TRUE))[,rep(1, length(idx.short.na))]


X_nmiss = length(idx.short.na)
d_nmiss = length(idx.short.na)

X_index = idx.short.na
d_index = idx.short.na

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
           n_knots = n_knots,
           X_nmiss = X_nmiss,
           d_nmiss = d_nmiss,
           X_index = X_index,
           d_index = d_index)

saveRDS(dat, 'scripts-stan/output/data_ecomem_basis_imp_test.RDS')

# 
# #######################################################################################
# ## inits
# #######################################################################################
# sigma_eta = 0.1;
# 
# library(mvtnorm)
# 
# eta = rnorm(6, 0, 1)
# 
# eta * rmvnorm(1, rep(0, n_basis), sigma_eta^2 * S_inv)
# 
# B_eta = eta * rmvnorm(1, rep(0, n_basis), sigma_eta^2 * S_inv)
# eta * rmvnorm(1, rep(0, n_basis), sigma_eta^2 * S_inv);
# 
# exp(B_eta)/sum(exp(B_eta))
# 
# tau_fire = 0.5
# sigma = abs(rcauchy(3, 0.01, 0.2))
# sigma_alpha = 0.2
# u = matrix(0, nrow=N_sites, ncol=N_years)
# beta = rep(0, N_covars)
# gamma0 = 0
# gamma1 = 0 
# gamma2 = 0
# 
# Sigma = diag(N_covars + 2)
# mu = matrix(0, nrow=N_sites, ncol=N_covars+2)
# 
# X1_imp = matrix(rep(rowMeans(X1), times = X_nmiss), nrow=N_sites, ncol=X_nmiss)
# X2_imp = matrix(rep(rowMeans(X2), times = X_nmiss), nrow=N_sites, ncol=X_nmiss)
# d_imp  = matrix(rep(rowMeans(d), times = X_nmiss), nrow=N_sites, ncol=X_nmiss) 
# 
# inits = list(list(tau_fire = tau_fire,
#                  sigma = sigma,
#                  sigma_alpha = sigma_alpha,
#                  u = u,
#                  beta = beta,
#                  gamma0 = gamma0,
#                  gamma1 = gamma1,
#                  gamma2 = gamma2,
#                  eta = eta,
#                  Sigma = Sigma,
#                  mu = mu,
#                  X1_imp = X1_imp,
#                  X2_imp = X2_imp,
#                  d_imp = d_imp))
    
#######################################################################################
## sampling
#######################################################################################

N_iter = 500

#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

sm<-stan_model("scripts-stan/ecomem_basis_imp_test.stan")

fit<-sampling(sm,
              data=dat,
              iter=N_iter,
              chains = 1, 
              cores = 1)#,
              #init = inits)#,control = list(adapt_delta=0.95))

saveRDS(fit, 'scripts-stan/output/fit_ecomem_basis_imp_test.RDS')

