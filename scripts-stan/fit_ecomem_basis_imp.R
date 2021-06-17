library("mgcv")
library("rstan")

#######################################################################################
## read in the data
#######################################################################################

# possible sites
# sites = c('BTP', 'CBS', 'TOW')
# fnames = c("data/BTP-BI.csv", "data/CBS-BI.csv", "data/TOW-BI.csv")

# indicate sites of interest and corresponding file names
sites = c('PSME')
fnames = c("data/WOLF_PSME_FINAL-std.csv")

#######################################################################################
## specify variables and lag
#######################################################################################

# continuous memory var
# currently must have exactly one
mem_var = 'ppt.pjja'
lag = 6

# covariates
# currently must have one or more
covars = c('tmax.paug', 'pdsi.poct')

#######################################################################################
## prepare data for stan model
#######################################################################################

N_sites = length(sites)

raw = list()
for (site in 1:N_sites){
  raw[[site]] = read.csv(fnames[site], header=TRUE)
}
names(raw) = sites

fire.raw = data.frame(year = seq(736, 2017, by=1), fire = rep(0))
fire.raw[which(fire.raw$year %in% c(957, 964, 1023, 1056, 1110, 1144, 1148, 1180, 1190, 1213, 1217, 1277, 1278, 1296, 1299, 1306, 1317, 1325, 1355, 1361, 1388, 1407, 1437, 1439, 1446, 1468, 1489, 1492, 1499, 1528, 1531, 1532, 1538, 1543, 1551, 1556, 1561, 1571, 1576, 1580, 1584, 1590, 1593, 1598, 1600, 1626, 1631, 1671, 1703, 1705, 1742, 1765, 1800, 1836, 1842, 1894, 1949, 1956, 1966, 1997)), 'fire'] = 1

## find the start year and end year
## for now require that chronologies and fire record all be the same length with no NA values
## climate data of shorter length NA values will be imputed
year_upper = NA
year_lower = NA
for (site in 1:N_sites){
  year_upper = min(year_upper, max(raw[[site]]$year, na.rm=TRUE), na.rm=TRUE)
  year_lower = max(year_lower, min(raw[[site]]$year, na.rm=TRUE), na.rm=TRUE)
}

year_upper = min(year_upper, max(fire.raw$year))
year_lower = max(year_lower, min(fire.raw$year))

years = seq(year_lower, year_upper)
N_years = length(years)

# subset site and fire data to selected years
for (site in 1:N_sites){ 
  raw[[site]] = raw[[site]][which(raw[[site]]$year %in% years),]
}

fire.raw = fire.raw[which(fire.raw$year %in% years),]

## define data objects

# chronologies
Y = t(matrix(unlist(lapply(raw, function(x) x$chron)), ncol=N_sites, byrow=FALSE))

# continuous memory var
d = t(matrix(unlist(lapply(raw, function(x) x[,mem_var])), ncol=N_sites, byrow=FALSE))

# binary memory var
fire.raw = as.matrix(fire.raw$fire)
fire = as.vector(fire.raw)

# covars

N_covars = length(covars)

X = array(NA, c(N_sites, N_years, N_covars))
for (i in 1:N_covars){
  X[,,i] = t(matrix(unlist(lapply(raw, function(x) x[,covars[i]])), ncol=N_sites, byrow=FALSE))
}

idx.short.na  = which(apply(X, 2, function(x) any(is.na(x))))
# idx.short.na.d  = which(apply(d, 1, function(x) any(is.na(x))))

for (i in 1:N_covars){
  if (N_sites != 1){
    X[,idx.short.na,i] =  matrix(rowMeans(X[,,i], na.rm=TRUE))[,rep(1, length(idx.short.na))]
    d[,idx.short.na]  = matrix(rowMeans(d, na.rm=TRUE))[,rep(1, length(idx.short.na))]
  } else if (N_sites == 1){
    X[,idx.short.na,i] =  rep(mean(X[,,i], na.rm=TRUE), length(idx.short.na))
    d[,idx.short.na]  = rep(mean(d, na.rm=TRUE), length(idx.short.na))
  }
}



X_nmiss = length(idx.short.na)
d_nmiss = length(idx.short.na)

X_index = idx.short.na
d_index = idx.short.na

#######################################################################################
## splines
#######################################################################################

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
## compile data as a list; save data as RDS object
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

saveRDS(dat, 'scripts-stan/output/data_ecomem_basis_imp.RDS')

# 
# #######################################################################################
# ## specify model initial conditions
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
## compile model and perform sampling
#######################################################################################

N_iter = 500

#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

# compile the stand model
# must be done every time a change is made to the .stan file
sm<-stan_model("scripts-stan/ecomem_basis_imp.stan")

# parameter estimation
fit<-sampling(sm,
              data=dat,
              iter=N_iter,
              chains = 1, 
              cores = 1)#,
#init = inits)#,control = list(adapt_delta=0.95))

# save stan fit object for subsequent analysis
saveRDS(fit, 'scripts-stan/output/fit_ecomem_basis_imp.RDS')

