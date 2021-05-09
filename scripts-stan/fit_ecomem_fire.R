library("mgcv")
library("rstan")

#######################################################################################
## read in the data
#######################################################################################

BTP = read.csv("data/BTP-BI.csv", header=TRUE)
CBS = read.csv("data/CBS-BI.csv", header=TRUE)
TOW = read.csv("data/TOW-BI.csv", header=TRUE)

# BTP.sub = BTP[apply(BTP[,c('chron', 'tmin.may', 'pcp.aug', 'pdsi.sep')], 1, function(x) all(!is.na(x))),]
# CBS.sub = CBS[apply(CBS[,c('chron', 'tmin.may', 'pcp.aug', 'pdsi.sep')], 1, function(x) all(!is.na(x))),]
# TOW.sub = TOW[apply(BTP[,c('chron', 'tmin.may', 'pcp.aug', 'pdsi.sep')], 1, function(x) all(!is.na(x))),]

# BTP.sub = BTP[apply(BTP[,c('chron')], 1, function(x) all(!is.na(x))),]
# CBS.sub = CBS[apply(CBS[,c('chron')], 1, function(x) all(!is.na(x))),]
# TOW.sub = TOW[apply(BTP[,c('chron', 'tmin.may', 'pcp.aug', 'pdsi.sep')], 1, function(x) all(!is.na(x))),]

year_upper = min(max(BTP$year, na.rm=TRUE), max(CBS$year, na.rm=TRUE), max(TOW$year, na.rm=TRUE))
year_lower = max(min(BTP$year, na.rm=TRUE), min(CBS$year, na.rm=TRUE), min(TOW$year, na.rm=TRUE))

BTP.sub = BTP[which(BTP$year %in% seq(year_lower, year_upper)),]
CBS.sub = CBS[which(CBS$year %in% seq(year_lower, year_upper)),]
TOW.sub = TOW[which(TOW$year %in% seq(year_lower, year_upper)),]

fire.raw = data.frame(year = seq(1650, 2017, by=1), fire = rep(0))
fire.raw[which(fire.raw$year %in% c(1664, 1804, 1900)), 'fire'] = 1
fire.raw = fire.raw[which(fire.raw$year %in% BTP.sub$year),]

Y  = as.matrix(data.frame(BTP.sub$chron, CBS.sub$chron, TOW.sub$chron))
# X1  = as.matrix(data.frame(BTP.sub$pcp.aug, CBS.sub$pcp.aug, TOW.sub$pcp.aug))
# X2  = as.matrix(data.frame(BTP.sub$pdsi.sep, CBS.sub$pdsi.sep, TOW.sub$pdsi.sep))
# d.raw  = as.matrix(data.frame(BTP.sub$tmin.may, CBS.sub$tmin.may, TOW.sub$tmin.may))
fire.raw = as.matrix(fire.raw$fire)

test.na = cbind(Y, fire.raw)
idx.na  = apply(test.na, 1, function(x) any(is.na(x)))

Y = Y[!idx.na,]
fire = fire.raw[!idx.na,]

Y = t(Y)

N_years = ncol(Y)
N_sites = nrow(Y)

lag = 6

#######################################################################################
## data list
#######################################################################################

dat = list(N_years = N_years,
           N_sites = N_sites,
           lag = lag,
           #N_covars = N_covars,
           Y = Y,
           #X = X,
           #d = d,
           fire = fire#,
           #B = B,
           #S_inv = S_inv,
           #n_basis = n_basis,
           #n_knots = n_knots
           )

#######################################################################################
## sampling
#######################################################################################

N_iter = 1000

#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

sm<-stan_model("scripts-stan/ecomem_fire.stan")

fit<-sampling(sm,
              data=dat,
              iter=N_iter,
              chains = 1, 
              cores = 1)#,control = list(adapt_delta=0.95))

saveRDS(fit, 'scripts-stan/output/fit_ecomem_fire.RDS')

