
#######################################################################################
## read in the data
#######################################################################################

# possible sites
sites = c('BTP', 'CBS', 'TOW')
fnames = c("data/BTP-BI.csv", "data/CBS-BI.csv", "data/TOW-BI.csv")
suffix = 'TEST1'

# indicate sites of interest and corresponding file names
#sites = c('BTP')
#fnames = c("data/BTP-BI.csv")

#######################################################################################
## specify variables and lag
#######################################################################################

# continuous memory var
# currently must have exactly one
mem_var = 'tmin.may'
lag = 6

# covariates
# currently must have one or more
covars = c('pcp.aug', 'pdsi.sep')

#######################################################################################
## prepare data for stan model
#######################################################################################

fire.raw = data.frame(year = seq(1600, 2017, by=1), fire = rep(0))
fire.raw[which(fire.raw$year %in% c(1664, 1804, 1900)), 'fire'] = 1


source("scripts-stan/common/fit_ecomem_basis_imp.R")