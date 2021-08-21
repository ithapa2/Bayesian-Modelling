
#######################################################################################
## read in the data
#######################################################################################

# possible sites
sites = c('BTP', 'CBS', 'TOW')
fnames = c("data/BTP-BI.csv", "data/CBS-BI.csv", "data/TOW-BI.csv")

# suffix = 'fire-stream-logy'
#suffix = 'insect-stream-logy'
# suffix = 'fire-stream-logy'
suffix = 'logy_clim_dat'

path_output = 'scripts-stan/three_climate/output'
path_figures = 'scripts-stan/three_climate/figures'

# indicate sites of interest and corresponding file names
#sites = c('BTP')
#fnames = c("data/BTP-BI.csv")

#######################################################################################
## specify variables and lag
#######################################################################################

# # MODEL: swe
# # continuous memory var
# # currently must have exactly one
# mem_var = 'swe.yel'
# # mem_var = 'ppt.jja'
# lag = 6
# 
# # covariates
# # currently must have one or more
# #covars = c('tmax.paug', 'pdsi.poct')
# covars = c('tmax.lsum', 'stream.yel')

# MODEL: stream
# continuous memory var
# currently must have exactly one
mem_var = 'tmin.may'
# mem_var = 'ppt.jja'
lag = 6

# covariates
# currently must have one or more
covars = c('pcp.aug', 'pdsi.sep')
#covars = c('tmax.paug', 'pdsi.poct')
#covars = c('tmax.lsum', 'swe.yel')

include_outbreak = 0
include_fire     = 0

#######################################################################################
## prepare data for stan model
#######################################################################################


fire.raw = data.frame(year = bi.years, fire = rep(0))
fire.raw[which(fire.raw$year %in% c(1664, 1804, 1900)), 'fire'] = 1
# 
# insect = read.csv('data/JPTD-budworm-outbreak.csv', header=TRUE)
# insect$bool = insect$percent
# insect$bool[which(insect$percent>60)] = 1
# insect$bool[which(insect$percent<=60)] = 0

if (!include_outbreak & !include_fire){
  model_name = 'ecomem_basis_imp_logy_0dmem.stan'
}

if (xor(include_outbreak, include_fire)){
  model_name = 'ecomem_basis_imp_logy.stan'
}

if ((include_outbreak)&(include_fire)){
  model_name = 'ecomem_basis_imp_logy_2dmem.stan'
}


#model_name = 'ecomem_basis_imp.stan'
# model_name = 'ecomem_basis_imp.stan'

source("scripts-stan/common/fit_ecomem_basis_imp_ndmem.R")
