
#######################################################################################
## read in the data
#######################################################################################

# possible sites
sites = c('PIFL')
fnames = c("data/WOLF_PIFL_FINAL-std-all.csv")

# suffix = 'fire-stream-logy'
#suffix = 'insect-stream-logy'
# suffix = 'fire-stream-logy'
suffix = 'insect-stream-logy_test'

path_output = 'scripts-stan/engelmann_spruce/output'
path_figures = 'scripts-stan/engelmann_spruce/figures'

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
mem_var = 'stream.yel'
# mem_var = 'ppt.jja'
lag = 6

# covariates
# currently must have one or more
#covars = c('tmax.paug', 'pdsi.poct')
covars = c('tmax.lsum', 'swe.yel')

include_outbreak = 1
include_fire     = 0

#######################################################################################
## prepare data for stan model
#######################################################################################

fire.raw = data.frame(year = seq(220, 2017, by=1), fire = rep(0))
fire.raw[which(fire.raw$year %in% c(388, 402, 432, 663, 690, 957, 964, 1023, 1056, 1110, 1144, 1148, 1180, 1190, 1213, 1217, 1277, 1278, 
                                    1296, 1299, 1306, 1317, 1325, 1355, 1361, 1388, 1407, 1437, 1439, 1446, 1468, 
                                    1489, 1492, 1499, 1528, 1531, 1532, 1538, 1543, 1551, 1556, 1561, 1571, 1576, 
                                    1580, 1584, 1590, 1593, 1598, 1600, 1626, 1631, 1671, 1703, 1705, 1742, 1765, 
                                    1800, 1836, 1842, 1894, 1949, 1956, 1966, 1997)), 'fire'] = 1

insect = read.csv('data/JPTD-budworm-outbreak.csv', header=TRUE)
insect$bool = insect$percent
insect$bool[which(insect$percent>60)] = 1
insect$bool[which(insect$percent<=60)] = 0

if (xor(include_outbreak, include_fire)){
  model_name = 'ecomem_basis_imp_logy.stan'
}

if ((include_outbreak)&(include_fire)){
  model_name = 'ecomem_basis_imp_logy_2dmem.stan'
}


#model_name = 'ecomem_basis_imp.stan'
# model_name = 'ecomem_basis_imp.stan'

source("scripts-stan/common/fit_ecomem_basis_imp_ndmem.R")
