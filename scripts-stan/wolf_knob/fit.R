
#######################################################################################
## read in the data
#######################################################################################

# possible sites
sites = c('PSME')
fnames = c("data/WOLF_PSME_FINAL-std-all.csv")

suffix = 'stream-insect'

path_output = 'scripts-stan/wolf_knob/output'
path_figures = 'scripts-stan/wolf_knob/figures'

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

#######################################################################################
## prepare data for stan model
#######################################################################################

fire.raw = data.frame(year = seq(736, 2017, by=1), fire = rep(0))
fire.raw[which(fire.raw$year %in% c(388, 402, 432, 663, 690, 957, 964, 1023, 1056, 1110, 1144, 1148, 1180, 1190, 1213, 1217, 1277, 1278, 
                                    1296, 1299, 1306, 1317, 1325, 1355, 1361, 1388, 1407, 1437, 1439, 1446, 1468, 
                                    1489, 1492, 1499, 1528, 1531, 1532, 1538, 1543, 1551, 1556, 1561, 1571, 1576, 
                                    1580, 1584, 1590, 1593, 1598, 1600, 1626, 1631, 1671, 1703, 1705, 1742, 1765, 
                                    1800, 1836, 1842, 1894, 1949, 1956, 1966, 1997)), 'fire'] = 1

include_outbreak = 0

insect = read.csv('data/JPTD-budworm-outbreak.csv', header=TRUE)
insect$bool = insect$percent
insect$bool[which(insect$percent>35)] = 1
insect$bool[which(insect$percent<=35)] = 0

#model_name = 'ecomem_basis_imp.stan'
model_name = 'ecomem_basis_imp_insect.stan'

source("scripts-stan/common/fit_ecomem_basis_imp.R")
