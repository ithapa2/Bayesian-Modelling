library(reshape2)

#######################################################################################
## read in and organize data data
#######################################################################################

BTP = read.csv("data/BTP-BI.csv", header=TRUE)
CBS = read.csv("data/CBS-BI.csv", header=TRUE)
TOW = read.csv("data/TOW-BI.csv", header=TRUE)

BTP.sub = BTP[apply(BTP[,c('chron', 'tmin.may', 'pcp.aug', 'pdsi.sep')], 1, function(x) all(!is.na(x))),]

CBS.sub = CBS[apply(CBS[,c('chron', 'tmin.may', 'pcp.aug', 'pdsi.sep')], 1, function(x) all(!is.na(x))),]

TOW.sub = TOW[apply(TOW[,c('chron', 'tmin.may', 'pcp.aug', 'pdsi.sep')], 1, function(x) all(!is.na(x))),]

fire.raw = data.frame(year = seq(1650, 2017, by=1), fire = rep(0))
fire.raw[which(fire.raw$year %in% c(1664, 1804, 1900)), 'fire'] = 1
fire.raw = fire.raw[which(fire.raw$year %in% BTP.sub$year),]


BTP.fire = merge(BTP.sub, fire.raw, by="year")
CBS.fire = merge(CBS.sub, fire.raw, by="year")
TOW.fire = merge(TOW.sub, fire.raw, by="year")


dat = rbind(data.frame(BTP.fire, site=rep("BTP", nrow(BTP.fire))),
            data.frame(CBS.fire, site=rep("CBS", nrow(BTP.fire))),
            data.frame(TOW.fire, site=rep("TOW", nrow(BTP.fire))))

write.csv(dat, 'data/BI-site-dat.csv', row.names=FALSE)


