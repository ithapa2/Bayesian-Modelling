library(corrplot)


# possible sites
sites = c('PSME')
fnames = c("data/WOLF_PSME_FINAL-std-all.csv")

raw = list()
for (site in 1:N_sites){
 raw[[site]] = read.csv(fnames[site], header=TRUE)
}
names(raw) = sites

res = cor(foo[,2:ncol(foo)], use='complete.obs')
corrplot(res)