library(ggplot2)
library(rstan)
#######################################################################################
## spherical decay function
#######################################################################################

spherical <- function(x, phi) {
  
  # if (x <= 0) {
  #   reject("x negative spherical function undefined; found x = ", x)
  # }
  if (x > phi) { g=0 } else {
    g = (1 - 3.0/2.0 * x/phi + 1.0/2.0 * (x/phi)^3) 
  }
  
  return(g)
}


#######################################################################################
## load and extract posterior samples
#######################################################################################
fit = readRDS('scripts-stan/output/fit_ecomem_fire.RDS')

plot(fit)

names(fit)
post = extract(fit)

names(post)

#######################################################################################
## plot output
#######################################################################################

## fire memory

phi_fire = 1/post$tau_fire
phi_fire.mean = mean(phi_fire)

x = seq(1e-6, lag, by=0.1)
mem.fire = unlist(lapply(x, function(x){spherical(x, phi_fire.mean)}))

plot(x, mem.fire, xlab="Year", ylab="Antecedent weight")

mem.fire.iter = matrix(NA, nrow=length(phi_fire), ncol=length(x))
for (i in 1:length(phi_fire)){
  print(i)
  mem.fire.iter[i,] = unlist(lapply(x, function(x){spherical(x, phi_fire[i])}))
}


mem.fire.quants = apply(mem.fire.iter, 2, function(x) quantile(x, c(0.05, 0.5, 0.95)))
mem.fire.dat = data.frame(x=x, t(mem.fire.quants))
colnames(mem.fire.dat) = c('x', 'q5', 'q50', 'q95')

ggplot(data=mem.fire.dat) + 
  geom_line(aes(x=x, y=q50), lwd=2) +
  geom_ribbon(aes(x=x, ymin=q5, ymax=q95), fill='dodgerblue', alpha=0.4) +
  theme_bw() +
  theme(text = element_text(size=22)) + 
 # xlab("Lag") +
  ylab("Fire Antecedent Weight") +
  scale_x_continuous(name="Lag", breaks=seq(0, 6))
ggsave("scripts-stan/figures/antecedent-weights-fire.png")

library(reshape2)
mem.fire.melt = melt(mem.fire.iter)
colnames(mem.fire.melt) = c('iter', 'lag', 'value')

## gamma parameters (memory magnitude parameters)

gamma = cbind(post$gamma0, post$gamma1, post$gamma2)
colnames(gamma) = paste0('gamma', seq(0,ncol(gamma)-1))
gamma.melt = melt(gamma)

ggplot(data=subset(gamma.melt, Var2 %in% c('gamma0', 'gamma1'))) +
  geom_density(aes(x=value), lwd=2, fill="dodgerblue", alpha=0.2) +
  geom_rug(aes(x =value, y = 0), position = position_jitter(height = 0)) +
  theme_bw() +
  theme(text = element_text(size=22)) +
  facet_wrap(Var2~., scales="free_x")
ggsave("scripts-stan/figures/gamma-small-post-density.png")


ggplot(data=subset(gamma.melt, Var2 %in% c('gamma0', 'gamma1'))) +
  geom_line(aes(x=Var1, y=value), lwd=2) +
  theme_bw() +
  theme(text = element_text(size=22)) +
  facet_grid(Var2~., scales="free_y") +
  xlab("Iteration")
ggsave("scripts-stan/figures/gamma-small-trace.png")


ggplot(data=subset(gamma.melt, Var2 %in% c('gamma2'))) +
  geom_density(aes(x=value), lwd=2, fill="dodgerblue", alpha=0.2) +
  geom_rug(aes(x =value, y = 0), position = position_jitter(height = 0)) +
  theme_bw() +
  theme(text = element_text(size=22)) +
  facet_wrap(Var2~., scales="free_x")
ggsave("scripts-stan/figures/gamma2-post-density.png")

#quantile(gamma1, c(0.10, 0.5, 0.90))
gamma.quants = apply(gamma, 2, quantile, c(0.05, 0.5, 0.95))
gamma.quants
coef.dat = data.frame(par=paste0('gamma', seq(0,ncol(gamma)-1)), t(gamma.quants))

beta = post$beta
beta.quants = apply(beta, 2, quantile, c(0.05, 0.5, 0.95))
beta.quants
coef.dat = rbind(coef.dat, data.frame(par=c('beta1', 'beta2'), t(beta.quants)))
coef.dat = data.frame(panel=c(1,1, 2, 3, 3), coef.dat)

coef.dat = coef.dat

ggplot(data=coef.dat) +
  geom_point(aes(y=par, x=X50.)) +
  geom_linerange(aes(y=par, xmin=X5., xmax=X95.)) +
  theme_bw() +
  theme(text = element_text(size=22))

ggplot(data=coef.dat[c(1,2,4),]) +
  geom_point(aes(y=par, x=X50.)) +
  geom_linerange(aes(y=par, xmin=X5., xmax=X95.)) +
  theme_bw() +
  theme(text = element_text(size=22)) +
  facet_grid(par~., scales="free")

ggplot(data=coef.dat) +
  geom_point(aes(y=par, x=X50.)) +
  geom_linerange(aes(y=par,xmin=X5., xmax=X95.)) +
  theme_bw() +
  theme(text = element_text(size=22)) +
  facet_grid(par~., scales="free")


ggplot(data=coef.dat[c(1,2),]) +
  geom_point(aes(y=par, x=X50.), size=2) +
  geom_linerange(aes(y=par,xmin=X5., xmax=X95.), lwd=1) +
  theme_bw() +
  theme(text = element_text(size=22),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank()) +
  facet_grid(par~., scales="free") + 
  ylab('') +
  geom_vline(xintercept=0, colour="red", linetype="dashed") #+
  #scale_x_continuous(name="Value", limits=c(-0.008, 0)) 
ggsave("scripts-stan/figures/gamma-credible.png")


ggplot(data=coef.dat[c(3),]) +
  geom_point(aes(y=par, x=X50.), size=2) +
  geom_linerange(aes(y=par,xmin=X5., xmax=X95.), lwd=1) +
  theme_bw() +
  theme(text = element_text(size=22),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank()) +
  facet_grid(par~., scales="free") + 
  ylab('') +
  geom_vline(xintercept=0, colour="red", linetype="dashed") #+
#scale_x_continuous(name="Value", limits=c(-0.008, 0)) 
ggsave("scripts-stan/figures/gamma-credible.png")


ggplot(data=coef.dat[c(4,5),]) +
  geom_point(aes(y=par, x=X50.), size=2) +
  geom_linerange(aes(y=par,xmin=X5., xmax=X95.), lwd=1) +
  theme_bw() +
  theme(text = element_text(size=22),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank()) +
  facet_grid(par~., scales="free") + 
  ylab('') +
  geom_vline(xintercept=0, colour="red", linetype="dashed") +
  scale_x_continuous(name="Value", limits=c(-0.008, 0)) 
ggsave("scripts-stan/figures/beta-credible.png")


quantile(post$sigma, c(0.10, 0.5, 0.90))

## Antecedent weight w parameters 

w = post$w
w.quants = apply(w, 2, quantile, c(0.05, 0.5, 0.95))
w.dat = data.frame(par=paste0('w', seq(0,n_basis)), t(w.quants))
colnames(w.dat) = c('par', 'q5', 'q50', 'q95')

ggplot(data=w.dat) +
  geom_point(aes(y=par, x=q50)) +
  geom_linerange(aes(y=par, xmin=q5, xmax=q95)) +
  theme_bw() +
  scale_y_discrete(limits=rev)

vals = seq(0, lag)
w.dat$vals = vals

ggplot(data=w.dat) + 
  geom_line(aes(x=vals, y=q50),) +
  geom_ribbon(aes(x=vals, ymin=q5, ymax=q95), fill='dodgerblue', alpha=0.5) +
  theme_bw() +
  theme(text = element_text(size=22)) +
  ylab("May Mininum Temperature \n Antecedent Weight") +
  scale_x_continuous(name="Lag", breaks=seq(0, 6))
ggsave("scripts-stan/figures/antecedent-weights-tmin-may.png")



w_mean = colMeans(post$w)
plot(seq(0, lag), w_mean)


hist(post$sigma)
hist(post$sigma_alpha)
