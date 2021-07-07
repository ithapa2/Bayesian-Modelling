library(ggplot2)
library(rstan)
library(reshape2)

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
# fit = readRDS('scripts-stan/output/fit_ecomem_basis_imp.RDS')
# 
# suff = 'imp'

fit = readRDS(paste0(path_output, '/fit_ecomem_basis_imp_', suffix, '.RDS'))

plot(fit)

names(fit)
post = extract(fit)

names(post)


#######################################################################################
## load data
#######################################################################################

dat = readRDS(paste0(path_output, '/data_ecomem_basis_imp_', suffix, '.RDS'))

n_basis = dat$n_basis
N_sites = dat$N_sites
lag = dat$lag

#######################################################################################
## plot output
#######################################################################################

if (!dir.exists(path_figures)){
  dir.create(path_figures)
}

## discrete memory

if ((include_outbreak==1)&(include_fire==0)){
  tau_df = post$tau_dmem 
  phi_df = 1/tau_df
  phi_df_mean = mean(phi_df)
  
  tau = data.frame(iter = seq(1, length(tau_df)), tau_df)
  tau_melt = melt(tau, id.vars='iter')
  
  phi = data.frame(iter = seq(1, length(phi_df)), phi_df)
  phi_melt = melt(phi, id.vars='iter')
  
  ggplot(data=tau_melt) + 
    geom_histogram(aes(x=value, y = (..count..)/sum(..count..))) + 
    facet_grid(variable~.) +
    xlab('phi') +
    ylab('frequency') +
    theme_bw()
  ggsave(paste0(path_figures, '/dmem_tau_hist_', suffix, '.png'))
  
  ggplot(data=tau_melt) + 
    geom_line(aes(x = iter, y = value)) + 
    facet_grid(variable~., scales = "free_y") +
    xlab('iter') +
    ylab('value') +
    theme_bw()
  ggsave(paste0(path_figures, '/dmem_tau_trace_', suffix, '.png'))
  
  
  x = seq(1e-6, lag, by=0.1)
  # mem_d = unlist(lapply(x, function(x){spherical(x, phi_d_mean)}))
  # 
  # plot(x, mem_dmem, xlab="Year", ylab="Antecedent weight")
  
  dfmem_iter = matrix(NA, nrow=length(phi_df), ncol=length(x))
  for (i in 1:length(phi_df)){
    print(i)
    dfmem_iter[i,] = unlist(lapply(x, function(x){spherical(x, phi_df[i])}))
  }
  
  dfmem_quants = apply(dfmem_iter, 2, function(x) quantile(x, c(0.05, 0.5, 0.95)))
  dfmem_dat = data.frame(type = rep('insect', length(x)), x=x, t(dfmem_quants))
  colnames(dfmem_dat) = c('type', 'x', 'q5', 'q50', 'q95')
  
  dmem_dat =dfmem_dat
  
  # ggplot(data=dimem_dat) +
  #   geom_point(aes(x=x, y=q50)) +
  #   # geom_linerange(aes(y=x, xmin=q5, xmax=q95)) +
  #   geom_linerange(aes(x=x, ymin=q5, ymax=q95)) +
  #   theme_bw() + 
  #   xlab#+
  #   #scale_y_discrete(limits=rev)
  
  # ggplot(data=dimem_dat) +
  #   geom_point(aes(x=x, y=q50)) +
  #   # geom_linerange(aes(y=x, xmin=q5, xmax=q95)) +
  #   geom_linerange(aes(x=x, ymin=q5, ymax=q95)) +
  #   theme_bw() +
  #   xlab#+
  #   #scale_y_discrete(limits=rev)
  
  ggplot(data=dmem_dat) + 
    geom_line(aes(x=x, y=q50, color=type), lwd=1) +
    geom_ribbon(aes(x=x, ymin=q5, ymax=q95, fill=type), alpha=0.4) +
    theme_bw() +
    theme(text = element_text(size=14)) + 
    # xlab("Lag") +
    ylab("Discrete Antecedent Weight") +
    scale_x_continuous(name="Lag", breaks=seq(0, 6))
  ggsave(paste0(path_figures, '/dmem_antecedent-weight_', suffix, '.png'))
  
}
if ((include_outbreak==0)&(include_fire==1)){
  
  tau_df = post$tau_dmem 
  phi_df = 1/tau_df
  phi_df_mean = mean(phi_df)
  
  tau = data.frame(iter = seq(1, length(tau_df)), tau_df)
  tau_melt = melt(tau, id.vars='iter')
  
  phi = data.frame(iter = seq(1, length(phi_df)), phi_df)
  phi_melt = melt(phi, id.vars='iter')
  
  ggplot(data=tau_melt) + 
    geom_histogram(aes(x=value, y = (..count..)/sum(..count..))) + 
    facet_grid(variable~.) +
    xlab('phi') +
    ylab('frequency') +
    theme_bw()
  ggsave(paste0(path_figures, '/dmem_tau_hist_', suffix, '.png'))
  
  ggplot(data=tau_melt) + 
    geom_line(aes(x = iter, y = value)) + 
    facet_grid(variable~., scales = "free_y") +
    xlab('iter') +
    ylab('value') +
    theme_bw()
  ggsave(paste0(path_figures, '/dmem_tau_trace_', suffix, '.png'))
  
  
  x = seq(1e-6, lag, by=0.1)
  # mem_d = unlist(lapply(x, function(x){spherical(x, phi_d_mean)}))
  # 
  # plot(x, mem_dmem, xlab="Year", ylab="Antecedent weight")
  
  dfmem_iter = matrix(NA, nrow=length(phi_df), ncol=length(x))
  for (i in 1:length(phi_df)){
    print(i)
    dfmem_iter[i,] = unlist(lapply(x, function(x){spherical(x, phi_df[i])}))
  }
  
  dfmem_quants = apply(dfmem_iter, 2, function(x) quantile(x, c(0.05, 0.5, 0.95)))
  dfmem_dat = data.frame(type = rep('fire', length(x)), x=x, t(dfmem_quants))
  colnames(dfmem_dat) = c('type', 'x', 'q5', 'q50', 'q95')
  
  dmem_dat =dfmem_dat
  
  # ggplot(data=dimem_dat) +
  #   geom_point(aes(x=x, y=q50)) +
  #   # geom_linerange(aes(y=x, xmin=q5, xmax=q95)) +
  #   geom_linerange(aes(x=x, ymin=q5, ymax=q95)) +
  #   theme_bw() + 
  #   xlab#+
  #   #scale_y_discrete(limits=rev)
  
  # ggplot(data=dimem_dat) +
  #   geom_point(aes(x=x, y=q50)) +
  #   # geom_linerange(aes(y=x, xmin=q5, xmax=q95)) +
  #   geom_linerange(aes(x=x, ymin=q5, ymax=q95)) +
  #   theme_bw() +
  #   xlab#+
  #   #scale_y_discrete(limits=rev)
  
  ggplot(data=dmem_dat) + 
    geom_line(aes(x=x, y=q50, color=type), lwd=1) +
    geom_ribbon(aes(x=x, ymin=q5, ymax=q95, fill=type), alpha=0.4) +
    theme_bw() +
    theme(text = element_text(size=14)) + 
    # xlab("Lag") +
    ylab("Discrete Antecedent Weight") +
    scale_x_continuous(name="Lag", breaks=seq(0, 6))
  ggsave(paste0(path_figures, '/dmem_antecedent-weight_', suffix, '.png'))

}

if ((include_outbreak)&(include_fire)){
 
  tau_di = post$tau_dimem 
  phi_di = 1/tau_di
  phi_di_mean = mean(phi_di)
  
  tau_df = post$tau_dfmem
  phi_df = 1/tau_df
  phi_df_mean = mean(phi_df)
  
  tau = data.frame(iter = seq(1, length(tau_df)), tau_di, tau_df)
  tau_melt = melt(tau, id.vars='iter')
  
  phi = data.frame(iter = seq(1, length(phi_df)), phi_di, phi_df)
  phi_melt = melt(phi, id.vars='iter')
  
  ggplot(data=tau_melt) + 
    geom_histogram(aes(x=value, y = (..count..)/sum(..count..))) + 
    facet_grid(variable~.) +
    xlab('phi') +
    ylab('frequency') +
    theme_bw()
  ggsave(paste0(path_figures, '/dmem_tau_hist_', suffix, '.png'))
  
  ggplot(data=tau_melt) + 
    geom_line(aes(x = iter, y = value)) + 
    facet_grid(variable~., scales = "free_y") +
    xlab('iter') +
    ylab('value') +
    theme_bw()
  ggsave(paste0(path_figures, '/dmem_tau_trace_', suffix, '.png'))
  
  
  x = seq(1e-6, lag, by=0.1)
  # mem_d = unlist(lapply(x, function(x){spherical(x, phi_d_mean)}))
  # 
  # plot(x, mem_dmem, xlab="Year", ylab="Antecedent weight")
  
  dimem_iter = matrix(NA, nrow=length(phi_di), ncol=length(x))
  for (i in 1:length(phi_di)){
    print(i)
    dimem_iter[i,] = unlist(lapply(x, function(x){spherical(x, phi_di[i])}))
  }
  
  dimem_quants = apply(dimem_iter, 2, function(x) quantile(x, c(0.05, 0.5, 0.95)))
  dimem_dat = data.frame(type = rep('insect', length(x)), x=x, t(dimem_quants))
  colnames(dimem_dat) = c('type', 'x', 'q5', 'q50', 'q95')
  
  dfmem_iter = matrix(NA, nrow=length(phi_df), ncol=length(x))
  for (i in 1:length(phi_df)){
    print(i)
    dfmem_iter[i,] = unlist(lapply(x, function(x){spherical(x, phi_df[i])}))
  }
  
  dfmem_quants = apply(dfmem_iter, 2, function(x) quantile(x, c(0.05, 0.5, 0.95)))
  dfmem_dat = data.frame(type = rep('fire', length(x)), x=x, t(dfmem_quants))
  colnames(dfmem_dat) = c('type', 'x', 'q5', 'q50', 'q95')
  
  dmem_dat = rbind(dimem_dat, dfmem_dat)
  
  # ggplot(data=dimem_dat) +
  #   geom_point(aes(x=x, y=q50)) +
  #   # geom_linerange(aes(y=x, xmin=q5, xmax=q95)) +
  #   geom_linerange(aes(x=x, ymin=q5, ymax=q95)) +
  #   theme_bw() + 
  #   xlab#+
  #   #scale_y_discrete(limits=rev)
  
  # ggplot(data=dimem_dat) +
  #   geom_point(aes(x=x, y=q50)) +
  #   # geom_linerange(aes(y=x, xmin=q5, xmax=q95)) +
  #   geom_linerange(aes(x=x, ymin=q5, ymax=q95)) +
  #   theme_bw() +
  #   xlab#+
  #   #scale_y_discrete(limits=rev)
  
  ggplot(data=dmem_dat) + 
    geom_line(aes(x=x, y=q50, color=type), lwd=1) +
    geom_ribbon(aes(x=x, ymin=q5, ymax=q95, fill=type), alpha=0.4) +
    theme_bw() +
    theme(text = element_text(size=14)) + 
    # xlab("Lag") +
    ylab("Discrete Antecedent Weight") +
    scale_x_continuous(name="Lag", breaks=seq(0, 6))
  ggsave(paste0(path_figures, '/dmem_antecedent-weight_', suffix, '.png'))
}

## continuous memory

w = post$w
w_quants = apply(w, 2, quantile, c(0.05, 0.5, 0.95))
w_dat = data.frame(par=paste0('w', seq(0,n_basis)), t(w_quants))
colnames(w_dat) = c('par', 'q5', 'q50', 'q95')

ggplot(data=w_dat) +
  geom_point(aes(x=par, y=q50)) +
  geom_linerange(aes(x=par, ymin=q5, ymax=q95)) +
  theme_bw() #+
  #scale_y_discrete(limits=rev)

vals = seq(0, lag)
w_dat$vals = vals

ggplot(data=w_dat) + 
  geom_line(aes(x=vals, y=q50),) +
  geom_ribbon(aes(x=vals, ymin=q5, ymax=q95), fill='dodgerblue', alpha=0.5) +
  theme_bw() +
  theme(text = element_text(size=16)) +
  ylab("Continuous \n Antecedent Weight") +
  scale_x_continuous(name="Lag", breaks=seq(0, 6))
ggsave(paste0(path_figures, '/cmem_antecedent-weight-_', suffix, '.png'))

## gamma
gamma_idx = which(substr(names(post), 1, 5) == "gamma")

gamma_dim = rep(NA, length(gamma_idx))
gamma_names = c()
for (i in 1:length(gamma_idx)){
  gamma_dim[i] = ncol(data.frame(post[gamma_idx[i]][[1]] ))
  gamma_name = names(post[gamma_idx[i]])
  if (gamma_dim[i] == 1) {
    gamma_name_rep = gamma_name
  } else {
    gamma_name_rep = paste0(gamma_name, '_', seq(1, gamma_dim[i]))
    #gamma_name_rep = rep(gamma_name, gamma_dim[i])
    
  }
  gamma_names = c(gamma_names, gamma_name_rep)
}

gamma = matrix(unlist(post[gamma_idx]), ncol = sum(gamma_dim), byrow = FALSE)
colnames(gamma) = gamma_names

gamma_melt = melt(gamma)
gamma_melt$type = substr(gamma_melt$Var2, 1, 6)
colnames(gamma_melt) = c('iter', 'variable', 'value', 'type')

ggplot(data=gamma_melt) + 
  geom_histogram(aes(x=value, y = (..count..)/sum(..count..))) + 
  geom_vline(xintercept=0, lty=2) +
  facet_wrap(~variable, scales='free_x') +
  xlab('gamma') +
  ylab('frequency') +
  theme_bw() #+ 
  #coord_flip()
ggsave(paste0(path_figures, '/gamma_hist_', suffix, '.png'))

ggplot(data=gamma_melt) + 
  geom_line(aes(x = iter, y = value)) + 
  facet_grid(variable~., scales = "free_y") +
  xlab('iter') +
  ylab('value') +
  theme_bw()
ggsave(paste0(path_figures, '/gamma_trace_', suffix, '.png'))


# #quantile(gamma1, c(0.10, 0.5, 0.90))
# gamma.quants = apply(gamma, 2, quantile, c(0.05, 0.5, 0.95))
# gamma.quants
# # coef.dat = data.frame(par=paste0('gamma', seq(0,ncol(gamma)-1)), t(gamma.quants))
# coef.dat = data.frame(par=colnames(gamma.quants), t(gamma.quants))

## beta

beta_idx = which(substr(names(post), 1, 4) == "beta")
beta_dim = ncol(data.frame(post[beta_idx][[1]]))

beta_names = paste0('beta', seq(0, beta_dim-1))

beta = matrix(unlist(post[beta_idx]), ncol = beta_dim, byrow = FALSE)
colnames(beta) = beta_names

beta_melt = melt(beta)
#beta_melt$type = substr(beta_melt$Var2, 1, 6)
colnames(beta_melt) = c('iter', 'variable', 'value')#, 'type')

ggplot(data=beta_melt) + 
  geom_histogram(aes(x=value, y = (..count..)/sum(..count..))) + 
  geom_vline(xintercept=0, lty=2) +
  facet_wrap(~variable, scales='free_x') +
  xlab('beta') +
  ylab('frequency') +
  theme_bw() #+ 
#coord_flip()
ggsave(paste0(path_figures, '/beta_hist_', suffix, '.png'))

ggplot(data=beta_melt) + 
  geom_line(aes(x = iter, y = value)) + 
  facet_grid(variable~., scales = "free_y") +
  xlab('iter') +
  ylab('value') +
  theme_bw()
ggsave(paste0(path_figures, '/beta_trace_', suffix, '.png'))


## sigma

sigma_idx = which(substr(names(post), 1, 5) == "sigma")

sigma_dim = rep(NA, length(sigma_idx))
sigma_names = c()
for (i in 1:length(sigma_idx)){
  sigma_dim[i] = ncol(data.frame(post[sigma_idx[i]][[1]] ))
  sigma_name = names(post[sigma_idx[i]])
  if (sigma_dim[i] == 1) {
    sigma_name_rep = sigma_name
  } else {
    sigma_name_rep = paste0(sigma_name, '_', seq(1, sigma_dim[i]))
    #gamma_name_rep = rep(gamma_name, gamma_dim[i])
    
  }
  sigma_names = c(sigma_names, sigma_name_rep)
}

sigma = matrix(unlist(post[sigma_idx]), ncol = sum(sigma_dim), byrow = FALSE)
colnames(sigma) = sigma_names

sigma_melt = melt(sigma)
#beta_melt$type = substr(beta_melt$Var2, 1, 6)
colnames(sigma_melt) = c('iter', 'variable', 'value')#, 'type')

ggplot(data=sigma_melt) + 
  geom_histogram(aes(x=value, y = (..count..)/sum(..count..))) + 
  #geom_vline(xintercept=0, lty=2) +
  facet_grid(variable~.) + #, scales='free_x') +
  xlab('sigma') +
  ylab('frequency') +
  theme_bw() #+ 
#coord_flip()
ggsave(paste0(path_figures, '/sigma_hist_', suffix, '.png'))

ggplot(data=sigma_melt) + 
  geom_line(aes(x = iter, y = value)) + 
  facet_grid(variable~., scales = "free_y") +
  xlab('iter') +
  ylab('value') +
  theme_bw()
ggsave(paste0(path_figures, '/sigma_trace_', suffix, '.png'))
