library(EcoMem)
library(reshape2)

#######################################################################################
## read in data
#######################################################################################

dat = read.csv("data/BI-site-dat.csv", header=TRUE)


#######################################################################################
## fit model 1: BI; memory: tmin.may
#######################################################################################

mod = ecomem(chron~tmin.may, data=dat, timeID='year', groupID='site', mem.vars='tmin.may', L=6)

# Convert ecomem output to mcmc object
post.samps = mem2mcmc(mod)

# Generate posterior sample traceplots
plot(post.samps[,1:3])

# Summarize marginal posterior distributions
coef.summ = memsum(mod, cred.int = 0.99, verbose = FALSE)

# Define parameter types (Beta, w1, w2)
# coef.summ$var.type = c(rep("Wts. v1", length(w1)))

# Plot posterior summary of model parameters relative to simulated values
ggplot(aes(x = var, y = mean, ymin = q0.005, ymax = q0.995),
       data = coef.summ) +
  geom_linerange() +
  geom_point(shape = 4, size = 2) +
  # geom_point(aes(y = sim.val, color = In),
  #            size = 2, alpha = 0.6) +
  #scale_color_manual(values = c("red","green")) +
  xlab("Parameter") + ylab("Posterior Summary") +
  # facet_wrap(~ var.type, ncol = 1, scale = "free") +
  coord_flip() +
  theme_bw() +
  labs(color = "Simulated value in\ncredible interval")

# Plot memory functions
p = plotmem(mod, cred.int = 0.99)
p
