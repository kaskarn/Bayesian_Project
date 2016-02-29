library(dplyr) #Used for data management
library(magrittr) #Pipe operators in R
library(MCMCpack) #inverse Wishart and inverse Gamma
library(reshape) #to graph things!
library(ggplot2) #to graph things nicely!

gdat <- cbind(melt(BETA))
n_vec <- character(H)
for(i in 1:H) n_vec[i] <- as.character(regions[[hospitals[[i]]$up]]$name)
gdat$reg <- rep(n_vec, rep(I, length(n_vec)))
for(i in 1:H) n_vec[i] <- as.character(regions[[hospitals[[i]]$up]]$p_death)
gdat$r_p <- rep(n_vec, rep(I, length(n_vec)))
for(i in 1:H) n_vec[i] <- hospitals[[i]]$p_death
gdat$h_p <- as.numeric(rep(n_vec, rep(I, length(n_vec))))
gdat$X2 <- factor(gdat$X2, labels = levels(docs$Hospital.Name))

g <-  ggplot(gdat[gdat$reg == "Bronx",], aes(x = X1, y=value)) + geom_line(aes(colour=X2)) +
  theme(legend.position="none") + ylab("Value") + xlab("Iteration") + facet_grid(X2 ~.)
g

reg.now <- "Manhattan"
g1 <-  ggplot(gdat %>% filter(reg == reg.now), aes(x = expit(value), fill=X2)) + geom_density(aes(colour=X2)) +
  theme(legend.position="none") + ylab("Value") + xlab("Iteration") +
  facet_grid(X2 ~.) + scale_x_continuous(limits=c(0,0.05),breaks=seq(0,0.05,0.01))+ 
  geom_vline(aes(xintercept = 0.016867)) + geom_vline(aes(xintercept = h_p), colour="red")
g1


mc.beta <- mcmc(dat[['gamma']])
sum.mcmc <- summary(mc.beta)
expit(sum.mcmc$quantiles)
autocorr.plot(mc.beta)

rejectionRate(mc.beta)
densplot(mc.beta)
