library(dplyr) #Used for data management
library(magrittr) #Pipe operators in R
library(MCMCpack) #inverse Wishart and inverse Gamma
library(reshape) #to graph things!
library(ggplot2) #to graph things nicely!

dat <- readRDS("results_mcmc.rds")
thin <- 1
burnin <- 1000
L <- I - burnin
for(i in 1:length(dat)) dat[[i]] <- dat[[i]][(!(1:nrow(dat[[i]]))%%thin),]
dat[['gamma']] <- dat[['gamma']][,-178]
# results <- list(BETA, THETA, SIGMA, GAMMA)

gdat.r <- melt(dat[[2]][-c(1:burnin),])
gdat.r$reg <- rep(unique(docs$Detailed.Region), rep(L, length(unique(docs$Detailed.Region))))

gdat.h <- melt(dat[[1]][-c(1:burnin),])
tm <- NULL
for(i in unique(as.numeric(docs$Hospital.Name))) tm <- cbind(tm, hos_list[[i]]$regs)
gdat.h$reg <- rep(tm, rep(L/thin, length(tm)))
gdat.h$reg <- factor(gdat.h$reg, labels = levels(docs$Detailed.Region))

tm <- NULL
for(i in unique(as.numeric(docs$Hospital.Name))) tm <- cbind(tm, hos_list[[i]]$e_mu)
gdat.h$e_mu <- rep(tm, rep(L/thin, length(tm)))

tm <- NULL
for(i in unique(as.numeric(docs$Hospital.Name))) tm <- cbind(tm, hos_list[[i]]$mu)
gdat.h$mu <- rep(tm, rep(L/thin, length(tm)))

tm <- NULL
for(i in unique(as.numeric(docs$Hospital.Name))) tm <- cbind(tm, hos_list[[i]]$mu.r)
gdat.h$mu.r <- rep(tm, rep(L/thin, length(tm)))

gdat.d <- melt(dat[["gamma"]][-c(1:burnin),])

reg.now <- "Manhattan"


g1 <-  ggplot(gdat.h %>% filter(reg == reg.now), aes(x = expit(value), fill=X2)) + geom_density(aes(colour=X2)) +
      theme(legend.position="none") + ylab("Value") + xlab("Iteration") +
      facet_grid(X2 ~.) + scale_x_continuous(limits=c(0,0.05),breaks=seq(0,0.05,0.01)) +
      geom_vline(aes(xintercept = mu.r)) + geom_vline(aes(xintercept = mu)) + geom_vline(aes(xintercept = mu.r)) 
g1

g2 <-  ggplot(gdat.r, aes(x = expit(value), fill=X2)) + geom_density(aes(colour=X2)) +
      theme(legend.position="none") + ylab("Value") + xlab("Iteration") +
      facet_grid(X2 ~.) + scale_x_continuous(limits=c(0,0.05),breaks=seq(0,0.05,0.01)) 
g2
expit <- function(theta){
  return(exp(theta)/(1+exp(theta)))
}

mc.beta <- mcmc(dat[['gamma']])
sum.mcmc <- summary(mc.beta)
expit(sum.mcmc$quantiles)
autocorr.plot(mc.beta)

rejectionRate(mc.beta)
densplot(mc.beta)
