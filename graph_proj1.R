library(dplyr) #Used for data management
library(magrittr) #Pipe operators in R
library(MCMCpack) #inverse Wishart and inverse Gamma
library(reshape) #to graph things!
library(ggplot2) #to graph things nicely!

dat <- readRDS("results_mcmc.rds")
thin <- 1
for(i in 1:length(dat)) dat[[i]] <- dat[[i]][(!(1:nrow(dat[[i]]))%%thin),]
dat[['gamma']] <- dat[['gamma']][,-178]
# results <- list(BETA, THETA, SIGMA, GAMMA)

gdat <- melt(dat[[1]])
tm <- NULL
for(i in unique(as.numeric(docs$Hospital.Name))) tm <- cbind(tm, hos_list[[i]]$regs)
gdat$reg <- rep(tm, rep(I/thin, length(tm)))
gdat$reg <- factor(gdat$reg, labels = levels(docs$Detailed.Region))

tm <- NULL
for(i in unique(as.numeric(docs$Physician.Name))) tm <- cbind(tm, doc_list[[i]]$hos)
gdat$reg <- rep(tm, rep(I/thin, length(tm)))
gdat$reg <- factor(gdat$reg, labels = levels(docs$Detailed.Region))

g <-  ggplot(gdat %>% filter(reg == "Bronx"), aes(x = X1, y=value)) + geom_line(aes(colour=X2)) +
  theme(legend.position="none") + ylab("Value") + xlab("Iteration") + facet_grid(X2 ~.)
g

g <-  ggplot(gdat %>% filter(reg == "Bronx"), aes(x = expit(value))) + geom_density(aes(colour=X2)) +
      theme(legend.position="none") + ylab("Value") + xlab("Iteration") + facet_grid(X2 ~.)
g


expit <- function(theta){
  return(exp(theta)/(1+exp(theta)))
}

mc.beta <- mcmc(dat[['gamma']])
sum.mcmc <- summary(mc.beta)
expit(sum.mcmc$quantiles)
autocorr.plot(mc.beta)

rejectionRate(mc.beta)
densplot(mc.beta)
