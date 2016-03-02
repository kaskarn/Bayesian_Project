library(dplyr) #Used for data management
library(magrittr) #Pipe operators in R
library(MCMCpack) #inverse Wishart and inverse Gamma
library(reshape) #to graph things!
library(ggplot2) #to graph things nicely!

res <- readRDS("results_mcmc_70000.rds")
burnin <- 2000; thin <- 10
IB <- nrow(res[[1]]) - burnin
keep <- (!((1:IB)%%thin))
size <- sum(keep)


tm <- (res[[3]][-(1:burnin),])[keep,]
res.regions <- melt(tm)
for(i in 1:R){
  tm.path <- paste("region_convergence_", i, ".png", sep = "")
  png(paste(tm.path))
  print( ggplot(res.regions[as.numeric(res.regions$X2) == i,], aes(x = X1, y = value)) + 
    geom_line(colour = rainbow(R, s = 0.7, v = 0.8)[i]) +
    xlab("Iteration") + ylab("logit(P)") + theme(legend.position = "none")
  )
}
graphics.off()


tm <- (res[[1]][-(1:burnin),])[keep,]
res.hospitals <- melt(tm)
n_vec <- tnam <- character(H)
for(i in 1:H) n_vec[i] <- as.character(regions[[hospitals[[i]]$up]]$name)
res.hospitals$reg <- rep(n_vec, rep(size, H))
for(i in 1:H) n_vec[i] <- regions[[hospitals[[i]]$up]]$p_death
res.hospitals$r_p <- as.numeric(rep(n_vec, rep(size, length(n_vec))))
for(i in 1:H){
  n_vec[i] <- hospitals[[i]]$p_death
  tnam[i] <- as.character(hospitals[[i]]$name)
}
res.hospitals$h_p <- as.numeric(rep(n_vec, rep(size, length(n_vec))))
res.hospitals$X2 <- factor(res.hospitals$X2, labels = tnam)

edges <- as.numeric(quantile(expit(res.doctors[res.doctors$reg == "Manhattan",]$value), probs = c(0.05, 0.90)))

res.hospitals %>% filter(reg == "Manhattan") %>% ggvis(x = ~ expit(value)) %>% 
  layer_densities(fillOpacity := 0.3) %>% 
  filter(X2 == eval(input_select(unique(as.character(res.hospitals[res.hospitals$reg == "Manhattan",]$X2))))) %>%
  layer_densities(fillOpacity := 0.3) %>%
  scale_numeric("x", domain = edges)

for(i in 1:H){
  tm.path <- paste("hospital_convergence_", i, ".png", sep = "")
  png(paste(tm.path))
  print( ggplot(res.hospitals[as.numeric(res.hospitals$X2) == i,], aes(x = X1, y = value)) + 
    geom_line(colour = rainbow(H, s = 0.7, v = 0.8)[i]) +
    xlab("Iteration") + ylab("logit(P)") + theme(legend.position = "none")
  )
}
graphics.off()

# newlabel <- docs %>% arrange(r.pdeath, h.pdeath) %$% as.character(unique(Hospital.Name))
# gdat_h$trylab <- factor(gdat_h$X2, levels = newlabel)
# g1 <-  ggplot(gdat_h, aes(x = expit(value), fill=X2)) + geom_density(colour="black", alpha = 0.1) + 
#   ylab("Density") + xlab("p(death)") + theme(legend.position = "none") +
#   facet_grid(reg ~.) + scale_x_continuous(limits=c(0,0.04),breaks=seq(0,0.10,0.01))
# g1

tm <- (res[["gamma"]][-(1:burnin),])[keep,]
res.doctors <- melt(tm)
n_vec <- character(D)
for(i in 1:D) n_vec[i] <- as.character(hospitals[[doctors[[i]]$up]]$name)
res.doctors$hos <- rep(n_vec, rep(size, D))
for(i in 1:D) n_vec[i] <- as.character(regions[[hospitals[[doctors[[i]]$up]]$up]]$name)
res.doctors$reg <- rep(n_vec, rep(size, D))

for(i in 1:D){
  tm.path <- paste("doctor_convergence_", i, ".png", sep = "")
  png(paste(tm.path))
  print( ggplot(res.doctors[as.numeric(res.doctors$X2) == i,], aes(x = X1, y = value)) + 
           geom_line(colour = rainbow(D, s = 0.7, v = 0.8)[i]) +
           xlab("Iteration") + ylab("logit(P)") + theme(legend.position = "none")
  )
  dev.off()
}
graphics.off()

g2 <- ggplot(res.doctors[res.doctors$reg == "Manhattan",], aes(x = expit(value))) + 
  geom_density(alpha=0.2, aes(fill = factor(X2)), position = "stack") +
  ylab("Density") + xlab("p(death)")  +  scale_y_continuous(breaks=NULL) +
  facet_grid(hos ~., scales = "free_y", drop = FALSE) +
  theme(legend.position = "none", strip.text.y = element_text(size = 8, angle = 0)) + 
  scale_x_continuous(limits=c(quantile(expit(res.doctors[res.doctors$reg == "Manhattan",]$value), probs = c(0.05, 0.90)))) 

#scale_x_continuous(limits=c(0,0.05),breaks=seq(0,0.05,0.01)) 
g2


  

tm <- res[[3]]
mc.beta <- mcmc(GAMMA)
sum.mcmc <- summary(mc.beta)
sum.mcmc$quantiles[]
autocorr.plot(mc.beta, lag.max = 50)
# 
rejectionRate(mc.beta)
# densplot(mc.beta)
