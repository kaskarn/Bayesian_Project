library(dplyr) #Used for data management
library(magrittr) #Pipe operators in R
library(MCMCpack) #inverse Wishart and inverse Gamma
library(reshape) #to graph things!
library(ggplot2) #to graph things nicely!

res <- readRDS("results_mcmc_70000.rds")
burnin <- 0; thin <- 1
keep <- !((1:nrow(res[[1]]))%%thin)
if(burnin > 0) keep <- keep[-(1:burnin)]
size <- sum(keep)

tm <- (res[[3]])[keep,]
res.regions <- melt(tm)
for(i in 1:R){
  tm.path <- paste("./convergence/regions/region_convergence_", i, ".png", sep = "")
  png(paste(tm.path), height = 200, width = 200)
  print( ggplot(res.regions[as.numeric(res.regions$X2) == i,], aes(x = X1, y = value)) + 
    geom_line(colour = rainbow(R, s = 0.7, v = 0.8)[i]) +
    xlab("Iteration") + ylab("logit(P)") + theme(legend.position = "none")
  )
  dev.off()
}
graphics.off()

tm <- (res[[1]])[keep,]
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

for(i in 1:H){
  tm.path <- paste("./convergence/hospitals/hospitals_convergence_", i, ".png", sep = "")
  png(paste(tm.path), height = 200, width = 200)
  print( ggplot(res.hospitals[as.numeric(res.hospitals$X2) == i,], aes(x = X1, y = value)) + 
           geom_line(colour = rainbow(H, s = 0.7, v = 0.8)[i]) +
           xlab("Iteration") + ylab("logit(P)") + theme(legend.position = "none")
  )
  dev.off()
}
graphics.off()

tm <- (res[["gamma"]])[keep,]
res.doctors <- melt(tm)
n_vec <- character(D)
for(i in 1:D) n_vec[i] <- as.character(hospitals[[doctors[[i]]$up]]$name)
res.doctors$hos <- rep(n_vec, rep(size, D))
for(i in 1:D) n_vec[i] <- as.character(regions[[hospitals[[doctors[[i]]$up]]$up]]$name)
res.doctors$reg <- rep(n_vec, rep(size, D))

for(i in 1:D){
  tm.path <- paste("./convergence/doctors/doctors_convergence_", i, ".png", sep = "")
  png(paste(tm.path), height = 200, width = 200)
  print( ggplot(res.doctors[as.numeric(res.doctors$X2) == i,], aes(x = X1, y = value)) + 
           geom_line(colour = rainbow(D, s = 0.7, v = 0.8)[i]) +
           xlab("Iteration") + ylab("logit(P)") + theme(legend.position = "none")
  )
  dev.off()
}
graphics.off()

edges <- as.numeric(quantile(expit(res.doctors[res.doctors$reg == "Manhattan",]$value), probs = c(0.05, 0.90)))

res.regions %>% ggvis(x = ~ expit(value)) %>% 
  layer_densities(fillOpacity := 0.7, fill := "#f2ede4") %>% 
  filter(X2 == "Manhattan") %>%
  layer_densities(fillOpacity := 0.7, fill := "#e4f2ed") %>%
  scale_numeric("x", domain = edges) %>% 
  set_options(width = "auto", height = "auto")






res.hospitals %>% filter(reg == "Manhattan") %>% ggvis(x = ~ expit(value)) %>% 
  layer_densities(fillOpacity := 0.3) %>% 
  filter(X2 == eval(input_select(unique(as.character(res.hospitals[res.hospitals$reg == "Manhattan",]$X2))))) %>%
  layer_densities(fillOpacity := 0.3) %>%
  scale_numeric("x", domain = edges)



# newlabel <- docs %>% arrange(r.pdeath, h.pdeath) %$% as.character(unique(Hospital.Name))
# gdat_h$trylab <- factor(gdat_h$X2, levels = newlabel)
# g1 <-  ggplot(gdat_h, aes(x = expit(value), fill=X2)) + geom_density(colour="black", alpha = 0.1) + 
#   ylab("Density") + xlab("p(death)") + theme(legend.position = "none") +
#   facet_grid(reg ~.) + scale_x_continuous(limits=c(0,0.04),breaks=seq(0,0.10,0.01))
# g1



g2 <- ggplot(res.doctors[res.doctors$reg == "Manhattan",], aes(x = expit(value))) + 
  geom_density(alpha=0.2, aes(fill = factor(X2)), position = "stack") +
  ylab("Density") + xlab("p(death)")  +  scale_y_continuous(breaks=NULL) +
  facet_grid(hos ~., scales = "free_y", drop = FALSE) +
  theme(legend.position = "none", strip.text.y = element_text(size = 8, angle = 0)) + 
  scale_x_continuous(limits=c(quantile(expit(res.doctors[res.doctors$reg == "Manhattan",]$value), probs = c(0.05, 0.90)))) 

#scale_x_continuous(limits=c(0,0.05),breaks=seq(0,0.05,0.01)) 
g2

res.doctors[res.doctors$hos == hospitals[[doctors[[1]]$up]]$name,] %>% ggvis(x = ~ expit(value)) %>%
  layer_densities(fillOpacity := 0.7, fill := "#f2ede4") %>% 
  filter(nam == doctors[[hospitals[[doctors[[1]]$up]]$down[1]]]$name) %>%
  layer_densities(fillOpacity := 0.7, fill := "#e4f2ed") %>%
  scale_numeric("x", domain = edges) %>% 
  set_options(width = "auto", height = "auto")
  
for(j in 1:floor(R/4)) for( i in ((j-1)*(floor(R/floor(R/4))) + 1) : (j*floor((R/floor(R/4)))) ) print(c(i,j))
y <- ceiling(H/6)
x <- 6
for(j in 1:6) for(i in ((j-1)*y+1):(j*y) ) {if(H >= i) print(min(i,H))}
((j-1)*ceiling(H/6)+1):(j*ceiling(H/6))

res.regions %>% ggvis(x = ~ expit(value)) %>% 
  layer_densities(fillOpacity := 0.7, fill := "#f2ede4") %>% 
  layer_densities(data = res.regions %>% filter(X2 == "Bronx"), fillOpacity := 0.7, fill := "#f2ede4") %>%
  set_options(width = "auto", height = "auto")

