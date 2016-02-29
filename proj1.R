#The libraries we need
library(xlsx) #Read the excel file
library(foreach) #Used to parallelize
library(dplyr) #Used for data management
library(magrittr) #Pipe operators in R
library(mvtnorm) #multivariate normal 
library(MCMCpack) #inverse Wishart and inverse Gamma
library(reshape) #to graph things!
library(ggplot2) #to graph things nicely!

#Read the data and switch from % to probability of death
docs <- read.xlsx("cardiac.xls",1,header=TRUE) %>% filter(Procedure == "CABG") %>% 
  mutate(pdeath = Number.of.Deaths/Number.of.Cases, 
         odeath = ifelse(Number.of.Deaths == 0, 0.000001, Number.of.Deaths / (Number.of.Cases - Number.of.Deaths))) %>%
  mutate(edeaths = Expected.Mortality.Rate/100 * Number.of.Cases) %>%
  group_by(Detailed.Region) %>% 
  mutate(r.cases = sum(Number.of.Cases), r.deaths = sum(Number.of.Deaths)) %>%
  mutate(r.pdeath = r.deaths/r.cases) %>% mutate(r.odeath = r.pdeath/(1-r.pdeath)) %>%
  mutate(r.e_pdeath = sum(edeaths)/r.cases) %>%
  group_by(Hospital.Name) %>% 
  mutate(h.cases = sum(Number.of.Cases), h.deaths = sum(Number.of.Deaths)) %>% 
  mutate(h.pdeath = h.deaths/h.cases, h.odeath = h.deaths/(h.cases - h.deaths)) %>% 
  mutate(h.e_pdeath = sum(edeaths)/h.cases) %>% ungroup() %>%
  mutate(Physician.Name = factor(as.character(Physician.Name), labels = unique(docs$Physician.Name))) %>%
  mutate(Hospital.Name = factor(as.character(Hospital.Name), labels = unique(docs$Hospital.Name))) %>%
  mutate(dnum = as.numeric(Physician.Name), hnum = as.numeric(Hospital.Name), rnum = as.numeric(Detailed.Region)) %>%
  arrange(Physician.Name)

D <- nrow(docs)
R <- length(unique(docs$Detailed.Region))
H <- length(unique(docs$Hospital.Name))
hos.N <-  docs %>% group_by(Detailed.Region) %>% summarise(tlen = length(unique(Hospital.Name))) %$% tlen
doc.N <- docs %>% group_by(Hospital.Name) %>% summarise(tlen = length(unique(Physician.Name))) %$% tlen

I <- 10000
#Create lists containing regional, hospital and doctor-level information
regions <- hospitals <- doctors <- list()
for(i in 1:R){
  tm <- docs[docs$rnum == i,]
  rname <- unique(tm$Detailed.Region)
  dl <- unique(tm$hnum)
  d <- unique(tm$r.deaths)
  pd <- unique(tm$r.pdeath)
  epd <- unique(tm$r.e_pdeath)
  regions[[i]] <- list(name = rname, down = dl, deaths = d, p_death = pd, e_pdeath = epd, post_e = numeric(I))
}
for(i in 1:H){
  tm <- docs[docs$hnum == i,]
  rname <- unique(tm$Hospital.Name)
  dl <- unique(tm$hnum)
  ul <- unique(tm$rnum)
  d <- unique(tm$r.deaths)
  pd <- unique(tm$r.pdeath)
  epd <- unique(tm$r.e_pdeath)
  hospitals[[i]] <- list(name = rname, down = dl, up = ul, deaths = d, p_death = pd, e_pdeath = epd, post_e = numeric(I))
}
for(i in 1:D){
  tm <- docs[docs$dnum == i,]
  rname <- unique(tm$Detailed.Region)
  ul <- unique(tm$hnum)
  d <- unique(tm$r.deaths)
  pd <- unique(tm$r.pdeath)
  epd <- unique(tm$r.e_pdeath)
  doctors[[i]] <- list(name = rname, up = ul, deaths = d, p_death = pd, e_pdeath = epd, post_e = numeric(I))
}
#priors
r.odeath <- docs[,c("Detailed.Region","r.pdeath")] %>% distinct() %>% mutate(odeath = r.pdeath/(1-r.pdeath)) %$% odeath
h.odeath <- docs[,c("Hospital.Name","h.pdeath")] %>% distinct() %>% mutate(odeath = h.pdeath/(1-h.pdeath)) %$% odeath
d.odeath <- docs[,c("Physician.Name","Number.of.Cases","Number.of.Deaths")] %>% group_by(Physician.Name) %>% 
  mutate(pdeath = sum(Number.of.Deaths)/sum(Number.of.Cases)) %>% distinct() %>% mutate(odeath = pdeath/(1-pdeath)) %$% odeath

theta_all <- theta.prior.E <- docs %>% mutate(tm = log(sum(Number.of.Deaths)/sum(Number.of.Cases))) %>% mutate(tm2 = tm/(1-tm)) %$% unique(tm2)
theta.prior.V <- 100
sigma.prior.df <- 1
sigma.now <- sigma.prior.s <- diag(var(log(r.odeath)),R)
delta.prior.df <- 1
delta.now <- delta.prior.s <- diag(var(log(h.odeath)),H)
beta.now <- log(h.odeath)
tm <- ifelse(d.odeath == 0, -10, log(d.odeath)) #Few docs haven't offed one yet! Way to go, but let's be real.
gamma.now <- (where_docs %*% diag(tm))

#Number of iterations and bookkeeping
I <- 10000
THETA <- SIGMA <- matrix(nrow=I, ncol=R)
colnames(THETA) <- colnames(SIGMA) <- unique(docs$Detailed.Region)
BETA <- DELTA <- matrix(nrow=I, ncol=H)
colnames(BETA) <- colnames(DELTA) <- unique(docs$Hospital.Name)
GAMMA <- matrix(nrow = I, ncol = D)
colnames(GAMMA) <- unique(docs$Physician.Name)
D.h <- length(docs$Physician.Name)
GAMMA <- matrix(nrow = I, ncol = D.h)
pb <- txtProgressBar(style=3)

sdocs <- as.numeric(apply(where_docs, 2, sum)-1)
test <- 0
for(i in 1:length(sdocs)) test <- cbind(test, test[i-1]+sdocs[i-1])
test <- test + 1:length(test)

for(t in 1:I){
  #GIBBS
  #Update Theta
  b_mean <- (where_hosp%*%beta.now)/hos.N
  V = solve(solve(theta.prior.V) + H*solve(sigma.now))
  E = V %*% t(theta.prior.E%*%solve(theta.prior.V) + t(H*solve(sigma.now)%*%b_mean))
  theta.now <- rmvnorm(1, E, V)
  
  #update sigma
  SS <- matrix(0, nrow=R, ncol=R)
  for(i in 1:H) {
    tb <- where_hosp[,i]*(beta.now[i]-theta.now)
    SS <- SS + t(tb)%*%tb
  }
  sigma.now <- riwish(sigma.prior.df + H, solve(sigma.prior.s + SS*H))
  
  #METROPOLIS
  #Hospitals
  for(i in 1:H){
    ind <- unlist(c(hos_list[[i]][1]))
    V <- sigma.now[ind,ind]
    beta.prop <- rnorm(1, beta.now[i], sqrt(V/2))
    p.beta.prop <- dnorm(beta.prop, theta.now[ind], sqrt(V), log = TRUE)
    p.beta.now <- dnorm(beta.now[i], theta.now[ind], sqrt(V), log = TRUE)
    p.dat.prop <- dbinom(hos_list[[i]]$dat[2], hos_list[[i]]$dat[1], prob=exp(beta.prop)/(1+exp(beta.prop)), log=T)
    p.dat.now <- dbinom(hos_list[[i]]$dat[2], hos_list[[i]]$dat[1], prob=exp(beta.now[i])/(1+exp(beta.now[i])), log=T)
    r = p.beta.prop + p.dat.prop - p.beta.now - p.dat.now
    u = log(runif(1))
    if(r > u) beta.now[i] = beta.prop
  }
  #GIBBS
  #update Delta
  SS <- matrix(0, nrow=H, ncol=H)
  for(i in 1:D) {
    tb <- where_docs[,i]*(gamma.now[i]-beta.now)
    SS <- SS + tb%*%t(tb)
  }
  delta.now <- riwish(delta.prior.df + D, solve(delta.prior.s + SS*D))
  
  #METROPOLIS
  #Doctors
  for(i in 1:D){
    ind <- unlist(c(doc_list[[i]]$hos))
    tm.gam <- unlist(gamma.now[ind,i])
    V <- delta.now[ind,ind]; Vp <- V/2
    if(length(ind) == 1){
      funR <- rnorm; V <- sqrt(V); Vp <- sqrt(Vp); funD <- dnorm
    }else
    {
      funR <- rmvnorm; funD <- dmvnorm
    }
    gamma.prop <- funR(1, tm.gam, Vp)
    p.gamma.prop <- funD(gamma.prop, beta.now[ind], V, log=TRUE)
    p.gamma.now <- funD(tm.gam, beta.now[ind], V, log=TRUE)
    p.dat.prop <- sum(dbinom(doc_list[[i]]$dat[,2], doc_list[[i]]$dat[,2], prob=exp(gamma.prop)/(1+exp(gamma.prop)), log=TRUE))
    p.dat.now <- sum(dbinom(doc_list[[i]]$dat[,2], doc_list[[i]]$dat[,2], prob=exp(tm.gam)/(1+exp(tm.gam)), log=TRUE))
    r = p.gamma.prop + p.dat.prop - p.gamma.now - p.dat.now
    u = log(runif(1))
    if(r > u) gamma.now[ind,i] <- gamma.prop
    GAMMA[t,(test[i]):(test[i]+length(ind)-1)] <- gamma.now[ind,i]
  }
  #Store results
  BETA[t,] <- beta.now
  THETA[t,] <- theta.now
  SIGMA[t,] <- diag(sigma.now)
  DELTA[t,] <- diag(delta.now)
  setTxtProgressBar(pb, t/I)
}
gdat <- cbind(melt(BETA))
gdat$reg <- rep(unique(docs$Detailed.Region), hos.N*I)
g <-  ggplot(gdat %>% filter(reg == "Bronx"), aes(x = X1, y=value)) + geom_line(aes(colour=X2)) +
  theme(legend.position="none") + ylab("Value") + xlab("Iteration") + facet_grid(X2 ~.)
g

expit <- function(theta){
  return(exp(theta)/(1+exp(theta)))
}
results <- list(beta = BETA, theta = THETA, sigma = SIGMA, gamma = GAMMA, delta = DELTA)
#saveRDS(results, file="results_mcmc.rds")
