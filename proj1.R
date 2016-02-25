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
docs <- read.xlsx("cardiac.xls",1,header=TRUE) 
docs$pdeath <- docs$Observed.Mortality.Rate/100
docs <- docs %>% group_by(Detailed.Region, Procedure) %>% mutate(r.cases = sum(Number.of.Cases), r.death = sum(Number.of.Deaths))
docs <- docs %>% group_by(Hospital.Name, Procedure) %>% mutate(h.cases = sum(Number.of.Cases), h.death = sum(Number.of.Deaths))

docs <- docs %>% filter(Procedure == "CABG")
docs$Physician.Name <- factor(as.character(docs$Physician.Name), labels = unique(docs$Physician.Name))
docs$Hospital.Name <- factor(as.character(docs$Hospital.Name), labels = unique(docs$Hospital.Name))

docs$r.pdeath <- with(docs, r.death/r.cases)
docs$r.ldeath <- with(docs, r.pdeath/(1+r.pdeath))
docs$h.pdeath <- with(docs, h.death/h.cases)
docs$h.ldeath <- with(docs, h.pdeath/(1+h.pdeath))

#docs %>% group_by(Hospital.Name) %>% summarise(n = length(unique(Detailed.Region))) %$% table(n)
#doc_hosps <- docs %>% group_by(Physician.Name) %>% summarise(n = length(unique(Hospital.Name))) %$% n
#docs %>% group_by(Physician.Name) %>% mutate(n = length(unique(Detailed.Region))) %$% table(n)

R <- length(unique(docs$Detailed.Region))
H <- length(unique(docs$Hospital.Name))
D <- length(unique(docs$Physician.Name))
hos.N <-  docs %>% group_by(Detailed.Region) %>% summarise(tlen = length(unique(Hospital.Name))) %$% tlen
doc.N <- docs %>% group_by(Hospital.Name) %>% summarise(tlen = length(unique(Physician.Name))) %$% tlen

#where_Hosp times a beta vector returns sum of betas by region
where_hosp <- matrix(nrow = R, ncol = H)
colnames(where_hosp) <- unique(docs$Hospital.Name)
rownames(where_hosp) <- unique(docs$Detailed.Region)
tm.dat <- docs[,c("Detailed.Region", "Hospital.Name")] %>% distinct()
for(i in as.numeric(unique(tm.dat$Detailed.Region))) where_hosp[i,] <- (as.numeric(tm.dat$Detailed.Region) == i)

where_docs <- matrix(FALSE, nrow = H, ncol = D)
colnames(where_docs) <- unique(docs$Physician.Name)
rownames(where_docs) <- unique(docs$Hospital.Name)
for(i in as.numeric(unique(docs$Physician.Name))) {
  temp <- as.numeric(unique( docs[ as.numeric(docs$Physician.Name) == i, ]$Hospital.Name )) 
  where_docs[temp, i] = TRUE
}


#Create list containing hospital-level and doctor-level information
tm.dat <- docs[, c("Hospital.Name", "Number.of.Cases", "Number.of.Deaths")] %>% group_by(Hospital.Name) %>% mutate(all.deaths = sum(Number.of.Deaths), all.cases = sum(Number.of.Cases)) 
tm.dat <- tm.dat[,c(1,4,5)] %>% distinct()
doc_list <- hos_list <- list()
for(i in 1:H){
  li<- as.numeric(unique(docs[as.numeric(docs$Hospital.Name) == i,]$Detailed.Region))
  tm <- with(tm.dat[as.numeric(tm.dat$Hospital.Name) == i,], cbind(all.cases, all.deaths))
  hos_list[[i]] <- list(index = li, dat = tm)
}

for(i in 1:D){
  li <- as.numeric(unique(docs[as.numeric(docs$Physician.Name) == i,]$Hospital.Name))
  tm <- with(docs[as.numeric(docs$Physician.Name) == i,], cbind(Number.of.Cases, Number.of.Deaths))
  doc_list[[i]] <- list(index = li, dat = tm)
}

#PRIORS
#Theta, the fixed regions effect, is R-variate-normal distributed with mean mu.prior.E and variance mu.prior.V
#Sigma, between-hospital variance is R-Inverse-Wishart distributed with df=R and simple shape matrix (for now) 
#Gamma, between-doctor variance is H-Inverse-Wishart distributed with df=H and simple shape matrix (for now) 
r.odeath <- docs[,c("Detailed.Region","r.pdeath")] %>% distinct() %>% mutate(odeath = r.pdeath/(1-r.pdeath)) %$% odeath
h.odeath <- docs[,c("Hospital.Name","h.pdeath")] %>% distinct() %>% mutate(odeath = h.pdeath/(1-h.pdeath)) %$% odeath
d.odeath <- docs[,c("Physician.Name","Number.of.Cases","Number.of.Deaths")] %>% group_by(Physician.Name) %>% 
            mutate(pdeath = sum(Number.of.Deaths)/sum(Number.of.Cases)) %>% distinct() %>% mutate(odeath = pdeath/(1-pdeath)) %$% odeath

theta.now <- theta.prior.E <- log(r.odeath)
theta.prior.V <- diag(var(log(r.odeath)),R)
sigma.prior.df <- 1
sigma.now <- sigma.prior.s <- diag(var(log(r.odeath)),R)
delta.prior.df <- 1
delta.now <- delta.prior.s <- diag(var(log(h.odeath)),H)
beta.now <- log(h.odeath)
tm <- ifelse(d.odeath == 0, -10, log(d.odeath)) #Few docs haven't offed one yet! Way to go, but let's be real.
gamma.now <- (where_docs %*% diag(tm))

#Number of iterations and bookkeeping
I <- 10000; a.gamma <- rep(0,D) ; a.beta <- rep(0, H)
THETA <- SIGMA <- matrix(nrow=I, ncol=R)
colnames(THETA) <- colnames(SIGMA) <- unique(docs$Detailed.Region)
BETA <- DELTA <- matrix(nrow=I, ncol=H)
colnames(BETA) <- colnames(DELTA) <- unique(docs$Hospital.Name)

RBETA <- matrix(nrow = I, ncol = H)
RGAMMA <- matrix(nrow = I, ncol = D)
D.h <- length(docs$Physician.Name)
GAMMA <- matrix(nrow = I, ncol = D.h)
# GAMMA <- matrix(nrow = I, ncol = D)
colnames(GAMMA) <- unique(docs$Physician.Name)
pb <- txtProgressBar(style=3)

sdocs <- as.numeric(apply(where_docs, 2, sum)-1)
test <- 0
for(i in 2:length(sdocs)) test <- cbind(test, test[i-1]+sdocs[i-1])
test <- test + 1:length(test)
# test <- c(1,test[1:length(test)])
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
  sigma.now <- riwish(sigma.prior.df + 1, solve(sigma.prior.s + SS))
  
  #METROPOLIS
  r.beta <- rep(0, H)
  #Hospitals
  for(i in 1:H){
    ind <- unlist(c(hos_list[[i]][1]))
    V <- sigma.now[ind,ind]
    beta.prop <- rnorm(1, beta.now[i], sqrt(V/2))
    p.beta.prop <- dnorm(beta.prop, theta.now[ind], sqrt(V), log = TRUE)
    p.beta.now <- dnorm(beta.now[i], theta.now[ind], sqrt(V), log = TRUE)
    p.dat.prop <- dbinom(hos_list[[i]]$dat[2], hos_list[[i]]$dat[1], prob=exp(beta.prop)/(1+exp(beta.prop)), log=T)
    p.dat.now <- dbinom(hos_list[[i]]$dat[2], hos_list[[i]]$dat[1], prob=exp(beta.now[i])/(1+exp(beta.now[i])), log=T)
    r.beta[i] = p.beta.prop + p.dat.prop - p.beta.now - p.dat.now
    u = log(runif(1))
    if(r.beta[i] > u){
      a.beta[i] <- a.beta[i] + 1
      beta.now[i] = beta.prop
    } 
  }
  #GIBBS
  #update Delta
  SS <- matrix(0, nrow=H, ncol=H)
  for(i in 1:D) {
    tb <- where_docs[,i]*(gamma.now[i]-beta.now)
    SS <- SS + tb%*%t(tb)
  }
  delta.now <- riwish(delta.prior.df + 1, solve(delta.prior.s + SS))

  #METROPOLIS
  r.gamma <- rep(0, D)
  #Doctors
  for(i in 1:D){
    ind <- unlist(c(doc_list[[i]]$index))
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
    r.gamma[i] = p.gamma.prop + p.dat.prop - p.gamma.now - p.dat.now
    u = log(runif(1))
    if(r.gamma[i] > u){
      gamma.now[ind,i] <- gamma.prop
      a.gamma[i] <- a.gamma[i] + 1
    }
    GAMMA[t,(test[i]):(test[i]+length(ind))] <- gamma.now[ind,i]
  }
  #Store results
  BETA[t,] <- beta.now
  THETA[t,] <- theta.now
  SIGMA[t,] <- diag(sigma.now)
  GAMMA[t,] <- apply(gamma.now, 2, mean)
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
results <- list(beta=BETA, theta=THETA, sigma=SIGMA, gamma=GAMMA, delta=DELTA, rbeta=RBETA, rgamma=RGAMMA)
saveRDS(results, file="results_mcmc2.rds")
saveRDS(docs, file="docs.rds")
saveRDS(where_docs, file="wheredocs.rds")
saveRDS(where_hosp, file="wheredocs.rds")

