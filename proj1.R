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
base <- read.xlsx("cardiac.xls",1,header=TRUE) 
docs <- base[base$Procedure == "CABG",]
docs <- docs %>% group_by(Detailed.Region) %>% mutate(r.cases = sum(Number.of.Cases), r.death = sum(Number.of.Deaths)) %>%
        mutate(r.pdeath = r.death/r.cases) %>% mutate(r.odeath = r.pdeath/(1-r.pdeath))
docs <- docs %>% group_by(Hospital.Name) %>% mutate(h.cases = sum(Number.of.Cases), h.death = sum(Number.of.Deaths)) %>% 
        mutate(h.pdeath = h.death/h.cases) %>% mutate(h.odeath = h.pdeath/(1-h.pdeath))
docs <- docs %>% mutate(pdeath = Number.of.Deaths/Number.of.Cases) %>% mutate(odeath = pdeath/(1-pdeath))

docs$edeaths <- docs$Expected.Mortality.Rate * docs$Number.of.Cases
docs$odeath <- ifelse(docs$odeath == 0, -10, log(docs$odeath))

docs$Physician.Name <- factor(as.character(docs$Physician.Name), labels = unique(docs$Physician.Name))
docs$Hospital.Name <- factor(as.character(docs$Hospital.Name), labels = unique(docs$Hospital.Name))

#docs %>% group_by(Hospital.Name) %>% summarise(n = length(unique(Detailed.Region))) %$% table(n)
#doc_hosps <- docs %>% group_by(Physician.Name) %>% summarise(n = length(unique(Hospital.Name))) %$% n
#docs %>% group_by(Physician.Name) %>% mutate(n = length(unique(Detailed.Region))) %$% table(n)

R <- length(unique(docs$Detailed.Region))
H <- length(unique(docs$Hospital.Name))
D <- length(unique(docs$Physician.Name))
hos.N <-  docs %>% group_by(Detailed.Region) %>% summarise(tlen = length(unique(Hospital.Name))) %$% tlen
doc.N <- docs %>% group_by(Hospital.Name) %>% summarise(tlen = length(unique(Physician.Name))) %$% tlen

#where_Hosp times a beta vector returns sum of betas by region
where_hosp <- matrix(FALSE, nrow = R, ncol = H)
colnames(where_hosp) <- unique(docs$Hospital.Name)
rownames(where_hosp) <- unique(docs$Detailed.Region)
for(i in as.numeric(unique(docs$Hospital.Name))) where_hosp[unique(as.numeric(docs[as.numeric(docs$Hospital.Name) == i,]$Detailed.Region)), i] <- TRUE

where_docs <- matrix(FALSE, nrow = H, ncol = D)
colnames(where_docs) <- unique(docs$Physician.Name)
rownames(where_docs) <- unique(docs$Hospital.Name)
for(i in as.numeric(unique(docs$Physician.Name))) {
  temp <- as.numeric(unique( docs[ as.numeric(docs$Physician.Name) == i, ]$Hospital.Name )) 
  where_docs[temp, i] <- TRUE
}

#Create list containing regional, hospital and doctor-level information
odeath <- docs[,c("odeath","Physician.Name")]  %>% distinct() %$% odeath
odeath[odeath==0] <- 0.0001
tm.dat <- docs[, c("Hospital.Name", "Number.of.Cases", "Number.of.Deaths")] %>% group_by(Hospital.Name) %>% mutate(all.deaths = sum(Number.of.Deaths), all.cases = sum(Number.of.Cases)) 
tm.dat <- tm.dat[,c(1,4,5)] %>% distinct()
reg_list <- doc_list <- hos_list <- list()
for(i in as.numeric(unique(docs$Detailed.Region))){
  bool <- (as.numeric(docs$Detailed.Region) == i)
  hos <- as.numeric(unique(docs[bool,]$Hospital.Name))
  doc <- as.numeric(unique(docs[bool,]$Physician.Name))
  e_deathrate <-  docs %>% filter(as.numeric(Detailed.Region) == i) %>% group_by(Detailed.Region) %>% 
                  mutate(mu = sum(edeaths/100)/sum(Number.of.Cases)) %$% unique(mu)
  deathrate <-  docs %>% filter(as.numeric(Detailed.Region) == i)  %>% group_by(Detailed.Region) %>%
                mutate(mu = sum(Number.of.Deaths)/sum(Number.of.Cases)) %$% unique(mu)
  reg_list[[i]] <- list(hos = hos, doc = doc, name = as.character(factor(docs$Detailed.Region)[i]), e_mu = e_deathrate, mu = deathrate)
}
for(i in as.numeric(unique(docs$Hospital.Name))){
  bool <- (as.numeric(docs$Hospital.Name) == i)
  li<- as.numeric(unique(docs[bool,]$Detailed.Region))
  do<- as.numeric(unique(docs[bool,]$Physician.Name))
  #delta.prior <- list(df = 1, SS = diag(var(log(docs[bool,]$odeath)),sum(bool)))
  tm <- with(tm.dat[as.numeric(tm.dat$Hospital.Name) == i,], cbind(all.cases, all.deaths))
  e_deathrate <- docs %>% filter(as.numeric(Hospital.Name) == i) %>% group_by(Hospital.Name) %>% 
                  mutate(mu = sum(edeaths/100)/sum(Number.of.Cases)) %$% unique(mu)
  deathrate <-  as.numeric(docs[bool,"h.pdeath"][1,1]) 
  r.deathrate <-  reg_list[[li]]$mu
  hos_list[[i]] <- list(regs = li, dat = tm, docs = do, e_mu = e_deathrate, mu = deathrate, mu.r = r.deathrate)
}
for(i in as.numeric(unique(docs$Physician.Name))){
  bool <- (as.numeric(docs$Physician.Name) == i)
  li <- as.numeric(unique(docs[bool,]$Hospital.Name))
  re <- as.numeric(unique(docs[bool,]$Detailed.Region))
  tm <- with(docs[bool,], cbind(Number.of.Cases, Number.of.Deaths))
  e_deathrate <- docs %>% filter(as.numeric(Physician.Name) == i) %>% group_by(Physician.Name) %>% 
    mutate(mu = sum(edeaths/100)/sum(Number.of.Cases)) %$% unique(mu)
  deathrate <-  docs %>% filter(as.numeric(Physician.Name) == i)  %>% group_by(Physician.Name) %>%
    mutate(mu = sum(Number.of.Deaths)/sum(Number.of.Cases)) %$% unique(mu)
  doc_list[[i]] <- list(hos = li, dat = tm, name = as.character(factor(docs$Hospital.Name)[i]), regs = re)
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
