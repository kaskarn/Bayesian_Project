#The libraries we need
library(xlsx) #Read the excel file
library(dplyr) #Used for data management
library(magrittr) #Pipe operators in R
library(mvtnorm) #multivariate normal 
library(MCMCpack) #inverse Wishart and inverse Gamma
library(reshape) #to graph things!
library(ggplot2) #to graph things nicely!

#inverse logit function
expit <- function(theta){
  return(exp(theta)/(1+exp(theta)))
}
#logit function
logit <- function(theta){
  return(log(theta/(1-theta)))
}
#Iterations
I <- 150000
ptm <- proc.time()

#Read the data and switch from % to probability of death
docs <- read.xlsx("cardiac.xls",1,header=TRUE) %>% filter(Procedure == "CABG") %>%
  mutate( pdeath = Number.of.Deaths/Number.of.Cases, 
          odeath = ifelse(Number.of.Deaths == 0, 0.000001, 
                          Number.of.Deaths / (Number.of.Cases - Number.of.Deaths)), 
          edeaths = Expected.Mortality.Rate/100 * Number.of.Cases) %>%
  group_by(Detailed.Region) %>% mutate(
    r.cases = sum(Number.of.Cases), 
    r.deaths = sum(Number.of.Deaths), 
    r.pdeath = r.deaths/r.cases, 
    r.odeath = r.pdeath/(1-r.pdeath), 
    r.e_pdeath = sum(edeaths)/r.cases) %>%
  group_by(Hospital.Name) %>% mutate(
    h.cases = sum(Number.of.Cases), 
    h.deaths = sum(Number.of.Deaths), 
    h.pdeath = h.deaths/h.cases,
    h.odeath = h.deaths/(h.cases - h.deaths), 
    h.e_pdeath = sum(edeaths)/h.cases) %>% 
  ungroup() %>% mutate(
    Physician.Name = factor(as.character(Physician.Name), labels = unique(Physician.Name)), 
    Hospital.Name = factor(as.character(Hospital.Name), labels = unique(Hospital.Name)),
    dnum = row_number(), 
    hnum = as.numeric(Hospital.Name), 
    rnum = as.numeric(Detailed.Region)) %>%
  arrange(Physician.Name, Hospital.Name)
saveRDS(docs, file = "docs.rds")

D <- nrow(docs)
R <- length(unique(docs$Detailed.Region))
H <- length(unique(docs$Hospital.Name))

#Create lists containing regional, hospital and doctor-level information, including priors!
regions <- hospitals <- doctors <- list()
theta.now <-  docs[,c("rnum", "r.odeath")] %>% distinct() %>% arrange(rnum) %$% log(r.odeath)
theta.prior.v <- 0.1
sigma.prior.df <- 10
delta.prior.df <- 20

sigma.now <- numeric(R)
for(i in 1:R){
  tm <- docs[docs$rnum == i,]
  rname <- unique(tm$Detailed.Region)
  dl <- unique(tm$hnum)
  d <- unique(tm$r.deaths)
  pd <- unique(tm$r.pdeath)
  epd <- unique(tm$r.e_pdeath)
  s.p.v <- tm[tm$rnum == i,c("rnum", "h.odeath")] %>% distinct() %$% var(log(h.odeath))
  sigma.now[i] <- 10
  regions[[i]] <- list(name = rname, down = dl, deaths = d, p_death = pd, e_pdeath = epd, sigma.prior.v = sigma.now[i])
}

delta.now <- beta.now <- numeric(H)
for(i in 1:H){
  tm <- docs[docs$hnum == i,]
  hname <- unique(tm$Hospital.Name)
  dl <- unique(tm$dnum)
  ul <- unique(tm$rnum)
  d <- unique(tm$h.deaths)
  ca <- unique(tm$h.cases)
  pd <- unique(tm$h.pdeath)
  epd <- unique(tm$h.e_pdeath)
  d.p <- tm[,c("dnum", "odeath")] %>% distinct() %$% var(log(odeath))
  delta.now[i] <- 10
  hospitals[[i]] <- list(name = hname, up = ul,  down = dl, deaths = d, cases = ca, p_death = pd, e_pdeath = epd,
                         delta.prior.v = delta.now[i])
  beta.now[i] <- tm[,c("hnum", "h.odeath")] %>% distinct() %$% log(h.odeath)
}
gamma.now <- numeric(D)
for(i in 1:D){
  tm <- docs[docs$dnum == i,]
  dname <- unique(tm$Physician.Name)
  ul <- unique(tm$hnum)
  d <- unique(tm$Number.of.Deaths)
  ca <- unique(tm$Number.of.Cases)
  pd <- unique(tm$pdeath)
  epd <- unique(tm$Expected.Mortality.Rate/100)
  doctors[[i]] <- list(name = dname, up = ul, deaths = d, cases = ca, p_death = pd, e_pdeath = epd)
  gamma.now[i] <- tm[,c("dnum", "odeath")] %>% distinct() %$% log(odeath)
}

theta.prior.e <- numeric(R)
for(i in 1:R) theta.prior.e[i] <- logit(regions[[i]]$e_pdeath)

#Bookkeeping
pb <- txtProgressBar(style=3)

n_vec <- numeric(R)
THETA <- SIGMA <- matrix(nrow = I, ncol = R)
for(i in 1:R) n_vec[i] <- as.numeric(regions[[i]]$name)
colnames(SIGMA) <- colnames(THETA) <- n_vec

n_vec <- numeric(H)
for(i in 1:H) n_vec[i] <- as.numeric(hospitals[[i]]$name)
BETA <- DELTA <- matrix(nrow = I, ncol = H)
colnames(BETA) <- colnames(DELTA) <- n_vec

n_vec <- numeric(D)
GAMMA <- matrix(nrow = I, ncol = D)
for(i in 1:D) n_vec[i] <- i
colnames(GAMMA) <- n_vec

#Launch algorithm!
for(t in 1:I){
  #GIBBS
  #Update Theta and sigma
  for(i in 1:R){
    H.now <- length(regions[[i]]$down)
    if(H.now >1){
      b_mean <- mean(beta.now[regions[[i]]$down])
      v <- solve(solve(theta.prior.v) + H.now*solve(sigma.now[i]))
      e <- v*(theta.prior.e*solve(theta.prior.v) + H.now*solve(sigma.now[i])*b_mean)
      theta.now[i] <- rnorm(1, e, sqrt(v))
    }else theta.now[i] <- rnorm(1, (beta.now[i] + theta.prior.e)/2, theta.prior.v)
  }
  for(i in 1:R){
    df <- sigma.prior.df + H.now
    ss <- regions[[i]]$sigma.prior.v*sigma.prior.df + sum((beta.now[regions[[i]]$down]-theta.now[i])^2)
    sigma.now[i] <- rinvgamma(1, df/2, ss/2)
  }
  
  #METROPOLIS
  #Hospitals
  for(i in 1:H) {
    if(length(regions[[hospitals[[i]]$up]]$down) == 1){ 
      beta.now[i] = theta.now[hospitals[[i]]$up]
    }else{
      beta.prop <- rnorm(1, beta.now[i], sqrt(sigma.now[hospitals[[i]]$up]/2))
      if(is.nan(expit(beta.prop))) next
      p.beta.prop <- dnorm(beta.prop, theta.now[hospitals[[i]]$up], sqrt(sigma.now[hospitals[[i]]$up]), log = TRUE)
      p.beta.now <- dnorm(beta.now[i], theta.now[hospitals[[i]]$up], sqrt(sigma.now[hospitals[[i]]$up]), log = TRUE)
      p.dat.prop <- dbinom(hospitals[[i]]$deaths, hospitals[[i]]$cases, prob=expit(beta.prop), log = TRUE)
      p.dat.now <- dbinom(hospitals[[i]]$deaths, hospitals[[i]]$cases, prob=expit(beta.now[i]), log = TRUE)
      r = p.beta.prop + p.dat.prop - p.beta.now - p.dat.now
      if(log(runif(1)) < r) beta.now[i] <- beta.prop
    }
  }
  for(i in 1:H){
    D.now <- length(hospitals[[i]]$down)
    df <- delta.prior.df + D.now
    ss <- hospitals[[i]]$delta.prior.v*delta.prior.df + sum((gamma.now[hospitals[[i]]$down] - beta.now[i])^2)
    delta.now[i] <- rinvgamma(1, df/2, ss/2)
  }
  for(i in 1:D){
    if(length(hospitals[[doctors[[i]]$up]]$down) == 1){ 
      gamma.now[i] = beta.now[doctors[[i]]$up]
    }else{ 
      gamma.prop <- rnorm(1, gamma.now[i], sqrt(delta.now[doctors[[i]]$up]/2))
      if(is.nan(expit(gamma.prop))) next
      p.gamma.prop <- dnorm(gamma.prop, beta.now[doctors[[i]]$up], sqrt(delta.now[doctors[[i]]$up]), log = TRUE)
      p.gamma.now <- dnorm(gamma.now[i], beta.now[doctors[[i]]$up], sqrt(delta.now[doctors[[i]]$up]), log = TRUE)
      p.dat.prop <- dbinom(doctors[[i]]$deaths, doctors[[i]]$cases, prob=expit(gamma.prop), log= TRUE)
      p.dat.now <- dbinom(doctors[[i]]$deaths, doctors[[i]]$cases, prob=expit(gamma.now[i]), log=TRUE)
      r = p.gamma.prop + p.dat.prop - p.gamma.now - p.dat.now
      if(log(runif(1)) < r) gamma.now[i] <- gamma.prop
    }
  }
  #Store results
  BETA[t,] <- beta.now
  THETA[t,] <- theta.now
  SIGMA[t,] <- sigma.now
  GAMMA[t,] <- gamma.now
  DELTA[t,] <- delta.now
  setTxtProgressBar(pb, t/I)
  
}
proc.time() - ptm

res <- list(beta = BETA, sigma = SIGMA, theta = THETA, gamma = GAMMA, delta = DELTA, hospitals = hospitals, regions = regions, doctors = doctors)
IB <- nrow(res[[1]]) - burnin; thin <- 50
keep <- (!((1:IB)%%thin))
size <- sum(keep)

tm <- (res[[3]][-(1:burnin),])[keep,]
res.regions <- melt(tm)

n_vec <- character(R)
for(i in 1:R) n_vec[i] <- as.character(regions[[i]]$name)
res.regions$X2 <- factor(res.regions$X2, labels = n_vec)
for(i in 1:R) n_vec[i] <- regions[[i]]$e_pdeath
res.regions$rat <- as.numeric(rep(n_vec, rep(size, R)))
res.regions$rat <- expit(res.regions$value) / res.regions$rat

tm <- (res[[1]][-(1:burnin),])[keep,]
res.hospitals <- melt(tm)

n_vec <- tnam <- character(H)
for(i in 1:H) n_vec[i] <- as.character(regions[[hospitals[[i]]$up]]$name)
res.hospitals$reg <- rep(n_vec, rep(size, H))
for(i in 1:H) n_vec[i] <- regions[[hospitals[[i]]$up]]$p_death
res.hospitals$r_p <- as.numeric(rep(n_vec, rep(size, length(n_vec))))

n_vec <- numeric(H)
for(i in 1:H){
  n_vec[i] <- hospitals[[i]]$p_death
  tnam[i] <- as.character(hospitals[[i]]$name)
}
res.hospitals$h_p <- as.numeric(rep(n_vec, rep(size, length(n_vec))))
res.hospitals$X2 <- factor(res.hospitals$X2, labels = tnam)
for(i in 1:H) n_vec[i] <- as.numeric(hospitals[[i]]$e_pdeath)
res.hospitals$rat <- as.numeric(rep(n_vec), rep(size, H))
res.hospitals$rat <- expit(res.hospitals$value)/res.hospitals$rat

tm <- (res[["gamma"]][-(1:burnin),])[keep,]
res.doctors <- melt(tm)

n_vec <- character(D)
for(i in 1:D) n_vec[i] <- as.character(hospitals[[doctors[[i]]$up]]$name)
res.doctors$hos <- rep(n_vec, rep(size, D))
for(i in 1:D) n_vec[i] <- as.character(regions[[hospitals[[doctors[[i]]$up]]$up]]$name)
res.doctors$reg <- rep(n_vec, rep(size, D))
for(i in 1:D) n_vec[i] <- as.character(doctors[[i]]$name)
res.doctors$nam <- rep(n_vec, rep(size, D))

n_vec <- numeric(D)
for(i in 1:D) n_vec[i] <- doctors[[i]]$e_pdeath
res.doctors$rat <- as.numeric(rep(n_vec, rep(size, D)))
res.doctors$rat <- expit(res.doctors$value) / res.doctors$rat
for(i in 1:D) n_vec[i] <- doctors[[i]]$p_death
res.doctors$p_death <- rep(n_vec, rep(size, D))


loc_table <- docs[,c("Physician.Name", "Hospital.Name", "Detailed.Region")] %>% distinct()

saveRDS(loc_table, "loc_table.rds")
saveRDS(hospitals, file="hospitals.rds")
saveRDS(doctors, file="doctors.rds")
saveRDS(regions, file="regions.rds")
saveRDS(res.regions, file = "res_regions.rds")
saveRDS(res.doctors, file = "res_doctors.rds")
saveRDS(res.hospitals, file = "res_hospitals.rds")

###### TESTING ######
res.hospitals <- res.hospitals %>% group_by(X2) %>% mutate(pmean = mean(expit(value)))
newlabel <- res.hospitals %>% arrange(pmean) %$% as.character(unique(X2))
res.hospitals$trylab <- factor(as.character(res.hospitals$X2), labels = newlabel)

ggplot(res.hospitals, aes(x = factor(trylab), y = expit(value))) + geom_boxplot() + coord_flip() 
