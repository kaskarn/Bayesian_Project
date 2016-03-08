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
I <- 1500
ptm <- proc.time()

##### Read Data #####
docs <- read.xlsx("cardiac.xls",1,header=TRUE) %>%
  group_by(Procedure) %>% mutate( 
    pdeath = Number.of.Deaths/Number.of.Cases, 
    odeath = ifelse(Number.of.Deaths == 0, 0.000001, 
                    Number.of.Deaths / (Number.of.Cases - Number.of.Deaths)), 
    edeaths = Expected.Mortality.Rate/100 * Number.of.Cases) %>%
  group_by(Detailed.Region, Procedure) %>% mutate(
    r.cases = sum(Number.of.Cases), 
    r.deaths = sum(Number.of.Deaths), 
    r.pdeath = r.deaths/r.cases, 
    r.odeath = r.pdeath/(1-r.pdeath), 
    r.e_pdeath = sum(edeaths)/r.cases) %>%
  group_by(Hospital.Name, Procedure) %>% mutate(
    h.cases = sum(Number.of.Cases), 
    h.deaths = sum(Number.of.Deaths), 
    h.pdeath = h.deaths/h.cases,
    h.odeath = h.deaths/(h.cases - h.deaths), 
    h.e_pdeath = sum(edeaths)/h.cases) %>% 
  ungroup() %>% mutate(
    Physician.Name = factor(as.character(Physician.Name), labels = unique(Physician.Name)), 
    Hospital.Name = factor(as.character(Hospital.Name), labels = unique(Hospital.Name)),
    dnum = as.numeric(Physician.Name), 
    hnum = as.numeric(Hospital.Name), 
    rnum = as.numeric(Detailed.Region),
    pnum = as.numeric(Procedure)) %>%
  arrange(Physician.Name, Hospital.Name)
saveRDS(docs, file = "docs.rds")

D <- length(unique(docs$Physician.Name))
R <- length(unique(docs$Detailed.Region))
H <- length(unique(docs$Hospital.Name))

##### lookup-tables #####
regions <- hospitals <- doctors <- list()
for(i in 1:R){
  tm <- docs[docs$rnum == i,]
  name <- unique(tm$Detailed.Region)
  dl <- rbind(unique(tm[tm$pnum == 1,]$hnum), unique(tm[tm$pnum == 2,]$hnum))
  deaths <- c(unique(tm[tm$pnum == 1,]$r.deaths), unique(tm[tm$pnum == 2,]$r.deaths))
  cases <- c(unique(tm[tm$pnum == 1,]$r.cases), unique(tm[tm$pnum == 2,]$r.cases))
  p_death <- c(unique(tm[tm$pnum == 1,]$r.pdeath), unique(tm[tm$pnum == 2,]$r.pdeath))
  e_pdeath <- c(unique(tm[tm$pnum == 1,]$r.e_pdeath), unique(tm[tm$pnum == 2,]$r.e_pdeath))
  regions[[i]] <- list(name = name, down = dl, cases = cases,
                       deaths = deaths, p_death = p_death, e_pdeath = e_pdeath)
}
beta.now <- numeric(H)
for(i in 1:H){
  tm <- docs[docs$hnum == i,]
  name <- unique(tm$Hospital.Name)
  dl <- rbind(unique(tm[tm$dnum == 1,]$rnum), unique(tm[tm$dnum == 2,]$rnum))
  ul <- rbind(unique(tm[tm$pnum == 1,]$rnum), unique(tm[tm$pnum == 2,]$rnum))
  deaths <- c(unique(tm[tm$pnum == 1,]$h.deaths), unique(tm[tm$pnum == 2,]$h.deaths))
  cases <- c(unique(tm[tm$pnum == 1,]$h.cases), unique(tm[tm$pnum == 2,]$h.cases))
  p_death <- c(unique(tm[tm$pnum == 1,]$h.pdeath), unique(tm[tm$pnum == 2,]$h.pdeath))
  e_pdeath <- c(unique(tm[tm$pnum == 1,]$h.e_pdeath), unique(tm[tm$pnum == 2,]$h.e_pdeath))
  hospitals[[i]] <- list(name = name, up = ul,  down = dl, deaths = deaths, cases = cases, 
                         p_death = p_death, e_pdeath = e_pdeath)
}
gamma.now <- numeric(D)
### I don't wanna talk about this dumb step. Maybe lapply could streamline? ###
for(i in 1:D){
  tm <- docs[docs$dnum == i,]
  name <- unique(tm$Physician.Name)
  ul <- deaths <- cases <- matrix(nrow = 2, 
    ncol = max(length(unique(tm[tm$pnum == 1,]$hnum)), length(unique(tm[tm$pnum == 2,]$hnum))))
  for(j in 1:ncol(ul)){
    if(is.na(unique(tm[tm$pnum == 1,]$hnum)[j]) == 0){
      ul[1,j] <- unique(tm[tm$pnum == 1,]$hnum)[j]
      cases[1,j] <- tm[tm$pnum == 1 & tm$hnum == unique(tm[tm$pnum == 2,]$hnum)[j],]$Number.of.Cases
      deaths[1,j] <- tm[tm$pnum == 1 & tm$hnum == unique(tm[tm$pnum == 2,]$hnum)[j],]$Number.of.Deaths
    }
    if(is.na(unique(tm[tm$pnum == 2,]$hnum)[j]) == 0){
      ul[2,j] <- unique(tm[tm$pnum == 2,]$hnum)[j]
      cases[2,j] <- tm[tm$pnum == 2 & tm$hnum == unique(tm[tm$pnum == 2,]$hnum)[j],]$Number.of.Cases
      deaths[2,j] <- tm[tm$pnum == 2 & tm$hnum == unique(tm[tm$pnum == 2,]$hnum)[j],]$Number.of.Deaths
    }
  } 
  p_death <- c(unique(tm[tm$pnum == 1,]$pdeath), unique(tm[tm$pnum == 2,]$pdeath))
  e_pdeath <- c(unique(tm[tm$pnum == 1,]$Expected.Mortality.Rate/100), unique(tm[tm$pnum == 2,]$Expected.Mortality.Rate/100))
  doctors[[i]] <- list(name = name, up = ul, deaths = deaths, cases = cases,
                       p_death = p_death, e_pdeath = e_pdeath)
}

##### priors #####
sigma.prior.df <- 10; sigma.prior.s <- sigma.now <- rep(1, R)
delta.prior.df <- 50; delta.prior.s <- delta.now <- rep(1, H)
theta.prior.v <- 0.1
theta.prior.e <- theta.now <- matrix(nrow=2, ncol=R)
for(i in 1:R) theta.prior.e[,i] <- theta.now[,i] <- c(logit(regions[[i]]$e_pdeath[1]), 
                                                      logit(regions[[i]]$e_pdeath[2]))

##### Bookkeeping #####
pb <- txtProgressBar(style=3)

reg.names <- lapply(1:R, function(i) regions[[i]]$name)
hos.names <- lapply(1:H, function(i) hospitals[[i]]$name)
doc.names <- lapply(1:D, function(i) doctors[[i]]$name)

THETA <- array(NA, dim=c(2, I, R))
colnames(THETA[1,,]) <- colnames(THETA[2,,]) <- unlist(reg.names)
BETA <- matrix(nrow = I, ncol = H)
colnames(BETA) <- unlist(hos.names)
GAMMA <- matrix(nrow = I, ncol = D)
colnames(GAMMA) <- unlist(doc.names)

beta.now <- rep(0,H)
gamma.now <- rep(0,D)

#Launch algorithm!
for(t in 1:I){
  #GIBBS
  #Update Theta and sigma
  for(j in 1:2){
    for(i in 1:R){
      H.now <- length(regions[[i]]$down[j,])
      b_mean <- mean(beta.now[regions[[i]]$down[j,]]+theta.now[j,i])
      v <- solve(solve(theta.prior.v) + H.now*solve(sigma.now[i]))
      e <- v*(theta.prior.e[j,i]*solve(theta.prior.v) + H.now*solve(sigma.now[i])*b_mean)
      theta.now[j,i] <- rnorm(1, e, sqrt(v))
    }
  }
  for(i in 1:R){
    H.now <- length(regions[[i]]$down)
    df <- sigma.prior.df + H.now
    ss <- sigma.prior.df*sigma.prior.s + sum(unlist(lapply(1:R, function (i) {beta.now[regions[[i]]$down]^2} )))
    sigma.now[i] <- rinvgamma(1, df/2, ss/2)
  }
  
#   attempt to have same hospital variance accross all regions. Keeping it around for a rainy day
#   df <- sigma.prior.df + H
#   ss <- sigma.prior.s * sigma.prior.df + sum(unlist(lapply(1:R, function (i) sum(beta.now[regions[[i]]$down]^2))))
#   sigma.now <- rinvgamma(1, df/2, ss/2)
  
  #METROPOLIS
  #Hospitals
  for(i in 1:H) {
    l_temp <- unlist(hospitals[[i]]$up)[1]
    beta.prop <- rnorm(1, beta.now[i], sigma.now[l_temp]/2)
    if(is.nan(expit(beta.prop))) next
    r <- sum(na.rm = TRUE,
      dnorm(beta.prop, 0, sqrt(sigma.now[l_temp]), log = TRUE),
      -dnorm(beta.now[i], 0, sqrt(sigma.now[l_temp]), log = TRUE),
      dbinom(hospitals[[i]]$deaths, hospitals[[i]]$cases, prob=expit(beta.prop + theta.now[,hospitals[[i]]$up[1]]), log = TRUE),
      -dbinom(hospitals[[i]]$deaths, hospitals[[i]]$cases, prob=expit(beta.now[i] + theta.now[,hospitals[[i]]$up[1]]), log = TRUE)
    )
    if(log(runif(1)) < r) beta.now[i] <- beta.prop
    
    D.now <- length(hospitals[[i]]$down)
    df <- delta.prior.df + D.now
    ss <- delta.prior.df*delta.prior.s + sum(unlist(lapply(1:H, function (j) {gamma.now[hospitals[[j]]$down]^2})))
    delta.now[i] <- rinvgamma(1, df/2, ss/2)
  }
  
#   attempt to have same doctor variance accross all hospitals Keeping it around for a rainy day
#   df <- delta.prior.df + D
#   ss <- delta.prior.s * sigma.prior.df + sum(unlist(lapply(1:H, function (i) sum(gamma.now[hospitals[[i]]$down]^2))))
#   delta.now <- rinvgamma(1, df/2, ss/2)

  for(i in 1:D){
    l_temp <- unlist(doctors[[i]]$up)[1]
    gamma.prop <- rnorm(1, gamma.now[i], sqrt(delta.now[l_temp]/2))
    #if(is.nan(expit(gamma.prop))) print(i) 
    r <- 0
    #for(j in 1:2) for(k in 1:ncol(doctors[[i]]$up)) {print(theta.now [hospitals[[doctors[[i]]$up[j,k]]]$up[j,]] )}
    prob_now <- do.call(rbind,
      lapply(1:2, function (j) lapply( 1:ncol(doctors[[i]]$up), function (k) {
        theta.now [hospitals[[doctors[[i]]$up[j,k]]]$up[j]] +
        beta.now [doctors[[i]]$up[j,k]] +
        gamma.now [i]
      }))
    )
    prob_p <- do.call(rbind,
      lapply(1:2, function (j) lapply( 1:ncol(doctors[[i]]$up), function (k) {
        theta.now [hospitals[[doctors[[i]]$up[j,k]]]$up[j]] +
        beta.now [doctors[[i]]$up[j,k]] +
        gamma.prop
      }))
    )
    r <- sum(na.rm = TRUE, 
      dbinom(as.vector(doctors[[i]]$deaths), as.vector(doctors[[i]]$cases), prob=expit(unlist(prob_p)), log = TRUE),
      -dbinom(as.vector(doctors[[i]]$deaths), as.vector(doctors[[i]]$cases), prob=expit(unlist(prob_now)), log = TRUE),
      dnorm(gamma.prop, beta.now[j], sqrt(delta.now[j]), log = TRUE),
      -dnorm(gamma.now[i], beta.now[j], sqrt(delta.now[j]), log = TRUE)
    )
    if(log(runif(1)) < r) gamma.now[i] <- gamma.prop
  }
  #Store results
  BETA[t,] <- beta.now
  THETA[,t,] <- theta.now
  #SIGMA[t,] <- sigma.now
  GAMMA[t,] <- gamma.now
  #DELTA[t,] <- delta.now
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
