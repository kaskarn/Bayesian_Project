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
I <- 150
ptm <- proc.time()

##### Read Data #####
docs <- read.xlsx("cardiac.xls",1,header=TRUE) %>%
  group_by(Physician.Name, Procedure) %>% mutate( 
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
  dl <- unique(tm$hnum)
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
  dl <- unique(tm$dnum)
  ul <- unique(tm[tm$pnum == 1,]$rnum) 
  deaths <- c(unique(tm[tm$pnum == 1,]$h.deaths), unique(tm[tm$pnum == 2,]$h.deaths))
  cases <- c(unique(tm[tm$pnum == 1,]$h.cases), unique(tm[tm$pnum == 2,]$h.cases))
  p_death <- c(unique(tm[tm$pnum == 1,]$h.pdeath), unique(tm[tm$pnum == 2,]$h.pdeath))
  e_pdeath <- c(unique(tm[tm$pnum == 1,]$h.e_pdeath), unique(tm[tm$pnum == 2,]$h.e_pdeath))
  p_death_ov <- sum(deaths, na.rm = TRUE) / sum(cases, na.rm = TRUE)
  e_pdeath_ov <- sum(e_pdeath*tm$edeaths)
  hospitals[[i]] <- list(name = name, up = ul,  down = dl, deaths = deaths, cases = cases, 
                         p_death = p_death, e_pdeath = e_pdeath, p_death_ov = p_death_ov)
}
gamma.now <- numeric(D)
### Probably a dumb step. Maybe lapply could streamline? ###
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
  p_death_ov <- sum(deaths, na.rm = TRUE) / sum(cases, na.rm = TRUE)
  doctors[[i]] <- list(name = name, up = ul, deaths = deaths, cases = cases,
                       p_death = p_death, e_pdeath = e_pdeath, p_death_ov = p_death_ov)
}

##### priors #####
sigma.prior.df <- 10; sigma.prior.s <- sigma.now <- 1
sigma.prior.s <- 1
delta.prior.df <- 50; delta.prior.s <- delta.now <- 1
delta.prior.s <- 1
theta.prior.v <- 0.1
theta.prior.e <- theta.now <- matrix(nrow=2, ncol=R)
for(i in 1:R) theta.prior.e[,i] <- theta.now[,i] <- c(logit(regions[[i]]$e_pdeath[1]), 
                                                      logit(regions[[i]]$e_pdeath[2]))

##### Bookkeeping #####
pb <- txtProgressBar(style=3)

reg.names <- lapply(1:R, function(i) regions[[i]]$name)
hos.names <- lapply(1:H, function(i) hospitals[[i]]$name)
doc.names <- lapply(1:D, function(i) doctors[[i]]$name)

THETA <- array(NA, dim=c(I, R, 2)); SIGMA <- matrix(nrow = I, ncol = R)
colnames(THETA[,,1]) <- colnames(THETA[,,2]) <- colnames(SIGMA) <- unlist(reg.names)
BETA <- DELTA <- matrix(nrow = I, ncol = H)
XI <- array(NA, dim=c(I,H,2))
colnames(BETA) <- colnames(DELTA) <- colnames(XI[,,1]) <- colnames(XI[,,2]) <- unlist(hos.names)
GAMMA <- matrix(nrow = I, ncol = D)
PSI <- array(NA, dim=c(I,D,2))
colnames(GAMMA) <- colnames(PSI[,,1]) <- colnames(PSI[,,2]) <- unlist(doc.names)

beta.now <- rep(0,H)
gamma.now <-  rep(0,D)
xi.now <- do.call(rbind, lapply(1:2, function (j) {
  unlist(lapply(1:H, function (i) regions[[hospitals[[i]]$up]]$p_death[j]))
  }))
psi.now <- do.call(rbind, lapply(1:2, function (j) {
  unlist(lapply(1:D, function (i) regions[[hospitals[[max(doctors[[i]]$up, na.rm = TRUE)]]$up]]$p_death[j]))
}))
  
#Launch algorithm!
for(t in 1:I){
  #GIBBS
  #Update Theta and sigma
  for(j in 1:2){
    for(i in 1:R){
      H.now <- length(regions[[i]]$down)
      b_mean <- mean(beta.now[regions[[i]]$down]+theta.now[j,i])
      v <- solve(solve(theta.prior.v) + H.now*solve(sigma.now))
      e <- v*(theta.prior.e[j,i]*solve(theta.prior.v) + H.now*solve(sigma.now)*b_mean)
      theta.now[j,i] <- rnorm(1, e, sqrt(v))
    }
  }
#   Different hospital variance by regions. Keeping it around for a rainy day
#   for(i in 1:R){
#     H.now <- length(regions[[i]]$down)
#     df <- sigma.prior.df + H.now
#     ss <- sigma.prior.df*sigma.prior.s[i] + sum(unlist(lapply(1:R, function (i) {beta.now[regions[[i]]$down]^2} )))
#     sigma.now[i] <- rinvgamma(1, df/2, ss/2)
#   }
  
  df <- sigma.prior.df + H
  ss <- sigma.prior.s * sigma.prior.df + sum(unlist(lapply(1:R, function (i) sum(beta.now[regions[[i]]$down]^2))))
  sigma.now <- rinvgamma(1, df/2, ss/2)
  
  #METROPOLIS
  #Hospitals
  for(i in 1:H) {
    beta.prop <- rnorm(1, beta.now[i], sigma.now/2)
    if(is.nan(expit(beta.prop))) next
    r <- sum(na.rm = TRUE,
      dnorm(beta.prop, 0, sqrt(sigma.now), log = TRUE),
      -dnorm(beta.now[i], 0, sqrt(sigma.now), log = TRUE),
      dbinom(hospitals[[i]]$deaths, hospitals[[i]]$cases, prob=expit(beta.prop + theta.now[,hospitals[[i]]$up]), log = TRUE),
      -dbinom(hospitals[[i]]$deaths, hospitals[[i]]$cases, prob=expit(beta.now[i] + theta.now[,hospitals[[i]]$up]), log = TRUE)
    )
    if(log(runif(1)) < r){
      beta.now[i] <- beta.prop
      xi.now[i] <- sum(expit( beta.now[i] + theta.now[,hospitals[[i]]$up] ) * hospitals[[i]]$cases)/sum(hospitals[[i]]$cases)
    }
#     Different doctor variance by hospital Keeping it around for a rainy day   
#     D.now <- length(hospitals[[i]]$down)
#     df <- delta.prior.df + D.now
#     ss <- delta.prior.df*delta.prior.s[i] + sum(unlist(lapply(1:H, function (j) {gamma.now[hospitals[[j]]$down]^2})))
#     delta.now[i] <- rinvgamma(1, df/2, ss/2)
  }
  
  # attempt to have same doctor variance accross all hospitals Keeping it around for a rainy day
  df <- delta.prior.df + D
  ss <- delta.prior.s * sigma.prior.df + sum(unlist(lapply(1:H, function (i) sum(gamma.now[hospitals[[i]]$down]^2))))
  delta.now <- rinvgamma(1, df/2, ss/2)

  for(i in 1:D){
    gamma.prop <- rnorm(1, gamma.now[i], sqrt(delta.now/2)) 
    prob_now <- do.call(rbind,
      lapply(1:2, function (j) lapply( 1:ncol(doctors[[i]]$up), function (k) {
        theta.now [hospitals[[doctors[[i]]$up[j,k]]]$up] +
        beta.now [doctors[[i]]$up[j,k]] +
        gamma.now [i]
      }))
    )
    prob_p <- do.call(rbind,
      lapply(1:2, function (j) lapply( 1:ncol(doctors[[i]]$up), function (k) {
        theta.now [hospitals[[doctors[[i]]$up[j,k]]]$up] +
        beta.now [doctors[[i]]$up[j,k]] +
        gamma.prop
      }))
    )
    r <- sum(na.rm = TRUE, 
      dbinom(as.vector(doctors[[i]]$deaths), as.vector(doctors[[i]]$cases), prob=expit(unlist(prob_p)), log = TRUE),
      -dbinom(as.vector(doctors[[i]]$deaths), as.vector(doctors[[i]]$cases), prob=expit(unlist(prob_now)), log = TRUE),
      dnorm(gamma.prop, beta.now[j], sqrt(delta.now), log = TRUE),
      -dnorm(gamma.now[i], beta.now[j], sqrt(delta.now), log = TRUE)
    )
    #if(is.na(psi.now[i])) print(i)
    if(log(runif(1)) < r){
      gamma.now[i] <- gamma.prop
      psi.now[i] <- sum(expit(unlist(prob_p)) * as.vector(doctors[[i]]$cases), na.rm = TRUE)/
                    sum(as.vector(doctors[[i]]$cases), na.rm = TRUE)
    }
  }
  #Store results
  THETA[t,,] <- theta.now
  PSI[t,,] <- psi.now
  XI[t,,] <- xi.now
  
  BETA[t,] <- beta.now
  SIGMA[t,] <- sigma.now
  GAMMA[t,] <- gamma.now
  DELTA[t,] <- delta.now
  setTxtProgressBar(pb, t/I)
  
}
proc.time() - ptm

thin <- 5; burnin <- 10; IB <- I - burnin
keep <- c(rep(0,burnin),(!((1:IB)%%thin)))
size <- sum(keep)

res.regions <- melt( THETA[keep,,] )
res.regions$X2 <- unlist(lapply(1:R, function (i) rep(regions[[i]]$name, size*2)))
res.regions$rat <- unlist(lapply(1:R, function (i) rep(regions[[i]]$e_pdeath / regions[[i]]$p_death, size)))

res.hospitals <- melt( XI[keep,] )
res.hospitals$X2 <- unlist(lapply(1:H, function (i) rep(hospitals[[i]]$name, size)))
res.hospitals$region <- unlist(lapply(1:H, function (i) rep(regions[[hospitals[[i]]$up]]$name, size)))
res.hospitals$pdeath <- unlist(lapply(1:H, function (i) rep(regions[[hospitals[[i]]$up]]$p_death, size)))

res.doctors <- melt( GAMMA[keep,] )
res.doctors$X2 <- unlist(lapply(1:D, function (i) rep(doctors[[i]]$name, size)))

res.doc_world <- melt( PSI[keep,,])
res.doc_world$X2 <- unlist(lapply(1:D, function (i) rep(doctors[[i]]$name, size*2)))
res.doc_world$pdeath <- unlist(lapply(1:D, function (i) rep(doctors[[i]]$p_death_ov, size)))

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
