#The libraries we need
library(xlsx) #Read the excel file
library(foreach) #Used to parallelize
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

#Iterations
I <- 3000

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

D <- nrow(docs)
R <- length(unique(docs$Detailed.Region))
H <- length(unique(docs$Hospital.Name))

#Create lists containing regional, hospital and doctor-level information, including priors!
regions <- hospitals <- doctors <- list()
theta.prior.e <- mean(log(unique(docs$r.odeath)))
theta.now <-  docs[,c("rnum", "r.odeath")] %>% distinct() %>% arrange(rnum) %$% log(r.odeath)
theta.prior.v <- 0.1
sigma.prior.df <- delta.prior.df <- 3

sigma.now <- numeric(R)
for(i in 1:R){
  tm <- docs[docs$rnum == i,]
  rname <- unique(tm$Detailed.Region)
  dl <- unique(tm$hnum)
  d <- unique(tm$r.deaths)
  pd <- unique(tm$r.pdeath)
  epd <- unique(tm$r.e_pdeath)
  s.p.v <- tm[tm$rnum == i,c("rnum", "h.odeath")] %>% distinct() %$% var(log(h.odeath))
  sigma.now[i] <- ifelse(is.na(s.p.v), 0.01, s.p.v)
  regions[[i]] <- list(name = rname, down = dl, deaths = d, p_death = pd, e_pdeath = epd, sigma.prior.v = sigma.now[i])
}

delta.now <- beta.now <- numeric(H)
for(i in 1:H){
  tm <- docs[docs$hnum == i,]
  hname <- unique(tm$Hospital.Name)
  dl <- unique(tm$hnum)
  ul <- unique(tm$rnum)
  d <- unique(tm$h.deaths)
  ca <- unique(tm$h.cases)
  pd <- unique(tm$h.pdeath)
  epd <- unique(tm$h.e_pdeath)
  d.p <- tm[,c("dnum", "odeath")] %>% distinct() %$% var(log(odeath))
  delta.now[i] <- ifelse(is.na(d.p), 0.01, d.p)
  hospitals[[i]] <- list(name = hname, up = ul,  down = dl, deaths = d, cases = ca, p_death = pd, e_pdeath = epd,
                         delta.prior.v = delta.now[i])
  beta.now[i] <- tm[,c("hnum", "h.odeath")] %>% distinct() %$% log(h.odeath)
}
gamma.now <- numeric(D)
for(i in 1:D){
  tm <- docs[docs$dnum == i,]
  dname <- unique(tm$Detailed.Region)
  ul <- unique(tm$hnum)
  d <- unique(tm$Number.of.Deaths)
  ca <- unique(tm$Number.of.Cases)
  pd <- unique(tm$pdeath)
  epd <- unique(tm$e_pdeath)
  doctors[[i]] <- list(name = dname, up = ul, deaths = d, cases = ca, p_death = pd, e_pdeath = epd, post_e = numeric(I))
  gamma.now[i] <- tm[,c("dnum", "odeath")] %>% distinct() %$% log(odeath)
}
#Bookkeeping
pb <- txtProgressBar(style=3)
n_vec <- numeric(H)
THETA <- SIGMA <- matrix(nrow = I, ncol = R)
colnames(SIGMA) <- colnames(THETA) <- unique(docs$rnum)
for(i in 1:H) n_vec[i] <- as.numeric(hospitals[[i]]$name)
BETA <- DELTA <- matrix(nrow = I, ncol = H)
colnames(BETA) <- colnames(DELTA) <- n_vec
 matrix(nrow = I, ncol = R)

GAMMA <- matrix(nrow = I, ncol = D)
for(i in 1:D) n_vec[i] <- as.numeric(doctors[[i]]$name)
colnames(GAMMA) <- n_vec

#Launch algorithm!
for(t in 1:I){
  #GIBBS
  #Update Theta and sigma
  for(i in 1:R){
    b_mean <- mean(beta.now[regions[[i]]$down])
    H.now <- length(regions[[i]]$down)
    v <- solve(solve(theta.prior.v) + H.now*solve(sigma.now[i]))
    e <- v*(theta.prior.e*solve(theta.prior.v) + H.now*solve(sigma.now[i])*b_mean)
    theta.now[i] <- rnorm(1, e, sqrt(v))
  }
  for(i in 1:R){
    df <- sigma.prior.df + H.now
    ss <- regions[[i]]$sigma.prior.v*sigma.prior.df + sum((beta.now[regions[[i]]$down]-theta.now[i])^2)
    sigma.now[i] <- rinvgamma(1, df/2, ss/2)
  }
  
  #METROPOLIS
  #Hospitals
  for(i in 1:H){
    beta.prop <- rnorm(1, beta.now[i], sqrt(sigma.now[hospitals[[i]]$up]/2))
    p.beta.prop <- dnorm(beta.prop, theta.now[hospitals[[i]]$up], sqrt(sigma.now[hospitals[[i]]$up]), log = TRUE)
    p.beta.now <- dnorm(beta.now[i], theta.now[hospitals[[i]]$up], sqrt(sigma.now[hospitals[[i]]$up]), log = TRUE)
    p.dat.prop <- dbinom(hospitals[[i]]$deaths, hospitals[[i]]$cases, prob=expit(beta.prop), log=T)
    p.dat.now <- dbinom(hospitals[[i]]$deaths, hospitals[[i]]$cases, prob=expit(beta.now), log=T)
    r = p.beta.prop + p.dat.prop - p.beta.now - p.dat.now
    if(log(runif(1)) < r) beta.now[i] = beta.prop
  }
  for(i in 1:H){
    D.now <- length(hospitals[[i]]$down)
    df <- delta.prior.df + D.now
    ss <- hospitals[[i]]$delta.prior.v + sum((gamma.now[hospitals[[i]]$down] - beta.now[i])^2)
    delta.now[i] <- rinvgamma(1, df/2, ss/2)
  }
  for(i in 1:D){
    gamma.prop <- rnorm(1, gamma.now[i], sqrt(delta.now[doctors[[i]]$up]/2))
    p.gamma.prop <- dnorm(gamma.prop, beta.now[doctors[[i]]$up], sqrt(delta.now[doctors[[i]]$up]), log = TRUE)
    p.gamma.now <- dnorm(gamma.now[i], beta.now[doctors[[i]]$up], sqrt(delta.now[doctors[[i]]$up]), log = TRUE)
    p.dat.prop <- dbinom(doctors[[i]]$deaths, doctors[[i]]$cases, prob=expit(gamma.prop), log=T)
    p.dat.now <- dbinom(doctors[[i]]$deaths, doctors[[i]]$cases, prob=expit(gamma.now), log=T)
    r = p.gamma.prop + p.dat.prop - p.gamma.now - p.dat.now
    if(log(runif(1)) < r) gamma.now[i] = gamma.now
  }
  #Store results
  BETA[t,] <- beta.now
  THETA[t,] <- theta.now
  SIGMA[t,] <- sigma.now
  GAMMA[t,] <- gamma.now
  DELTA[t,] <- delta.now
  setTxtProgressBar(pb, t/I)

}
# gdat <- cbind(melt(BETA))
# n_vec <- character(H)
# for(i in 1:H) n_vec[i] <- as.character(regions[[hospitals[[i]]$up]]$name)
# gdat$reg <- rep(n_vec, rep(I, length(n_vec)))
# for(i in 1:H) n_vec[i] <- as.character(regions[[hospitals[[i]]$up]]$p_death)
# gdat$r_p <- rep(n_vec, rep(I, length(n_vec)))
# for(i in 1:H) n_vec[i] <- hospitals[[i]]$p_death
# gdat$h_p <- as.numeric(rep(n_vec, rep(I, length(n_vec))))
# gdat$X2 <- factor(gdat$X2, labels = levels(docs$Hospital.Name))
# 
# g <-  ggplot(gdat[gdat$reg == "Bronx",], aes(x = X1, y=value)) + geom_line(aes(colour=X2)) +
#   theme(legend.position="none") + ylab("Value") + xlab("Iteration") + facet_grid(X2 ~.)
# g
# 
# reg.now <- "Manhattan"
# g1 <-  ggplot(gdat %>% filter(reg == reg.now), aes(x = expit(value), fill=X2)) + geom_density(aes(colour=X2)) +
#   theme(legend.position="none") + ylab("Value") + xlab("Iteration") +
#   facet_grid(X2 ~.) + scale_x_continuous(limits=c(0,0.05),breaks=seq(0,0.05,0.01))+ 
#   geom_vline(aes(xintercept = 0.016867)) + geom_vline(aes(xintercept = h_p), colour="red")
# g1
saveRDS(c(BETA, SIGMA, THETA), file="results_mcmc_new.rds")
