#The libraries we need
library(xlsx) #Read the excel file
library(foreach) #Used to parallelize
library(dplyr) #Used for data management
library(magrittr) #Pipe operators in R
library(mvtnorm) #multivariate normal 
library(MCMCpack) #inverse Wishart and inverse Gamma

#Read the data and switch from % to probability of death
docs <- read.xlsx("cardiac.xls",1,header=TRUE) 
docs$pdeath <- docs$Observed.Mortality.Rate/100
docs <- docs %>% group_by(Detailed.Region) %>% mutate(r.cases = sum(Number.of.Cases), r.death = sum(Number.of.Deaths))
docs <- docs %>% group_by(Hospital.Name) %>% mutate(h.cases = sum(Number.of.Cases), h.death = sum(Number.of.Deaths))

docs$r.pdeath <- with(docs, r.death/r.cases)
docs$r.ldeath <- with(docs, r.pdeath/(1+r.pdeath))
docs$h.pdeath <- with(docs, h.death/h.cases)
docs$h.ldeath <- with(docs, h.pdeath/(1+h.pdeath))

docs %>% group_by(Hospital.Name) %>% summarise(n = length(unique(Detailed.Region))) %$% table(n)
doc_hosps <- docs %>% group_by(Physician.Name) %>% summarise(n = length(unique(Hospital.Name))) %$% n
docs %>% group_by(Physician.Name) %>% mutate(n = length(unique(Detailed.Region))) %$% table(n)

hosps <- docs[,c("Detailed.Region", "Hospital.Name")] %>% distinct()

R <- length(unique(docs$Detailed.Region))
H <- length(unique(docs$Hospital.Name))
D <- length(unique(docs$Physician.Name))
hos.N <- docs %>% group_by(Detailed.Region) %>% summarise(tlen = length(unique(Hospital.Name))) %$% tlen
doc.N <- docs %>% group_by(Hospital.Name) %>% summarise(tlen = length(unique(Physician.Name))) %$% tlen

#where_Hosp times a beta vector returns sum of betas by region
where_hosp <- matrix(nrow = R, ncol = H)
colnames(where_hosp) <- unique(docs$Hospital.Name)
rownames(where_hosp) <- unique(docs$Detailed.Region)
for(i in as.numeric(unique(hosps$Detailed.Region))) where_hosp[i,] <- (as.numeric(hosps$Detailed.Region) == i)

where_docs <- matrix(FALSE, nrow = H, ncol = D)
colnames(where_docs) <- unique(docs$Physician.Name)
rownames(where_docs) <- unique(docs$Hospital.Name)
for(i in as.numeric(unique(docs$Physician.Name))) {
  temp <- as.numeric(unique( docs[ as.numeric(docs$Physician.Name) == i, ]$Hospital.Name )) 
  where_docs[temp, i] = TRUE
}

beta.try <- rep(1,H)
test <- cbind(where_hosp%*%beta.try)

#Theta, the fixed regions effect, is R-variate-normal distributed with mean mu.prior.E and variance mu.prior.V
#Sigma, between-hospital variance is R-Inverse-Wishart distributed with df=R and simple shape matrix (for now) 
#Gamma, between-doctor variance is H-Inverse-Wishart distributed with df=H and simple shape matrix (for now) 
r.odeath <- docs[,c("Detailed.Region","r.pdeath")] %>% distinct() %>% mutate(odeath = r.pdeath/(1-r.pdeath)) %$% odeath
h.odeath <- docs[,c("Hospital.Name","h.pdeath")] %>% distinct() %>% mutate(odeath = h.pdeath/(1-h.pdeath)) %$% odeath

theta.now <- theta.prior.E <- log(r.odeath)
theta.prior.V <- diag(var(log(r.odeath)),R)
sigma.prior.df <- R
sigma.now <- sigma.prior.s <- diag(var(log(r.odeath)),R)
delta.prior.df <- H
delta.prior.s <- diag(var(log(h.odeath)),H)
beta.now <- log(h.odeath)

#Number of iterations
I <- 100
THETA <- SIGMA <- matrix(nrow=I, ncol=R)
BETA <- DELTA <- matrix(nrow=I, ncol=H)

#Gibbs step

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
sigma.now <- riwish(sigma.prior.df + H, solve(sigma.prior.s + SS))
  
#Metropolis step
beta.prop <- rmvnorm(1, (!where_hosp[,1]*theta.now)+beta.now[1], 0.5*sigma.now)
tm_p <- exp(beta.prop%*%where_hosp[,1])
num <- log(dbinom(100,110,prob=(tm_p/(1+tm_p)))) 
p_bet <- rmvnorm(as.vector(beta.prop[1,]), mean=theta.now, sigma=sigma.now)

unlist(beta.prop[1,])

sigma <- matrix(c(4,2,2,3), ncol=2)
x <- rmvnorm(n=500, mean=c(1,2), sigma=sigma)

S <- rep(0, n.regions)
for(i in 1:n.hospitals) S <- S + sum(beta.now[i,] - theta.now)%*%t(beta.now[i,] - theta.now)

return(sigma)






# THINGS WE'LL STORE, AND THEIR STARTING VALUES (PRIOR)
  #fixed effects and between-region variance
  THETA <- SIGMA <- matrix(nrow=I, ncol=R)
  THETA[1,] <- mu.prior.E
  #fixed region effect + hospital random effects and between-hospital variance
  BETA <- DELTA <- matrix(nrow=I, ncol=H)
  BETA[1,] <- rep(mu.prior.E, H)
  #fixed region effect + hospital random effects + doctor random effects
  GAMMA  <- matrix(nrow=I, ncol=D)
  GAMMA[1,] <- rep(mu.prior.E, D)

# OUR KEY FUNCTIONS
  #Sample Theta from full conditional posterior distribution
  get_theta <- function(n.hospitals, n.regions, theta.prior.V, sigma.now, theta.prior.E, beta.now){
    b_mean <- rep(0, n.regions)
    for(i in 1:n.regions) b_mean <- b_mean + sum(beta.now[i,])/n.regions
    v = solve(solve(theta.prior.V) + n.hospitals*sigma.now)
    E = V *(theta.prior.V*theta.prior.E + n.hospitals*sigma.now%*%b_mean)
    theta <- rmvnorm(1, E, V)
    return(theta)
  }
  #Sample between-hospital variance from full conditional posterior distribution
  get_sigma <- function(n.hospitals, n.regions, sigma.prior.nu, sigma.prior.S, theta.now, beta.now){
    S <- rep(0, n.regions)
    for(i in 1:n.hospitals) S <- S + sum(beta.now[i,] - theta.now)%*%t(beta.now[i,] - theta.now)
    sigma <- riwish(sigma.prior.nu + n.hospitals, solve(sigma.prior.S + S))
    return(sigma)
  }
  #Sample between-doctor variance from full conditional posterior distribution
  get_delta <- function(n.doctors, n.hospitals, delta.prior.nu, delta.prior.S, beta.now, gamma.now){
    S <- rep(0, n.hospitals)
    for(i in 1:n.doctors) S <- S + sum(gamma.now[i,] - beta.now)%*%t(gamma.now[i,] - beta.now)
    sigma <- riwish(sigma.prior.nu + n.hospitals, solve(sigma.prior.S + S))
    return(delta)
  }
# DESCRIPTION OF ALGORITHM:
  #A. Sample fixed effects -- Gibbs
  #B. Hospital step:
    #Sample fixed + hospital variance -- Gibbs
    #Sample fixed + hospital estimate -- Metropolis
  #C. Doctor step:
    #Sample fixed + hospital + doctor variance -- Gibbs
    #Sample fixed + hospital + doctor estimate -- Metropolis

