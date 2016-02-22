#The libraries we need
library(xlsx) #Read the excel file
library(foreach) #Used to parallelize
library(dplyr) #Used for data management
library(magrittr) #Pipe operators in R
library(mvtnorm) #multivariate normal 
library(MCMCpack) #inverse Wishart and inverse Gamma

#Read the data and switch from % to probability of death
docs<- read.xlsx("cardiac.xls",1,header=TRUE) 
docs <- docs[docs$Procedure == "CABG",]
  
docs$pdeath <- docs$Observed.Mortality.Rate/100
#shortcuts we'll need:
  #total number of regions, hospitals and doctors:
  R <- length(unique(docs$Detailed.Region))
  H <- length(unique(docs$Hospital.Name))
  D <- length(unique(docs$Physician.Name))
  #number of hospitals in each region
  hos.N <- docs %>% group_by(Detailed.Region) %>% summarise(tlen = length(unique(Hospital.Name))) %$% tlen
  #variance of regional prior (sd 0.25%)
  reg.var <- (0.0025)^2
  #variance of hospital-level prior (sd 0.5%)
  hos.var <- (0.005)^2
  #variance of doctor-level prior (sd 1%)
  doc.var <- (0.01)^2
  #Expected prior death probability for the whole of NY
  death.prior <- 0.01 
  #Number of iterations
  I <- 100
  
# PRIORS
  #Theta, the fixed regions effect, is multivariate-normal distributed with mean mu.prior.E and variance mu.prior.V
  #Gamma.i, the random hospital effect, is multivariate-normal distributed with mean 0 and variance matrix hos.i.prior.V
  #Detla.i.j, the random doctor effect, is multivariate-normal distributed with mean 0 and variance matrix doc.i.j.prior.V
  
  mu.prior.E = rep(death.prior, R)
  mu.prior.V = diag(R)*reg.var
  for(i in 1:R)
  {
    tnam <- paste(eval(sprintf("hos.%d.prior.V", i)))
    assign(tnam, diag(hos.N[i])*hos.var)
    docs.N <- docs %>% filter(as.numeric(Detailed.Region) == i) %>% group_by(Hospital.Name) %>% summarise(tlen = length(unique(Physician.Name))) %$% tlen
    for(j in 1:hos.N[i]){
      tnam <- paste(eval(sprintf("doc.%d.%d.prior.V", i, j)))
      assign(tnam, diag(docs.N[j])*doc.var)
    }
  }

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

