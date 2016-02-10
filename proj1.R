library(xlsx)
library(foreach)

docs<-read.xlsx("cardiac.xls",1,header=TRUE)
docs$pdeath <- docs$Observed.Mortality.Rate/100
docs$ldeath <- with(docs, log(pdeath/(1-pdeath)))

nReg <- length(unique(docs$Detailed.Region))
nHos <- doclist <- NULL 
nDoc <- matrix(nrow = nReg, ncol = 8)
for(i in 1:nReg){
  hoslist <- unique(factor(docs[as.numeric(docs$Detailed.Region)==i,]$Hospital.Name))
  nHos <- cbind(nHos, length(hoslist))
  for(j in 1:nHos[i]){
    doclist <- unique(factor(docs[docs$Hospital.Name == as.character(hoslist[j]),]$Physician.Name))
    nDoc[i,j] <- length(doclist)
  }
}



#########################################################
## WEAKLY INFORMATIVE PRIORS  
########################################################
#### priors for grand mean ###
mu0 <- 0.0325 #from http://www.state.nj.us/health/healthcarequality/cabgs98/cabgs98t.htm
g20 <- 0.001  #since mu0 comes from large literature

### inverse-gamma prior for inter-region variance Tau ###
e20 <- 1
t20 <- 0.01
### inverse-gamma prior for inter-hopsital variance Xi ###
n20 <- 1
x20 <- 0.01
### inverse-gamma prior for inter-doctor variance Psi ###
k20 <- 1
p20 <- 0.01

##########################################################

