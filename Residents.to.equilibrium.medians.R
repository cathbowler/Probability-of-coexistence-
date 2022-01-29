##############################################################################################
# This code simulates population growth to an equilibrium for each focal species, 
# using the median value of each parameter posterior 
##############################################################################################

rm(list=ls())
# import data

library(coda)
library(rstan)
library(reshape2)
library(ggplot2)


sp <- 9

sp_list <- c("dagl", "hygl", "plde", 
             "poca", "trcy", "vero", "arca", "medi", "peai")

source("Bayes.data.R")
# -----------------------------------------------------------------------------------------

# # single species equilibrium calculations NO LOTTERY MODEL
# do.abundance.single.species <- function(N_t, surv, g, lambda, alpha_intra) {
#   N_tp1 <- N_t*surv*(1-g) + N_t*g*lambda/(1+alpha_intra*N_t*g)
#   return(N_tp1)
# }

# lottery model
do.new.seeds <- function(lived, lambda, alpha_intra) {
  N_tp1 <- lived*lambda/(1+alpha_intra*lived)
  return(N_tp1)
}


time <- 200
runs <- 1000
#ARCA--------------------------------------------------
seed.abundances <- matrix(data=NA,ncol=time)
seed.abundances[,1] <- 2 # start seed abundances at 2 for each species

germ=median(germination$arca)
surv=median(survival$arca)
lambdas=median(arca$lambda)
intras=median(arca$alpha_intra)

# lottery model
for (t in 1:(time-1)) {
  #x1 <- sample(seq(1, 4500),1) # germ
  #x2 <- sample(seq(1, 4500),1) # surv
  #x3 <- sample(seq(1, 4500),1) # lambdas, intras
  #for (t in 1:(time-1)) { # place this line here to not change posterior sample value each time step
  # calculate how many seeds remain in the seedbank
  seedbank <- (1-germ)*surv*seed.abundances[t]
  
  # calculate germination from the seedbank
  germinated <- seed.abundances[t]*germ
  
  # calculate how many germinants live that year
  #lived <- min(germinated, ag_carrying$x[y])
  lived=germinated
  
  # calculate how many seeds are produced by living germinants 
  new.seeds <- do.new.seeds(lived=lived, lambda=lambdas, alpha_intra=intras)
  
  # add above ground seeds produced to those that survived and didn't germinate from the seedbank
  seed.abundances[t+1] <- seedbank + new.seeds
}

# plot the seed abundance trajectory through time (200 years)
library(ggplot2)
library(reshape2)
data.melt.arca <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.arca <- as.data.frame(spline(data.melt.arca$time, data.melt.arca$abundance))

arca.equ.median <- median(seed.abundances[,200])
# save as rdata
abundances.arca <- seed.abundances
save(abundances.arca, file="median.abundance.arca2.rdata")
remove(seed.abundances)

#DAGL-----------------------------------------------
seed.abundances <- matrix(data=NA,ncol=time)
seed.abundances[,1] <- 2

germ=median(germination$dagl)
surv=median(survival$dagl)
lambdas=median(dagl$lambda)
intras=median(dagl$alpha_intra)
#intras <- ifelse(intras > 0, 0, intras)

# lottery model
for (t in 1:(time-1)) {
  # x1 <- sample(seq(1, 4500),1) # germ
  # x2 <- sample(seq(1, 4500),1) # surv
  # x3 <- sample(seq(1, 4500),1) # lambdas, intras
  
  # calculate how many seeds remain in the seedbank
  seedbank <- (1-germ)*surv*seed.abundances[t]
  
  # calculate germination from the seedbank
  germinated <- seed.abundances[t]*germ
  
  # calculate how many germinants live that year
  #lived <- min(germinated, ag_carrying$x[y])
  lived = germinated
  
  # calculate how many seeds are produced by living germinants 
  new.seeds <- do.new.seeds(lived=lived, lambda=lambdas, alpha_intra=intras)
  
  # add above ground seeds produced to those that survived and didn't germinate from the seedbank
  seed.abundances[t+1] <- seedbank + new.seeds
}


data.melt.dagl <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.dagl <- as.data.frame(spline(data.melt.dagl$time, data.melt.dagl$abundance))

dagl.equ.median <- median(seed.abundances[,200])

abundances.dagl <- seed.abundances
save(abundances.dagl, file="median.abundance.dagl2.rdata")
remove(seed.abundances)

# HYGL ---------------------------------------------
seed.abundances <- matrix(data=NA,ncol=time)
seed.abundances[,1] <- 2

germ=median(germination$hygl)
surv=median(survival$hygl)
lambdas=median(hygl$lambda)
intras=median(hygl$alpha_intra)
#intras <- ifelse(intras > 0, 0, intras)

# lottery model
for (t in 1:(time-1)) {
  #x1 <- sample(seq(1, 4500),1) # germ
  #x2 <- sample(seq(1, 4500),1) # surv
  #x3 <- sample(seq(1, 4500),1) # lambdas, intras
  
  # calculate how many seeds remain in the seedbank
  seedbank <- (1-germ)*surv*seed.abundances[t]
  
  # calculate germination from the seedbank
  germinated <- seed.abundances[t]*germ
  
  # calculate how many germinants live that year
  #lived <- min(germinated, ag_carrying$x[y])
  lived = germinated
  
  # calculate how many seeds are produced by living germinants 
  new.seeds <- do.new.seeds(lived=lived, lambda=lambdas, alpha_intra=intras)
  
  # add above ground seeds produced to those that survived and didn't germinate from the seedbank
  seed.abundances[t+1] <- seedbank + new.seeds
}

data.melt.hygl <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.hygl <- as.data.frame(spline(data.melt.hygl$time, data.melt.hygl$abundance))

hygl.equ.median <- median(seed.abundances[,200])

abundances.hygl <- seed.abundances
save(abundances.hygl, file="median.abundance.hygl2.rdata")
remove(seed.abundances)

# PLDE -------------------------------------
seed.abundances <- matrix(data=NA,ncol=time)
seed.abundances[,1] <- 2

germ=median(germination$plde)
surv=median(survival$plde)
lambdas=median(plde$lambda)
intras=median(plde$alpha_intra)

# lottery model
for (t in 1:(time-1)) {
  #x1 <- sample(seq(1, 4500),1) # germ
  #x2 <- sample(seq(1, 4500),1) # surv
  #x3 <- sample(seq(1, 4500),1) # lambdas, intras
  
  # calculate how many seeds remain in the seedbank
  seedbank <- (1-germ)*surv*seed.abundances[t]
  
  # calculate germination from the seedbank
  germinated <- seed.abundances[t]*germ
  
  # calculate how many germinants live that year
  #lived <- min(germinated, ag_carrying$x[y])
  lived = germinated
  
  # calculate how many seeds are produced by living germinants 
  new.seeds <- do.new.seeds(lived=lived, lambda=lambdas, alpha_intra=intras)
  
  # add above ground seeds produced to those that survived and didn't germinate from the seedbank
  seed.abundances[t+1] <- seedbank + new.seeds
}

data.melt.plde <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.plde <- as.data.frame(spline(data.melt.plde$time, data.melt.plde$abundance))

plde.equ.median <- median(seed.abundances[,200])

abundances.plde <- seed.abundances
save(abundances.plde, file="median.abundance.plde2.rdata")
remove(seed.abundances)

# POCA ---------------------------------------
seed.abundances <- matrix(data=NA,ncol=time)
seed.abundances[,1] <- 2

germ=median(germination$poca)
surv=median(survival$poca)
lambdas=median(poca$lambda)
intras=median(poca$alpha_intra)
#intras <- ifelse(intras > 0, 0, intras)

# lottery model
for (t in 1:(time-1)) {
  #x1 <- sample(seq(1, 4500),1) # germ
  #x2 <- sample(seq(1, 4500),1) # surv
  #x3 <- sample(seq(1, 4500),1) # lambdas, intras
  
  # calculate how many seeds remain in the seedbank
  seedbank <- (1-germ)*surv*seed.abundances[t]
  
  # calculate germination from the seedbank
  germinated <- seed.abundances[t]*germ
  
  # calculate how many germinants live that year
  #lived <- min(germinated, ag_carrying$x[y])
  lived = germinated
  
  # calculate how many seeds are produced by living germinants 
  new.seeds <- do.new.seeds(lived=lived, lambda=lambdas, alpha_intra=intras)
  
  # add above ground seeds produced to those that survived and didn't germinate from the seedbank
  seed.abundances[t+1] <- seedbank + new.seeds
}

data.melt.poca <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.poca <- as.data.frame(spline(data.melt.poca$time, data.melt.poca$abundance))

poca.equ.median <- median(seed.abundances[,200])

abundances.poca <- seed.abundances
save(abundances.poca, file="median.abundance.poca2.rdata")
remove(seed.abundances)

# TRCY -------------------------------------------------
seed.abundances <- matrix(data=NA,ncol=time)
seed.abundances[,1] <- 2

germ=median(germination$trcy)
surv=median(survival$trcy)
lambdas=median(trcy$lambda)
intras=median(trcy$alpha_intra)
#intras <- ifelse(intras > 0, 0, intras)

# lottery model
for (t in 1:(time-1)) {
  #x1 <- sample(seq(1, 4500),1) # germ
  #x2 <- sample(seq(1, 4500),1) # surv
  #x3 <- sample(seq(1, 4500),1) # lambdas, intras
  
  # calculate how many seeds remain in the seedbank
  seedbank <- (1-germ)*surv*seed.abundances[t]
  
  # calculate germination from the seedbank
  germinated <- seed.abundances[t]*germ
  
  # calculate how many germinants live that year
  #lived <- min(germinated, ag_carrying$x[y])
  lived = germinated
  
  # calculate how many seeds are produced by living germinants 
  new.seeds <- do.new.seeds(lived=lived, lambda=lambdas, alpha_intra=intras)
  
  # add above ground seeds produced to those that survived and didn't germinate from the seedbank
  seed.abundances[t+1] <- seedbank + new.seeds
}

data.melt.trcy <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.trcy <- as.data.frame(spline(data.melt.trcy$time, data.melt.trcy$abundance))

trcy.equ.median <- median(seed.abundances[,200])

abundances.trcy <- seed.abundances
save(abundances.trcy, file="median.abundance.trcy2.rdata")
remove(seed.abundances)

# VERO ---------------------------------------
seed.abundances <- matrix(data=NA,ncol=time)
seed.abundances[,1] <- 2

germ=median(germination$vero)
surv=median(survival$vero)
lambdas=median(vero$lambda)
intras=median(vero$alpha_intra)

# lottery model
for (t in 1:(time-1)) {
  #x1 <- sample(seq(1, 4500),1) # germ
  #x2 <- sample(seq(1, 4500),1) # surv
  #x3 <- sample(seq(1, 4500),1) # lambdas, intras
  
  # calculate how many seeds remain in the seedbank
  seedbank <- (1-germ)*surv*seed.abundances[t]
  
  # calculate germination from the seedbank
  germinated <- seed.abundances[t]*germ
  
  # calculate how many germinants live that year
  #lived <- min(germinated, ag_carrying$x[y])
  lived = germinated
  
  # calculate how many seeds are produced by living germinants 
  new.seeds <- do.new.seeds(lived=lived, lambda=lambdas, alpha_intra=intras)
  
  # add above ground seeds produced to those that survived and didn't germinate from the seedbank
  seed.abundances[t+1] <- seedbank + new.seeds
}

data.melt.vero <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.vero <- as.data.frame(spline(data.melt.vero$time, data.melt.vero$abundance))


vero.equ.median <- median(seed.abundances[,200])

abundances.vero <- seed.abundances
save(abundances.vero, file="median.abundance.vero2.rdata")
remove(seed.abundances)

# MEDI ----------------------------------
seed.abundances <- matrix(data=NA,ncol=time)
seed.abundances[,1] <- 2

germ=median(germination$medi)
surv=median(survival$medi)
lambdas=median(medi$lambda)
intras=median(medi$alpha_intra)

# lottery model
for (t in 1:(time-1)) {
  #x1 <- sample(seq(1, 4500),1) # germ
  #x2 <- sample(seq(1, 4500),1) # surv
  #x3 <- sample(seq(1, 4500),1) # lambdas, intras
  
  # calculate how many seeds remain in the seedbank
  seedbank <- (1-germ)*surv*seed.abundances[t]
  
  # calculate germination from the seedbank
  germinated <- seed.abundances[t]*germ
  
  # calculate how many germinants live that year
  #lived <- min(germinated, ag_carrying$x[y])
  lived = germinated
  
  # calculate how many seeds are produced by living germinants 
  new.seeds <- do.new.seeds(lived=lived, lambda=lambdas, alpha_intra=intras)
  
  # add above ground seeds produced to those that survived and didn't germinate from the seedbank
  seed.abundances[t+1] <- seedbank + new.seeds
}

data.melt.medi <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.medi <- as.data.frame(spline(data.melt.medi$time, data.melt.medi$abundance))

medi.equ.median <- median(seed.abundances[,200])

abundances.medi <- seed.abundances
save(abundances.medi, file="median.abundance.medi2.rdata")
remove(seed.abundances)

# PEAI -------------------------------------------
seed.abundances <- matrix(data=NA,ncol=time)
seed.abundances[,1] <- 2

germ=germination$peai
surv=survival$peai
lambdas=peai$lambda
intras=median(peai$alpha_intra)

#x1 <- sample(seq(1, 4500),1) # germ
#x2 <- sample(seq(1, 4500),1) # surv
#x3 <- sample(seq(1, 4500),1) # lambdas, intras

# lottery model
for (t in 1:(time-1)) {
  x <- sample(seq(1, 4500),1)
  #y <- sample(dim(ag_carrying)[1], 1)
  
  # calculate how many seeds remain in the seedbank
  seedbank <- (1-germ)*surv*seed.abundances[t]
  
  # calculate germination from the seedbank
  germinated <- seed.abundances[t]*germ
  
  # calculate how many germinants live that year 
  #lived <- min(germinated, ag_carrying$x[y])
  lived = germinated
  
  # calculate how many seeds are produced by living germinants 
  new.seeds <- do.new.seeds(lived=lived, lambda=lambdas, alpha_intra=intras)
  
  # add above ground seeds produced to those that survived and didn't germinate from the seedbank
  seed.abundances[t+1] <- seedbank + new.seeds
}

data.melt.peai <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.peai <- as.data.frame(spline(data.melt.peai$time, data.melt.peai$abundance))

peai.equ.median <- median(seed.abundances[,200])

abundances.peai <- seed.abundances
save(abundances.peai, file="median.abundance.peai2.rdata")
remove(seed.abundances)
