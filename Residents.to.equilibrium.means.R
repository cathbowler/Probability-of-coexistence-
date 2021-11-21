##############################################################################################
# This code simulates population growth to an equilibrium for each focal species, 
# using the mean value of each parameter posterior 
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

germ=mean(germination$arca)
surv=mean(survival$arca)
lambdas=mean(arca$lambda)
intras=mean(arca$alpha_intra)

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

arca.equ.mean <- mean(seed.abundances[,200])
# save as rdata
abundances.arca <- seed.abundances
save(abundances.arca, file="mean.abundance.arca2.rdata")
remove(seed.abundances)

#DAGL-----------------------------------------------
seed.abundances <- matrix(data=NA,ncol=time)
seed.abundances[,1] <- 2

germ=mean(germination$dagl)
surv=mean(survival$dagl)
lambdas=mean(dagl$lambda)
intras=mean(dagl$alpha_intra)
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

dagl.equ.mean <- mean(seed.abundances[,200])

abundances.dagl <- seed.abundances
save(abundances.dagl, file="mean.abundance.dagl2.rdata")
remove(seed.abundances)

# HYGL ---------------------------------------------
seed.abundances <- matrix(data=NA,ncol=time)
seed.abundances[,1] <- 2

germ=mean(germination$hygl)
surv=mean(survival$hygl)
lambdas=mean(hygl$lambda)
intras=mean(hygl$alpha_intra)
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

hygl.equ.mean <- mean(seed.abundances[,200])

abundances.hygl <- seed.abundances
save(abundances.hygl, file="mean.abundance.hygl2.rdata")
remove(seed.abundances)

# PLDE -------------------------------------
seed.abundances <- matrix(data=NA,ncol=time)
seed.abundances[,1] <- 2

germ=mean(germination$plde)
surv=mean(survival$plde)
lambdas=mean(plde$lambda)
intras=mean(plde$alpha_intra)

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

plde.equ.mean <- mean(seed.abundances[,200])

abundances.plde <- seed.abundances
save(abundances.plde, file="mean.abundance.plde2.rdata")
remove(seed.abundances)

# POCA ---------------------------------------
seed.abundances <- matrix(data=NA,ncol=time)
seed.abundances[,1] <- 2

germ=mean(germination$poca)
surv=mean(survival$poca)
lambdas=mean(poca$lambda)
intras=mean(poca$alpha_intra)
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

poca.equ.mean <- mean(seed.abundances[,200])

abundances.poca <- seed.abundances
save(abundances.poca, file="mean.abundance.poca2.rdata")
remove(seed.abundances)

# TRCY -------------------------------------------------
seed.abundances <- matrix(data=NA,ncol=time)
seed.abundances[,1] <- 2

germ=mean(germination$trcy)
surv=mean(survival$trcy)
lambdas=mean(trcy$lambda)
intras=mean(trcy$alpha_intra)
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

trcy.equ.mean <- mean(seed.abundances[,200])

abundances.trcy <- seed.abundances
save(abundances.trcy, file="mean.abundance.trcy2.rdata")
remove(seed.abundances)

# VERO ---------------------------------------
seed.abundances <- matrix(data=NA,ncol=time)
seed.abundances[,1] <- 2

germ=mean(germination$vero)
surv=mean(survival$vero)
lambdas=mean(vero$lambda)
intras=mean(vero$alpha_intra)

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


vero.equ.mean <- mean(seed.abundances[,200])

abundances.vero <- seed.abundances
save(abundances.vero, file="mean.abundance.vero2.rdata")
remove(seed.abundances)

# MEDI ----------------------------------
seed.abundances <- matrix(data=NA,ncol=time)
seed.abundances[,1] <- 2

germ=mean(germination$medi)
surv=mean(survival$medi)
lambdas=mean(medi$lambda)
intras=mean(medi$alpha_intra)

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

medi.equ.mean <- mean(seed.abundances[,200])

abundances.medi <- seed.abundances
save(abundances.medi, file="mean.abundance.medi2.rdata")
remove(seed.abundances)

# PEAI -------------------------------------------
seed.abundances <- matrix(data=NA,ncol=time)
seed.abundances[,1] <- 2

germ=germination$peai
surv=survival$peai
lambdas=peai$lambda
intras=mean(peai$alpha_intra)

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

peai.equ.mean <- mean(seed.abundances[,200])

abundances.peai <- seed.abundances
save(abundances.peai, file="mean.abundance.peai2.rdata")
remove(seed.abundances)

### SUPP PLOT ###
pdf("Figures/supp.single.sp.equi.pdf")
par(mfrow=c(3,3), mar=c(2,2,2,2), oma=c(4,4,1,1))
plot(spline.arca, main=expression(italic("A. calendula")))
mtext(side=2, line=3, "abundance")
plot(spline.medi, main=expression(italic("M. minima")))
plot(spline.peai, main=expression(italic("P. airoides")))
plot(spline.dagl, main=expression(italic("D. glochidiatus")))
mtext(side=2, line=3, "abundance")
plot(spline.hygl, main=expression(italic("H. glutinosum")))
plot(spline.plde, main=expression(italic("P. debilis")))
plot(spline.poca, main=expression(italic("P. canescens")))
mtext(side=1, line=3, "time (years)")
mtext(side=2, line=3, "abundance")
plot(spline.trcy, main=expression(italic("T. cyanopetala")))
mtext(side=1, line=3, "time (years)")
plot(spline.vero, main=expression(italic("G. rosea")))
mtext(side=1, line=3, "time (years)")
dev.off()

pdf("Figures/supp.single.sp.200.pdf")

#####
pdf("Figures/single.sp.t=200.pdf")
par(mfrow=c(3,3), mar=c(2,2,2,2), oma=c(4,4,1,1))
plot(density(abundances.arca[,200]), main=expression(italic("A. calendula")))
mtext(side=2, line=3, "density")
plot(density(abundances.medi[,200]), main=expression(italic("M. minima")))
plot(density(abundances.peai[,200]), main=expression(italic("P. airoides")))
plot(density(abundances.dagl[,200]), main=expression(italic("D. glochidiatus")))
mtext(side=2, line=3, "density")
plot(density(abundances.hygl[,200]), main=expression(italic("H. glutinosum")))
plot(density(abundances.plde[,200]), main=expression(italic("P. debilis")))
plot(density(abundances.poca[,200]), main=expression(italic("P. canescens")))
mtext(side=2, line=3, "density")
plot(density(abundances.trcy[,200]), main=expression(italic("T. cyanopetala")))
mtext(side=1, line=3, "steady state population size")
plot(density(abundances.vero[,200]), main=expression(italic("G. rosea")))
dev.off()

# save the t200 mean values 
t200.mean.vals <- c(dagl.equ.mean, hygl.equ.mean, plde.equ.mean, poca.equ.mean,
                    trcy.equ.mean, vero.equ.mean, arca.equ.mean, medi.equ.mean,
                    peai.equ.mean)
save(t200.mean.vals, file="t200.mean.values")

