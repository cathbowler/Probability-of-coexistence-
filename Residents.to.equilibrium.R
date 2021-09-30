### Code pairwise coexistence #### This should run fine once I get Trcy working to get species to equilibrium 
# then will need to use separate scripts to pull out array values for each species to invade into (with the 
# different neighbours present) And can add carrying capacity back in if numbers are huge. 
rm(list=ls())
# import data

library(coda)
library(rstan)
library(reshape2)
library(ggplot2)


sp <- 9
# removing MOMO for now

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
do.lottery.new.seeds <- function(lived, lambda, alpha_intra) {
  N_tp1 <- lived*lambda/(1+alpha_intra*lived)
  return(N_tp1)
}


time <- 200
runs <- 4500
#ARCA--------------------------------------------------
seed.abundances <- matrix(data=NA,nrow=runs,ncol=time)
seed.abundances[,1] <- 2 # start seed abundances at 2 for each species

germ=germination$arca
surv=survival$arca
lambdas=arca$lambda
intras=arca$alpha_intra

# lottery model
#for (r in 1:runs){
  for (t in 1:(time-1)) {
    p <- sample(seq(1, 4500),1) # germ
    #for (t in 1:(time-1)) { 
    for (r in 1:runs){ 
    # calculate how many seeds remain in the seedbank
    seedbank <- (1-germ[p])*surv[p]*seed.abundances[r,t]
    
    # calculate germination from the seedbank
    germinated <- seed.abundances[r,t]*germ[p]
    
    # calculate how many germinants live that year
    #lived <- min(germinated, ag_carrying$x[y])
    lived=germinated
    
    # calculate how many seeds are produced by living germinants 
    new.seeds <- do.lottery.new.seeds(lived=lived, lambda=lambdas[p], alpha_intra=intras[p])
    
    # add above ground seeds produced to those that survived and didn't germinate from the seedbank
    seed.abundances[r, t+1] <- seedbank + new.seeds
  }
}

# plot the seed abundance trajectory through time (200 years)
library(ggplot2)
library(reshape2)
data.melt.arca <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.arca <- as.data.frame(spline(data.melt.arca$time, data.melt.arca$abundance))

arca.equ.mean <- mean(seed.abundances[,200])
# save as rdata
abundances.arca <- seed.abundances
save(abundances.arca, file="nolottery.abundance.arca.rdata")
remove(seed.abundances)

#DAGL-----------------------------------------------
seed.abundances <- matrix(data=NA,nrow=runs,ncol=time)
seed.abundances[,1] <- 2

germ=germination$dagl
surv=survival$dagl
lambdas=dagl$lambda
intras=dagl$alpha_intra
#intras <- ifelse(intras > 0, 0, intras)

# lottery model
#for (r in 1:runs){
  for (t in 1:(time-1)) {
    p <- sample(seq(1, 4500),1) # germ
    #for (t in 1:(time-1)) {
    for (r in 1:runs){
    # calculate how many seeds remain in the seedbank
    seedbank <- (1-germ[p])*surv[p]*seed.abundances[r,t]
    
    # calculate germination from the seedbank
    germinated <- seed.abundances[r,t]*germ[p]
    
    # calculate how many germinants live that year
    #lived <- min(germinated, ag_carrying$x[y])
    lived = germinated
    
    # calculate how many seeds are produced by living germinants 
    new.seeds <- do.lottery.new.seeds(lived=lived, lambda=lambdas[p], alpha_intra=intras[p])
    
    # add above ground seeds produced to those that survived and didn't germinate from the seedbank
    seed.abundances[r, t+1] <- seedbank + new.seeds
  }
}

data.melt.dagl <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.dagl <- as.data.frame(spline(data.melt.dagl$time, data.melt.dagl$abundance))

dagl.equ.mean <- mean(seed.abundances[,200])

abundances.dagl <- seed.abundances
save(abundances.dagl, file="nolottery.abundance.dagl.rdata")
remove(seed.abundances)

# HYGL ---------------------------------------------
seed.abundances <- matrix(data=NA,nrow=runs,ncol=time)
seed.abundances[,1] <- 2

germ=germination$hygl
surv=survival$hygl
lambdas=hygl$lambda
intras=hygl$alpha_intra
#intras <- ifelse(intras > 0, 0, intras)
# set x to r (r being 4500), then don't use sample, just us x everywhere 
# lottery model
#for (r in 1:runs){
 for (t in 1:(time-1)) {
    p <- sample(seq(1, 4500),1) # germ

    #for (t in 1:(time-1)) {
    for (r in 1:runs){
    # calculate how many seeds remain in the seedbank
    seedbank <- (1-germ[p])*surv[p]*seed.abundances[r,t]
    
    # calculate germination from the seedbank
    germinated <- seed.abundances[r,t]*germ[p]
    
    # calculate how many germinants live that year
    #lived <- min(germinated, ag_carrying$x[y])
    lived = germinated
    
    # calculate how many seeds are produced by living germinants 
    new.seeds <- do.lottery.new.seeds(lived=lived, lambda=lambdas[p], alpha_intra=intras[p])
    
    # add above ground seeds produced to those that survived and didn't germinate from the seedbank
    seed.abundances[r, t+1] <- seedbank + new.seeds
  }
}

data.melt.hygl <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.hygl <- as.data.frame(spline(data.melt.hygl$time, data.melt.hygl$abundance))

hygl.equ.mean <- mean(seed.abundances[,200])

abundances.hygl <- seed.abundances
save(abundances.hygl, file="nolottery.abundance.hygl.rdata")
remove(seed.abundances)

# PLDE -------------------------------------
seed.abundances <- matrix(data=NA,nrow=runs,ncol=time)
seed.abundances[,1] <- 2

germ=germination$plde
surv=survival$plde
lambdas=plde$lambda
intras=plde$alpha_intra

# lottery model
#for (r in 1:runs){
  for (t in 1:(time-1)) {
    p <- sample(seq(1, 4500),1) # germ

    #for (t in 1:(time-1)) {
    for (r in 1:runs){
    # calculate how many seeds remain in the seedbank
    seedbank <- (1-germ[p])*surv[p]*seed.abundances[r,t]
    
    # calculate germination from the seedbank
    germinated <- seed.abundances[r,t]*germ[p]
    
    # calculate how many germinants live that year
    #lived <- min(germinated, ag_carrying$x[y])
    lived = germinated
    
    # calculate how many seeds are produced by living germinants 
    new.seeds <- do.lottery.new.seeds(lived=lived, lambda=lambdas[p], alpha_intra=intras[p])
    
    # add above ground seeds produced to those that survived and didn't germinate from the seedbank
    seed.abundances[r, t+1] <- seedbank + new.seeds
  }
}

data.melt.plde <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.plde <- as.data.frame(spline(data.melt.plde$time, data.melt.plde$abundance))

plde.equ.mean <- mean(seed.abundances[,200])

abundances.plde <- seed.abundances
save(abundances.plde, file="nolottery.abundance.plde.rdata")
remove(seed.abundances)

# POCA ---------------------------------------
seed.abundances <- matrix(data=NA,nrow=runs,ncol=time)
seed.abundances[,1] <- 2

germ=germination$poca
surv=survival$poca
lambdas=poca$lambda
intras=poca$alpha_intra
#intras <- ifelse(intras > 0, 0, intras)

# lottery model
#for (r in 1:runs){
  for (t in 1:(time-1)) {
    x1 <- sample(seq(1, 4500),1) # germ
    x2 <- sample(seq(1, 4500),1) # surv
    x3 <- sample(seq(1, 4500),1) # lambdas, intras
    #for (t in 1:(time-1)) {
    for (r in 1:runs){
    # calculate how many seeds remain in the seedbank
    seedbank <- (1-germ[x1])*surv[x2]*seed.abundances[r,t]
    
    # calculate germination from the seedbank
    germinated <- seed.abundances[r,t]*germ[x1]
    
    # calculate how many germinants live that year
    #lived <- min(germinated, ag_carrying$x[y])
    lived = germinated
    
    # calculate how many seeds are produced by living germinants 
    new.seeds <- do.lottery.new.seeds(lived=lived, lambda=lambdas[x3], alpha_intra=intras[x3])
    
    # add above ground seeds produced to those that survived and didn't germinate from the seedbank
    seed.abundances[r, t+1] <- seedbank + new.seeds
  }
}

data.melt.poca <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.poca <- as.data.frame(spline(data.melt.poca$time, data.melt.poca$abundance))

poca.equ.mean <- mean(seed.abundances[,200])

abundances.poca <- seed.abundances
save(abundances.poca, file="nolottery.abundance.poca.rdata")
remove(seed.abundances)

# TRCY -------------------------------------------------
seed.abundances <- matrix(data=NA,nrow=runs,ncol=time)
seed.abundances[,1] <- 2

germ=germination$trcy
surv=survival$trcy
lambdas=trcy$lambda
intras=trcy$alpha_intra
#intras <- ifelse(intras > 0, 0, intras)

# lottery model
#for (r in 1:runs){
  for (t in 1:(time-1)) {
    p <- sample(seq(1, 4500),1) # germ

    #for (t in 1:(time-1)) {
    for (r in 1:runs){
    # calculate how many seeds remain in the seedbank
    seedbank <- (1-germ[p])*surv[p]*seed.abundances[r,t]
    
    # calculate germination from the seedbank
    germinated <- seed.abundances[r,t]*germ[p]
    
    # calculate how many germinants live that year
    #lived <- min(germinated, ag_carrying$x[y])
    lived = germinated
    
    # calculate how many seeds are produced by living germinants 
    new.seeds <- do.lottery.new.seeds(lived=lived, lambda=lambdas[p], alpha_intra=intras[p])
    
    # add above ground seeds produced to those that survived and didn't germinate from the seedbank
    seed.abundances[r, t+1] <- seedbank + new.seeds
  }
}

data.melt.trcy <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.trcy <- as.data.frame(spline(data.melt.trcy$time, data.melt.trcy$abundance))

trcy.equ.mean <- mean(seed.abundances[,200])

abundances.trcy <- seed.abundances
save(abundances.trcy, file="nolottery.abundance.trcy.rdata")
remove(seed.abundances)
# VERO ---------------------------------------
seed.abundances <- matrix(data=NA,nrow=runs,ncol=time)
seed.abundances[,1] <- 2

germ=germination$vero
surv=survival$vero
lambdas=vero$lambda
intras=vero$alpha_intra

# lottery model
#for (r in 1:runs){
  for (t in 1:(time-1)) {
    p <- sample(seq(1, 4500),1) # germ
    #for (t in 1:(time-1)) { 
    for (r in 1:runs){
    # calculate how many seeds remain in the seedbank
    seedbank <- (1-germ[p])*surv[p]*seed.abundances[r,t]
    
    # calculate germination from the seedbank
    germinated <- seed.abundances[r,t]*germ[p]
    
    # calculate how many germinants live that year
    #lived <- min(germinated, ag_carrying$x[y])
    lived = germinated
    
    # calculate how many seeds are produced by living germinants 
    new.seeds <- do.lottery.new.seeds(lived=lived, lambda=lambdas[p], alpha_intra=intras[p])
    
    # add above ground seeds produced to those that survived and didn't germinate from the seedbank
    seed.abundances[r, t+1] <- seedbank + new.seeds
  }
}

data.melt.vero <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.vero <- as.data.frame(spline(data.melt.vero$time, data.melt.vero$abundance))


vero.equ.mean <- mean(seed.abundances[,200])

abundances.vero <- seed.abundances
save(abundances.vero, file="nolottery.abundance.vero.rdata")
remove(seed.abundances)

# MEDI ----------------------------------
seed.abundances <- matrix(data=NA,nrow=runs,ncol=time)
seed.abundances[,1] <- 2

germ=germination$medi
surv=survival$medi
lambdas=medi$lambda
intras=medi$alpha_intra

# lottery model
#for (r in 1:runs){
  for (t in 1:(time-1)) {
    p <- sample(seq(1, 4500),1) # germ

    #for (t in 1:(time-1)) { 
    for (r in 1:runs){
    # calculate how many seeds remain in the seedbank
    seedbank <- (1-germ[p])*surv[p]*seed.abundances[r,t]
    
    # calculate germination from the seedbank
    germinated <- seed.abundances[r,t]*germ[p]
    
    # calculate how many germinants live that year
    #lived <- min(germinated, ag_carrying$x[y])
    lived = germinated
    
    # calculate how many seeds are produced by living germinants 
    new.seeds <- do.lottery.new.seeds(lived=lived, lambda=lambdas[p], alpha_intra=intras[p])
    
    # add above ground seeds produced to those that survived and didn't germinate from the seedbank
    seed.abundances[r, t+1] <- seedbank + new.seeds
  }
}

data.melt.medi <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.medi <- as.data.frame(spline(data.melt.medi$time, data.melt.medi$abundance))

medi.equ.mean <- mean(seed.abundances[,200])

abundances.medi <- seed.abundances
save(abundances.medi, file="nolottery.abundance.medi.rdata")
remove(seed.abundances)

# PEAI -------------------------------------------
seed.abundances <- matrix(data=NA,nrow=runs,ncol=time)
seed.abundances[,1] <- 2

germ=germination$peai
surv=survival$peai
lambdas=peai$lambda
intras=peai$alpha_intra

# lottery model
#for (r in 1:runs){
  for (t in 1:(time-1)) {
    p <- sample(seq(1, 4500),1)
    #for (t in 1:(time-1)) {
    for (r in 1:runs){
    
    # calculate how many seeds remain in the seedbank
    seedbank <- (1-germ[p])*surv[p]*seed.abundances[r,t]
    
    # calculate germination from the seedbank
    germinated <- seed.abundances[r,t]*germ[p]
    
    # calculate how many germinants live that year 
    #lived <- min(germinated, ag_carrying$x[y])
    lived = germinated
    
    # calculate how many seeds are produced by living germinants 
    new.seeds <- do.lottery.new.seeds(lived=lived, lambda=lambdas[p], alpha_intra=intras[p])
    
    # add above ground seeds produced to those that survived and didn't germinate from the seedbank
    seed.abundances[r, t+1] <- seedbank + new.seeds
  }
}

data.melt.peai <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.peai <- as.data.frame(spline(data.melt.peai$time, data.melt.peai$abundance))

peai.equ.mean <- mean(seed.abundances[,200])

abundances.peai <- seed.abundances
save(abundances.peai, file="nolottery.abundance.peai.rdata")
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

