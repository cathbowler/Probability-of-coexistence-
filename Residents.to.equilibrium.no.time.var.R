##############################################################################################
# This code simulates population growth to an equilibrium for each focal species, 
# NOT sampling the posterior values at each time step 
##############################################################################################

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

# lottery model
do.new.seeds <- function(lived, lambda, alpha_intra) {
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
for (r in 1:runs){ # then make runs=4500 to systematically go through each posterior 
   p <- r
  for (t in 1:(time-1)) {
  #  p <- sample(seq(1, 4500),1) 
    #for (t in 1:(time-1)) { 
    # calculate how many seeds remain in the seedbank
    seedbank <- (1-germ[p])*surv[p]*seed.abundances[r,t]
    
    # calculate germination from the seedbank
    germinated <- seed.abundances[r,t]*germ[p]
    
    # calculate how many germinants live that year
    #lived <- min(germinated, ag_carrying$x[y])
    lived=germinated
    
    # calculate how many seeds are produced by living germinants 
    new.seeds <- do.new.seeds(lived=lived, lambda=lambdas[p], alpha_intra=intras[p])
    
    # add above ground seeds produced to those that survived and didn't germinate from the seedbank
    seed.abundances[r, t+1] <- seedbank + new.seeds
  }
}

# plot the seed abundance trajectory through time (200 years)
library(ggplot2)
library(reshape2)
plot(seed.abundances[1,])
data.melt.arca <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.arca <- as.data.frame(spline(data.melt.arca$time, data.melt.arca$abundance))

arca.equ.mean <- mean(seed.abundances[,200])
# save as rdata
abundances.arca.novar <- seed.abundances
save(abundances.arca.novar, file="abundance.arca.novar.rdata")
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
for (r in 1:runs){ # then make runs=4500 to systematically go through each posterior 
   p <- r
  for (t in 1:(time-1)) {
  #  p <- sample(seq(1, 4500),1) 
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

abundances.hygl.novar <- seed.abundances
save(abundances.hygl.novar, file="abundance.hygl.novar.rdata")
remove(seed.abundances)

# PLDE -------------------------------------
seed.abundances <- matrix(data=NA,nrow=runs,ncol=time)
seed.abundances[,1] <- 2

germ=germination$plde
surv=survival$plde
lambdas=plde$lambda
intras=plde$alpha_intra

# lottery model
for (r in 1:runs){ # then make runs=4500 to systematically go through each posterior 
   p <- r
  for (t in 1:(time-1)) {
 #   p <- sample(seq(1, 4500),1) 
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

abundances.plde.novar <- seed.abundances
save(abundances.plde.novar, file="abundance.plde.novar.rdata")
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
for (r in 1:runs){ # then make runs=4500 to systematically go through each posterior 
   p <- r
  for (t in 1:(time-1)) {
   # p <- sample(seq(1, 4500),1) 
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

data.melt.poca <- melt(seed.abundances, varnames = c("run", "time"), value.name = "abundance")
spline.poca <- as.data.frame(spline(data.melt.poca$time, data.melt.poca$abundance))

poca.equ.mean <- mean(seed.abundances[,200])

abundances.poca.novar <- seed.abundances
save(abundances.poca.novar, file="abundance.poca.novar.rdata")
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
  for (r in 1:runs){ # then make runs=4500 to systematically go through each posterior 
   for (t in 1:(time-1)) {
     p <- r
    #p <- sample(seq(1, 4500),1) 
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

abundances.trcy.novar <- seed.abundances
save(abundances.trcy.novar, file="abundance.trcy.novar.rdata")
remove(seed.abundances)

# VERO ---------------------------------------
seed.abundances <- matrix(data=NA,nrow=runs,ncol=time)
seed.abundances[,1] <- 2

germ=germination$vero
surv=survival$vero
lambdas=vero$lambda
intras=vero$alpha_intra

# lottery model
for (r in 1:runs){ # then make runs=4500 to systematically go through each posterior 
  # p <- r
  for (t in 1:(time-1)) {
    p <- sample(seq(1, 4500),1) 
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

abundances.vero.novar <- seed.abundances
save(abundances.vero.novar, file="abundance.vero.novar.rdata")
remove(seed.abundances)

# MEDI ----------------------------------
seed.abundances <- matrix(data=NA,nrow=runs,ncol=time)
seed.abundances[,1] <- 2

germ=germination$medi
surv=survival$medi
lambdas=medi$lambda
intras=medi$alpha_intra

# lottery model
for (r in 1:runs){ # then make runs=4500 to systematically go through each posterior 
  p <- r
  for (t in 1:(time-1)) {
   # p <- sample(seq(1, 4500),1) 
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

abundances.medi.novar <- seed.abundances
save(abundances.medi.novar, file="abundance.medi.novar.rdata")
remove(seed.abundances)

# PEAI -------------------------------------------
seed.abundances <- matrix(data=NA,nrow=runs,ncol=time)
seed.abundances[,1] <- 2

germ=germination$peai
surv=survival$peai
lambdas=peai$lambda
intras=peai$alpha_intra

# lottery model
for (r in 1:runs){ # then make runs=4500 to systematically go through each posterior 
   p <- r
  for (t in 1:(time-1)) {
   # p <- sample(seq(1, 4500),1) 
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

abundances.peai.novar <- seed.abundances
save(abundances.peai.novar, file="abundance.peai.novar.rdata")
remove(seed.abundances)

### SUPP PLOT ###
pdf("Figures/supp.single.sp.equi_notime.pdf")
par(mfrow=c(3,3), mar=c(2,2,2,2), oma=c(4,4,1,1))
plot(spline.arca, main=expression(italic("A. calendula")))
mtext(side=2, line=3, "abundance")
plot(spline.medi, main=expression(italic("M. minima")))
plot(spline.peai, main=expression(italic("P. airoides")))
plot(spline.hygl, main=expression(italic("H. glutinosum")))
mtext(side=2, line=3, "abundance")
plot(spline.plde, main=expression(italic("P. debilis")))
plot(spline.poca, main=expression(italic("P. canescens")))
mtext(side=1, line=3, "time (years)")
plot(spline.trcy, main=expression(italic("T. cyanopetala")))
mtext(side=2, line=3, "abundance")
mtext(side=1, line=3, "time (years)")
plot(spline.vero, main=expression(italic("G. rosea")))
mtext(side=1, line=3, "time (years)")
dev.off()

#####
pdf("Figures/single.sp.t=200_notime.pdf")
par(mfrow=c(3,3), mar=c(2,2,2,2), oma=c(4,4,1,1))
plot(density(abundances.arca.novar[,200][which(abundances.arca.novar[,200]<HPDinterval(as.mcmc(abundances.arca.novar[,200]))[2])]), main=expression(italic("A. calendula")))
mtext(side=2, line=3, "density")
plot(density(abundances.medi.novar[,200][which(abundances.medi.novar[,200]<HPDinterval(as.mcmc(abundances.medi.novar[,200]))[2])]), main=expression(italic("M. minima")))
plot(density(abundances.peai.novar[,200][which(abundances.peai.novar[,200]<HPDinterval(as.mcmc(abundances.peai.novar[,200]))[2])]), main=expression(italic("P. airoides")))
plot(density(abundances.hygl.novar[,200][which(abundances.hygl.novar[,200]<HPDinterval(as.mcmc(abundances.hygl.novar[,200]))[2])]), main=expression(italic("H. glutinosum")))
mtext(side=2, line=3, "density")
plot(density(abundances.plde.novar[,200][which(abundances.plde.novar[,200]<HPDinterval(as.mcmc(abundances.plde.novar[,200]))[2])]), main=expression(italic("P. debilis")))
plot(density(abundances.poca.novar[,200][which(abundances.poca.novar[,200]<HPDinterval(as.mcmc(abundances.poca.novar[,200]))[2])]), main=expression(italic("P. canescens")))
plot(density(abundances.trcy.novar[,200][which(abundances.trcy.novar[,200]<HPDinterval(as.mcmc(abundances.trcy.novar[,200]))[2])]), 
     main=expression(italic("T. cyanopetala")))
mtext(side=1, line=3, "steady state population size")
mtext(side=2, line=3, "density")
plot(density(abundances.vero.novar[,200][which(abundances.vero.novar[,200]<HPDinterval(as.mcmc(abundances.vero.novar[,200]))[2])]), main=expression(italic("G. rosea")))
mtext(side=1, line=3, "steady state population size")
dev.off()

# save the t200 mean values 
t200.mean.vals <- c(hygl.equ.mean, plde.equ.mean, poca.equ.mean,
                    trcy.equ.mean, vero.equ.mean, arca.equ.mean, medi.equ.mean,
                    peai.equ.mean)
save(t200.mean.vals, file="t200.mean.values")


