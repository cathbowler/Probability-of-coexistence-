sp <- 9
library(coda)
library(rstan)
source("BH.equalizing.stabilizing.R")
# removing MOMO for now

#sp_list <- c(dagl, gite, hygl, plde, 
#             poca, trcy, vero, arca, medi, momo, peai)

#setwd("/Users/lash1937/Dropbox/Shared/Bayesian data - Cath Collaboration/Updated Coexistence code")
#setwd("/Users/catherinebowler/Dropbox/Bayesian data - Cath Collaboration/Coexistence code")
load("BH stan models/Arca_posteriors.rdata")
arca <- extract(PrelimFit)
remove(PrelimFit)

load("BH stan models/Dagl_posteriors.rdata")
dagl <- extract(PrelimFit)
remove(PrelimFit)

load("BH stan models//Hygl_posteriors.rdata")
hygl <- extract(PrelimFit)
remove(PrelimFit)

load("BH stan models//Medi_posteriors.rdata")
medi <- extract(PrelimFit)
remove(PrelimFit)

#load("Momo/Momo_posteriors.rdata")
#momo <- extract(PrelimFit)
#remove(PrelimFit)

load("BH stan models//Peai_posteriors.rdata")
peai <- extract(PrelimFit)
remove(PrelimFit)

load("BH stan models//Plde_posteriors.rdata")
plde <- extract(PrelimFit)
remove(PrelimFit)

load("BH stan models//Poca_posteriors.rdata")
poca <- extract(PrelimFit)
remove(PrelimFit)

load("BH stan models//Trcy_posteriors.rdata")
trcy <- extract(PrelimFit)
remove(PrelimFit)

load("BH stan models//Vero_posteriors.rdata")
vero <- extract(PrelimFit)
remove(PrelimFit)

load("SurvivalAndGermination/Germination.rdata")
germ <- extract(PrelimFit)
germ<-as.data.frame(germ)
remove(PrelimFit)

load("SurvivalAndGermination/Survival.rdata")
surv <- extract(PrelimFit)
surv<-as.data.frame(surv)
remove(PrelimFit)

germination <- data.frame(dagl=germ$p.1, hygl=germ$p.2, plde=germ$p.3, poca=germ$p.4, trcy=germ$p.5, 
                          vero=germ$p.6, arca=germ$p.7, medi=germ$p.8, peai=germ$p.9)

dagl[["germination"]] <- germ$p.1
hygl[["germination"]] <- germ$p.2
plde[["germination"]] <- germ$p.3
poca[["germination"]] <- germ$p.4
trcy[["germination"]] <- germ$p.5
vero[["germination"]] <- germ$p.6
arca[["germination"]] <- germ$p.7
medi[["germination"]] <- germ$p.8
peai[["germination"]] <- germ$p.9

survival <- data.frame(dagl=surv$p.1, hygl=surv$p.2, plde=surv$p.3, poca=surv$p.4, trcy=surv$p.5, 
                       vero=surv$p.6, arca=surv$p.7, medi=surv$p.8, peai=surv$p.9)

dagl[["survival"]] <-surv$p.1
hygl[["survival"]] <- surv$p.2
plde[["survival"]] <-surv$p.3
poca[["survival"]] <- surv$p.4
trcy[["survival"]] <- surv$p.5
vero[["survival"]] <- surv$p.6
arca[["survival"]] <- surv$p.7
medi[["survival"]] <-surv$p.8
peai[["survival"]] <-surv$p.9

#-----------------------------------------------------------------------
# calculate stabilizing

do.stabilization <- function(alpha.ii, alpha.jj, alpha.ij, alpha.ji) {
  rho <- sqrt((alpha.ij*alpha.ji)/(alpha.jj*alpha.ii))
  stabilizing <- 1- rho
  return(stabilizing)
}

# calculate average fitness differences
do.kj.over.ki <- function(lambda.i, lambda.j, g.i, g.j, s.i, s.j, alpha.ii, alpha.jj, alpha.ij, alpha.ji) {
  n.i <- lambda.i*g.i/(1-(1-g.i)*s.i)
  n.j <- lambda.j*g.j/(1-(1-g.j)*s.j)
  kj.over.ki <- ((n.j-1)/(n.i-1))*sqrt((alpha.ij*alpha.ii)/(alpha.jj*alpha.ji))
  #kj.over.ki <- sqrt((alpha.ij*alpha.ii)/(alpha.jj*alpha.ji))
  return(kj.over.ki)
}

# calculate coexistence
do.coexistence <- function(stabilizing, kj.over.ki) {
  rho <- 1 - stabilizing
  ki.over.kj <- 1/kj.over.ki
  if (rho < ki.over.kj) {
    coexists <- 1
  } else {
    coexists <- 0
  }
  return(coexists)
}


runs <- 4500
# -----------------------------------------------------------------------------------------
# Make things generic

sp_list <- c("dagl", "hygl", "plde", 
             "poca", "trcy", "vero", "arca", "medi", "peai")

total_combos <- 14

s1_all <- list(arca, arca, hygl, hygl, hygl, hygl, medi, peai, peai, peai, peai, plde, poca, poca)
s2_all <- list(peai, poca, medi, peai, plde, vero, peai, plde, poca, trcy, vero, poca, trcy, vero)
alpha_s1_all <- c("alpha_arca", "alpha_arca", "alpha_hygl", "alpha_hygl", "alpha_hygl", "alpha_hygl", "alpha_medi", 
                  "alpha_peai", "alpha_peai", "alpha_peai", "alpha_peai", "alpha_plde", "alpha_poca", "alpha_poca")
alpha_s2_all <- c("alpha_peai", "alpha_poca", "alpha_medi", "alpha_peai", "alpha_plde", "alpha_vero", "alpha_peai",
                  "alpha_plde", "alpha_poca", "alpha_trcy", "alpha_vero", "alpha_poca", "alpha_trcy", "alpha_vero")


stabilizing <- matrix(NA, nrow = runs, ncol = total_combos)
equalizing <- matrix(NA, nrow = runs, ncol = total_combos)
coexists <- matrix(NA, nrow = runs, ncol = total_combos)
#p=r # set x to number of runs which equals number of posterior values 
# -----------------------------------------------------------------------------------------
# Run through all possible combos
for (r in 1:runs){
for (x in 1:total_combos) {
 # for (p in 1:runs){
  p <- sample(seq(1, 4500),1) 
#  y2 <- sample(seq(1, 4500),1) # for .j
#  y3 <- sample(seq(1, 4500),1) # for germ i
#  y4 <- sample(seq(1, 4500),1) # for germ j
#  y5 <- sample(seq(1, 4500),1) # for surv i
#  y6 <- sample(seq(1, 4500),1) # for surv j
  
  # assign s1, s2, alpha_s1, alpha_s2
  s1 <- s1_all[[x]]
  s2 <- s2_all[[x]]
  alpha_s1 <- alpha_s1_all[x]
  alpha_s2 <- alpha_s2_all[x]
  
  germ.j=s1$germination
  surv.j=s1$survival
  lambdas.j=s1$lambda
  intras.j=s1$alpha_intra
  
  germ.i=s2$germination
  surv.i=s2$survival
  lambdas.i=s2$lambda
  intras.i=s2$alpha_intra
  
  inter.ij <- s2[[alpha_s1]]
  inter.ji <- s1[[alpha_s2]]
  
  stabilizing[r,x] <- do.stabilization(alpha.ii=intras.i[p], alpha.jj=intras.j[p], 
                                          alpha.ij=inter.ij[p], alpha.ji=inter.ji[p])
  
  equalizing[r,x] <- do.kj.over.ki(lambda.i=lambdas.i[p], lambda.j=lambdas.j[p], 
                                      g.i=germ.i[p], g.j=germ.j[p], s.i=surv.i[p], 
                                      s.j=surv.j[p], alpha.ii=intras.i[p], alpha.jj=intras.j[p], 
                                      alpha.ij=inter.ij[p], alpha.ji=inter.ji[p])
  if(equalizing[r,x] < 1) {
    equalizing[r,x] <- 1/equalizing[r,x]
  }
  coexists[r,x] <- do.coexistence(stabilizing=stabilizing[r,x], kj.over.ki=equalizing[r,x]) # negs in both of these
}
  }
#}
corrected.stabilizing <- ifelse(stabilizing< 0, 0, stabilizing)

# calculate HPD invtervals - something like this
stab.int <- (as.numeric(HPDinterval(mcmc(corrected.stabilizing))[,2]) - as.numeric(HPDinterval(mcmc(corrected.stabilizing))[,1]))/2
eq.int <- (as.numeric(HPDinterval(mcmc(equalizing))[,2]) - as.numeric(HPDinterval(mcmc(equalizing))[,1]))/2

# -----------------------------------------------------------------------------------------
# spit to a plot of the ones you want in text and then the others in supps 
do.transparent <- function(InCol, NewAlpha){
  rgb_vals <- col2rgb(InCol)[,1]
  OutCol <- rgb(red = rgb_vals[1], green = rgb_vals[2], blue = rgb_vals[3],
                alpha = NewAlpha, maxColorValue = 255)
  return(OutCol)
}

pairs <- c("A. calendula & P. airoides", "A.calendula & P.canescens", "H.glutinosum & M.minima", 
           "H. glutinosum & P. airoides", "H. glutinosum & P.deibilis", "H. glutinosum & G. rosea", 
           "M.minima & P. airoides", "P.airoides & P. debilis", "P. airoides & P. canescens", 
           "P. airoides & T. cyanopetala", "P. airoides & G. rosea", "P. debilis & P. canescens", 
           "P.canescens & T. cyanopetala", "P. canescens & G. rosea")

pdf("Figures/EqualizingStabilizingDistributions.pdf")
par(mfrow=c(5,3))
for (i in 1:14) {
  par(mar=c(4,4,2,2))
  plot(corrected.stabilizing[,i], log(equalizing[,i]), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
  points(mean.corrected.stabilizing[i], log(mean.equalizing[i]), pch=19, col="black")
  lines(xseq, log(yseq), col="grey", lwd=2)
  #arrows(mean.corrected.stabilizing[i], mean.equalizing[i]+eq.int[i], mean.corrected.stabilizing[i], mean.equalizing[i]-eq.int[i], length=0.05, angle=90, code=3)
  #arrows(mean.corrected.stabilizing[i]+stab.int[i], mean.equalizing[i], mean.corrected.stabilizing[i]-stab.int[i], mean.equalizing[i], length=0.05, angle=90, code=3)
  polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
  mtext(pairs[i], side=3, cex = 0.7)
  
}
dev.off()

par(mar=c(4,4,2,2))
plot(corrected.stabilizing[,6], log(equalizing[,6]), pch=16, col="skyblue", xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(mean.corrected.stabilizing[6], log(mean.equalizing[6]), pch=19, col="black")
lines(xseq, log(yseq), col="grey", lwd=2)
abline(v=HPDinterval(mcmc(corrected.stabilizing[,6])), h=log(HPDinterval(mcmc(equalizing[,6]))))
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)


# matrix plot ####

par(mar=c(4,4,2,2))
plot(corrected.stabilizing[,1], log(equalizing[,1]), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(mean.corrected.stabilizing[1], log(mean.equalizing[1]), pch=19, col="black")
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
mtext(pairs[1], side=3, cex = 0.7)

plot(corrected.stabilizing[,1], log(equalizing[,1]), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(mean.corrected.stabilizing[1], log(mean.equalizing[1]), pch=19, col="black")
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
mtext(pairs[1], side=3, cex = 0.7)
