# load in Stan model fits
# removing MOMO for now

sp_list <- c("dagl", "hygl", "plde", 
             "poca", "trcy", "vero", "arca", "medi", "peai")

source("Bayes.data.R")
source("Residents.to.equilibrium.no.time.var.R")

sp <- 8
library(coda)
library(rstan)

source("BH.equalizing.stabilizing.R")
#dagl[["germination"]] <- germ$p.1
hygl[["germination"]] <- germ$p.2
plde[["germination"]] <- germ$p.3
poca[["germination"]] <- germ$p.4
trcy[["germination"]] <- germ$p.5
vero[["germination"]] <- germ$p.6
arca[["germination"]] <- germ$p.7
medi[["germination"]] <- germ$p.8
peai[["germination"]] <- germ$p.9


#dagl[["survival"]] <-surv$p.1
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


# -----------------------------------------------------------------------------------------

runs <- 4500
stabilizing <- matrix(NA, nrow = runs, ncol = total_combos)
equalizing <- matrix(NA, nrow = runs, ncol = total_combos)
coexists <- matrix(NA, nrow = runs, ncol = total_combos)
hygl.equilibrum <- abundances.hygl.novar[,200]
peai.equilibrum <- abundances.peai.novar[,200]
invader.abund <- 1

# create all of our LDGR matrices
nolot.medi.into.hygl <- matrix(NA, nrow=runs)
nolot.peai.into.medi <- matrix(NA, nrow=runs)

# Run through all possible combos
for (r in 1:runs){
  #p <- sample(seq(1, 4500),1)
  p <- r
  Nj.hygl <- hygl.equilibrum[r]
  germ_j_hygl <- germination$hygl[p]
  lived.hygl<-germ_j_hygl*Nj.hygl
  
  Nj.peai <- peai.equilibrum[r]
  germ_j_peai <- germination$peai[p]
  lived.peai<-germ_j_peai*Nj.peai
  
  # invade MEDI
  peai_tp1 <- (1-germination$peai[p])*survival$peai[p]*invader.abund + 
    invader.abund*germination$peai[p]*peai$lambda[p]/(1+peai$alpha_hygl[p]*lived.hygl)
  # calculate LDGR of medi
  nolot.peai.into.hygl[r] <- log(peai_tp1/invader.abund)
  
  hygl_tp1 <- (1-germination$hygl[p])*survival$hygl[p]*invader.abund + 
    invader.abund*germination$hygl[p]*hygl$lambda[p]/(1+hygl$alpha_peai[p]*lived.peai)
  # calculate LDGR of medi
  nolot.hygl.into.peai[r] <- log(hygl_tp1/invader.abund)
  
 # for (x in 1:total_combos) {
    # assign s1, s2, alpha_s1, alpha_s2
    s1 <- s1_all[[4]]
    s2 <- s2_all[[4]]
    alpha_s1 <- alpha_s1_all[4]
    alpha_s2 <- alpha_s2_all[4]
    
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
    
    stabilizing[r,4] <- do.stabilization(alpha.ii=intras.i[p], alpha.jj=intras.j[p], 
                                         alpha.ij=inter.ij[p], alpha.ji=inter.ji[p])
    
    equalizing[r,4] <- do.kj.over.ki(lambda.i=lambdas.i[p], lambda.j=lambdas.j[p], 
                                     g.i=germ.i[p], g.j=germ.j[p], s.i=surv.i[p], 
                                     s.j=surv.j[p], alpha.ii=intras.i[p], alpha.jj=intras.j[p], 
                                     alpha.ij=inter.ij[p], alpha.ji=inter.ji[p])
    if(equalizing[r,4] < 1) {
      equalizing[r,4] <- 1/equalizing[r,4]
    }
    coexists[r,4] <- do.coexistence(stabilizing=stabilizing[r,4], kj.over.ki=equalizing[r,4]) # negs in both of these
  }

#}
corrected.stabilizing <- ifelse(stabilizing< 0, 0, stabilizing)

# plots ####
# hygl and medi
sum(nolot.hygl.into.peai[which(nolot.hygl.into.peai>0)])
sum(nolot.peai.into.hygl[which(nolot.peai.into.hygl>0)])

auc <- ecdf(nolot.hygl.into.peai)
1-auc(0)
auc <- ecdf(nolot.peai.into.hygl)
1-auc(0)

dat <- cbind(nolot.hygl.into.peai,nolot.peai.into.hygl)
coex <- which(dat[,1]>0&dat[,2]>0) #43 coexist (without variance allowed at all)
sub.stab <- corrected.stabilizing[,4]
sub.stab <- sub.stab[-coex]
sub.ex <- equalizing[,4]
sub.ex <- sub.ex[-coex]

plot(corrected.stabilizing[,4][coex], log(equalizing[,4][coex]), pch=16, col=alpha("orange", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(sub.stab, log(sub.ex), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
#points(mean.corrected.stabilizing[3], log(mean.equalizing[3]), pch=19, col="black")
#abline(v=HPDinterval(as.mcmc(corrected.stabilizing[,3]))[2])
#abline(v=HPDinterval(as.mcmc(corrected.stabilizing[,3]))[1])
#abline(h=(HPDinterval(as.mcmc(log(equalizing[,3])))[2]))
#abline(h=(HPDinterval(as.mcmc(log(equalizing[,3])))[1]))
lines(xseq, log(yseq), col="grey", lwd=2)
#polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
mtext(pairs[3], side=3, cex = 0.7)

# calculate probability of coexistence 
coex.thresh <- vector()
stab.val<-vector()
for (i in 1:4500){
  coex.thresh[i] <- 1/(1-corrected.stabilizing[,3][i])
  if (log(equalizing[,3][i]) < log(coex.thresh[i])){
    coex.thresh[i] <- 1
  } else {
    coex.thresh[i] <- 0
  }
}

sum(coex.thresh)/4500



############ adding back in time variance for invasion approach ######
source("Residents.to.equilibrium.R")
sp <- 8
library(coda)
library(rstan)

source("BH.equalizing.stabilizing.R")
#dagl[["germination"]] <- germ$p.1
hygl[["germination"]] <- germ$p.2
plde[["germination"]] <- germ$p.3
poca[["germination"]] <- germ$p.4
trcy[["germination"]] <- germ$p.5
vero[["germination"]] <- germ$p.6
arca[["germination"]] <- germ$p.7
medi[["germination"]] <- germ$p.8
peai[["germination"]] <- germ$p.9


#dagl[["survival"]] <-surv$p.1
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


# -----------------------------------------------------------------------------------------

runs <- 4500
stabilizing <- matrix(NA, nrow = runs, ncol = total_combos)
equalizing <- matrix(NA, nrow = runs, ncol = total_combos)
coexists <- matrix(NA, nrow = runs, ncol = total_combos)
hygl.equilibrum <- abundances.hygl[,200]
medi.equilibrum <- abundances.medi[,200]
invader.abund <- 1

# create all of our LDGR matrices
nolot.medi.into.hygl <- matrix(NA, nrow=runs)
nolot.hygl.into.medi <- matrix(NA, nrow=runs)

# Run through all possible combos
for (r in 1:runs){
  p <- r
  Nj.hygl <- hygl.equilibrum[r]
  germ_j_hygl <- germination$hygl[p]
  lived.hygl<-germ_j_hygl*Nj.hygl
  
  Nj.medi <- medi.equilibrum[r]
  germ_j_medi <- germination$medi[p]
  lived.medi<-germ_j_medi*Nj.medi
  
  # invade MEDI
  medi_tp1 <- (1-germination$medi[p])*survival$medi[p]*invader.abund + 
    invader.abund*germination$medi[p]*medi$lambda[p]/(1+medi$alpha_hygl[p]*lived.hygl)
  # calculate LDGR of medi
  nolot.medi.into.hygl[r] <- log(medi_tp1/invader.abund)
  
  hygl_tp1 <- (1-germination$hygl[p])*survival$hygl[p]*invader.abund + 
    invader.abund*germination$hygl[p]*hygl$lambda[p]/(1+hygl$alpha_medi[p]*lived.medi)
  # calculate LDGR of medi
  nolot.hygl.into.medi[r] <- log(hygl_tp1/invader.abund)
  
  # for (x in 1:total_combos) {
  # assign s1, s2, alpha_s1, alpha_s2
  s1 <- s1_all[[3]]
  s2 <- s2_all[[3]]
  alpha_s1 <- alpha_s1_all[3]
  alpha_s2 <- alpha_s2_all[3]
  
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
  
  stabilizing[r,3] <- do.stabilization(alpha.ii=intras.i[p], alpha.jj=intras.j[p], 
                                       alpha.ij=inter.ij[p], alpha.ji=inter.ji[p])
  
  equalizing[r,3] <- do.kj.over.ki(lambda.i=lambdas.i[p], lambda.j=lambdas.j[p], 
                                   g.i=germ.i[p], g.j=germ.j[p], s.i=surv.i[p], 
                                   s.j=surv.j[p], alpha.ii=intras.i[p], alpha.jj=intras.j[p], 
                                   alpha.ij=inter.ij[p], alpha.ji=inter.ji[p])
  if(equalizing[r,3] < 1) {
    equalizing[r,3] <- 1/equalizing[r,3]
  }
  coexists[r,3] <- do.coexistence(stabilizing=stabilizing[r,3], kj.over.ki=equalizing[r,3]) # negs in both of these
}

#}
corrected.stabilizing <- ifelse(stabilizing< 0, 0, stabilizing)

# plots ####
# hygl and medi
HPDinterval(as.mcmc(nolot.medi.into.hygl))
HPDinterval(as.mcmc(nolot.hygl.into.medi))

dat <- cbind(nolot.hygl.into.medi,nolot.medi.into.hygl)
coex <- which(dat[,1]>0&dat[,2]>0) #48 coex

sub.stab <- corrected.stabilizing[,3] 
sub.stab <- sub.stab[-coex]
sub.ex <- equalizing[,3]
sub.ex <- sub.ex[-coex]

plot(corrected.stabilizing[,3][coex], log(equalizing[,3][coex]), pch=16, col=alpha("orange", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(sub.stab, log(sub.ex), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
#points(mean.corrected.stabilizing[3], log(mean.equalizing[3]), pch=19, col="black")
abline(v=HPDinterval(as.mcmc(corrected.stabilizing[,3]))[2])
abline(v=HPDinterval(as.mcmc(corrected.stabilizing[,3]))[1])
abline(h=(HPDinterval(as.mcmc(log(equalizing[,3])))[2]))
abline(h=(HPDinterval(as.mcmc(log(equalizing[,3])))[1]))
lines(xseq, log(yseq), col="grey", lwd=2)
#polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
mtext(pairs[3], side=3, cex = 0.7)

