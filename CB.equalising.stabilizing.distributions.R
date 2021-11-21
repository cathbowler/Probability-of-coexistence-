############################################################################################################
# This code calculate pair-wise low-density growth rate calculations using the full posterior probabilities
# for each parameter
############################################################################################################


sp <- 8
library(coda)
library(rstan)
source("Bayes.data.R")
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
  #p <- sample(seq(1, 4500),1) 
p <- r
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
plot(corrected.stabilizing[,3], log(equalizing[,3]), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(mean.corrected.stabilizing[3], log(mean.equalizing[3]), pch=19, col="black")
abline(v=HPDinterval(as.mcmc(corrected.stabilizing[,3]))[2])
abline(v=HPDinterval(as.mcmc(corrected.stabilizing[,3]))[1])
abline(h=log(HPDinterval(as.mcmc(equalizing[,3]))[2]))
abline(h=log(HPDinterval(as.mcmc(equalizing[,3]))[1]))
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
mtext(pairs[5], side=3, cex = 0.7)

# -----------------------------------------------------------------------------------------
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

# pdf("Figures/EqualizingStabilizingDistributions.pdf")
# par(mfrow=c(5,3))
# for (i in 1:14) {
#   par(mar=c(4,4,2,2))
#   plot(corrected.stabilizing[,i], log(equalizing[,i]), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
#   points(mean.corrected.stabilizing[i], log(mean.equalizing[i]), pch=19, col="black")
#   lines(xseq, log(yseq), col="grey", lwd=2)
#   #arrows(mean.corrected.stabilizing[i], mean.equalizing[i]+eq.int[i], mean.corrected.stabilizing[i], mean.equalizing[i]-eq.int[i], length=0.05, angle=90, code=3)
#   #arrows(mean.corrected.stabilizing[i]+stab.int[i], mean.equalizing[i], mean.corrected.stabilizing[i]-stab.int[i], mean.equalizing[i], length=0.05, angle=90, code=3)
#   polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
#   mtext(pairs[i], side=3, cex = 0.7)
#   
# }


# matrix plot ####
pdf("Figures/EqualizingStabilizingDistributions_revisions_order.pdf")
par(mar=c(4,4,2,2), mfrow=c(5,3))

plot(corrected.stabilizing[,5], log(equalizing[,5]), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(mean.corrected.stabilizing[5], log(mean.equalizing[5]), pch=19, col="black")
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
mtext(pairs[5], side=3, cex = 0.7)

# calculate the probability of coexistence 
coex.thresh <- vector()
stab.val<-vector()
for (i in 1:4500){
  coex.thresh[i] <- 1/(1-corrected.stabilizing[,5][i])
  if (log(equalizing[,5][i]) < log(coex.thresh[i])){
    coex.thresh[i] <- 1
  } else {
    coex.thresh[i] <- 0}}
tx <- sum(coex.thresh)/4500*100
text(.9, 11, round(tx,2), cex=0.8)

plot(corrected.stabilizing[,8], log(equalizing[,8]), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(mean.corrected.stabilizing[8], log(mean.equalizing[8]), pch=19, col="black")
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
mtext(pairs[8], side=3, cex = 0.7)

coex.thresh <- vector()
stab.val<-vector()
for (i in 1:4500){
  coex.thresh[i] <- 1/(1-corrected.stabilizing[,8][i])
  if (log(equalizing[,8][i]) < log(coex.thresh[i])){
    coex.thresh[i] <- 1
  } else {
    coex.thresh[i] <- 0}}
tx <- sum(coex.thresh)/4500*100
text(.9, 11, round(tx,2), cex=0.8)

plot(corrected.stabilizing[,4], log(equalizing[,4]), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(mean.corrected.stabilizing[4], log(mean.equalizing[4]), pch=19, col="black")
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
mtext(pairs[4], side=3, cex = 0.7)

coex.thresh <- vector()
stab.val<-vector()
for (i in 1:4500){
  coex.thresh[i] <- 1/(1-corrected.stabilizing[,4][i])
  if (log(equalizing[,4][i]) < log(coex.thresh[i])){
    coex.thresh[i] <- 1
  } else {
    coex.thresh[i] <- 0}}
tx <- sum(coex.thresh)/4500*100
text(.9, 11, round(tx,2), cex=0.8)

plot(corrected.stabilizing[,9], log(equalizing[,9]), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(mean.corrected.stabilizing[9], log(mean.equalizing[9]), pch=19, col="black")
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
mtext(pairs[9], side=3, cex = 0.7)

coex.thresh <- vector()
stab.val<-vector()
for (i in 1:4500){
  coex.thresh[i] <- 1/(1-corrected.stabilizing[,9][i])
  if (log(equalizing[,9][i]) < log(coex.thresh[i])){
    coex.thresh[i] <- 1
  } else {
    coex.thresh[i] <- 0}}
tx <- sum(coex.thresh)/4500*100
text(.9, 11, round(tx,2), cex=0.8)

plot(corrected.stabilizing[,1], log(equalizing[,1]), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(mean.corrected.stabilizing[1], log(mean.equalizing[1]), pch=19, col="black")
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
mtext(pairs[1], side=3, cex = 0.7)

coex.thresh <- vector()
stab.val<-vector()
for (i in 1:4500){
  coex.thresh[i] <- 1/(1-corrected.stabilizing[,1][i])
  if (log(equalizing[,1][i]) < log(coex.thresh[i])){
    coex.thresh[i] <- 1
  } else {
    coex.thresh[i] <- 0}}
tx <- sum(coex.thresh)/4500*100
text(.9, 11, round(tx,2), cex=0.8)

plot(corrected.stabilizing[,3], log(equalizing[,3]), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(mean.corrected.stabilizing[3], log(mean.equalizing[3]), pch=19, col="black")
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
mtext(pairs[3], side=3, cex = 0.7)

coex.thresh <- vector()
stab.val<-vector()
for (i in 1:4500){
  coex.thresh[i] <- 1/(1-corrected.stabilizing[,3][i])
  if (log(equalizing[,3][i]) < log(coex.thresh[i])){
    coex.thresh[i] <- 1
  } else {
    coex.thresh[i] <- 0}}
tx <- sum(coex.thresh)/4500*100
text(.9, 11, round(tx,2), cex=0.8)

plot(corrected.stabilizing[,12], log(equalizing[,12]), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(mean.corrected.stabilizing[12], log(mean.equalizing[12]), pch=19, col="black")
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
mtext(pairs[12], side=3, cex = 0.7)

coex.thresh <- vector()
stab.val<-vector()
for (i in 1:4500){
  coex.thresh[i] <- 1/(1-corrected.stabilizing[,12][i])
  if (log(equalizing[,12][i]) < log(coex.thresh[i])){
    coex.thresh[i] <- 1
  } else {
    coex.thresh[i] <- 0}}
tx <- sum(coex.thresh)/4500*100
text(.9, 11, round(tx,2), cex=0.8)

plot(corrected.stabilizing[,2], log(equalizing[,2]), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(mean.corrected.stabilizing[2], log(mean.equalizing[2]), pch=19, col="black")
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
mtext(pairs[2], side=3, cex = 0.7)

coex.thresh <- vector()
stab.val<-vector()
for (i in 1:4500){
  coex.thresh[i] <- 1/(1-corrected.stabilizing[,2][i])
  if (log(equalizing[,2][i]) < log(coex.thresh[i])){
    coex.thresh[i] <- 1
  } else {
    coex.thresh[i] <- 0}}
tx <- sum(coex.thresh)/4500*100
text(.9, 11, round(tx,2), cex=0.8)

plot(corrected.stabilizing[,6], log(equalizing[,6]), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(mean.corrected.stabilizing[6], log(mean.equalizing[6]), pch=19, col="black")
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
mtext(pairs[6], side=3, cex = 0.7)

coex.thresh <- vector()
stab.val<-vector()
for (i in 1:4500){
  coex.thresh[i] <- 1/(1-corrected.stabilizing[,6][i])
  corrected.equalizing <- ifelse(equalizing< 0, 100, equalizing) # just making 100 to ignore it coz coex not possible anyway
  if (log(corrected.equalizing[,6][i]) < log(coex.thresh[i])){
    coex.thresh[i] <- 1
  } else {
    coex.thresh[i] <- 0}}
tx <- sum(coex.thresh)/4500*100
text(.9, 11, round(tx,2), cex=0.8)

plot(corrected.stabilizing[,7], log(equalizing[,7]), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(mean.corrected.stabilizing[7], log(mean.equalizing[7]), pch=19, col="black")
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
mtext(pairs[7], side=3, cex = 0.7)

coex.thresh <- vector()
stab.val<-vector()
for (i in 1:4500){
  coex.thresh[i] <- 1/(1-corrected.stabilizing[,7][i])
  if (log(equalizing[,7][i]) < log(coex.thresh[i])){
    coex.thresh[i] <- 1
  } else {
    coex.thresh[i] <- 0}}
tx <- sum(coex.thresh)/4500*100
text(.9, 11, round(tx,2), cex=0.8)

plot(corrected.stabilizing[,10], log(equalizing[,10]), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(mean.corrected.stabilizing[10], log(mean.equalizing[10]), pch=19, col="black")
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
mtext(pairs[10], side=3, cex = 0.7)

coex.thresh <- vector()
stab.val<-vector()
for (i in 1:4500){
  coex.thresh[i] <- 1/(1-corrected.stabilizing[,10][i])
  corrected.equalizing <- ifelse(equalizing< 0, 100, equalizing) # just making 100 to ignore it coz coex not possible anyway
  if (log(corrected.equalizing[,10][i]) < log(coex.thresh[i])){
    coex.thresh[i] <- 1
  } else {
    coex.thresh[i] <- 0}}
tx <- sum(coex.thresh)/4500*100
text(.9, 11, round(tx,2), cex=0.8)

plot(corrected.stabilizing[,11], log(equalizing[,11]), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(mean.corrected.stabilizing[11], log(mean.equalizing[11]), pch=19, col="black")
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
mtext(pairs[11], side=3, cex = 0.7)

coex.thresh <- vector()
stab.val<-vector()
for (i in 1:4500){
  coex.thresh[i] <- 1/(1-corrected.stabilizing[,11][i])
  corrected.equalizing <- ifelse(equalizing< 0, 100, equalizing) # just making 100 to ignore it coz coex not possible anyway
  if (log(corrected.equalizing[,11][i]) < log(coex.thresh[i])){
    coex.thresh[i] <- 1
  } else {
    coex.thresh[i] <- 0}}
tx <- sum(coex.thresh)/4500*100
text(.9, 11, round(tx,2), cex=0.8)

plot(corrected.stabilizing[,13], log(equalizing[,13]), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(mean.corrected.stabilizing[13], log(mean.equalizing[13]), pch=19, col="black")
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
mtext(pairs[13], side=3, cex = 0.7)

coex.thresh <- vector()
stab.val<-vector()
for (i in 1:4500){
  coex.thresh[i] <- 1/(1-corrected.stabilizing[,13][i])
  corrected.equalizing <- ifelse(equalizing< 0, 100, equalizing) # just making 100 to ignore it coz coex not possible anyway
  if (log(corrected.equalizing[,13][i]) < log(coex.thresh[i])){
    coex.thresh[i] <- 1
  } else {
    coex.thresh[i] <- 0}}
tx <- sum(coex.thresh)/4500*100
text(.9, 11, round(tx,2), cex=0.8)

plot(corrected.stabilizing[,14], log(equalizing[,14]), pch=16, col=alpha("skyblue", 0.5), xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(mean.corrected.stabilizing[14], log(mean.equalizing[14]), pch=19, col="black")
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
mtext(pairs[14], side=3, cex = 0.7)

coex.thresh <- vector()
stab.val<-vector()
for (i in 1:4500){
  coex.thresh[i] <- 1/(1-corrected.stabilizing[,14][i])
  corrected.equalizing <- ifelse(equalizing< 0, 100, equalizing) # just making 100 to ignore it coz coex not possible anyway
  if (log(corrected.equalizing[,14][i]) < log(coex.thresh[i])){
    coex.thresh[i] <- 1
  } else {
    coex.thresh[i] <- 0}}
tx <- sum(coex.thresh)/4500*100
text(.9, 11, round(tx,2), cex=0.8)

dev.off()

