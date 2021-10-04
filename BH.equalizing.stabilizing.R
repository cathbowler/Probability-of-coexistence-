# Calculate pair-wise low-density growth rate calculations
# load in Stan model fits
sp <- 9
library(coda)
library(rstan)
source("Bayes.data.R")

#dagl[["germination"]] <- germ$p.1
hygl[["germination"]] <- germ$p.2
plde[["germination"]] <- germ$p.3
poca[["germination"]] <- germ$p.4
trcy[["germination"]] <- germ$p.5
vero[["germination"]] <- germ$p.6
arca[["germination"]] <- germ$p.7
medi[["germination"]] <- germ$p.8
peai[["germination"]] <- germ$p.9

survival <- data.frame(hygl=surv$p.2, plde=surv$p.3, poca=surv$p.4, trcy=surv$p.5, 
                       vero=surv$p.6, arca=surv$p.7, medi=surv$p.8, peai=surv$p.9)

#dagl[["survival"]] <-surv$p.1
hygl[["survival"]] <- surv$p.2
plde[["survival"]] <-surv$p.3
poca[["survival"]] <- surv$p.4
trcy[["survival"]] <- surv$p.5
vero[["survival"]] <- surv$p.6
arca[["survival"]] <- surv$p.7
medi[["survival"]] <-surv$p.8
peai[["survival"]] <-surv$p.9

# -----------------------------------------------------------------------------------------
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


runs <- 10000
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

mean.stabilizing <- rep(NA, total_combos)
mean.equalizing <- rep(NA, total_combos)
coexists <- rep(NA, total_combos)
# -----------------------------------------------------------------------------------------
# Run through all possible combos

for (x in 1:total_combos) {
  
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
  
  mean.stabilizing[x] <- do.stabilization(alpha.ii=mean(intras.i), alpha.jj=mean(intras.j), 
                                                 alpha.ij=mean(inter.ij), alpha.ji=mean(inter.ji))
  
  mean.equalizing[x] <- do.kj.over.ki(lambda.i=mean(lambdas.i), lambda.j=mean(lambdas.j), 
                                             g.i=mean(germ.i), g.j=mean(germ.j), s.i=mean(surv.i), 
                                             s.j=mean(surv.j), alpha.ii=mean(intras.i), alpha.jj=mean(intras.j), 
                                             alpha.ij=mean(inter.ij), alpha.ji=mean(inter.ji))
  if(mean.equalizing[x] < 1) {
    mean.equalizing[x] <- 1/mean.equalizing[x]
  }
  
  coexists[x] <- do.coexistence(stabilizing=mean.stabilizing[x], kj.over.ki=mean.equalizing[x])
}

mean.corrected.stabilizing <- ifelse(mean.stabilizing < 0, 0, mean.stabilizing)
# -----------------------------------------------------------------------------------------
# make coexistence line
xseq <- seq(from=0, to=.99, by=.01)
yseq <- 1/(1-xseq)

#colors <- ifelse(coexists == 0, "grey", "black")
cols = c("turquoise","turquoise3", "paleturquoise", "darkblue", "black",
         "paleturquoise4", "blue", "lightseagreen", "palevioletred",
         "palevioletred4", "goldenrod", "plum", "orange","red", "gold")

pairs <- c("ARCA & PEAI", "ARCA & POCA", "HYGL & MEDI", "HYGL & PEAI", "HYGL & PLDE",
           "HYGL & VERO", "MEDI & PEAI", "PEAI & PLDE", "PEAI & POCA", "PEAI & TRCY", "PEAI & VERO", "PLDE & POCA",
           "POCA & TRCY", "POCA & VERO")

do.transparent <- function(InCol, NewAlpha){
  rgb_vals <- col2rgb(InCol)[,1]
  OutCol <- rgb(red = rgb_vals[1], green = rgb_vals[2], blue = rgb_vals[3],
                alpha = NewAlpha, maxColorValue = 255)
  return(OutCol)
}

tgreen <- do.transparent("seagreen4", 100)
# make graph
pdf("Figures/EqualizingStabilizing.pdf")
plot(mean.corrected.stabilizing, log(mean.equalizing), cex=1.5, pch=c(1:14), col=cols, xlab="Stabilizing Niche Differences", ylab="Log Fitness Differences",  xlim=c(0,1), ylim=c(0,12))
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
text("coexistence", x=0.8, y=1)
text("comp. exclusion", x=0.5, y=3)
legend("topright", legend=pairs, col=cols, ncol=2, bty="n", pch=c(1:14), border=FALSE, cex=1)

dev.off()

