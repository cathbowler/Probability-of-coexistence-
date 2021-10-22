##############################################################################################
# This code is used to load the Bayesian stan model fit posteriors 
# this script can be called at the start of each subsequent script that uses this data 
##############################################################################################

library(coda)
library(rstan)

sp <- 8
# removing MOMO for now

sp_list <- c("hygl", "plde", 
             "poca", "trcy", "vero", "arca", "medi", "peai")
#sp_list <- c(dagl, gite, hygl, plde, 
#             poca, trcy, vero, arca, medi, momo, peai)

#setwd("/Users/lash1937/Dropbox/Shared/Bayesian data - Cath Collaboration/Updated Coexistence code")
#setwd("/Users/catherinebowler/Dropbox/Bayesian data - Cath Collaboration/Coexistence code")
load("BH stan models//Arca_posteriors.rdata")
arca <- rstan::extract(PrelimFit)
remove(PrelimFit)
load("BH stan models//Arca_posteriors_5.rdata")
arca_5 <- rstan::extract(PrelimFit)
remove(PrelimFit)
load("BH stan models//Arca_posteriors_10.rdata")
arca_10 <- rstan::extract(PrelimFit)
remove(PrelimFit)

load("BH stan models//Hygl_posteriors.rdata")
hygl <- rstan::extract(PrelimFit)
remove(PrelimFit)
load("BH stan models//Hygl_posteriors_5.rdata")
hygl_5 <- rstan::extract(PrelimFit)
remove(PrelimFit)
load("BH stan models//Hygl_posteriors_10.rdata")
hygl_10 <- rstan::extract(PrelimFit)
remove(PrelimFit)

load("BH stan models//Medi_posteriors.rdata")
medi <- rstan::extract(PrelimFit)
remove(PrelimFit)
load("BH stan models//Medi_posteriors_5.rdata")
medi_5 <- rstan::extract(PrelimFit)
remove(PrelimFit)
load("BH stan models//Medi_posteriors_10.rdata")
medi_10 <- rstan::extract(PrelimFit)
remove(PrelimFit)

load("BH stan models//Peai_posteriors.rdata")
peai <- rstan::extract(PrelimFit)
remove(PrelimFit)
load("BH stan models//Peai_posteriors_5.rdata")
peai_5 <- rstan::extract(PrelimFit)
remove(PrelimFit)
load("BH stan models//Peai_posteriors_10.rdata")
peai_10 <- rstan::extract(PrelimFit)
remove(PrelimFit)

load("BH stan models//Plde_posteriors.rdata")
plde <- rstan::extract(PrelimFit)
remove(PrelimFit)
load("BH stan models//Plde_posteriors_5.rdata")
plde_5 <- rstan::extract(PrelimFit)
remove(PrelimFit)
load("BH stan models//Plde_posteriors_10.rdata")
plde_10 <- rstan::extract(PrelimFit)
remove(PrelimFit)

load("BH stan models//Poca_posteriors.rdata")
poca <- rstan::extract(PrelimFit)
remove(PrelimFit)
load("BH stan models//Poca_posteriors_5.rdata")
poca_5 <- rstan::extract(PrelimFit)
remove(PrelimFit)
load("BH stan models//Poca_posteriors_10.rdata")
poca_10 <- rstan::extract(PrelimFit)
remove(PrelimFit)

load("BH stan models//Trcy_posteriors.rdata")
trcy <- rstan::extract(PrelimFit)
remove(PrelimFit)
load("BH stan models//Trcy_posteriors_5.rdata")
trcy_5 <- rstan::extract(PrelimFit)
remove(PrelimFit)
load("BH stan models//Trcy_posteriors_10.rdata")
trcy_10 <- rstan::extract(PrelimFit)
remove(PrelimFit)

load("BH stan models//Vero_posteriors.rdata")
vero <- rstan::extract(PrelimFit)
remove(PrelimFit)
load("BH stan models//Vero_posteriors_5.rdata")
vero_5 <- rstan::extract(PrelimFit)
remove(PrelimFit)
load("BH stan models//Vero_posteriors_10.rdata")
vero_10 <- rstan::extract(PrelimFit)
remove(PrelimFit)

load("SurvivalAndGermination/Germination.rdata")
germ <- rstan::extract(PrelimFit)
germ<-as.data.frame(germ)
remove(PrelimFit)

load("SurvivalAndGermination/Survival.rdata")
surv <- rstan::extract(PrelimFit)
surv<-as.data.frame(surv)
remove(PrelimFit)

germination <- data.frame(hygl=germ$p.2, plde=germ$p.3, poca=germ$p.4, trcy=germ$p.5, 
                          vero=germ$p.6, arca=germ$p.7, medi=germ$p.8, peai=germ$p.9)


survival <- data.frame(hygl=surv$p.2, plde=surv$p.3, poca=surv$p.4, trcy=surv$p.5, 
                       vero=surv$p.6, arca=surv$p.7, medi=surv$p.8, peai=surv$p.9)

# plotting difference between priors 
cols = c("black", "#F8766D","#E7861B", "#95A900", "#39B600", "#00C19C",
         "#0088D8", "#00A5FF", "#7997FF", "grey")
species.names = c("Intra", "H. glutinosum", "P. debilis", "P. canescens", 
                  "T. cyanopetala", "G. rosea", "A. calendula", "M. minima", "P. airoides", "Other")

# dev.off()
pdf("Figures/supp.priors_5.pdf")
 par(mfrow=c(4,2), mar=c(4,4,2,2))
 #arca
 plot(density(abs(arca$alpha_intra)-abs(arca_5$alpha_intra)), col="black", xlim=c(-3,1), main=expression(italic("A. calendula")), xlab="")
 lines(density(abs(arca$alpha_poca)-abs(arca_5$alpha_poca)), col=cols[3])
 lines(density(abs(arca$alpha_plde)-abs(arca_5$alpha_plde)), col=cols[2])
 lines(density(abs(arca$alpha_trcy)-abs(arca_5$alpha_trcy)), col=cols[4])
 lines(density(abs(arca$alpha_vero)-abs(arca_5$alpha_vero)), col=cols[5]) 
 lines(density(abs(arca$alpha_peai)-abs(arca_5$alpha_peai)), col=cols[8])
 lines(density(abs(arca$alpha_other)-abs(arca_5$alpha_other)), col=cols[9])
 
 lines(density(abs(arca$alpha_intra)-abs(arca_10$alpha_intra)), col="black", lty=2)
 lines(density(abs(arca$alpha_poca)-abs(arca_10$alpha_poca)), col=cols[3], lty=2)
 lines(density(abs(arca$alpha_plde)-abs(arca_10$alpha_plde)), col=cols[2], lty=2)
 lines(density(abs(arca$alpha_trcy)-abs(arca_10$alpha_trcy)), col=cols[4], lty=2)
 lines(density(abs(arca$alpha_vero)-abs(arca_10$alpha_vero)), col=cols[5], lty=2) 
 lines(density(abs(arca$alpha_peai)-abs(arca_10$alpha_peai)), col=cols[8], lty=2)
 lines(density(abs(arca$alpha_other)-abs(arca_10$alpha_other)), col=cols[9], lty=2)
 legend("topleft", legend = species.names, col = cols, lty=1, ncol=2, bty="n")

#hygl
 plot(density(abs(hygl$alpha_intra)-abs(hygl_5$alpha_intra)), col="black", xlim=c(-.2,.2), main=expression(italic("H. glutinosum")), xlab="")
 lines(density(abs(hygl$alpha_plde)-abs(hygl_5$alpha_plde)), col=cols[2])
 lines(density(abs(hygl$alpha_poca)-abs(hygl_5$alpha_poca)), col=cols[3])
 lines(density(abs(hygl$alpha_medi)-abs(hygl_5$alpha_medi)), col=cols[7])
 lines(density(abs(hygl$alpha_trcy)-abs(hygl_5$alpha_trcy)), col=cols[4])
 lines(density(abs(hygl$alpha_vero)-abs(hygl_5$alpha_vero)), col=cols[5])
 lines(density(abs(hygl$alpha_peai)-abs(hygl_5$alpha_peai)), col=cols[8])
 lines(density(abs(hygl$alpha_other)-abs(hygl_5$alpha_other)), col=cols[9])

 lines(density(abs(hygl$alpha_intra)-abs(hygl_10$alpha_intra)), col="black", lty=2)
 lines(density(abs(hygl$alpha_plde)-abs(hygl_10$alpha_plde)), col=cols[2], lty=2)
 lines(density(abs(hygl$alpha_poca)-abs(hygl_10$alpha_poca)), col=cols[3], lty=2)
 lines(density(abs(hygl$alpha_medi)-abs(hygl_10$alpha_medi)), col=cols[7], lty=2)
 lines(density(abs(hygl$alpha_trcy)-abs(hygl_10$alpha_trcy)), col=cols[4], lty=2)
 lines(density(abs(hygl$alpha_vero)-abs(hygl_10$alpha_vero)), col=cols[5], lty=2)
 lines(density(abs(hygl$alpha_peai)-abs(hygl_10$alpha_peai)), col=cols[8], lty=2)
 lines(density(abs(hygl$alpha_other)-abs(hygl_10$alpha_other)), col=cols[9], lty=2)

 #plde 
 plot(density(abs(plde$alpha_intra)-abs(plde_5$alpha_intra)), col="black", main=expression(italic("P. debilis")), xlim =c(-.2,.2), xlab="")
 lines(density(abs(plde$alpha_hygl)-abs(plde_5$alpha_hygl)), col=cols[2])
 lines(density(abs(plde$alpha_poca)-abs(plde_5$alpha_poca)), col=cols[3])
 lines(density(abs(plde$alpha_vero)-abs(plde_5$alpha_vero)), col=cols[5])
 lines(density(abs(plde$alpha_peai)-abs(plde_5$alpha_peai)), col=cols[8])
 lines(density(abs(plde$alpha_other)-abs(plde_5$alpha_other)), col=cols[9])
 
 lines(density(abs(plde$alpha_intra)-abs(plde_10$alpha_intra)), col="black", lty=2)
 lines(density(abs(plde$alpha_hygl)-abs(plde_10$alpha_hygl)), col=cols[2], lty=2)
 lines(density(abs(plde$alpha_poca)-abs(plde_10$alpha_poca)), col=cols[3], lty=2)
 lines(density(abs(plde$alpha_vero)-abs(plde_10$alpha_vero)), col=cols[5], lty=2)
 lines(density(abs(plde$alpha_peai)-abs(plde_10$alpha_peai)), col=cols[8], lty=2)
 lines(density(abs(plde$alpha_other)-abs(plde_10$alpha_other)), col=cols[9], lty=2)

 #poca
 plot(density(abs(poca$alpha_intra)-abs(poca_5$alpha_intra)), col="black", main=expression(italic("P. canescens")), xlim =c(-.2,.2), xlab="")
 lines(density(abs(poca$alpha_plde)-abs(poca_5$alpha_plde)), col=cols[2])
 lines(density(abs(poca$alpha_arca)-abs(poca_5$alpha_arca)), col=cols[6])
 lines(density(abs(poca$alpha_trcy)-abs(poca_5$alpha_trcy)), col=cols[4])
 lines(density(abs(poca$alpha_vero)-abs(poca_5$alpha_vero)), col=cols[5])
 lines(density(abs(poca$alpha_peai)-abs(poca_5$alpha_peai)), col=cols[8])
 lines(density(abs(poca$alpha_other)-abs(poca_5$alpha_other)), col=cols[9])
 
 lines(density(abs(poca$alpha_intra)-abs(poca_10$alpha_intra)), col="black", lty=2)
 lines(density(abs(poca$alpha_plde)-abs(poca_10$alpha_plde)), col=cols[2], lty=2)
 lines(density(abs(poca$alpha_arca)-abs(poca_10$alpha_arca)), col=cols[6], lty=2)
 lines(density(abs(poca$alpha_trcy)-abs(poca_10$alpha_trcy)), col=cols[4], lty=2)
 lines(density(abs(poca$alpha_vero)-abs(poca_10$alpha_vero)), col=cols[5], lty=2)
 lines(density(abs(poca$alpha_peai)-abs(poca_10$alpha_peai)), col=cols[8], lty=2)
 lines(density(abs(poca$alpha_other)-abs(poca_10$alpha_other)), col=cols[9], lty=2) 

 #peai
 plot(density(abs(peai$alpha_intra)-abs(peai_5$alpha_intra)), col="black", main=expression(italic("P. airoides")), xlim =c(-.2,.2), xlab="")
 lines(density(abs(peai$alpha_arca)-abs(peai_5$alpha_arca)), col=cols[6])
 lines(density(abs(peai$alpha_plde)-abs(peai_5$alpha_plde)), col=cols[2])
  lines(density(abs(peai$alpha_arca)-abs(peai_5$alpha_arca)), col=cols[3])
 lines(density(abs(peai$alpha_trcy)-abs(peai_5$alpha_trcy)), col=cols[4])
 lines(density(abs(peai$alpha_vero)-abs(peai_5$alpha_vero)), col=cols[5])
 lines(density(abs(peai$alpha_medi)-abs(peai_5$alpha_medi)), col=cols[7])
 lines(density(abs(peai$alpha_other)-abs(peai_5$alpha_other)), col=cols[9])
 
 lines(density(abs(peai$alpha_intra)-abs(peai_10$alpha_intra)), col="black", lty=2)
 lines(density(abs(peai$alpha_arca)-abs(peai_10$alpha_arca)), col=cols[6], lty=2)
 lines(density(abs(peai$alpha_plde)-abs(peai_10$alpha_plde)), col=cols[2], lty=2)
 lines(density(abs(peai$alpha_arca)-abs(peai_10$alpha_arca)), col=cols[3], lty=2)
 lines(density(abs(peai$alpha_trcy)-abs(peai_10$alpha_trcy)), col=cols[4], lty=2)
 lines(density(abs(peai$alpha_vero)-abs(peai_10$alpha_vero)), col=cols[5], lty=2)
 lines(density(abs(peai$alpha_medi)-abs(peai_10$alpha_medi)), col=cols[7], lty=2)
 lines(density(abs(peai$alpha_other)-abs(peai_10$alpha_other)), col=cols[9], lty=2)

 #medi
 plot(density(abs(medi$alpha_intra)-abs(medi_5$alpha_intra)), col="black", main=expression(italic("M. minima")), xlim =c(-.2,.2), xlab="")
 lines(density(abs(medi$alpha_hygl)-abs(medi_5$alpha_hygl)), col=cols[2])
 lines(density(abs(medi$alpha_poca)-abs(medi_5$alpha_poca)), col=cols[3])
 lines(density(abs(medi$alpha_vero)-abs(medi_5$alpha_vero)), col=cols[5])
 lines(density(abs(medi$alpha_peai)-abs(medi_5$alpha_peai)), col=cols[8])
 lines(density(abs(medi$alpha_other)-abs(medi_5$alpha_other)), col=cols[9])
 
 lines(density(abs(medi$alpha_intra)-abs(medi_10$alpha_intra)), col="black", lty=2)
 lines(density(abs(medi$alpha_hygl)-abs(medi_10$alpha_hygl)), col=cols[2], lty=2)
 lines(density(abs(medi$alpha_poca)-abs(medi_10$alpha_poca)), col=cols[3], lty=2)
 lines(density(abs(medi$alpha_vero)-abs(medi_10$alpha_vero)), col=cols[5], lty=2)
 lines(density(abs(medi$alpha_peai)-abs(medi_10$alpha_peai)), col=cols[8], lty=2)
 lines(density(abs(medi$alpha_other)-abs(medi_10$alpha_other)), col=cols[9], lty=2)

 #trcy
 plot(density(abs(trcy$alpha_intra)-abs(trcy_10$alpha_intra)), col="black", main=expression(italic("T. cyanopetala")), xlim =c(-.2,.2), xlab="absolute difference")
 lines(density(abs(trcy$alpha_plde)-abs(trcy_10$alpha_plde)), col=cols[2])
 lines(density(abs(trcy$alpha_poca)-abs(trcy_10$alpha_poca)), col=cols[3])
 lines(density(abs(trcy$alpha_vero)-abs(trcy_10$alpha_vero)), col=cols[5])
 lines(density(abs(trcy$alpha_peai)-abs(trcy_10$alpha_peai)), col=cols[8])
 lines(density(abs(trcy$alpha_other)-abs(trcy_10$alpha_other)), col=cols[9])
 
 lines(density(abs(trcy$alpha_intra)-abs(trcy_5$alpha_intra)), col="black", lty=2)
 lines(density(abs(trcy$alpha_plde)-abs(trcy_5$alpha_plde)), col=cols[2], lty=2)
 lines(density(abs(trcy$alpha_poca)-abs(trcy_5$alpha_poca)), col=cols[3], lty=2)
 lines(density(abs(trcy$alpha_vero)-abs(trcy_5$alpha_vero)), col=cols[5], lty=2)
 lines(density(abs(trcy$alpha_peai)-abs(trcy_5$alpha_peai)), col=cols[8], lty=2)
 lines(density(abs(trcy$alpha_other)-abs(trcy_5$alpha_other)), col=cols[9], lty=2)

 #vero
 plot(density(abs(vero$alpha_intra)-abs(vero_5$alpha_intra)), col="black", main=expression(italic("G. rosea")), xlim =c(-.2,.2), xlab="absolute difference")
 lines(density(abs(vero$alpha_hygl)-abs(vero_5$alpha_hygl)), col=cols[2])
 lines(density(abs(vero$alpha_poca)-abs(vero_5$alpha_poca)), col=cols[3])
 lines(density(abs(vero$alpha_trcy)-abs(vero_5$alpha_trcy)), col=cols[4])
 lines(density(abs(vero$alpha_peai)-abs(vero_5$alpha_peai)), col=cols[8])
 lines(density(abs(vero$alpha_other)-abs(vero_5$alpha_other)), col=cols[9])

 lines(density(abs(vero$alpha_intra)-abs(vero_10$alpha_intra)), col=cols[1], lty=2)
 lines(density(abs(vero$alpha_hygl)-abs(vero_10$alpha_hygl)), col=cols[2], lty=2)
 lines(density(abs(vero$alpha_poca)-abs(vero_10$alpha_poca)), col=cols[3], lty=2)
 lines(density(abs(vero$alpha_trcy)-abs(vero_10$alpha_trcy)), col=cols[4], lty=2)
 lines(density(abs(vero$alpha_peai)-abs(vero_10$alpha_peai)), col=cols[8], lty=2)
 lines(density(abs(vero$alpha_other)-abs(vero_10$alpha_other)), col=cols[9], lty=2)
dev.off() 

