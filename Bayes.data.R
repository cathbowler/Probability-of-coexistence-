# import data

library(coda)
library(rstan)

sp <- 9
# removing MOMO for now

sp_list <- c("dagl", "hygl", "plde", 
             "poca", "trcy", "vero", "arca", "medi", "peai")
#sp_list <- c(dagl, gite, hygl, plde, 
#             poca, trcy, vero, arca, medi, momo, peai)

#setwd("/Users/lash1937/Dropbox/Shared/Bayesian data - Cath Collaboration/Updated Coexistence code")
#setwd("/Users/catherinebowler/Dropbox/Bayesian data - Cath Collaboration/Coexistence code")
load("BH stan models//Arca_posteriors.rdata")
arca <- extract(PrelimFit)
remove(PrelimFit)

load("BH stan models//Dagl_posteriors.rdata")
dagl <- extract(PrelimFit)
remove(PrelimFit)

load("BH stan models//Hygl_posteriors.rdata")
hygl <- extract(PrelimFit)
remove(PrelimFit)

load("BH stan models//Medi_posteriors.rdata")
medi <- extract(PrelimFit)
remove(PrelimFit)

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


survival <- data.frame(dagl=surv$p.1, hygl=surv$p.2, plde=surv$p.3, poca=surv$p.4, trcy=surv$p.5, 
                       vero=surv$p.6, arca=surv$p.7, medi=surv$p.8, peai=surv$p.9)





