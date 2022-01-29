# Calculate pair-wise low-density growth rate calculations
library(rstan)

source("Bayes.data.R")

# LDGR into ARCA as resident ---------------------------
load("median.abundance.arca2.rdata")
arca.equilibrum <- abundances.arca[200]
invader.abund <- 1

# hygl, peai, plde, poca
# create all of our LDGR matrices
median.hygl.into.arca <- rep(NA, runs)
median.peai.into.arca <- rep(NA, runs)
median.plde.into.arca <- rep(NA, runs)
median.poca.into.arca <- rep(NA, runs)

#y <- sample(dim(ag_carrying)[1], 1)

# calculate resident abundance
Nj <- arca.equilibrum
germ_j <- median(germination$arca)
lived.arca <- germ_j*Nj
#lived.arca <- min(germ_j*Nj, ag_carrying$x[y])

# invade HYGL
hygl_tp1 <- (1-median(germination$hygl))*median(survival$hygl)*invader.abund + 
  invader.abund*median(germination$hygl)*median(hygl$lambda)/(1+median(hygl$alpha_arca)*lived.arca)
# calculate LDGR of HYGL
median.hygl.into.arca <- log(hygl_tp1/invader.abund)

# invade PEAI
peai_tp1 <- (1-median(germination$peai))*median(survival$peai)*invader.abund + 
  invader.abund*median(germination$peai)*median(peai$lambda)/(1+median(peai$alpha_arca)*lived.arca)
# calculate LDGR of PEAI
median.peai.into.arca <- log(peai_tp1/invader.abund)

# invade POCA
poca_tp1 <- (1-median(germination$poca))*median(survival$poca)*invader.abund + 
  invader.abund*median(germination$poca)*median(poca$lambda)/(1+median(poca$alpha_arca)*lived.arca)
# calculate LDGR of PEAI
median.poca.into.arca <- log(poca_tp1/invader.abund)

# LDGR INTO HYGL AS RESIDENT------------------------------------------------
load("median.abundance.hygl2.rdata")
hygl.equilibrum <- abundances.hygl[200]
invader.abund <- 1

# medi, peai, plde, vero
#posteriors <- 100

# create all of our LDGR matrices
median.medi.into.hygl <- rep(NA, runs)
median.peai.into.hygl <- rep(NA, runs)
median.plde.into.hygl <- rep(NA, runs)
median.vero.into.hygl <- rep(NA, runs)

#y <- sample(dim(ag_carrying)[1], 1)

# calculate resident abundance
Nj <- hygl.equilibrum
germ_j <- median(germination$hygl)
lived.hygl <- germ_j*Nj
#lived.hygl<- min(germ_j*Nj, ag_carrying$x[y])

# invade MEDI
medi_tp1 <- (1-median(germination$medi))*median(survival$medi)*invader.abund + 
  invader.abund*median(germination$medi)*median(medi$lambda)/(1+median(medi$alpha_hygl)*lived.hygl)
# calculate LDGR of medi
median.medi.into.hygl <- log(medi_tp1/invader.abund)

# invade PEAI
peai_tp1 <- (1-median(germination$peai))*median(survival$peai)*invader.abund + 
  invader.abund*median(germination$peai)*median(peai$lambda)/(1+median(peai$alpha_hygl)*lived.hygl)
# calculate LDGR of peai
median.peai.into.hygl <- log(peai_tp1/invader.abund)

# invade PLDE 
plde_tp1 <- (1-median(germination$plde))*median(survival$plde)*invader.abund + 
  invader.abund*median(germination$plde)*median(plde$lambda)/(1+median(plde$alpha_hygl)*lived.hygl)
# calculate LDGR of plde
median.plde.into.hygl <- log(plde_tp1/invader.abund)

# invade VERO
vero_tp1 <- (1-median(germination$vero))*median(survival$vero)*invader.abund + 
  invader.abund*median(germination$vero)*median(vero$lambda)/(1+median(vero$alpha_hygl)*lived.hygl)
# calculate LDGR of vero
median.vero.into.hygl <- log(vero_tp1/invader.abund)

#}

# LDGR INTO MEDI AS RESIDENT ---------------------------------------
load("median.abundance.medi2.rdata")
medi.equilibrum <- abundances.medi[200]
invader.abund <- 1

# dagl, hygl, peai 
#posteriors <- 100
# create all of our LDGR matrices
median.dagl.into.medi <- rep(NA, runs)
median.hygl.into.medi <- rep(NA, runs)
median.peai.into.medi <- rep(NA, runs)

#y <- sample(dim(ag_carrying)[1], 1)

# calculate resident abundance
Nj <- medi.equilibrum
germ_j <- median(germination$medi)
lived.medi <- germ_j*Nj
#lived.medi <- min(germ_j*Nj, ag_carrying$x[y])

# invade dagl
#dagl_tp1 <- (1-median(germination$dagl))*median(survival$dagl)*invader.abund + 
#  invader.abund*median(germination$dagl)*median(dagl$lambda)/(1+median(dagl$alpha_medi)*lived.medi)
# calculate LDGR of dagl
#median.dagl.into.medi <- log(dagl_tp1/invader.abund)

# invade hygl
hygl_tp1 <- (1-median(germination$hygl))*median(survival$hygl)*invader.abund + 
  invader.abund*median(germination$hygl)*median(hygl$lambda)/(1+median(hygl$alpha_medi)*lived.medi)
# calculate LDGR of hygl
median.hygl.into.medi <- log(hygl_tp1/invader.abund)

# invade peai
peai_tp1 <- (1-median(germination$peai))*median(survival$peai)*invader.abund + 
  invader.abund*median(germination$peai)*median(peai$lambda)/(1+median(peai$alpha_medi)*lived.medi)
# calculate LDGR of peai
median.peai.into.medi <- log(peai_tp1/invader.abund)

# LDGR INTO PEAI AS RESIDENT ---------------------------
load("median.abundance.peai2.rdata")
peai.equilibrum <- abundances.peai[200]
invader.abund <- 1

# dagl, hygl, peai 

# create all of our LDGR matrices
#arca, dagl, hygl, medi, plde, poca, trcy, vero
median.arca.into.peai <- rep(NA, runs)
#median.dagl.into.peai <- rep(NA, runs)
median.hygl.into.peai <- rep(NA, runs)
median.medi.into.peai <- rep(NA, runs)
median.plde.into.peai <- rep(NA, runs)
median.poca.into.peai <- rep(NA, runs)
median.trcy.into.peai <- rep(NA, runs)
median.vero.into.peai <- rep(NA, runs)

#y <- sample(dim(ag_carrying)[1], 1)

# calculate resident abundance
Nj <- peai.equilibrum
germ_j <- median(germination$peai)
#lived.peai <- min(germ_j*Nj, ag_carrying$x[y])
lived.peai <- germ_j*Nj

# invade arca
arca_tp1 <- (1-median(germination$arca))*median(survival$arca)*invader.abund + 
  invader.abund*median(germination$arca)*median(arca$lambda)/(1+median(arca$alpha_peai)*lived.peai)
# calculate LDGR of arca
median.arca.into.peai <- log(arca_tp1/invader.abund)

# invade hygl
hygl_tp1 <- (1-median(germination$hygl))*median(survival$hygl)*invader.abund + 
  invader.abund*median(germination$hygl)*median(hygl$lambda)/(1+median(hygl$alpha_peai)*lived.peai)
# calculate LDGR of hygl
median.hygl.into.peai <- log(hygl_tp1/invader.abund)

# invade medi
medi_tp1 <- (1-median(germination$medi))*median(survival$medi)*invader.abund + 
  invader.abund*median(germination$medi)*median(medi$lambda)/(1+median(medi$alpha_peai)*lived.peai)
# calculate LDGR of medi
median.medi.into.peai <- log(medi_tp1/invader.abund)

# invade plde
plde_tp1 <- (1-median(germination$plde))*median(survival$plde)*invader.abund + 
  invader.abund*median(germination$plde)*median(plde$lambda)/(1+median(plde$alpha_peai)*lived.peai)
# calculate LDGR of plde
median.plde.into.peai <- log(plde_tp1/invader.abund)

# invade poca
poca_tp1 <- (1-median(germination$poca))*median(survival$poca)*invader.abund + 
  invader.abund*median(germination$poca)*median(poca$lambda)/(1+median(poca$alpha_peai)*lived.peai)
# calculate LDGR of poca
median.poca.into.peai <- log(poca_tp1/invader.abund)

# invade trcy
trcy_tp1 <- (1-median(germination$trcy))*median(survival$trcy)*invader.abund + 
  invader.abund*median(germination$trcy)*median(trcy$lambda)/(1+median(trcy$alpha_peai)*lived.peai)
# calculate LDGR of trcy
median.trcy.into.peai <- log(trcy_tp1/invader.abund)

# invade vero
vero_tp1 <- (1-median(germination$vero))*median(survival$vero)*invader.abund + 
  invader.abund*median(germination$vero)*median(vero$lambda)/(1+median(vero$alpha_peai)*lived.peai)
# calculate LDGR of vero
median.vero.into.peai <- log(vero_tp1/invader.abund)


# LDGR INTO PLDE AS RESIDENT ----------------------------
load("median.abundance.plde2.rdata")
plde.equilibrum <- abundances.plde[200]
invader.abund <- 1

# dagl, hygl, peai 

# create all of our LDGR matrices
#arca, dagl, medi, peai, poca, trcy
median.hygl.into.plde <- rep(NA, runs)
median.peai.into.plde <- rep(NA, runs)
median.poca.into.plde <- rep(NA, runs)
median.trcy.into.plde <- rep(NA, runs)

#y <- sample(dim(ag_carrying)[1], 1)

# calculate resident abundance
Nj <- plde.equilibrum
germ_j <- median(germination$plde)
#lived.plde <- min(germ_j*Nj, ag_carrying$x[y])
lived.plde <- germ_j*Nj

# invade arca
#arca_tp1 <- (1-median(germination$arca))*median(survival$arca)*invader.abund + 
#  invader.abund*median(germination$arca)*median(arca$lambda)/(1+median(arca$alpha_plde)*lived.plde)
# calculate LDGR of arca
#median.arca.into.plde <- log(arca_tp1/invader.abund)

# invade hygl
hygl_tp1 <- (1-median(germination$hygl))*median(survival$hygl)*invader.abund + 
  invader.abund*median(germination$hygl)*median(hygl$lambda)/(1+median(hygl$alpha_plde)*lived.plde)
# calculate LDGR of hygl
median.hygl.into.plde <- log(hygl_tp1/invader.abund)

# invade peai
peai_tp1 <- (1-median(germination$peai))*median(survival$peai)*invader.abund + 
  invader.abund*median(germination$peai)*median(peai$lambda)/(1+median(peai$alpha_plde)*lived.plde)
# calculate LDGR of peai
median.peai.into.plde <- log(peai_tp1/invader.abund)

# invade poca
poca_tp1 <- (1-median(germination$poca))*median(survival$poca)*invader.abund + 
  invader.abund*median(germination$poca)*median(poca$lambda)/(1+median(poca$alpha_plde)*lived.plde)
# calculate LDGR of poca
median.poca.into.plde <- log(poca_tp1/invader.abund)

# invade trcy
trcy_tp1 <- (1-median(germination$trcy))*median(survival$trcy)*invader.abund + 
  invader.abund*median(germination$trcy)*median(trcy$lambda)/(1+median(trcy$alpha_plde)*lived.plde)
# calculate LDGR of trcy
median.trcy.into.plde <- log(trcy_tp1/invader.abund)

# LDGR INTO POCA RESIDENT -----------------------------
load("median.abundance.poca2.rdata")
poca.equilibrum <- abundances.poca[200]
invader.abund <- 1

# create all of our LDGR matrices
# arca, dagl, medi, peai, plde, trcy, vero
median.arca.into.poca <- rep(NA, runs)
median.medi.into.poca <- rep(NA, runs)
median.peai.into.poca <- rep(NA, runs)
median.plde.into.poca <- rep(NA, runs)
median.trcy.into.poca <- rep(NA, runs)
median.vero.into.poca <- rep(NA, runs)

#y <- sample(dim(ag_carrying)[1], 1)

# calculate resident abundance
Nj <- poca.equilibrum
germ_j <- median(germination$poca)
#lived.poca <- min(germ_j*Nj, ag_carrying$x[y])
lived.poca <- germ_j*Nj

# invade arca
arca_tp1 <- (1-median(germination$arca))*median(survival$arca)*invader.abund + 
  invader.abund*median(germination$arca)*median(arca$lambda)/(1+median(arca$alpha_poca)*lived.poca)
# calculate LDGR of arca
median.arca.into.poca <- log(arca_tp1/invader.abund)

# invade medi
medi_tp1 <- (1-median(germination$medi))*median(survival$medi)*invader.abund + 
  invader.abund*median(germination$medi)*median(medi$lambda)/(1+median(medi$alpha_poca)*lived.poca)
# calculate LDGR of medi
median.medi.into.poca <- log(medi_tp1/invader.abund)

# invade peai
peai_tp1 <- (1-median(germination$peai))*median(survival$peai)*invader.abund + 
  invader.abund*median(germination$peai)*median(peai$lambda)/(1+median(peai$alpha_poca)*lived.poca)
# calculate LDGR of peai
median.peai.into.poca <- log(peai_tp1/invader.abund)

# invade plde
plde_tp1 <- (1-median(germination$plde))*median(survival$plde)*invader.abund + 
  invader.abund*median(germination$plde)*median(plde$lambda)/(1+median(plde$alpha_poca)*lived.poca)
# calculate LDGR of plde
median.plde.into.poca <- log(plde_tp1/invader.abund)

# invade trcy
trcy_tp1 <- (1-median(germination$trcy))*median(survival$trcy)*invader.abund + 
  invader.abund*median(germination$trcy)*median(trcy$lambda)/(1+median(trcy$alpha_poca)*lived.poca)
# calculate LDGR of trcy
median.trcy.into.poca <- log(trcy_tp1/invader.abund)

# invade vero
vero_tp1 <- (1-median(germination$vero))*median(survival$vero)*invader.abund + 
  invader.abund*median(germination$vero)*median(vero$lambda)/(1+median(vero$alpha_poca)*lived.poca)
# calculate LDGR of vero
median.vero.into.poca <- log(vero_tp1/invader.abund)


# LDGR INTO TRCY AS RESIDENT --------------------------
load("median.abundance.trcy2.rdata")
trcy.equilibrum <- abundances.trcy[200]
invader.abund <- 1

# create all of our LDGR matrices
# arca, dagl, peai, plde, poca
median.arca.into.trcy <- rep(NA, runs)
median.peai.into.trcy <- rep(NA, runs)
median.poca.into.trcy <- rep(NA, runs)

#y <- sample(dim(ag_carrying)[1], 1)

# calculate resident abundance
Nj <- trcy.equilibrum
germ_j <- median(germination$trcy)
#lived.trcy <- min(germ_j*Nj, ag_carrying$x[y])
lived.trcy <- germ_j*Nj

# invade arca
arca_tp1 <- (1-median(germination$arca))*median(survival$arca)*invader.abund + 
  invader.abund*median(germination$arca)*median(arca$lambda)/(1+median(arca$alpha_trcy)*lived.trcy)
# calculate LDGR of arca
median.arca.into.trcy <- log(arca_tp1/invader.abund)

# invade peai
peai_tp1 <- (1-median(germination$peai))*median(survival$peai)*invader.abund + 
  invader.abund*median(germination$peai)*median(peai$lambda)/(1+median(peai$alpha_trcy)*lived.trcy)
# calculate LDGR of peai
median.peai.into.trcy <- log(peai_tp1/invader.abund)

# invade poca
poca_tp1 <- (1-median(germination$poca))*median(survival$poca)*invader.abund + 
  invader.abund*median(germination$poca)*median(poca$lambda)/(1+median(poca$alpha_trcy)*lived.trcy)
# calculate LDGR of poca
median.poca.into.trcy <- log(poca_tp1/invader.abund)

# LDGR INTO VERO AS RESIDENT --------------------------
load("median.abundance.vero2.rdata")
vero.equilibrum <- abundances.vero[200]
invader.abund <- 1


# create all of our LDGR matrices
# arca, dagl, hygl, medi, peai, plde, poca, trcy, vero
median.arca.into.vero <- rep(NA, runs)
median.hygl.into.vero <- rep(NA, runs)
median.medi.into.vero <- rep(NA, runs)
median.peai.into.vero <- rep(NA, runs)
median.plde.into.vero <- rep(NA, runs)
median.poca.into.vero <- rep(NA, runs)
median.trcy.into.vero <- rep(NA, runs)

#y <- sample(dim(ag_carrying)[1], 1)

# calculate resident abundance
Nj <- vero.equilibrum
germ_j <- median(germination$vero)
#lived.vero <- min(germ_j*Nj, ag_carrying$x[y])
lived.vero <- germ_j*Nj

# invade arca
arca_tp1 <- (1-median(germination$arca))*median(survival$arca)*invader.abund + 
  invader.abund*median(germination$arca)*median(arca$lambda)/(1+median(arca$alpha_vero)*lived.vero)
# calculate LDGR of arca
median.arca.into.vero <- log(arca_tp1/invader.abund)

# invade hygl
hygl_tp1 <- (1-median(germination$hygl))*median(survival$hygl)*invader.abund + 
  invader.abund*median(germination$hygl)*median(hygl$lambda)/(1+median(hygl$alpha_vero)*lived.vero)
# calculate LDGR of hygl
median.hygl.into.vero <- log(hygl_tp1/invader.abund)

# invade peai
peai_tp1 <- (1-median(germination$peai))*median(survival$peai)*invader.abund + 
  invader.abund*median(germination$peai)*median(peai$lambda)/(1+median(peai$alpha_vero)*lived.vero)
# calculate LDGR of peai
median.peai.into.vero <- log(peai_tp1/invader.abund)

# invade plde
plde_tp1 <- (1-median(germination$plde))*median(survival$plde)*invader.abund + 
  invader.abund*median(germination$plde)*median(plde$lambda)/(1+median(plde$alpha_vero)*lived.vero)
# calculate LDGR of plde
median.plde.into.vero <- log(plde_tp1/invader.abund)

# invade poca
poca_tp1 <- (1-median(germination$poca))*median(survival$poca)*invader.abund + 
  invader.abund*median(germination$poca)*median(poca$lambda)/(1+median(poca$alpha_vero)*lived.vero)
# calculate LDGR of poca
median.poca.into.vero <- log(poca_tp1/invader.abund)

# invade medi
medi_tp1 <- (1-median(germination$medi))*median(survival$medi)*invader.abund + 
  invader.abund*median(germination$medi)*median(medi$lambda)/(1+median(medi$alpha_vero)*lived.vero)
# calculate LDGR of medi
median.medi.into.vero <- log(medi_tp1/invader.abund)

# invade trcy
trcy_tp1 <- (1-median(germination$trcy))*median(survival$trcy)*invader.abund + 
  invader.abund*median(germination$trcy)*median(trcy$lambda)/(1+median(trcy$alpha_vero)*lived.vero)
# calculate LDGR of trcy
median.trcy.into.vero <- log(trcy_tp1/invader.abund)


