# Calculate pair-wise low-density growth rate calculations
# load in Stan model fits
rm(list=ls())
sp <- 9
# removing MOMO for now

sp_list <- c("dagl", "hygl", "plde", 
             "poca", "trcy", "vero", "arca", "medi", "peai")
#sp_list <- c(dagl, gite, hygl, plde, 
#             poca, trcy, vero, arca, medi, momo, peai)

#setwd("/Users/lash1937/Dropbox/Shared/Bayesian data - Cath Collaboration/Updated Coexistence code")

source("Bayes.data.R")

# LDGR into ARCA as resident ---------------------------
load("nolottery.abundance.arca.rdata")
arca.equilibrum <- abundances.arca[,200]
invader.abund <- 1 # start invader abundance at 1 (seed)

# hygl, peai, plde, poca can all invade into arca (have enough data)

runs=4500
# create all of our LDGR matrices
nolot.hygl.into.arca <- matrix(NA, nrow=runs)
nolot.peai.into.arca <- matrix(NA, nrow=runs)
nolot.poca.into.arca <- matrix(NA, nrow=runs)

for (r in 1:runs) {
  p <- sample(seq(1, 4500),1) 
    # calculate resident abundance
    Nj <- arca.equilibrum[r]
    germ_j <- germination$arca[p]
    #lived.arca <- min(germ_j*Nj, ag_carrying$x[y])
    lived.arca <- germ_j*Nj
    
    # invade HYGL
    hygl_tp1 <- (1-germination$hygl[p])*survival$hygl[p]*invader.abund + 
      invader.abund*germination$hygl[p]*hygl$lambda[p]/(1+hygl$alpha_arca[p]*lived.arca)
    # calculate LDGR of HYGL
    nolot.hygl.into.arca[r] <- log(hygl_tp1/invader.abund)
    
    # invade PEAI
    peai_tp1 <- (1-germination$peai[p])*survival$peai[p]*invader.abund + 
      invader.abund*germination$peai[p]*peai$lambda[p]/(1+peai$alpha_arca[p]*lived.arca)
    # calculate LDGR of PEAI
    nolot.peai.into.arca[r] <- log(peai_tp1/invader.abund)
    
    # invade PLDE
    #plde_tp1 <- (1-germination$plde[x])*survival$plde[x]*invader.abund + 
    #  invader.abund*germination$plde[x]*plde$lambda[x]/(1+plde$alpha_arca[x]*lived.arca)
    # calculate LDGR of PLDE
    #nolot.plde.into.arca[r,p] <- log(plde_tp1/invader.abund)
    
    # invade POCA
    poca_tp1 <- (1-germination$poca[p])*survival$poca[p]*invader.abund + 
      invader.abund*germination$poca[p]*poca$lambda[p]/(1+poca$alpha_arca[p]*lived.arca)
    # calculate LDGR of POCA
    nolot.poca.into.arca[r] <- log(poca_tp1/invader.abund)
    
  }


# LDGR INTO HYGL AS RESIDENT------------------------------------------------
load("nolottery.abundance.hygl.rdata")
hygl.equilibrum <- abundances.hygl[,200]
invader.abund <- 1

# create all of our LDGR matrices
nolot.medi.into.hygl <- matrix(NA, nrow=runs)
nolot.peai.into.hygl <- matrix(NA, nrow=runs)
nolot.plde.into.hygl <- matrix(NA, nrow=runs)
nolot.vero.into.hygl <- matrix(NA, nrow=runs)

for (r in 1:runs) {
  p <- sample(seq(1, 4500),1) 
    # calculate resident abundance
    Nj <- hygl.equilibrum[r]
    germ_j <- germination$hygl[p]
    #lived.hygl<- min(germ_j*Nj, ag_carrying$x[y])
    lived.hygl<-germ_j*Nj
    
    # invade MEDI
    medi_tp1 <- (1-germination$medi[p])*survival$medi[p]*invader.abund + 
      invader.abund*germination$medi[p]*medi$lambda[p]/(1+medi$alpha_hygl[p]*lived.hygl)
    # calculate LDGR of medi
    nolot.medi.into.hygl[r] <- log(medi_tp1/invader.abund)
    
    # invade PEAI
    peai_tp1 <- (1-germination$peai[p])*survival$peai[p]*invader.abund + 
      invader.abund*germination$peai[p]*peai$lambda[p]/(1+peai$alpha_hygl[p]*lived.hygl)
    # calculate LDGR of peai
    nolot.peai.into.hygl[r] <- log(peai_tp1/invader.abund)
    
    # invade PLDE 
    plde_tp1 <- (1-germination$plde[p])*survival$plde[p]*invader.abund + 
      invader.abund*germination$plde[p]*plde$lambda[p]/(1+plde$alpha_hygl[p]*lived.hygl)
    # calculate LDGR of plde
    nolot.plde.into.hygl[r] <- log(plde_tp1/invader.abund)
    
    # invade VERO
    vero_tp1 <- (1-germination$vero[p])*survival$vero[p]*invader.abund + 
      invader.abund*germination$vero[p]*vero$lambda[p]/(1+vero$alpha_hygl[p]*lived.hygl)
    # calculate LDGR of vero
    nolot.vero.into.hygl[r] <- log(vero_tp1/invader.abund)
  }


# LDGR INTO MEDI AS RESIDENT ---------------------------------------
load("nolottery.abundance.medi.rdata")
medi.equilibrum <- abundances.medi[,200]
invader.abund <- 1


# create all of our LDGR matrices
nolot.dagl.into.medi <- matrix(NA, nrow=runs)
nolot.hygl.into.medi <- matrix(NA, nrow=runs)
nolot.peai.into.medi <- matrix(NA, nrow=runs)

for (r in 1:runs) {
  p <- sample(seq(1, 4500),1) 
    # calculate resident abundance
    Nj <- medi.equilibrum[r]
    germ_j <- germination$medi[p]
    #lived.medi <- min(germ_j*Nj, ag_carrying$x[y])
    lived.medi <- germ_j*Nj
    
    # invade dagl
    dagl_tp1 <- (1-germination$dagl[p])*survival$dagl[p]*invader.abund + 
      invader.abund*germination$dagl[p]*dagl$lambda[p]/(1+dagl$alpha_medi[p]*lived.medi)
    # calculate LDGR of dagl
    nolot.dagl.into.medi[r] <- log(dagl_tp1/invader.abund)
    
    # invade hygl
    hygl_tp1 <- (1-germination$hygl[p])*survival$hygl[p]*invader.abund + 
      invader.abund*germination$hygl[p]*hygl$lambda[p]/(1+hygl$alpha_medi[p]*lived.medi)
    # calculate LDGR of HYGL
    nolot.hygl.into.medi[r] <- log(hygl_tp1/invader.abund)
    
    # invade peai
    peai_tp1 <- (1-germination$peai[p])*survival$peai[p]*invader.abund + 
      invader.abund*germination$peai[p]*peai$lambda[p]/(1+peai$alpha_medi[p]*lived.medi)
    # calculate LDGR of peai
    nolot.peai.into.medi[r] <- log(peai_tp1/invader.abund)
    
  }

# LDGR INTO PEAI AS RESIDENT ---------------------------
load("nolottery.abundance.peai.rdata")
peai.equilibrum <- abundances.peai[,200]
invader.abund <- 1

# create all of our LDGR matrices
#arca, dagl, hygl, medi, plde, poca, trcy, vero
nolot.arca.into.peai <- matrix(NA, nrow=runs)
nolot.dagl.into.peai <- matrix(NA, nrow=runs)
nolot.hygl.into.peai <- matrix(NA, nrow=runs)
nolot.medi.into.peai <- matrix(NA, nrow=runs)
nolot.plde.into.peai <- matrix(NA, nrow=runs)
nolot.poca.into.peai <- matrix(NA, nrow=runs)
nolot.trcy.into.peai <- matrix(NA, nrow=runs)
nolot.vero.into.peai <- matrix(NA, nrow=runs)


for (r in 1:runs) {
  p <- sample(seq(1, 4500),1) 
    
    # calculate resident abundance
    Nj <- peai.equilibrum[r]
    germ_j <- germination$peai[p]
    #lived.peai <- min(germ_j*Nj, ag_carrying$x[y])
    lived.peai <- germ_j*Nj
    
    # invade arca
    arca_tp1 <- (1-germination$arca[p])*survival$arca[p]*invader.abund + 
      invader.abund*germination$arca[p]*arca$lambda[p]/(1+arca$alpha_peai[p]*lived.peai)
    # calculate LDGR of arca
    nolot.arca.into.peai[r] <- log(arca_tp1/invader.abund)
    
    # invade dagl
    dagl_tp1 <- (1-germination$dagl[p])*survival$dagl[p]*invader.abund + 
      invader.abund*germination$dagl[p]*dagl$lambda[p]/(1+dagl$alpha_peai[p]*lived.peai)
    # calculate LDGR of dagl
    nolot.dagl.into.peai[r] <- log(dagl_tp1/invader.abund)
    
    # invade hygl
    hygl_tp1 <- (1-germination$hygl[p])*survival$hygl[p]*invader.abund + 
      invader.abund*germination$hygl[p]*hygl$lambda[p]/(1+hygl$alpha_peai[p]*lived.peai)
    # calculate LDGR of hygl
    nolot.hygl.into.peai[r] <- log(hygl_tp1/invader.abund)
    
    # invade medi
    medi_tp1 <- (1-germination$medi[p])*survival$medi[p]*invader.abund + 
      invader.abund*germination$medi[p]*medi$lambda[p]/(1+medi$alpha_peai[p]*lived.peai)
    # calculate LDGR of medi
    nolot.medi.into.peai[r] <- log(medi_tp1/invader.abund)
    
    # invade plde
    plde_tp1 <- (1-germination$plde[p])*survival$plde[p]*invader.abund + 
      invader.abund*germination$plde[p]*plde$lambda[p]/(1+plde$alpha_peai[p]*lived.peai)
    # calculate LDGR of plde
    nolot.plde.into.peai[r] <- log(plde_tp1/invader.abund)
    
    # invade poca
    poca_tp1 <- (1-germination$poca[p])*survival$poca[p]*invader.abund + 
      invader.abund*germination$poca[p]*poca$lambda[p]/(1+poca$alpha_peai[p]*lived.peai)
    # calculate LDGR of poca
    nolot.poca.into.peai[r] <- log(poca_tp1/invader.abund)
    
    # invade trcy
    trcy_tp1 <- (1-germination$trcy[p])*survival$trcy[p]*invader.abund + 
      invader.abund*germination$trcy[p]*trcy$lambda[p]/(1+trcy$alpha_peai[p]*lived.peai)
    # calculate LDGR of trcy
    nolot.trcy.into.peai[r] <- log(trcy_tp1/invader.abund)
    
    # invade vero
    vero_tp1 <- (1-germination$vero[p])*survival$vero[p]*invader.abund + 
      invader.abund*germination$vero[p]*vero$lambda[p]/(1+vero$alpha_peai[p]*lived.peai)
    # calculate LDGR of vero
    nolot.vero.into.peai[r] <- log(vero_tp1/invader.abund)
  }


# LDGR INTO PLDE AS RESIDENT ----------------------------
load("nolottery.abundance.plde.rdata")
plde.equilibrum <- abundances.plde[,200]
invader.abund <- 1


# create all of our LDGR matrices
#arca, dagl, medi, peai, poca, trcy
#nolot.arca.into.plde <- matrix(NA, nrow=runs, ncol=posteriors)
nolot.arca.into.plde <- matrix(NA, nrow=runs)
nolot.dagl.into.plde <- matrix(NA, nrow=runs)
nolot.hygl.into.plde <- matrix(NA, nrow=runs)
nolot.peai.into.plde <- matrix(NA, nrow=runs)
nolot.poca.into.plde <- matrix(NA, nrow=runs)
nolot.trcy.into.plde <- matrix(NA, nrow=runs)

for (r in 1:runs) {
  p <- sample(seq(1, 4500),1) 

    # calculate resident abundance
    Nj <- plde.equilibrum[r]
    germ_j <- germination$plde[p]
    #lived.plde <- min(germ_j*Nj, ag_carrying$x[y])
    lived.plde <- germ_j*Nj
    
    # invade arca
    arca_tp1 <- (1-germination$arca[p])*survival$arca[p]*invader.abund + 
      invader.abund*germination$arca[p]*arca$lambda[p]/(1+arca$alpha_plde[p]*lived.plde)
    # calculate LDGR of arca
    nolot.arca.into.plde[r] <- log(arca_tp1/invader.abund)
    
    # invade dagl
    dagl_tp1 <- (1-germination$dagl[p])*survival$dagl[p]*invader.abund + 
      invader.abund*germination$dagl[p]*dagl$lambda[p]/(1+dagl$alpha_plde[p]*lived.plde)
    # calculate LDGR of dagl
    nolot.dagl.into.plde[r] <- log(dagl_tp1/invader.abund)
    
    # invade hygl
    hygl_tp1 <- (1-germination$hygl[p])*survival$hygl[p]*invader.abund + 
      invader.abund*germination$hygl[p]*hygl$lambda[p]/(1+hygl$alpha_plde[p]*lived.plde)
    # calculate LDGR of hygl
    nolot.hygl.into.plde[r] <- log(hygl_tp1/invader.abund)
    
    # invade peai
    peai_tp1 <- (1-germination$peai[p])*survival$peai[p]*invader.abund + 
      invader.abund*germination$peai[p]*peai$lambda[p]/(1+peai$alpha_plde[p]*lived.plde)
    # calculate LDGR of dagl
    nolot.peai.into.plde[r] <- log(peai_tp1/invader.abund)
    
    # invade poca
    poca_tp1 <- (1-germination$poca[p])*survival$poca[p]*invader.abund + 
      invader.abund*germination$poca[p]*poca$lambda[p]/(1+poca$alpha_plde[p]*lived.plde)
    # calculate LDGR of poca
    nolot.poca.into.plde[r] <- log(poca_tp1/invader.abund)
    
    # invade trcy
    trcy_tp1 <- (1-germination$trcy[p])*survival$trcy[p]*invader.abund + 
      invader.abund*germination$trcy[p]*trcy$lambda[p]/(1+trcy$alpha_plde[p]*lived.plde)
    # calculate LDGR of trcy
    nolot.trcy.into.plde[r] <- log(trcy_tp1/invader.abund)
  }

# LDGR INTO POCA RESIDENT -----------------------------
load("nolottery.abundance.poca.rdata")
poca.equilibrum <- abundances.poca[,200]
invader.abund <- 1

# create all of our LDGR matrices
# arca, dagl, medi, peai, plde, trcy, vero
nolot.arca.into.poca <- matrix(NA, nrow=runs)
nolot.dagl.into.poca <- matrix(NA, nrow=runs)
nolot.medi.into.poca <- matrix(NA, nrow=runs)
nolot.peai.into.poca <- matrix(NA, nrow=runs)
nolot.plde.into.poca <- matrix(NA, nrow=runs)
nolot.trcy.into.poca <- matrix(NA, nrow=runs)
nolot.vero.into.poca <- matrix(NA, nrow=runs)


for (r in 1:runs) {
  p <- sample(seq(1, 4500),1) 
    
    # calculate resident abundance
    Nj <- poca.equilibrum[r]
    germ_j <- germination$poca[p]
    #lived.poca <- min(germ_j*Nj, ag_carrying$x[y])
    lived.poca <- germ_j*Nj
    
    # invade arca
    arca_tp1 <- (1-germination$arca[p])*survival$arca[p]*invader.abund + 
      invader.abund*germination$arca[p]*arca$lambda[p]/(1+arca$alpha_poca[p]*lived.poca)
    # calculate LDGR of arca
    nolot.arca.into.poca[r] <- log(arca_tp1/invader.abund)
    
    # invade dagl
    dagl_tp1 <- (1-germination$dagl[p])*survival$dagl[p]*invader.abund + 
      invader.abund*germination$dagl[p]*dagl$lambda[p]/(1+dagl$alpha_poca[p]*lived.poca)
    # calculate LDGR of dagl
    nolot.dagl.into.poca[r] <- log(dagl_tp1/invader.abund)
    
    # invade medi
    medi_tp1 <- (1-germination$medi[p])*survival$medi[p]*invader.abund + 
      invader.abund*germination$medi[p]*medi$lambda[p]/(1+medi$alpha_poca[p]*lived.poca)
    # calculate LDGR of medi
    nolot.medi.into.poca[r] <- log(medi_tp1/invader.abund)
    
    # invade peai
    peai_tp1 <- (1-germination$peai[p])*survival$peai[p]*invader.abund + 
      invader.abund*germination$peai[p]*peai$lambda[p]/(1+peai$alpha_poca[p]*lived.poca)
    # calculate LDGR of dagl
    nolot.peai.into.poca[r] <- log(peai_tp1/invader.abund)
    
    # invade plde
    plde_tp1 <- (1-germination$plde[p])*survival$plde[p]*invader.abund + 
      invader.abund*germination$plde[p]*plde$lambda[p]/(1+plde$alpha_poca[p]*lived.poca)
    # calculate LDGR of plde
    nolot.plde.into.poca[r] <- log(plde_tp1/invader.abund)
    
    # invade trcy
    trcy_tp1 <- (1-germination$trcy[p])*survival$trcy[p]*invader.abund + 
      invader.abund*germination$trcy[p]*trcy$lambda[p]/(1+trcy$alpha_poca[p]*lived.poca)
    # calculate LDGR of trcy
    nolot.trcy.into.poca[r] <- log(trcy_tp1/invader.abund)
    
    # invade vero
    vero_tp1 <- (1-germination$vero[p])*survival$vero[p]*invader.abund + 
      invader.abund*germination$vero[p]*vero$lambda[p]/(1+vero$alpha_poca[p]*lived.poca)
    # calculate LDGR of vero
    nolot.vero.into.poca[r] <- log(vero_tp1/invader.abund)
  }

# LDGR INTO TRCY AS RESIDENT --------------------------
load("nolottery.abundance.trcy.rdata")
trcy.equilibrum <- abundances.trcy[,200]
invader.abund <- 1

# create all of our LDGR matrices
# arca, dagl, peai, plde, poca
nolot.arca.into.trcy <- matrix(NA, nrow=runs)
nolot.dagl.into.trcy <- matrix(NA, nrow=runs)
nolot.peai.into.trcy <- matrix(NA, nrow=runs)
nolot.poca.into.trcy <- matrix(NA, nrow=runs)

for (r in 1:runs) {
  p <- sample(seq(1, 4500),1) 
    # calculate resident abundance
    Nj <- trcy.equilibrum[r]
    germ_j <- germination$trcy[p]
    #lived.trcy <- min(germ_j*Nj, ag_carrying$x[y])
    lived.trcy <- germ_j*Nj
    
    # invade arca
    arca_tp1 <- (1-germination$arca[p])*survival$arca[p]*invader.abund + 
      invader.abund*germination$arca[p]*arca$lambda[p]/(1+arca$alpha_trcy[p]*lived.trcy)
    # calculate LDGR of arca
    nolot.arca.into.trcy[r] <- log(arca_tp1/invader.abund)
    
    # invade dagl
    dagl_tp1 <- (1-germination$dagl[p])*survival$dagl[p]*invader.abund + 
      invader.abund*germination$dagl[p]*dagl$lambda[p]/(1+dagl$alpha_trcy[p]*lived.trcy)
    # calculate LDGR of dagl
    nolot.dagl.into.trcy[r] <- log(dagl_tp1/invader.abund)
    
    # invade peai
    peai_tp1 <- (1-germination$peai[p])*survival$peai[p]*invader.abund + 
      invader.abund*germination$peai[p]*peai$lambda[p]/(1+peai$alpha_trcy[p]*lived.trcy)
    # calculate LDGR of dagl
    nolot.peai.into.trcy[r] <- log(peai_tp1/invader.abund)
    
    # invade poca
    poca_tp1 <- (1-germination$poca[p])*survival$poca[p]*invader.abund + 
      invader.abund*germination$poca[p]*poca$lambda[p]/(1+poca$alpha_trcy[p]*lived.trcy)
    # calculate LDGR of poca
    nolot.poca.into.trcy[r] <- log(poca_tp1/invader.abund)
  }

# LDGR INTO VERO AS RESIDENT --------------------------
load("nolottery.abundance.vero.rdata")
vero.equilibrum <- abundances.vero[,200]
invader.abund <- 1

# create all of our LDGR matrices
# arca, dagl, hygl, medi, peai, plde, poca, trcy, vero
nolot.arca.into.vero <- matrix(NA, nrow=runs)
nolot.dagl.into.vero <- matrix(NA, nrow=runs)
nolot.hygl.into.vero <- matrix(NA, nrow=runs)
nolot.medi.into.vero <- matrix(NA, nrow=runs)
nolot.peai.into.vero <- matrix(NA, nrow=runs)
nolot.plde.into.vero <- matrix(NA, nrow=runs)
nolot.poca.into.vero <- matrix(NA, nrow=runs)
nolot.trcy.into.vero <- matrix(NA, nrow=runs)

for (r in 1:runs) {
  p <- sample(seq(1, 4500),1) 
  
 # calculate resident abundance
    Nj <- vero.equilibrum[r]
    germ_j <- germination$vero[p]
    #lived.vero <- min(germ_j*Nj, ag_carrying$x[y])
    lived.vero <- germ_j*Nj
    
    # invade arca
    arca_tp1 <- (1-germination$arca[p])*survival$arca[p]*invader.abund + 
      invader.abund*germination$arca[p]*arca$lambda[p]/(1+arca$alpha_vero[p]*lived.vero)
    # calculate LDGR of arca
    nolot.arca.into.vero[r] <- log(arca_tp1/invader.abund)
    
    # invade dagl
    dagl_tp1 <- (1-germination$dagl[p])*survival$dagl[p]*invader.abund + 
      invader.abund*germination$dagl[p]*dagl$lambda[p]/(1+dagl$alpha_vero[p]*lived.vero)
    # calculate LDGR of dagl
    nolot.dagl.into.vero[r] <- log(dagl_tp1/invader.abund)
    
    # invade hygl
    hygl_tp1 <- (1-germination$hygl[p])*survival$hygl[p]*invader.abund + 
      invader.abund*germination$hygl[p]*hygl$lambda[p]/(1+hygl$alpha_vero[p]*lived.vero)
    # calculate LDGR of hygl
    nolot.hygl.into.vero[r] <- log(hygl_tp1/invader.abund)
    
    # invade peai
    peai_tp1 <- (1-germination$peai[p])*survival$peai[p]*invader.abund + 
      invader.abund*germination$peai[p]*peai$lambda[p]/(1+peai$alpha_vero[p]*lived.vero)
    # calculate LDGR of dagl
    nolot.peai.into.vero[r] <- log(peai_tp1/invader.abund)
    
    # invade plde
    plde_tp1 <- (1-germination$plde[p])*survival$plde[p]*invader.abund + 
      invader.abund*germination$plde[p]*plde$lambda[p]/(1+plde$alpha_vero[p]*lived.vero)
    # calculate LDGR of plde
    nolot.plde.into.vero[r] <- log(plde_tp1/invader.abund)
    
    # invade poca
    poca_tp1 <- (1-germination$poca[p])*survival$poca[p]*invader.abund + 
      invader.abund*germination$poca[p]*poca$lambda[p]/(1+poca$alpha_vero[p]*lived.vero)
    # calculate LDGR of poca
    nolot.poca.into.vero[r] <- log(poca_tp1/invader.abund)
    
    # invade medi
    medi_tp1 <- (1-germination$medi[p])*survival$medi[p]*invader.abund + 
      invader.abund*germination$medi[p]*medi$lambda[p]/(1+medi$alpha_vero[p]*lived.vero)
    # calculate LDGR of medi
    nolot.medi.into.vero[r] <- log(medi_tp1/invader.abund)
    
    # invade trcy
    trcy_tp1 <- (1-germination$trcy[p])*survival$trcy[p]*invader.abund + 
      invader.abund*germination$trcy[p]*trcy$lambda[p]/(1+trcy$alpha_vero[p]*lived.vero)
    # calculate LDGR of trcy
    nolot.trcy.into.vero[r] <- log(trcy_tp1/invader.abund)
  }

