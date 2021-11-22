##############################################################################################
# This Code simulates pair-wise low-density growth rates where the posterior values have not 
# been sampled through time 
##############################################################################################

sp <- 9

sp_list <- c("dagl", "hygl", "plde", 
             "poca", "trcy", "vero", "arca", "medi", "peai")

# load in Stan model fits
source("Bayes.data.R")
# load in resident species population where the posterior values have not been sampled through time
source("Residents.to.equilibrium.no.time.var.R")

# LDGR into ARCA as resident ---------------------------
arca.equilibrum <- abundances.arca.novar[,200]
invader.abund <- 1 # start invader abundance at 1 (seed)

# hygl, peai, plde, poca can all invade into arca (have enough data)

runs=4500

# create all of our LDGR matrices
hygl.into.arca <- matrix(NA, nrow=runs)
peai.into.arca <- matrix(NA, nrow=runs)
poca.into.arca <- matrix(NA, nrow=runs)


for (r in 1:runs) {
  #p <- sample(seq(1, 4500),1)
  p <- r
  # calculate resident abundance
  Nj <- arca.equilibrum[r]
  germ_j <- germination$arca[p]
  #lived.arca <- min(germ_j*Nj, ag_carrying$x[y])
  lived.arca <- germ_j*Nj
  
  # invade HYGL
  hygl_tp1 <- (1-germination$hygl[p])*survival$hygl[p]*invader.abund + 
    invader.abund*germination$hygl[p]*hygl$lambda[p]/(1+hygl$alpha_arca[p]*lived.arca)
  # calculate LDGR of HYGL
  hygl.into.arca[r] <- log(hygl_tp1/invader.abund)
  
  # invade PEAI
  peai_tp1 <- (1-germination$peai[p])*survival$peai[p]*invader.abund + 
    invader.abund*germination$peai[p]*peai$lambda[p]/(1+peai$alpha_arca[p]*lived.arca)
  # calculate LDGR of PEAI
  peai.into.arca[r] <- log(peai_tp1/invader.abund)
  
  # invade PLDE
  #plde_tp1 <- (1-germination$plde[x])*survival$plde[x]*invader.abund + 
  #  invader.abund*germination$plde[x]*plde$lambda[x]/(1+plde$alpha_arca[x]*lived.arca)
  # calculate LDGR of PLDE
  #plde.into.arca[r,p] <- log(plde_tp1/invader.abund)
  
  # invade POCA
  poca_tp1 <- (1-germination$poca[p])*survival$poca[p]*invader.abund + 
    invader.abund*germination$poca[p]*poca$lambda[p]/(1+poca$alpha_arca[p]*lived.arca)
  # calculate LDGR of POCA
  poca.into.arca[r] <- log(poca_tp1/invader.abund)
  
  # ratio approach 
  
}


# LDGR INTO HYGL AS RESIDENT------------------------------------------------
hygl.equilibrum <- abundances.hygl.novar[,200]
invader.abund <- 1

# create all of our LDGR matrices
medi.into.hygl <- matrix(NA, nrow=runs)
peai.into.hygl <- matrix(NA, nrow=runs)
plde.into.hygl <- matrix(NA, nrow=runs)
vero.into.hygl <- matrix(NA, nrow=runs)

for (r in 1:runs) {
  #p <- sample(seq(1, 4500),1) 
  p<-r
  # calculate resident abundance
  Nj <- hygl.equilibrum[r]
  germ_j <- germination$hygl[p]
  #lived.hygl<- min(germ_j*Nj, ag_carrying$x[y])
  lived.hygl<-germ_j*Nj
  
  # invade MEDI
  medi_tp1 <- (1-germination$medi[p])*survival$medi[p]*invader.abund + 
    invader.abund*germination$medi[p]*medi$lambda[p]/(1+medi$alpha_hygl[p]*lived.hygl)
  # calculate LDGR of medi
  medi.into.hygl[r] <- log(medi_tp1/invader.abund)
  
  # invade PEAI
  peai_tp1 <- (1-germination$peai[p])*survival$peai[p]*invader.abund + 
    invader.abund*germination$peai[p]*peai$lambda[p]/(1+peai$alpha_hygl[p]*lived.hygl)
  # calculate LDGR of peai
  peai.into.hygl[r] <- log(peai_tp1/invader.abund)
  
  # invade PLDE 
  plde_tp1 <- (1-germination$plde[p])*survival$plde[p]*invader.abund + 
    invader.abund*germination$plde[p]*plde$lambda[p]/(1+plde$alpha_hygl[p]*lived.hygl)
  # calculate LDGR of plde
  plde.into.hygl[r] <- log(plde_tp1/invader.abund)
  
  # invade VERO
  vero_tp1 <- (1-germination$vero[p])*survival$vero[p]*invader.abund + 
    invader.abund*germination$vero[p]*vero$lambda[p]/(1+vero$alpha_hygl[p]*lived.hygl)
  # calculate LDGR of vero
  vero.into.hygl[r] <- log(vero_tp1/invader.abund)
}


# LDGR INTO MEDI AS RESIDENT ---------------------------------------
medi.equilibrum <- abundances.medi.novar[,200]
invader.abund <- 1


# create all of our LDGR matrices
#dagl.into.medi <- matrix(NA, nrow=runs)
hygl.into.medi <- matrix(NA, nrow=runs)
peai.into.medi <- matrix(NA, nrow=runs)

for (r in 1:runs) {
  #p <- sample(seq(1, 4500),1) 
  p <- r
  # calculate resident abundance
  Nj <- medi.equilibrum[r]
  germ_j <- germination$medi[p]
  #lived.medi <- min(germ_j*Nj, ag_carrying$x[y])
  lived.medi <- germ_j*Nj
  
  # invade dagl
  #dagl_tp1 <- (1-germination$dagl[p])*survival$dagl[p]*invader.abund + 
  #  invader.abund*germination$dagl[p]*dagl$lambda[p]/(1+dagl$alpha_medi[p]*lived.medi)
  # calculate LDGR of dagl
  #dagl.into.medi[r] <- log(dagl_tp1/invader.abund)
  
  # invade hygl
  hygl_tp1 <- (1-germination$hygl[p])*survival$hygl[p]*invader.abund + 
    invader.abund*germination$hygl[p]*hygl$lambda[p]/(1+hygl$alpha_medi[p]*lived.medi)
  # calculate LDGR of HYGL
  hygl.into.medi[r] <- log(hygl_tp1/invader.abund)
  
  # invade peai
  peai_tp1 <- (1-germination$peai[p])*survival$peai[p]*invader.abund + 
    invader.abund*germination$peai[p]*peai$lambda[p]/(1+peai$alpha_medi[p]*lived.medi)
  # calculate LDGR of peai
  peai.into.medi[r] <- log(peai_tp1/invader.abund)
  
}

# LDGR INTO PEAI AS RESIDENT ---------------------------
peai.equilibrum <- abundances.peai.novar[,200]
invader.abund <- 1

# create all of our LDGR matrices
#arca, dagl, hygl, medi, plde, poca, trcy, vero
arca.into.peai <- matrix(NA, nrow=runs)
#dagl.into.peai <- matrix(NA, nrow=runs)
hygl.into.peai <- matrix(NA, nrow=runs)
medi.into.peai <- matrix(NA, nrow=runs)
plde.into.peai <- matrix(NA, nrow=runs)
poca.into.peai <- matrix(NA, nrow=runs)
trcy.into.peai <- matrix(NA, nrow=runs)
vero.into.peai <- matrix(NA, nrow=runs)


for (r in 1:runs) {
  #p <- sample(seq(1, 4500),1) 
  p <- r
  # calculate resident abundance
  Nj <- peai.equilibrum[r]
  germ_j <- germination$peai[p]
  #lived.peai <- min(germ_j*Nj, ag_carrying$x[y])
  lived.peai <- germ_j*Nj
  
  # invade arca
  arca_tp1 <- (1-germination$arca[p])*survival$arca[p]*invader.abund + 
    invader.abund*germination$arca[p]*arca$lambda[p]/(1+arca$alpha_peai[p]*lived.peai)
  # calculate LDGR of arca
  arca.into.peai[r] <- log(arca_tp1/invader.abund)
  
  # invade dagl
  #dagl_tp1 <- (1-germination$dagl[p])*survival$dagl[p]*invader.abund + 
  #  invader.abund*germination$dagl[p]*dagl$lambda[p]/(1+dagl$alpha_peai[p]*lived.peai)
  # calculate LDGR of dagl
  #dagl.into.peai[r] <- log(dagl_tp1/invader.abund)
  
  # invade hygl
  hygl_tp1 <- (1-germination$hygl[p])*survival$hygl[p]*invader.abund + 
    invader.abund*germination$hygl[p]*hygl$lambda[p]/(1+hygl$alpha_peai[p]*lived.peai)
  # calculate LDGR of hygl
  hygl.into.peai[r] <- log(hygl_tp1/invader.abund)
  
  # invade medi
  medi_tp1 <- (1-germination$medi[p])*survival$medi[p]*invader.abund + 
    invader.abund*germination$medi[p]*medi$lambda[p]/(1+medi$alpha_peai[p]*lived.peai)
  # calculate LDGR of medi
  medi.into.peai[r] <- log(medi_tp1/invader.abund)
  
  # invade plde
  plde_tp1 <- (1-germination$plde[p])*survival$plde[p]*invader.abund + 
    invader.abund*germination$plde[p]*plde$lambda[p]/(1+plde$alpha_peai[p]*lived.peai)
  # calculate LDGR of plde
  plde.into.peai[r] <- log(plde_tp1/invader.abund)
  
  # invade poca
  poca_tp1 <- (1-germination$poca[p])*survival$poca[p]*invader.abund + 
    invader.abund*germination$poca[p]*poca$lambda[p]/(1+poca$alpha_peai[p]*lived.peai)
  # calculate LDGR of poca
  poca.into.peai[r] <- log(poca_tp1/invader.abund)
  
  # invade trcy
  trcy_tp1 <- (1-germination$trcy[p])*survival$trcy[p]*invader.abund + 
    invader.abund*germination$trcy[p]*trcy$lambda[p]/(1+trcy$alpha_peai[p]*lived.peai)
  # calculate LDGR of trcy
  trcy.into.peai[r] <- log(trcy_tp1/invader.abund)
  
  # invade vero
  vero_tp1 <- (1-germination$vero[p])*survival$vero[p]*invader.abund + 
    invader.abund*germination$vero[p]*vero$lambda[p]/(1+vero$alpha_peai[p]*lived.peai)
  # calculate LDGR of vero
  vero.into.peai[r] <- log(vero_tp1/invader.abund)
}


# LDGR INTO PLDE AS RESIDENT ----------------------------
plde.equilibrum <- abundances.plde.novar[,200]
invader.abund <- 1


# create all of our LDGR matrices
#arca, dagl, medi, peai, poca, trcy
#arca.into.plde <- matrix(NA, nrow=runs, ncol=posteriors)
arca.into.plde <- matrix(NA, nrow=runs)
#dagl.into.plde <- matrix(NA, nrow=runs)
hygl.into.plde <- matrix(NA, nrow=runs)
peai.into.plde <- matrix(NA, nrow=runs)
poca.into.plde <- matrix(NA, nrow=runs)
trcy.into.plde <- matrix(NA, nrow=runs)

for (r in 1:runs) {
  #p <- sample(seq(1, 4500),1) 
  p <- r
  # calculate resident abundance
  Nj <- plde.equilibrum[r]
  germ_j <- germination$plde[p]
  #lived.plde <- min(germ_j*Nj, ag_carrying$x[y])
  lived.plde <- germ_j*Nj
  
  # invade arca
  arca_tp1 <- (1-germination$arca[p])*survival$arca[p]*invader.abund + 
    invader.abund*germination$arca[p]*arca$lambda[p]/(1+arca$alpha_plde[p]*lived.plde)
  # calculate LDGR of arca
  arca.into.plde[r] <- log(arca_tp1/invader.abund)
  
  # invade dagl
  #dagl_tp1 <- (1-germination$dagl[p])*survival$dagl[p]*invader.abund + 
  #  invader.abund*germination$dagl[p]*dagl$lambda[p]/(1+dagl$alpha_plde[p]*lived.plde)
  # calculate LDGR of dagl
  #dagl.into.plde[r] <- log(dagl_tp1/invader.abund)
  
  # invade hygl
  hygl_tp1 <- (1-germination$hygl[p])*survival$hygl[p]*invader.abund + 
    invader.abund*germination$hygl[p]*hygl$lambda[p]/(1+hygl$alpha_plde[p]*lived.plde)
  # calculate LDGR of hygl
  hygl.into.plde[r] <- log(hygl_tp1/invader.abund)
  
  # invade peai
  peai_tp1 <- (1-germination$peai[p])*survival$peai[p]*invader.abund + 
    invader.abund*germination$peai[p]*peai$lambda[p]/(1+peai$alpha_plde[p]*lived.plde)
  # calculate LDGR of dagl
  peai.into.plde[r] <- log(peai_tp1/invader.abund)
  
  # invade poca
  poca_tp1 <- (1-germination$poca[p])*survival$poca[p]*invader.abund + 
    invader.abund*germination$poca[p]*poca$lambda[p]/(1+poca$alpha_plde[p]*lived.plde)
  # calculate LDGR of poca
  poca.into.plde[r] <- log(poca_tp1/invader.abund)
  
  # invade trcy
  trcy_tp1 <- (1-germination$trcy[p])*survival$trcy[p]*invader.abund + 
    invader.abund*germination$trcy[p]*trcy$lambda[p]/(1+trcy$alpha_plde[p]*lived.plde)
  # calculate LDGR of trcy
  trcy.into.plde[r] <- log(trcy_tp1/invader.abund)
}

# LDGR INTO POCA RESIDENT -----------------------------
poca.equilibrum <- abundances.poca.novar[,200]
invader.abund <- 1

# create all of our LDGR matrices
# arca, dagl, medi, peai, plde, trcy, vero
arca.into.poca <- matrix(NA, nrow=runs)
#dagl.into.poca <- matrix(NA, nrow=runs)
medi.into.poca <- matrix(NA, nrow=runs)
peai.into.poca <- matrix(NA, nrow=runs)
plde.into.poca <- matrix(NA, nrow=runs)
trcy.into.poca <- matrix(NA, nrow=runs)
vero.into.poca <- matrix(NA, nrow=runs)


for (r in 1:runs) {
  #p <- sample(seq(1, 4500),1) 
  p <- r
  # calculate resident abundance
  Nj <- poca.equilibrum[r]
  germ_j <- germination$poca[p]
  #lived.poca <- min(germ_j*Nj, ag_carrying$x[y])
  lived.poca <- germ_j*Nj
  
  # invade arca
  arca_tp1 <- (1-germination$arca[p])*survival$arca[p]*invader.abund + 
    invader.abund*germination$arca[p]*arca$lambda[p]/(1+arca$alpha_poca[p]*lived.poca)
  # calculate LDGR of arca
  arca.into.poca[r] <- log(arca_tp1/invader.abund)
  
  # invade dagl
  #dagl_tp1 <- (1-germination$dagl[p])*survival$dagl[p]*invader.abund + 
  #  invader.abund*germination$dagl[p]*dagl$lambda[p]/(1+dagl$alpha_poca[p]*lived.poca)
  # calculate LDGR of dagl
  #dagl.into.poca[r] <- log(dagl_tp1/invader.abund)
  
  # invade medi
  medi_tp1 <- (1-germination$medi[p])*survival$medi[p]*invader.abund + 
    invader.abund*germination$medi[p]*medi$lambda[p]/(1+medi$alpha_poca[p]*lived.poca)
  # calculate LDGR of medi
  medi.into.poca[r] <- log(medi_tp1/invader.abund)
  
  # invade peai
  peai_tp1 <- (1-germination$peai[p])*survival$peai[p]*invader.abund + 
    invader.abund*germination$peai[p]*peai$lambda[p]/(1+peai$alpha_poca[p]*lived.poca)
  # calculate LDGR of dagl
  peai.into.poca[r] <- log(peai_tp1/invader.abund)
  
  # invade plde
  plde_tp1 <- (1-germination$plde[p])*survival$plde[p]*invader.abund + 
    invader.abund*germination$plde[p]*plde$lambda[p]/(1+plde$alpha_poca[p]*lived.poca)
  # calculate LDGR of plde
  plde.into.poca[r] <- log(plde_tp1/invader.abund)
  
  # invade trcy
  trcy_tp1 <- (1-germination$trcy[p])*survival$trcy[p]*invader.abund + 
    invader.abund*germination$trcy[p]*trcy$lambda[p]/(1+trcy$alpha_poca[p]*lived.poca)
  # calculate LDGR of trcy
  trcy.into.poca[r] <- log(trcy_tp1/invader.abund)
  
  # invade vero
  vero_tp1 <- (1-germination$vero[p])*survival$vero[p]*invader.abund + 
    invader.abund*germination$vero[p]*vero$lambda[p]/(1+vero$alpha_poca[p]*lived.poca)
  # calculate LDGR of vero
  vero.into.poca[r] <- log(vero_tp1/invader.abund)
}

# LDGR INTO TRCY AS RESIDENT --------------------------
trcy.equilibrum <- abundances.trcy.novar[,200]
invader.abund <- 1

# create all of our LDGR matrices
# arca, dagl, peai, plde, poca
arca.into.trcy <- matrix(NA, nrow=runs)
#dagl.into.trcy <- matrix(NA, nrow=runs)
peai.into.trcy <- matrix(NA, nrow=runs)
poca.into.trcy <- matrix(NA, nrow=runs)

for (r in 1:runs) {
  #p <- sample(seq(1, 4500),1) 
  p <- r
  # calculate resident abundance
  Nj <- trcy.equilibrum[r]
  germ_j <- germination$trcy[p]
  #lived.trcy <- min(germ_j*Nj, ag_carrying$x[y])
  lived.trcy <- germ_j*Nj
  
  # invade arca
  arca_tp1 <- (1-germination$arca[p])*survival$arca[p]*invader.abund + 
    invader.abund*germination$arca[p]*arca$lambda[p]/(1+arca$alpha_trcy[p]*lived.trcy)
  # calculate LDGR of arca
  arca.into.trcy[r] <- log(arca_tp1/invader.abund)
  
  # invade dagl
  #dagl_tp1 <- (1-germination$dagl[p])*survival$dagl[p]*invader.abund + 
  #  invader.abund*germination$dagl[p]*dagl$lambda[p]/(1+dagl$alpha_trcy[p]*lived.trcy)
  # calculate LDGR of dagl
  #dagl.into.trcy[r] <- log(dagl_tp1/invader.abund)
  
  # invade peai
  peai_tp1 <- (1-germination$peai[p])*survival$peai[p]*invader.abund + 
    invader.abund*germination$peai[p]*peai$lambda[p]/(1+peai$alpha_trcy[p]*lived.trcy)
  # calculate LDGR of dagl
  peai.into.trcy[r] <- log(peai_tp1/invader.abund)
  
  # invade poca
  poca_tp1 <- (1-germination$poca[p])*survival$poca[p]*invader.abund + 
    invader.abund*germination$poca[p]*poca$lambda[p]/(1+poca$alpha_trcy[p]*lived.trcy)
  # calculate LDGR of poca
  poca.into.trcy[r] <- log(poca_tp1/invader.abund)
}

# LDGR INTO VERO AS RESIDENT --------------------------
vero.equilibrum <- abundances.vero.novar[,200]
invader.abund <- 1

# create all of our LDGR matrices
# arca, dagl, hygl, medi, peai, plde, poca, trcy, vero
arca.into.vero <- matrix(NA, nrow=runs)
#dagl.into.vero <- matrix(NA, nrow=runs)
hygl.into.vero <- matrix(NA, nrow=runs)
medi.into.vero <- matrix(NA, nrow=runs)
peai.into.vero <- matrix(NA, nrow=runs)
plde.into.vero <- matrix(NA, nrow=runs)
poca.into.vero <- matrix(NA, nrow=runs)
trcy.into.vero <- matrix(NA, nrow=runs)

for (r in 1:runs) {
  #p <- sample(seq(1, 4500),1) 
  p <- r
  # calculate resident abundance
  Nj <- vero.equilibrum[r]
  germ_j <- germination$vero[p]
  #lived.vero <- min(germ_j*Nj, ag_carrying$x[y])
  lived.vero <- germ_j*Nj
  
  # invade arca
  arca_tp1 <- (1-germination$arca[p])*survival$arca[p]*invader.abund + 
    invader.abund*germination$arca[p]*arca$lambda[p]/(1+arca$alpha_vero[p]*lived.vero)
  # calculate LDGR of arca
  arca.into.vero[r] <- log(arca_tp1/invader.abund)
  
  # invade dagl
  #dagl_tp1 <- (1-germination$dagl[p])*survival$dagl[p]*invader.abund + 
  #  invader.abund*germination$dagl[p]*dagl$lambda[p]/(1+dagl$alpha_vero[p]*lived.vero)
  # calculate LDGR of dagl
  #dagl.into.vero[r] <- log(dagl_tp1/invader.abund)
  
  # invade hygl
  hygl_tp1 <- (1-germination$hygl[p])*survival$hygl[p]*invader.abund + 
    invader.abund*germination$hygl[p]*hygl$lambda[p]/(1+hygl$alpha_vero[p]*lived.vero)
  # calculate LDGR of hygl
  hygl.into.vero[r] <- log(hygl_tp1/invader.abund)
  
  # invade peai
  peai_tp1 <- (1-germination$peai[p])*survival$peai[p]*invader.abund + 
    invader.abund*germination$peai[p]*peai$lambda[p]/(1+peai$alpha_vero[p]*lived.vero)
  # calculate LDGR of dagl
  peai.into.vero[r] <- log(peai_tp1/invader.abund)
  
  # invade plde
  plde_tp1 <- (1-germination$plde[p])*survival$plde[p]*invader.abund + 
    invader.abund*germination$plde[p]*plde$lambda[p]/(1+plde$alpha_vero[p]*lived.vero)
  # calculate LDGR of plde
  plde.into.vero[r] <- log(plde_tp1/invader.abund)
  
  # invade poca
  poca_tp1 <- (1-germination$poca[p])*survival$poca[p]*invader.abund + 
    invader.abund*germination$poca[p]*poca$lambda[p]/(1+poca$alpha_vero[p]*lived.vero)
  # calculate LDGR of poca
  poca.into.vero[r] <- log(poca_tp1/invader.abund)
  
  # invade medi
  medi_tp1 <- (1-germination$medi[p])*survival$medi[p]*invader.abund + 
    invader.abund*germination$medi[p]*medi$lambda[p]/(1+medi$alpha_vero[p]*lived.vero)
  # calculate LDGR of medi
  medi.into.vero[r] <- log(medi_tp1/invader.abund)
  
  # invade trcy
  trcy_tp1 <- (1-germination$trcy[p])*survival$trcy[p]*invader.abund + 
    invader.abund*germination$trcy[p]*trcy$lambda[p]/(1+trcy$alpha_vero[p]*lived.vero)
  # calculate LDGR of trcy
  trcy.into.vero[r] <- log(trcy_tp1/invader.abund)
}

