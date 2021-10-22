# Calculate pair-wise low-density growth rate calculations
library(rstan)

source("Bayes.data.R")

# LDGR into ARCA as resident ---------------------------
load("mean.abundance.arca2.rdata")
arca.equilibrum <- abundances.arca[200]
invader.abund <- 1

# hygl, peai, plde, poca
# create all of our LDGR matrices
mean.hygl.into.arca <- rep(NA, runs)
mean.peai.into.arca <- rep(NA, runs)
mean.plde.into.arca <- rep(NA, runs)
mean.poca.into.arca <- rep(NA, runs)

#y <- sample(dim(ag_carrying)[1], 1)
    
    # calculate resident abundance
    Nj <- arca.equilibrum
    germ_j <- mean(germination$arca)
    lived.arca <- germ_j*Nj
    #lived.arca <- min(germ_j*Nj, ag_carrying$x[y])
    
    # invade HYGL
    hygl_tp1 <- (1-mean(germination$hygl))*mean(survival$hygl)*invader.abund + 
      invader.abund*mean(germination$hygl)*mean(hygl$lambda)/(1+mean(hygl$alpha_arca)*lived.arca)
    # calculate LDGR of HYGL
    mean.hygl.into.arca <- log(hygl_tp1/invader.abund)
    
    # invade PEAI
    peai_tp1 <- (1-mean(germination$peai))*mean(survival$peai)*invader.abund + 
      invader.abund*mean(germination$peai)*mean(peai$lambda)/(1+mean(peai$alpha_arca)*lived.arca)
    # calculate LDGR of PEAI
    
    # invade POCA
    poca_tp1 <- (1-mean(germination$poca))*mean(survival$poca)*invader.abund + 
      invader.abund*mean(germination$poca)*mean(poca$lambda)/(1+mean(poca$alpha_arca)*lived.arca)
    # calculate LDGR of PEAI
    mean.poca.into.arca <- log(poca_tp1/invader.abund)

# LDGR INTO HYGL AS RESIDENT------------------------------------------------
load("mean.abundance.hygl2.rdata")
hygl.equilibrum <- abundances.hygl[200]
invader.abund <- 1

# medi, peai, plde, vero
#posteriors <- 100

# create all of our LDGR matrices
mean.medi.into.hygl <- rep(NA, runs)
mean.peai.into.hygl <- rep(NA, runs)
mean.plde.into.hygl <- rep(NA, runs)
mean.vero.into.hygl <- rep(NA, runs)

#y <- sample(dim(ag_carrying)[1], 1)
    
    # calculate resident abundance
    Nj <- hygl.equilibrum
    germ_j <- mean(germination$hygl)
    lived.hygl <- germ_j*Nj
    #lived.hygl<- min(germ_j*Nj, ag_carrying$x[y])
    
    # invade MEDI
    medi_tp1 <- (1-mean(germination$medi))*mean(survival$medi)*invader.abund + 
      invader.abund*mean(germination$medi)*mean(medi$lambda)/(1+mean(medi$alpha_hygl)*lived.hygl)
    # calculate LDGR of medi
    mean.medi.into.hygl <- log(medi_tp1/invader.abund)
    
    # invade PEAI
    peai_tp1 <- (1-mean(germination$peai))*mean(survival$peai)*invader.abund + 
      invader.abund*mean(germination$peai)*mean(peai$lambda)/(1+mean(peai$alpha_hygl)*lived.hygl)
    # calculate LDGR of peai
    mean.peai.into.hygl <- log(peai_tp1/invader.abund)
    
    # invade PLDE 
    plde_tp1 <- (1-mean(germination$plde))*mean(survival$plde)*invader.abund + 
      invader.abund*mean(germination$plde)*mean(plde$lambda)/(1+mean(plde$alpha_hygl)*lived.hygl)
    # calculate LDGR of plde
    mean.plde.into.hygl <- log(plde_tp1/invader.abund)
    
    # invade VERO
    vero_tp1 <- (1-mean(germination$vero))*mean(survival$vero)*invader.abund + 
      invader.abund*mean(germination$vero)*mean(vero$lambda)/(1+mean(vero$alpha_hygl)*lived.hygl)
    # calculate LDGR of vero
    mean.vero.into.hygl <- log(vero_tp1/invader.abund)

#}

# LDGR INTO MEDI AS RESIDENT ---------------------------------------
load("mean.abundance.medi2.rdata")
medi.equilibrum <- abundances.medi[200]
invader.abund <- 1

# dagl, hygl, peai 
#posteriors <- 100
# create all of our LDGR matrices
mean.dagl.into.medi <- rep(NA, runs)
mean.hygl.into.medi <- rep(NA, runs)
mean.peai.into.medi <- rep(NA, runs)

#y <- sample(dim(ag_carrying)[1], 1)
    
    # calculate resident abundance
    Nj <- medi.equilibrum
    germ_j <- mean(germination$medi)
    lived.medi <- germ_j*Nj
    #lived.medi <- min(germ_j*Nj, ag_carrying$x[y])
    
    # invade dagl
    #dagl_tp1 <- (1-mean(germination$dagl))*mean(survival$dagl)*invader.abund + 
    #  invader.abund*mean(germination$dagl)*mean(dagl$lambda)/(1+mean(dagl$alpha_medi)*lived.medi)
    # calculate LDGR of dagl
    #mean.dagl.into.medi <- log(dagl_tp1/invader.abund)
    
    # invade hygl
    hygl_tp1 <- (1-mean(germination$hygl))*mean(survival$hygl)*invader.abund + 
      invader.abund*mean(germination$hygl)*mean(hygl$lambda)/(1+mean(hygl$alpha_medi)*lived.medi)
    # calculate LDGR of hygl
    mean.hygl.into.medi <- log(hygl_tp1/invader.abund)
    
    # invade peai
    peai_tp1 <- (1-mean(germination$peai))*mean(survival$peai)*invader.abund + 
      invader.abund*mean(germination$peai)*mean(peai$lambda)/(1+mean(peai$alpha_medi)*lived.medi)
    # calculate LDGR of peai
    mean.peai.into.medi <- log(peai_tp1/invader.abund)

# LDGR INTO PEAI AS RESIDENT ---------------------------
load("mean.abundance.peai2.rdata")
peai.equilibrum <- abundances.peai[200]
invader.abund <- 1

# dagl, hygl, peai 

# create all of our LDGR matrices
#arca, dagl, hygl, medi, plde, poca, trcy, vero
mean.arca.into.peai <- rep(NA, runs)
#mean.dagl.into.peai <- rep(NA, runs)
mean.hygl.into.peai <- rep(NA, runs)
mean.medi.into.peai <- rep(NA, runs)
mean.plde.into.peai <- rep(NA, runs)
mean.poca.into.peai <- rep(NA, runs)
mean.trcy.into.peai <- rep(NA, runs)
mean.vero.into.peai <- rep(NA, runs)

#y <- sample(dim(ag_carrying)[1], 1)
    
    # calculate resident abundance
    Nj <- peai.equilibrum
    germ_j <- mean(germination$peai)
    #lived.peai <- min(germ_j*Nj, ag_carrying$x[y])
    lived.peai <- germ_j*Nj
    
    # invade arca
    arca_tp1 <- (1-mean(germination$arca))*mean(survival$arca)*invader.abund + 
      invader.abund*mean(germination$arca)*mean(arca$lambda)/(1+mean(arca$alpha_peai)*lived.peai)
    # calculate LDGR of arca
    mean.arca.into.peai <- log(arca_tp1/invader.abund)
    
    # invade hygl
    hygl_tp1 <- (1-mean(germination$hygl))*mean(survival$hygl)*invader.abund + 
      invader.abund*mean(germination$hygl)*mean(hygl$lambda)/(1+mean(hygl$alpha_peai)*lived.peai)
    # calculate LDGR of hygl
    mean.hygl.into.peai <- log(hygl_tp1/invader.abund)
    
    # invade medi
    medi_tp1 <- (1-mean(germination$medi))*mean(survival$medi)*invader.abund + 
      invader.abund*mean(germination$medi)*mean(medi$lambda)/(1+mean(medi$alpha_peai)*lived.peai)
    # calculate LDGR of medi
    mean.medi.into.peai <- log(medi_tp1/invader.abund)
    
    # invade plde
    plde_tp1 <- (1-mean(germination$plde))*mean(survival$plde)*invader.abund + 
      invader.abund*mean(germination$plde)*mean(plde$lambda)/(1+mean(plde$alpha_peai)*lived.peai)
    # calculate LDGR of plde
    mean.plde.into.peai <- log(plde_tp1/invader.abund)
    
    # invade poca
    poca_tp1 <- (1-mean(germination$poca))*mean(survival$poca)*invader.abund + 
      invader.abund*mean(germination$poca)*mean(poca$lambda)/(1+mean(poca$alpha_peai)*lived.peai)
    # calculate LDGR of poca
    mean.poca.into.peai <- log(poca_tp1/invader.abund)
    
    # invade trcy
    trcy_tp1 <- (1-mean(germination$trcy))*mean(survival$trcy)*invader.abund + 
      invader.abund*mean(germination$trcy)*mean(trcy$lambda)/(1+mean(trcy$alpha_peai)*lived.peai)
    # calculate LDGR of trcy
    mean.trcy.into.peai <- log(trcy_tp1/invader.abund)
    
    # invade vero
    vero_tp1 <- (1-mean(germination$vero))*mean(survival$vero)*invader.abund + 
      invader.abund*mean(germination$vero)*mean(vero$lambda)/(1+mean(vero$alpha_peai)*lived.peai)
    # calculate LDGR of vero
    mean.vero.into.peai <- log(vero_tp1/invader.abund)


# LDGR INTO PLDE AS RESIDENT ----------------------------
load("mean.abundance.plde2.rdata")
plde.equilibrum <- abundances.plde[200]
invader.abund <- 1

# dagl, hygl, peai 

# create all of our LDGR matrices
#arca, dagl, medi, peai, poca, trcy
mean.hygl.into.plde <- rep(NA, runs)
mean.peai.into.plde <- rep(NA, runs)
mean.poca.into.plde <- rep(NA, runs)
mean.trcy.into.plde <- rep(NA, runs)

#y <- sample(dim(ag_carrying)[1], 1)
    
    # calculate resident abundance
    Nj <- plde.equilibrum
    germ_j <- mean(germination$plde)
    #lived.plde <- min(germ_j*Nj, ag_carrying$x[y])
    lived.plde <- germ_j*Nj
    
    # invade arca
    #arca_tp1 <- (1-mean(germination$arca))*mean(survival$arca)*invader.abund + 
    #  invader.abund*mean(germination$arca)*mean(arca$lambda)/(1+mean(arca$alpha_plde)*lived.plde)
    # calculate LDGR of arca
    #mean.arca.into.plde <- log(arca_tp1/invader.abund)
    
    # invade hygl
    hygl_tp1 <- (1-mean(germination$hygl))*mean(survival$hygl)*invader.abund + 
      invader.abund*mean(germination$hygl)*mean(hygl$lambda)/(1+mean(hygl$alpha_plde)*lived.plde)
    # calculate LDGR of hygl
    mean.hygl.into.plde <- log(hygl_tp1/invader.abund)
    
    # invade peai
    peai_tp1 <- (1-mean(germination$peai))*mean(survival$peai)*invader.abund + 
      invader.abund*mean(germination$peai)*mean(peai$lambda)/(1+mean(peai$alpha_plde)*lived.plde)
    # calculate LDGR of peai
    mean.peai.into.plde <- log(peai_tp1/invader.abund)
    
    # invade poca
    poca_tp1 <- (1-mean(germination$poca))*mean(survival$poca)*invader.abund + 
      invader.abund*mean(germination$poca)*mean(poca$lambda)/(1+mean(poca$alpha_plde)*lived.plde)
    # calculate LDGR of poca
    mean.poca.into.plde <- log(poca_tp1/invader.abund)
    
    # invade trcy
    trcy_tp1 <- (1-mean(germination$trcy))*mean(survival$trcy)*invader.abund + 
      invader.abund*mean(germination$trcy)*mean(trcy$lambda)/(1+mean(trcy$alpha_plde)*lived.plde)
    # calculate LDGR of trcy
    mean.trcy.into.plde <- log(trcy_tp1/invader.abund)

# LDGR INTO POCA RESIDENT -----------------------------
load("mean.abundance.poca2.rdata")
poca.equilibrum <- abundances.poca[200]
invader.abund <- 1

# create all of our LDGR matrices
# arca, dagl, medi, peai, plde, trcy, vero
mean.arca.into.poca <- rep(NA, runs)
mean.medi.into.poca <- rep(NA, runs)
mean.peai.into.poca <- rep(NA, runs)
mean.plde.into.poca <- rep(NA, runs)
mean.trcy.into.poca <- rep(NA, runs)
mean.vero.into.poca <- rep(NA, runs)

#y <- sample(dim(ag_carrying)[1], 1)
    
    # calculate resident abundance
    Nj <- poca.equilibrum
    germ_j <- mean(germination$poca)
    #lived.poca <- min(germ_j*Nj, ag_carrying$x[y])
    lived.poca <- germ_j*Nj
    
    # invade arca
    arca_tp1 <- (1-mean(germination$arca))*mean(survival$arca)*invader.abund + 
      invader.abund*mean(germination$arca)*mean(arca$lambda)/(1+mean(arca$alpha_poca)*lived.poca)
    # calculate LDGR of arca
    mean.arca.into.poca <- log(arca_tp1/invader.abund)
    
    # invade medi
    medi_tp1 <- (1-mean(germination$medi))*mean(survival$medi)*invader.abund + 
      invader.abund*mean(germination$medi)*mean(medi$lambda)/(1+mean(medi$alpha_poca)*lived.poca)
    # calculate LDGR of medi
    mean.medi.into.poca <- log(medi_tp1/invader.abund)
    
    # invade peai
    peai_tp1 <- (1-mean(germination$peai))*mean(survival$peai)*invader.abund + 
      invader.abund*mean(germination$peai)*mean(peai$lambda)/(1+mean(peai$alpha_poca)*lived.poca)
    # calculate LDGR of peai
    mean.peai.into.poca <- log(peai_tp1/invader.abund)
    
    # invade plde
    plde_tp1 <- (1-mean(germination$plde))*mean(survival$plde)*invader.abund + 
      invader.abund*mean(germination$plde)*mean(plde$lambda)/(1+mean(plde$alpha_poca)*lived.poca)
    # calculate LDGR of plde
    mean.plde.into.poca <- log(plde_tp1/invader.abund)
    
    # invade trcy
    trcy_tp1 <- (1-mean(germination$trcy))*mean(survival$trcy)*invader.abund + 
      invader.abund*mean(germination$trcy)*mean(trcy$lambda)/(1+mean(trcy$alpha_poca)*lived.poca)
    # calculate LDGR of trcy
    mean.trcy.into.poca <- log(trcy_tp1/invader.abund)
    
    # invade vero
    vero_tp1 <- (1-mean(germination$vero))*mean(survival$vero)*invader.abund + 
      invader.abund*mean(germination$vero)*mean(vero$lambda)/(1+mean(vero$alpha_poca)*lived.poca)
    # calculate LDGR of vero
    mean.vero.into.poca <- log(vero_tp1/invader.abund)


# LDGR INTO TRCY AS RESIDENT --------------------------
load("mean.abundance.trcy2.rdata")
trcy.equilibrum <- abundances.trcy[200]
invader.abund <- 1

# create all of our LDGR matrices
# arca, dagl, peai, plde, poca
mean.arca.into.trcy <- rep(NA, runs)
mean.peai.into.trcy <- rep(NA, runs)
mean.poca.into.trcy <- rep(NA, runs)

#y <- sample(dim(ag_carrying)[1], 1)
    
    # calculate resident abundance
    Nj <- trcy.equilibrum
    germ_j <- mean(germination$trcy)
    #lived.trcy <- min(germ_j*Nj, ag_carrying$x[y])
    lived.trcy <- germ_j*Nj
    
    # invade arca
    arca_tp1 <- (1-mean(germination$arca))*mean(survival$arca)*invader.abund + 
      invader.abund*mean(germination$arca)*mean(arca$lambda)/(1+mean(arca$alpha_trcy)*lived.trcy)
    # calculate LDGR of arca
    mean.arca.into.trcy <- log(arca_tp1/invader.abund)
    
    # invade peai
    peai_tp1 <- (1-mean(germination$peai))*mean(survival$peai)*invader.abund + 
      invader.abund*mean(germination$peai)*mean(peai$lambda)/(1+mean(peai$alpha_trcy)*lived.trcy)
    # calculate LDGR of peai
    mean.peai.into.trcy <- log(peai_tp1/invader.abund)
    
    # invade poca
    poca_tp1 <- (1-mean(germination$poca))*mean(survival$poca)*invader.abund + 
      invader.abund*mean(germination$poca)*mean(poca$lambda)/(1+mean(poca$alpha_trcy)*lived.trcy)
    # calculate LDGR of poca
    mean.poca.into.trcy <- log(poca_tp1/invader.abund)

# LDGR INTO VERO AS RESIDENT --------------------------
load("mean.abundance.vero2.rdata")
vero.equilibrum <- abundances.vero[200]
invader.abund <- 1


# create all of our LDGR matrices
# arca, dagl, hygl, medi, peai, plde, poca, trcy, vero
mean.arca.into.vero <- rep(NA, runs)
mean.hygl.into.vero <- rep(NA, runs)
mean.medi.into.vero <- rep(NA, runs)
mean.peai.into.vero <- rep(NA, runs)
mean.plde.into.vero <- rep(NA, runs)
mean.poca.into.vero <- rep(NA, runs)
mean.trcy.into.vero <- rep(NA, runs)

#y <- sample(dim(ag_carrying)[1], 1)
    
    # calculate resident abundance
    Nj <- vero.equilibrum
    germ_j <- mean(germination$vero)
    #lived.vero <- min(germ_j*Nj, ag_carrying$x[y])
    lived.vero <- germ_j*Nj
    
    # invade arca
    arca_tp1 <- (1-mean(germination$arca))*mean(survival$arca)*invader.abund + 
      invader.abund*mean(germination$arca)*mean(arca$lambda)/(1+mean(arca$alpha_vero)*lived.vero)
    # calculate LDGR of arca
    mean.arca.into.vero <- log(arca_tp1/invader.abund)
    
    # invade hygl
    hygl_tp1 <- (1-mean(germination$hygl))*mean(survival$hygl)*invader.abund + 
      invader.abund*mean(germination$hygl)*mean(hygl$lambda)/(1+mean(hygl$alpha_vero)*lived.vero)
    # calculate LDGR of hygl
    mean.hygl.into.vero <- log(hygl_tp1/invader.abund)
    
    # invade peai
    peai_tp1 <- (1-mean(germination$peai))*mean(survival$peai)*invader.abund + 
      invader.abund*mean(germination$peai)*mean(peai$lambda)/(1+mean(peai$alpha_vero)*lived.vero)
    # calculate LDGR of peai
    mean.peai.into.vero <- log(peai_tp1/invader.abund)
    
    # invade plde
    plde_tp1 <- (1-mean(germination$plde))*mean(survival$plde)*invader.abund + 
      invader.abund*mean(germination$plde)*mean(plde$lambda)/(1+mean(plde$alpha_vero)*lived.vero)
    # calculate LDGR of plde
    mean.plde.into.vero <- log(plde_tp1/invader.abund)
    
    # invade poca
    poca_tp1 <- (1-mean(germination$poca))*mean(survival$poca)*invader.abund + 
      invader.abund*mean(germination$poca)*mean(poca$lambda)/(1+mean(poca$alpha_vero)*lived.vero)
    # calculate LDGR of poca
    mean.poca.into.vero <- log(poca_tp1/invader.abund)
    
    # invade medi
    medi_tp1 <- (1-mean(germination$medi))*mean(survival$medi)*invader.abund + 
      invader.abund*mean(germination$medi)*mean(medi$lambda)/(1+mean(medi$alpha_vero)*lived.vero)
    # calculate LDGR of medi
    mean.medi.into.vero <- log(medi_tp1/invader.abund)
    
    # invade trcy
    trcy_tp1 <- (1-mean(germination$trcy))*mean(survival$trcy)*invader.abund + 
      invader.abund*mean(germination$trcy)*mean(trcy$lambda)/(1+mean(trcy$alpha_vero)*lived.vero)
    # calculate LDGR of trcy
    mean.trcy.into.vero <- log(trcy_tp1/invader.abund)


