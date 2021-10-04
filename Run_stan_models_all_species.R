################################################################################
# Code to run stan models to estimate lambda and alpha values for each species #
################################################################################
#install.packages("rstan")
library("rstan")
library("coda")
# This line of code detects the number of cores available to your machine so that
#    chains can be run in parallel. If you will be doing other things on your 
#    computer while the model is running, I suggest not running this line of code
#    as it can slow down other applications on your computer
#options(mc.cores = parallel::detectCores())

# This line allows a compiled model to be saved directly, so that stan does
#    not need to recompile it if it hasn't changed which can save time
#rstan_options(auto_write = TRUE)

# ARCA --------------------------------------------------------------------------
# input data and subset to respose data (total.fecundity)
SpData <- read.csv("Data/groups.ARCA.csv")
SpData <- subset(SpData, select = c(total.fecundity, plot, intra, Daucus.glochidiatus, Hyalosperma.glutinosum, 
                                    Plantago.debilis, Podolepis.canescens, Trachymene.cyanopetala, Velleia.rosea,
                                    Medicago.minima, Pentameris.airoides, other))
colSums(SpData[,-c(1:2)])
SpData <- na.omit(SpData)
# make variables for stan model
N <- as.integer(dim(SpData)[1]) # number of observations
P <- as.integer(16) # number of plots 
Fecundity <- as.integer(SpData$total.fecundity)
# check which species don't co-occur above the threshold of 3 obvservations (total, regardless of abundances)
sum(as.integer(SpData$Daucus.glochidiatus != 0))
sum(as.integer(SpData$Hyalosperma.glutinosum != 0))
sum(as.integer(SpData$Plantago.debilis != 0))
sum(as.integer(SpData$Podolepis.canescens != 0))
sum(as.integer(SpData$Trachymene.cyanopetala != 0))
sum(as.integer(SpData$Velleia.rosea != 0))
sum(as.integer(SpData$Medicago.minima != 0))
sum(as.integer(SpData$Pentameris.airoides != 0))
sum(as.integer(SpData$other != 0))
# and remove them
SpData$other <- SpData$other + SpData$Daucus.glochidiatus + SpData$Hyalosperma.glutinosum + SpData$Medicago.minima 
SpData <- subset(SpData, select = c(total.fecundity, plot, intra, Plantago.debilis, Podolepis.canescens, Trachymene.cyanopetala, Velleia.rosea, Pentameris.airoides,
                                    other))
# Create the model matrix of neigbour abundances
ModelMatrix <- subset(SpData, select = c(intra, Plantago.debilis, Podolepis.canescens,Trachymene.cyanopetala, Velleia.rosea, Pentameris.airoides,
                                         other))

ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(SpData, select = plot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}

#initials <- list(epsilon=rep(1,P), sigma = .01)
#initials1<- list(initials, initials, initials)
# run the stan model 
PrelimFit <- stan(file = 'BH stan models/ARCA.stan', data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 6000, chains = 3, thin = 2) 
PrelimFit
save(PrelimFit, file = "BH stan models/Arca_posteriors.rdata")

arca.posts <- PrelimFit
plot(arca.posts, pars=c("alpha_intra", "alpha_plde", "alpha_poca", "alpha_trcy", "alpha_vero",  "alpha_peai",
                        "alpha_other"), show_density = TRUE, ci_level = 0.5, outer_level = 0.95, fill_color = "grey") + ggtitle(expression(italic("A. calendula")))
plot(arca.posts, pars=c("lambda"), show_density = TRUE, ci_level = 0.5, outer_level = 0.95, fill_color = "grey") + ggtitle(expression(italic("A. calendula")))


# HYGL ----------------------------------------------
SpData <- read.csv("Data/groups.HYGL.csv")
SpData <- subset(SpData, select = c(total.fecundity, plot, intra, Arctotheca.calendula, Daucus.glochidiatus,
                                    Plantago.debilis, Podolepis.canescens, Trachymene.cyanopetala, Velleia.rosea,
                                    Medicago.minima, Pentameris.airoides, other))

SpData <- na.omit(SpData)
N <- as.integer(dim(SpData)[1])
P <- as.integer(16)
Fecundity <- as.integer(SpData$total.fecundity)

sum(as.integer(SpData$Daucus.glochidiatus != 0))
sum(as.integer(SpData$Arctotheca.calendula != 0))
sum(as.integer(SpData$Plantago.debilis != 0))
sum(as.integer(SpData$Podolepis.canescens != 0))
sum(as.integer(SpData$Trachymene.cyanopetala != 0))
sum(as.integer(SpData$Velleia.rosea != 0))
sum(as.integer(SpData$Medicago.minima != 0))
sum(as.integer(SpData$Pentameris.airoides != 0))
sum(as.integer(SpData$other != 0))

SpData$other <- SpData$other + SpData$Daucus.glochidiatus 

SpData <- subset(SpData, select = c(total.fecundity, plot, intra, Arctotheca.calendula,
                                    Plantago.debilis, Podolepis.canescens, Trachymene.cyanopetala, Velleia.rosea,
                                    Medicago.minima, Pentameris.airoides, other))

# Create the model matrix
ModelMatrix <- subset(SpData, select = c(intra, Arctotheca.calendula,
                                         Plantago.debilis, Podolepis.canescens, Trachymene.cyanopetala, Velleia.rosea,
                                         Medicago.minima, Pentameris.airoides, other))
ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(SpData, select = plot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}

initials <- list(epsilon=rep(1,P), sigma = .01)
initials1<- list(initials, initials, initials)

PrelimFit <- stan(file = 'BH stan models/HYGL.stan', data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 6000, chains = 3, thin = 2, init = initials1)
PrelimFit
save(PrelimFit, file = "BH stan models/Hygl_posteriors.rdata")

hygl.posts <- PrelimFit


# MEDI ------------------------------------------------
SpData <- read.csv("Data/groups.MEDI.csv")
SpData <- subset(SpData, select = c(total.fecundity, plot, intra, Arctotheca.calendula, Daucus.glochidiatus, Hyalosperma.glutinosum, 
                                    Plantago.debilis, Podolepis.canescens, Trachymene.cyanopetala, Velleia.rosea,
                                    Pentameris.airoides, other))
SpData <- na.omit(SpData)
N <- as.integer(dim(SpData)[1])
P <- as.integer(16)
Fecundity <- as.integer(SpData$total.fecundity)
# check which species don't co-occur above the threshold of 3 obvservations (total, regardless of abundances)
sum(as.integer(SpData$Daucus.glochidiatus != 0))
sum(as.integer(SpData$Hyalosperma.glutinosum != 0))
sum(as.integer(SpData$Plantago.debilis != 0))
sum(as.integer(SpData$Podolepis.canescens != 0))
sum(as.integer(SpData$Trachymene.cyanopetala != 0))
sum(as.integer(SpData$Velleia.rosea != 0))
sum(as.integer(SpData$Arctotheca.calendula != 0))
sum(as.integer(SpData$Pentameris.airoides != 0))
sum(as.integer(SpData$other != 0))

SpData$other <- SpData$other + SpData$Daucus.glochidiatus + SpData$Plantago.debilis+ SpData$Arctotheca.calendula + SpData$Trachymene.cyanopetala


SpData <- subset(SpData, select = c(total.fecundity, plot, intra, Hyalosperma.glutinosum, 
                                    Podolepis.canescens, Velleia.rosea,
                                    Pentameris.airoides, other))
# Create the model matrix
ModelMatrix <- subset(SpData, select = c(intra, Hyalosperma.glutinosum, 
                                         Podolepis.canescens, Velleia.rosea,
                                         Pentameris.airoides, other))
ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(SpData, select = plot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}

initials <- list(epsilon=rep(1,P), sigma = .01)
initials1<- list(initials, initials, initials)

PrelimFit <- stan(file = 'BH stan models/MEDI.stan', data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 6000, chains = 3, thin = 2, init = initials1)
PrelimFit
save(PrelimFit, file = "BH stan models/Medi_posteriors.rdata")

medi.posts <- PrelimFit

plot(medi.posts, pars=c("alpha_intra", "alpha_hygl", "alpha_poca", "alpha_vero", "alpha_peai",
                        "alpha_other"), show_density = TRUE, ci_level = 0.5, outer_level = 0.95, fill_color = "grey") + ggtitle(expression(italic("M. minima")))


# PEAI ------------------------------------------------
SpData <- read.csv("Data/groups.PEAI.csv")
SpData <- subset(SpData, select = c(total.fecundity, plot, intra, Arctotheca.calendula, Daucus.glochidiatus, Hyalosperma.glutinosum, 
                                    Plantago.debilis, Podolepis.canescens, Trachymene.cyanopetala, Velleia.rosea, Medicago.minima,
                                    other))
SpData <- na.omit(SpData)
N <- as.integer(dim(SpData)[1])
P <- as.integer(16)
Fecundity <- as.integer(SpData$total.fecundity)
# check which species don't co-occur above the threshold of 3 obvservations (total, regardless of abundances)
sum(as.integer(SpData$Daucus.glochidiatus != 0))
sum(as.integer(SpData$Hyalosperma.glutinosum != 0))
sum(as.integer(SpData$Plantago.debilis != 0))
sum(as.integer(SpData$Podolepis.canescens != 0))
sum(as.integer(SpData$Trachymene.cyanopetala != 0))
sum(as.integer(SpData$Velleia.rosea != 0))
sum(as.integer(SpData$Arctotheca.calendula != 0))
sum(as.integer(SpData$Medicago.minima != 0))
sum(as.integer(SpData$other != 0))

SpData$other <- SpData$other + SpData$Daucus.glochidiatus 


SpData <- subset(SpData, select = c(total.fecundity, plot, intra, Arctotheca.calendula, Daucus.glochidiatus, Hyalosperma.glutinosum, 
                                    Plantago.debilis, Podolepis.canescens, Trachymene.cyanopetala, Velleia.rosea, Medicago.minima,
                                    other))
SpData$other <- SpData$other + SpData$Daucus.glochidiatus
# subset for model matrix 
ModelMatrix <- subset(SpData, select = c(intra, Arctotheca.calendula, Hyalosperma.glutinosum, 
                                         Plantago.debilis, Podolepis.canescens, Trachymene.cyanopetala, Velleia.rosea, Medicago.minima,
                                         other))
ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(SpData, select = plot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}

initials <- list(epsilon=rep(1,P), sigma = .01)
initials1<- list(initials, initials, initials)

PrelimFit <- stan(file = 'BH stan models/PEAI.stan', data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 6000, chains = 3, thin = 2, init = initials1)
PrelimFit
save(PrelimFit, file = "BH stan models/Peai_posteriors.rdata")

peai.posts <- PrelimFit


# PLDE -----------------------------------------------
SpData <- read.csv("Data/groups.PLDE.csv")
SpData <- subset(SpData, select = c(total.fecundity, plot, intra, Arctotheca.calendula, Daucus.glochidiatus, Hyalosperma.glutinosum, 
                                    Podolepis.canescens, Trachymene.cyanopetala, Velleia.rosea,
                                    Medicago.minima, Pentameris.airoides, other
))
SpData <- na.omit(SpData)
N <- as.integer(dim(SpData)[1])
P <- as.integer(16)
Fecundity <- as.integer(SpData$total.fecundity)
# check which species don't co-occur above the threshold of 3 obvservations (total, regardless of abundances)
sum(as.integer(SpData$Daucus.glochidiatus != 0))
sum(as.integer(SpData$Arctotheca.calendula != 0))
sum(as.integer(SpData$Hyalosperma.glutinosum != 0))
sum(as.integer(SpData$Podolepis.canescens != 0))
sum(as.integer(SpData$Trachymene.cyanopetala != 0))
sum(as.integer(SpData$Velleia.rosea != 0))
sum(as.integer(SpData$Medicago.minima != 0))
sum(as.integer(SpData$Pentameris.airoides != 0))
sum(as.integer(SpData$other != 0))
# and remove them
SpData$other <- SpData$other + SpData$Daucus.glochidiatus + SpData$Arctotheca.calendula + SpData$Medicago.minima + SpData$Trachymene.cyanopetala

SpData <- subset(SpData, select = c(total.fecundity, plot, intra, Hyalosperma.glutinosum, 
                                    Podolepis.canescens, Velleia.rosea,
                                    Pentameris.airoides, other))


# Create the model matrix
ModelMatrix <- subset(SpData, select = c(intra, Hyalosperma.glutinosum, 
                                         Podolepis.canescens, Velleia.rosea,
                                         Pentameris.airoides, other))
ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(SpData, select = plot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}
dim(ModelMatrix)

initials <- list(epsilon=rep(1,P), sigma = .01)
initials1<- list(initials, initials, initials)

PrelimFit <- stan(file = 'BH stan models/PLDE.stan', data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 9000, chains = 3, thin = 3, init = initials1, control = list(adapt_delta=0.9999,max_treedepth=15))


PrelimFit
#traceplot(PrelimFit)
#stan_dens(PrelimFit, pars=c("alpha_intra", "alpha_hygl", "alpha_poca",
#                            "alpha_vero", "alpha_peai"))
save(PrelimFit, file = "BH stan models/Plde_posteriors.rdata")

plde.posts <- PrelimFit

# POCA -----------------------------------------------
SpData <- read.csv("Data/groups.POCA.csv")
SpData <- subset(SpData, select = c(total.fecundity, plot, intra, Arctotheca.calendula, Daucus.glochidiatus, Hyalosperma.glutinosum, 
                                    Plantago.debilis, Trachymene.cyanopetala, Velleia.rosea,
                                    Medicago.minima, Pentameris.airoides, other
))
SpData <- na.omit(SpData)
N <- as.integer(dim(SpData)[1])
P <- as.integer(16)
Fecundity <- as.integer(SpData$total.fecundity)
# check which species don't co-occur above the threshold of 3 obvservations (total, regardless of abundances)
sum(as.integer(SpData$Daucus.glochidiatus != 0))
sum(as.integer(SpData$Arctotheca.calendula != 0))
sum(as.integer(SpData$Hyalosperma.glutinosum != 0))
sum(as.integer(SpData$Plantago.debilis != 0))
#sum(as.integer(SpData$Podolepis.canescens != 0))
sum(as.integer(SpData$Trachymene.cyanopetala != 0))
sum(as.integer(SpData$Velleia.rosea != 0))
sum(as.integer(SpData$Medicago.minima != 0))
sum(as.integer(SpData$Pentameris.airoides != 0))
sum(as.integer(SpData$other != 0))
# and remove them
SpData$other <- SpData$other + SpData$Daucus.glochidiatus + SpData$Hyalosperma.glutinosum + SpData$Medicago.minima 

SpData <- subset(SpData, select = c(total.fecundity, plot, intra, Arctotheca.calendula,  
                                    Plantago.debilis, Trachymene.cyanopetala, Velleia.rosea,
                                    Pentameris.airoides, other
))


# Create the model matrix
ModelMatrix <- subset(SpData, select = c(intra, Arctotheca.calendula,  
                                         Plantago.debilis, Trachymene.cyanopetala, Velleia.rosea,
                                         Pentameris.airoides, other))
ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(SpData, select = plot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}
dim(ModelMatrix)
# Do a preliminary model fit to compile the stan model and check for autocorrelation
#    and convergence
initials <- list(epsilon=rep(1,P), sigma = .01)
initials1<- list(initials, initials, initials)

PrelimFit <- stan(file = 'BH stan models/POCA.stan', data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 6000, chains = 3, thin = 2, init = initials1)

PrelimFit
save(PrelimFit, file = "BH stan models/Poca_posteriors.rdata")

poca.posts <- PrelimFit


# TRCY -----------------------------------------------
SpData <- read.csv("Data/groups.TRCY.csv")
#recode plot column 
SpData <- subset(SpData, select = c(total.fecundity, plot, intra, Arctotheca.calendula, Daucus.glochidiatus, Hyalosperma.glutinosum, 
                                    Plantago.debilis, Podolepis.canescens, Velleia.rosea,
                                    Medicago.minima, Pentameris.airoides, other
))


plotcodes<-data.frame(cbind(sort(unique(SpData$plot)), rep(1:22,1)))
names(plotcodes)[names(plotcodes) == "X1"] <- "plot"
names(plotcodes)[names(plotcodes) == "X2"] <- "newplot"
SpData <- merge(SpData, plotcodes, by = "plot", all.x = TRUE)


SpData <- na.omit(SpData)
N <- as.integer(dim(SpData)[1])
P <- as.integer(22)
Fecundity <- as.integer(SpData$total.fecundity)
# check which species don't co-occur above the threshold of 3 obvservations (total, regardless of abundances)
sum(as.integer(SpData$Arctotheca.calendula != 0))
sum(as.integer(SpData$Daucus.glochidiatus != 0))
sum(as.integer(SpData$Hyalosperma.glutinosum != 0))
sum(as.integer(SpData$Plantago.debilis != 0))
sum(as.integer(SpData$Podolepis.canescens != 0))
#sum(as.integer(SpData$Trachymene.cyanopetala != 0))
sum(as.integer(SpData$Velleia.rosea != 0))
sum(as.integer(SpData$Medicago.minima != 0))
sum(as.integer(SpData$Pentameris.airoides != 0))
sum(as.integer(SpData$other != 0))
# and remove them
SpData$other <- SpData$other + SpData$Daucus.glochidiatus + SpData$Arctotheca.calendula + SpData$Hyalosperma.glutinosum

SpData <- subset(SpData, select = c(total.fecundity, newplot, intra, 
                                    Plantago.debilis, Podolepis.canescens, Velleia.rosea,
                                    Pentameris.airoides, other
))


# Create the model matrix
ModelMatrix <- subset(SpData, select = c(intra, 
                                         Plantago.debilis, Podolepis.canescens, Velleia.rosea,
                                         Pentameris.airoides, other))
ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(SpData, select = newplot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}
dim(ModelMatrix)
# Do a preliminary model fit to compile the stan model and check for autocorrelation
#    and convergence
initials <- list(epsilon=rep(1,P), sigma = .01)
initials1<- list(initials, initials, initials)

PrelimFit <- stan(file = 'BH stan models/TRCY.stan', data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 6000, chains = 3, thin = 2, init = initials1)

#control = list(adapt_delta = 0.9
PrelimFit
save(PrelimFit, file = "BH stan models/Trcy_posteriors.rdata")

trcy.posts <- PrelimFit

# VERO ----------------------------------------------
SpData <- read.csv("Data/groups.VERO.csv")
SpData <- subset(SpData, select = c(total.fecundity, plot, intra, Arctotheca.calendula, Daucus.glochidiatus, Hyalosperma.glutinosum, 
                                    Plantago.debilis, Podolepis.canescens, Trachymene.cyanopetala,
                                    Medicago.minima, Pentameris.airoides, other))
SpData <- na.omit(SpData)
N <- as.integer(dim(SpData)[1])
P <- as.integer(16)
Fecundity <- as.integer(SpData$total.fecundity)
# check which species don't co-occur above the threshold of 3 obvservations (total, regardless of abundances)
sum(as.integer(SpData$Arctotheca.calendula != 0))
sum(as.integer(SpData$Daucus.glochidiatus != 0))
sum(as.integer(SpData$Hyalosperma.glutinosum != 0))
sum(as.integer(SpData$Plantago.debilis != 0))
sum(as.integer(SpData$Podolepis.canescens != 0))
sum(as.integer(SpData$Trachymene.cyanopetala != 0))
sum(as.integer(SpData$Velleia.rosea != 0))
sum(as.integer(SpData$Medicago.minima != 0))
sum(as.integer(SpData$Pentameris.airoides != 0))
sum(as.integer(SpData$other != 0))
# and remove them
SpData$other <- SpData$other + SpData$Arctotheca.calendula + SpData$Daucus.glochidiatus + SpData$Plantago.debilis +SpData$Trachymene.cyanopetala + SpData$Medicago.minima

SpData <- subset(SpData, select = c(total.fecundity, plot, intra, Hyalosperma.glutinosum, 
                                    Podolepis.canescens, Trachymene.cyanopetala,
                                    Pentameris.airoides, other
))


# Create the model matrix
ModelMatrix <- subset(SpData, select = c(intra, Hyalosperma.glutinosum, 
                                         Podolepis.canescens, Trachymene.cyanopetala,
                                         Pentameris.airoides, other))
ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(SpData, select = plot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}
dim(ModelMatrix)

initials <- list(epsilon=rep(1,P), sigma = .01)
initials1<- list(initials, initials, initials)

PrelimFit <- stan(file = 'BH stan models/VERO.stan', data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 6000, chains = 3, thin = 2, init = initials1)
PrelimFit
save(PrelimFit, file = "BH stan models/Vero_posteriors.rdata")

vero.posts <- PrelimFit

##################
## make plot of lambda and alpha posteriors for supplementary information ####
library(ggplot2)
library(cowplot)
bb <- plot(hygl.posts, pars=c("lambda"), show_density = TRUE, ci_level = 0.5, outer_level = 0.95, fill_color = "grey") + ggtitle(expression(italic("H. glutinosum")))+coord_cartesian(xlim = c(0, 360))
cc <- plot(plde.posts, pars=c("lambda"), show_density = TRUE, ci_level = 0.5, outer_level = 0.95, fill_color = "grey") + ggtitle(expression(italic("P. debilis")))+coord_cartesian(xlim = c(0, 360))
dd <- plot(poca.posts, pars=c("lambda"), show_density = TRUE, ci_level = 0.5, outer_level = 0.95, fill_color = "grey") + ggtitle(expression(italic("*P. canescens")))
ee <- plot(trcy.posts, pars=c("lambda"), show_density = TRUE, ci_level = 0.5, outer_level = 0.95, fill_color = "grey") + ggtitle(expression(italic("T. cyanopetala")))+coord_cartesian(xlim = c(0, 360))
ff <- plot(vero.posts, pars=c("lambda"), show_density = TRUE, ci_level = 0.5, outer_level = 0.95, fill_color = "grey") + ggtitle(expression(italic("G. rosea")))+coord_cartesian(xlim = c(0, 360))
gg <- plot(arca.posts, pars=c("lambda"), show_density = TRUE, ci_level = 0.5, outer_level = 0.95, fill_color = "grey") + ggtitle(expression(italic("A. calendula")))+coord_cartesian(xlim = c(0, 360))
hh <- plot(medi.posts, pars=c("lambda"), show_density = TRUE, ci_level = 0.5, outer_level = 0.95, fill_color = "grey") + ggtitle(expression(italic("M. minima")))+coord_cartesian(xlim = c(0, 360))
ii <- plot(peai.posts, pars=c("lambda"), show_density = TRUE, ci_level = 0.5, outer_level = 0.95, fill_color = "grey") + ggtitle(expression(italic("P. airoides")))+coord_cartesian(xlim = c(0, 360))

posts1 <- plot_grid(bb, cc, dd, ee, ff, gg, hh, ii, nrow = 4)
pdf("Figures/lambda_posteriors_supp.pdf", width = 10, height=13)
posts1
dev.off()

# make the x axis the same
b <- plot(hygl.posts, pars=c("alpha_intra", "alpha_arca", "alpha_plde", "alpha_poca", "alpha_trcy", "alpha_vero", "alpha_medi", "alpha_peai",
                        "alpha_other"), show_density = TRUE, ci_level = 0.5, outer_level = 0.95, fill_color = "grey") + ggtitle(expression(italic("H. glutinosum")))+coord_cartesian(xlim = c(0, 1.2))
c <- plot(plde.posts, pars=c("alpha_intra", "alpha_hygl", "alpha_poca",  "alpha_vero", "alpha_peai",
                        "alpha_other"), show_density = TRUE, ci_level = 0.5, outer_level = 0.95, fill_color = "grey") + ggtitle(expression(italic("P. debilis")))+coord_cartesian(xlim = c(0, 1.2))
d <- plot(poca.posts, pars=c("alpha_intra", "alpha_arca", "alpha_plde", "alpha_trcy", "alpha_vero", "alpha_peai",
                        "alpha_other"), show_density = TRUE, ci_level = 0.5, outer_level = 0.95, fill_color = "grey") + ggtitle(expression(italic("P. canescens")))+coord_cartesian(xlim = c(0, 1.2))
e <- plot(trcy.posts, pars=c("alpha_intra", "alpha_plde", "alpha_poca", "alpha_vero", "alpha_peai",
                        "alpha_other"), show_density = TRUE, ci_level = 0.5, outer_level = 0.95, fill_color = "grey") + ggtitle(expression(italic("T. cyanopetala")))+coord_cartesian(xlim = c(0, 1.2))
f <- plot(vero.posts, pars=c("alpha_intra", "alpha_hygl", "alpha_poca", "alpha_trcy", "alpha_peai",
                        "alpha_other"), show_density = TRUE, ci_level = 0.5, outer_level = 0.95, fill_color = "grey") + ggtitle(expression(italic("G. rosea")))+coord_cartesian(xlim = c(0, 1.2))
g <- plot(arca.posts, pars=c("alpha_intra", "alpha_plde", "alpha_poca", "alpha_trcy", "alpha_vero",  "alpha_peai",
                        "alpha_other"), show_density = TRUE, ci_level = 0.5, outer_level = 0.95, fill_color = "grey") + ggtitle(expression(italic("*A. calendula")))+coord_cartesian(xlim = c(0, 7.5))
h <- plot(medi.posts, pars=c("alpha_intra", "alpha_hygl", "alpha_poca", "alpha_vero", "alpha_peai",
                        "alpha_other"), show_density = TRUE, ci_level = 0.5, outer_level = 0.95, fill_color = "grey") + ggtitle(expression(italic("M. minima")))+coord_cartesian(xlim = c(0, 1.2))
i <- plot(peai.posts, pars=c("alpha_intra", "alpha_arca", "alpha_hygl", "alpha_plde", "alpha_poca", "alpha_trcy", "alpha_vero", "alpha_medi",
                        "alpha_other"), show_density = TRUE, ci_level = 0.5, outer_level = 0.95, fill_color = "grey") + ggtitle(expression(italic("P. airoides")))+coord_cartesian(xlim = c(0, 1.2))

posts <- plot_grid(b,c,d,e,f,g,h,i, nrow = 4)
pdf("Figures/alpha_posteriors_supp.pdf", width = 10, height=13)
posts
dev.off()

