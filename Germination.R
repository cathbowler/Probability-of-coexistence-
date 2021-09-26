# this code gets posteriors for germination rates to be used in coexistence calculations 

setwd("/Users/catherinebowler/Dropbox/Bayesian data - Cath Collaboration")
library(here)
library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Load in the data
GermData <- read.csv(here("SurvivalAndGermination","AllSpecies.csv"))
GermData <- subset(GermData, GermData$treatment=="dry")
SpeciesToUse <- c("DAGL", "HYGL", "PLDE", 
                  "POCA", "TRCY", "VERO", "ARCA", "MEDI", "PEAI")
GermData$SpeciesCodeNew <- ifelse(GermData$species.code == "TROR", "TRCY", as.character(GermData$species.code))
GermData <- subset(GermData, SpeciesCodeNew %in% SpeciesToUse) 
GermData$Spent.total<-GermData$spent.1+ GermData$spent.2+ GermData$spent.3
GermData <- subset(GermData, select=c("SpeciesCodeNew", "present", "Spent.total"))
GermData<-na.omit(GermData)

N <- as.integer(dim(GermData)[1])
S <- as.integer(length(SpeciesToUse))
present <- as.integer(GermData$present)
spent <- as.integer(GermData$Spent.total)
species <- rep(NA, N)
for(i in 1:N){
     species[i] <- which(SpeciesToUse == GermData$SpeciesCodeNew[i])
}

# Do a preliminary model fit to compile the stan model and check for autocorrelation
#    and convergence
PrelimFit <- stan(file = here('SurvivalAndGermination','Germination.stan'), data = c("N", "S", "present", "spent", "species"),
                  iter = 6000, chains = 3, thin = 2)
PrelimFit

# Take a look at the different convergence of the different chains
traceplot(PrelimFit)

# Check out the posteriors for the different parameters
stan_dens(PrelimFit, pars = "p")

# Plot them together with credible intervals
plot(PrelimFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, fill_color = "salmon")

# Finally check the autocorrelation of the MCMC samples
quartz()
par(mfcol = c(3,3))
# Now extract the posteriors
posteriors <- extract(PrelimFit)
# Now plot each posterior
for(i in 1:9){
     acf(posteriors$p[,i])
}

save(PrelimFit, file = "SurvivalAndGermination/Germination.rdata")


