# this code gets posteriors for survival rates to be used in coexistence calculations 

library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Load in the data
SurvData <- read.csv("SurvivalAndGermination/AllSpecies.csv")
SurvData <- subset(SurvData, SurvData$treatment=="dry")

SpeciesToUse <- c("DAGL", "HYGL", "PLDE", 
                   "POCA", "TRCY", "VERO", "ARCA", "MEDI", "PEAI")
SurvData$SpeciesCodeNew <- ifelse(SurvData$species.code == "TROR", "TRCY", as.character(SurvData$species.code))
SurvData <- subset(SurvData, SpeciesCodeNew %in% SpeciesToUse)
SurvData$Viable.total<-SurvData$germ.1+ SurvData$germ.2+ SurvData$viable 
subset(SurvData, select=c("SpeciesCodeNew", "present", "Viable.total"))
SurvData<-na.omit(SurvData)

N <- as.integer(dim(SurvData)[1])
S <- as.integer(length(SpeciesToUse))
present <- as.integer(SurvData$present)
surv <- as.integer(SurvData$Viable.total)
species <- rep(NA, N)
for(i in 1:N){
  species[i] <- which(SpeciesToUse == SurvData$SpeciesCodeNew[i])
}

# Do a preliminary model fit to compile the stan model and check for autocorrelation
#    and convergence
PrelimFit <- stan(file = 'SurvivalAndGermination/Survival.stan', data = c("N", "S", "present", "surv", "species"),
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
extract(PrelimFit)
# Now plot each posterior
for(i in 1:9){
  acf(posteriors[[i]])
}

save(PrelimFit, file = "SurvivalAndGermination/Survival.rdata")
fit <- extract(PrelimFit)

