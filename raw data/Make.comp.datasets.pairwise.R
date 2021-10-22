# Code to convert the community composition data from long form to wide form 
#        and then join it with the fecundity data for each species 
#        and then make species specific datasets 
library(tidyverse)

rm(list=ls())
# the latest data on RDM
Fecundity <- read.csv("raw data/Fecundity_2018.csv") 
Community <- read.csv("raw data/Community_composition.csv") 

f <- subset(Fecundity, select = c(focal.species, plot, treatment, subplot, mature, immature))
comm <- subset(Community, select = c(focal.species, plot, treatment, subplot, neighbour.species,
                                     neighbour.abundance))
# summing the neighbor abundances
a <- subset(comm, comm$neighbour.species=="Arctotheca calendula")
a <- sum(a$neighbour.abundance)
h <- subset(comm, comm$neighbour.species=="Hyalosperma glutinosum")
h <- sum(h$neighbour.abundance)
p <- subset(comm, comm$neighbour.species=="Pentameris airoides")
p <- sum(p$neighbour.abundance, na.rm = T)
pl <- subset(comm, comm$neighbour.species=="Plantago debilis")
pl <- sum(pl$neighbour.abundance)
t <- subset(comm, comm$neighbour.species=="Trachymene cyanopetala")
t <- sum(t$neighbour.abundance)
m <- subset(comm, comm$neighbour.species=="Medicago minima")
m <- sum(m$neighbour.abundance)
v <- subset(comm, comm$neighbour.species=="Velleia rosea")
v <- sum(v$neighbour.abundance)


# remove any values of NA in mature and immature columns 
f[f=='NA'] <- NA # check there are real NA values for both datasets
comm[comm=="NA"] <- NA

# need to to remove rows where mature and immature both have NA's
library(data.table)
f=data.table(f)
setkeyv(f, c("mature", "immature"))
dontwant <- f[J(NA, NA)] 
f <- f[!J(NA, NA)]  # Pull everything mature=na and immature=na
f=as.data.frame(f)

# then change the NA when its only in one of these to a zero 
f$mature[is.na(f$mature)] <- 0
f$immature[is.na(f$immature)] <- 0

f <- tbl_df(f)
comm <- tbl_df(comm)

# change community comp data to wide format 
spread.comm <- comm
spread.comm <- spread.comm %>% 
  group_by_at(vars(-neighbour.abundance)) %>% # group by everything other than the number of neighbours
  mutate(row_id=1:n()) %>% ungroup() %>%  # build group index column
  spread(neighbour.species, neighbour.abundance, fill = 0) %>% 
  select(-row_id)  # drop the index

# combine with fecundity data 
#not.in.spread.data <- anti_join(f, spread.comm) # these should be the E plots, we want them read in with neighbour abundance NA, then we change it to 0 
spread.data<- right_join(spread.comm, f) 

# udentift the mismatch in row number between spread.data and f (they're duplicates so need to remove them)
spread.data <- spread.data %>% distinct(focal.species, plot, treatment, subplot, mature, immature, .keep_all = TRUE)

# check you have removed the duplicates (there shouldn't be any here now)
# f %>% filter(focal.species=="Gilberta tenuifolia" & treatment=="C" & plot==7) %>%
#   subset(select=c("focal.species", "plot", "treatment", "subplot", "mature", "immature"))
# spread.data %>% filter(focal.species=="Gilberta tenuifolia" & treatment=="C" & plot==7) %>%
#   subset(select=c("focal.species", "plot", "treatment", "subplot", "mature", "immature"))
# 
# f %>% filter(focal.species=="Pentameris airoides" & treatment=="C" & plot==2) %>%
#   subset(select=c("focal.species", "plot", "treatment", "subplot", "mature", "immature"))
# spread.data %>% filter(focal.species=="Pentameris airoides" & treatment=="C" & plot==2) %>%
#   subset(select=c("focal.species", "plot", "treatment", "subplot", "mature", "immature"))
# 
# f %>% filter(focal.species=="Daucus glochidiatus" & treatment=="C" & plot==2) %>%
#   subset(select=c("focal.species", "plot", "treatment", "subplot", "mature", "immature"))
# spread.data %>% filter(focal.species=="Daucus glochidiatus" & treatment=="C" & plot==2) %>%
#   subset(select=c("focal.species", "plot", "treatment", "subplot", "mature", "immature"))


# make a total fecundity column because we are treating immature same as mature - presumably would also be mature if plant harvested later
spread.data$total.fecundity <- spread.data$mature + spread.data$immature

# reorder the neighbour columns so you can add things more easily for grouping neighbours     
spread.data.ordered <- spread.data[,c(1:4, 89:91, 18, 29, 40, 41, 51, 86, 5, 33, 39, 6:17, 19:28, 30:32, 34:38, 42:50, 52:85, 87)]

# change the NA's in competitor columns to 0's 
spread.data.ordered[rowSums(is.na(spread.data.ordered)) > 0,]
spread.data.ordered[ , 8:90][is.na(spread.data.ordered[ , 8:90])] = 0 

#sort(colSums(spread.data.ordered[8:90]))

#################################################################################
######## NEED TO SUBSET BY FOCAL SPECIES SO INTRA COLUMN STAYS SEPARATE #########
arca <- subset(spread.data.ordered, spread.data.ordered$focal.species=="Arctotheca calendula")
# grouping neighbours 
grouped <- arca
grouped$intra <- grouped$`Arctotheca calendula`
grouped <- select(grouped, -contains("Arctotheca calendula"))
grouped <- grouped %>% mutate(other = rowSums(.[16:89]))
grouped <- grouped[,-c(16:89)]
write.csv(grouped, "Data/groups.ARCA.csv")
remove(grouped)

#################
medi <- subset(spread.data.ordered, spread.data.ordered$focal.species=="Medicago minima")
# grouping neighbours 
grouped <- medi
grouped$intra <- grouped$`Medicago minima`
grouped <- select(grouped, -contains("Medicago minima")) 
grouped <- grouped %>% mutate(other = rowSums(.[16:89]))
grouped <- grouped[,-c(16:89)]
write.csv(grouped, "Data/groups.MEDI.csv")
remove(grouped)


#################
peai <- subset(spread.data.ordered, spread.data.ordered$focal.species=="Pentameris airoides")
# grouping neighbours 
grouped <- peai
grouped$intra <- grouped$`Pentameris airoides`
grouped <- grouped %>% mutate(other = rowSums(.[16:89]))
grouped <- grouped[,-c(16:89)]
write.csv(grouped, "Data/groups.PEAI.csv")
remove(grouped)

#################
dagl <- subset(spread.data.ordered, spread.data.ordered$focal.species=="Daucus glochidiatus")
# grouping neighbours 
grouped <- dagl
grouped$intra <- grouped$`Daucus glochidiatus`
grouped <- select(grouped, -contains("Daucus glochidiatus")) 
grouped <- grouped %>% mutate(other = rowSums(.[16:89]))
grouped <- grouped[,-c(16:89)]
write.csv(grouped, "Data/groups.DAGL.csv")
remove(grouped)

#################
hygl <- subset(spread.data.ordered, spread.data.ordered$focal.species=="Hyalosperma glutinosum")
# grouping neighbours 
grouped <- hygl
grouped$intra <- grouped$`Hyalosperma glutinosum`
grouped <- select(grouped, -contains("Hyalosperma glutinosum")) 
grouped <- grouped %>% mutate(other = rowSums(.[16:89]))
grouped <- grouped[,-c(16:89)]
write.csv(grouped, "Data/groups.HYGL.csv")
remove(grouped)


#################
plde <- subset(spread.data.ordered, spread.data.ordered$focal.species=="Plantago debilis")
# grouping neighbours 
grouped <- plde
grouped$intra <- grouped$`Plantago debilis`
grouped <- select(grouped, -contains("Plantago debilis")) 
grouped <- grouped %>% mutate(other = rowSums(.[16:89]))
grouped <- grouped[,-c(16:89)]
write.csv(grouped, "Data/groups.PLDE.csv")
remove(grouped)

#################
poca <- subset(spread.data.ordered, spread.data.ordered$focal.species=="Podolepis canescens")
# grouping neighbours 
grouped <- poca
grouped$intra <- grouped$`Podolepis canescens`
grouped <- select(grouped, -contains("Podolepis canescens")) 
grouped <- grouped %>% mutate(other = rowSums(.[16:89]))
grouped <- grouped[,-c(16:89)]
write.csv(grouped, "Data/groups.POCA.csv")
remove(grouped)

##################
trcy <- subset(spread.data.ordered, spread.data.ordered$focal.species=="Trachymene cyanopetala")
# grouping neighbours 
grouped <- trcy
grouped$intra <- grouped$`Trachymene cyanopetala`
grouped <- select(grouped, -contains("Trachymene cyanopetala")) 
grouped <- grouped %>% mutate(other = rowSums(.[16:89]))
grouped <- grouped[,-c(16:89)]
write.csv(grouped, "Data/groups.TRCY.csv")
remove(grouped)


##################
vero <- subset(spread.data.ordered, spread.data.ordered$focal.species=="Velleia rosea")
# grouping neighbours 
grouped <- vero
grouped$intra <- grouped$`Velleia rosea`
grouped <- select(grouped, -contains("Velleia rosea")) 
grouped <- grouped %>% mutate(other = rowSums(.[16:89]))
grouped <- grouped[,-c(16:89)]
write.csv(grouped, "Data/groups.VERO.csv")
remove(grouped)


