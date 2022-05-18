#### McLaughlin HMM ####
rm(list = ls())
source("HMMs/HMM_Functions.R")

#### Load Libraries ####
library(ggplot2)
library(plyr)
library(tidyverse)

#### Prep Data ####
trait <- read.csv("/Users/Marina/Google Drive/02_McLaughlin_80Sites_Organized/Modified_CombinedFiles/McL_80SitesSpeciesTraits_012615.csv")
mcl <- read.csv("/Users/Marina/Google Drive/02_McLaughlin_80Sites_Organized/Data-Storage-Files/Core_Community_Data2019.csv")
serp <- read.csv("/Users/Marina/Google Drive/02_McLaughlin_80Sites_Organized/Data-Storage-Files/McL_Abiotic-Data.csv")

rm <- c("Logfia spp./Micropus californicus", "Aegilops triuncialis") # rm commonly mis-ID'd species and goatgrass which is treated at the reserve

mcl <- filter(mcl, Species_Name %in% trait[trait$Annual.Perennial == "Annual",]$Species_Name, !(Species_Name %in% rm)) 

mcl <- ddply(mcl, .(Site, Year, Species_Name), summarize, PA = 1)
mcl <- pivot_wider(mcl, names_from = "Year", values_from = "PA")
mcl <- as.data.frame(mcl)
mcl[is.na(mcl)] <- 0
mcl <- merge(serp[,1:2], mcl, by = "Site")

mcl.S <- filter(mcl, Serpentine == "S")
mcl.N <- filter(mcl, Serpentine == "N")

#### Create Species List ####

#Full
species.list <- list()

for(i in unique(mcl$Species_Name)) {
  tmp <- filter(mcl, Species_Name == i)
  tmp <- tmp[,-c(1:3)]
  tmp <- tmp[ , sort(colnames(tmp))]
  if(nrow(tmp) > 8) species.list[[i]] <- tmp # for now I got rid of species that were present in less than 8 sites throughout the entire time series
}

# Serp only
species.list <- list()

for(i in unique(mcl.S$Species_Name)) {
  tmp <- filter(mcl.S, Species_Name == i)
  tmp <- tmp[,-c(1:3)]
  tmp <- tmp[ , sort(colnames(tmp))]
  if(nrow(tmp) > 8) species.list[[i]] <- tmp # for now I got rid of species that were present in less than 8 sites throughout the entire time series
}

# Nonserp only
species.list <- list()
for(i in unique(mcl.N$Species_Name)) {
  tmp <- filter(mcl.N, Species_Name == i)
  tmp <- tmp[,-c(1:3)]
  tmp <- tmp[ , sort(colnames(tmp))]
  if(nrow(tmp) > 8) species.list[[i]] <- tmp # for now I got rid of species that were present in less than 8 sites throughout the entire time series
}

df <- data.frame(names = NA, obs = NA)

for(i in names(species.list)) {
  df <- data.frame(rbind(df, c(names = i, obs = nrow(species.list[[i]]))))
}

df$obs <- as.numeric(as.character(df$obs))

#saveRDS(species.list, "HMMs/McLaughlin/McL_Species-List.RDS")
#saveRDS(species.list, "HMMs/McLaughlin/McLS_Species-List.RDS")
#saveRDS(species.list, "HMMs/McLaughlin/McLN_Species-List.RDS")

#### Run HMM ####
# starting values. might not necessarily help much to change these, but we can try
trueParam = numeric(0) #object grouping together the parameters and other quantities which depend on them.
trueParam$c = 0.2      #colonization rate
trueParam$g = 0.5      #germination rate (sigma)
trueParam$r = 1        #reproductive rate (phi) keep this at 1 (if a seed germinates and survives to flower, we assume it produces a seed)
trueParam$s = 0.5      #seed survival
trueParam$p0 = 0.5     #initial state of the seed bank (probability that there were seeds in the soil the year before the first obs of existing flora)

trueParam = makeParametersCalculations(trueParam)

#species.list <- readRDS("McL_Species-List.RDS") # change depending on which one is being run; each object is a dataframe of PA data for a single species where each row is a site and each column is a year


mcl.df <- expand.grid(Species_Name = names(species.list), p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA) # empty df to fill with rates


# Took about 20 minutes or so
for(j in names(species.list)){ # for each species
  X = as.matrix(species.list[[j]]) # use their time series PA data
  n = 5 # why n = 5?
  p0Results = rep(0,n)
  gResults = rep(0,n)
  cResults = rep(0,n)
  sResults = rep(0,n)
  rResults = rep(0,n)
  llResults = rep(0,n)
  lltrueParam = rep(0,n)
  
  for (i in 1:n) {
    for (k in 1:5) print(i) # I don't understand what's going on here, why print it 5 times?
    print(logLikelihood(X, trueParam)) # print the log-likelihood given the "true params"
    EMresult = EMestimation(X, r = 1) # update the model params, stop when new log likelihood is lower than the old log likelihood plus some precision 
    print(logLikelihood(X, trueParam)) # print the log-likelihood given the "true params"; why is this on here twice?
    # fill in results with new params from the best log likelihood model
    p0Results[i] = EMresult$param$p0 # initial seed bank prob
    gResults[i] = EMresult$param$g # germ
    cResults[i] = EMresult$param$c #col
    sResults[i] = EMresult$param$s #surv
    rResults[i] = EMresult$param$r #prod
    llResults[i] = EMresult$ll # log-likelihood
    lltrueParam[i] = logLikelihood(X,trueParam) # log likelihood of x given "true" parameters
  }

  # here I assumed that the last estimated parameters, were what the model converged on but I'm not sure that's correct?
  mcl.df[mcl.df$Species_Name == j,]$p0 <- EMresult$param$p0
  mcl.df[mcl.df$Species_Name == j,]$g <- EMresult$param$g
  mcl.df[mcl.df$Species_Name == j,]$c <- EMresult$param$c
  mcl.df[mcl.df$Species_Name == j,]$s <- EMresult$param$s
  mcl.df[mcl.df$Species_Name == j,]$r <- EMresult$param$r
  mcl.df[mcl.df$Species_Name == j,]$iter <- which.max(EMresult$lllist)
  
}

# remove species that didn't converge in 100 steps
mcl.df <- filter(mcl.df, iter < 100)

#### Explore Output ####
trait <- read.csv("Data/Long-Term-Datasets/McLaughlin/McL_80SitesSpeciesTraits_012615.csv")
trait$FunGroup <- paste(trait$Native.Exotic, trait$Grass.Forb.Shrub)
mcl.df <- merge(mcl.df, trait[,c(3,24)], by = "Species_Name", all.y = F)

# mcl.df <- mcl.df[,c(1,8,9,2:7)]
mclN.df <- mcl.df[,c(1,8,7,2:6)]
#mclS.df <- mcl.df[,c(1,8,7,2:6)]

saveRDS(mclS.df, "Hmms/McLaughlin/McL-S_HMM.RDS")
saveRDS(mcl.df, "Hmms/McLaughlin/McL-N_HMM.RDS")
saveRDS(mcl.df, "Hmms/McLaughlin/McL_HMM.RDS")

mclS.df <- readRDS("Hmms/McLaughlin/McL-S_HMM.RDS")
mclN.df <- readRDS("Hmms/McLaughlin/McL-N_HMM.RDS")
mcl.df <- readRDS("Hmms/McLaughlin/McL_HMM.RDS")

# No soil types #
ggplot(mcl.df[mcl.df$FunGroup != "Native Grass",], aes(x = s, y = c, col = FunGroup, group = FunGroup)) +
  #geom_point() +
  geom_text(aes(label = Code)) +
  geom_smooth(method = "lm", formula = y ~ x, se = F)

# Soil types #
ggplot(mclN.df, aes(x = s, y = c, col = FunGroup, group = FunGroup)) +
  geom_point() +
  geom_text(aes(label = Code)) +
  geom_smooth(method = "lm", formula = y ~ x, se = F)

ggplot(mclS.df, aes(x = s, y = c, col = FunGroup, group = FunGroup)) +
  #geom_point() +
  geom_text(aes(label = Code)) +
  geom_smooth(method = "lm", formula = y ~ x, se = F)

mclS.df$Serp <- "S"
mclN.df$Serp <- "N"
mcl.df.Soil <- rbind(mclS.df, mclN.df)

ggplot(mcl.df.Soil, aes(x = Serp, y = s)) +
  geom_boxplot(aes(col = FunGroup))

ggplot(mcl.df.Soil, aes(x = Serp, y = c)) +
  geom_boxplot(aes(col = FunGroup))
