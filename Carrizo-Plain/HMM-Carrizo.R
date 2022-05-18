#### Carrizo HMM ####
rm(list=ls())
source("HMMs/HMM_Functions.R")

#### Load Libraries ####
library(ggplot2)
library(tidyverse)

#### Prep Data ####
cp <- read.csv("Data/Long-Term-Datasets/Carrizo-Plain/Carrizo_2007_2021_coverdata_long_format_with_cover_percentage.csv")
rm <- c("BARE", "LITTER","FRESHDIRT","BPHOLE","HOLE","ASTspp","COWPIE","VULPIA","ANT","ALLIUM","BURROW","GOPHER","unknown10", "UNK 10", "UNK 11", "MOSS", "DICCAP")
cp$Precinct <- ifelse(cp$Precinct == "n" | cp$Precinct == "N", "N", "P")
cp <- filter(cp, !(SpeciesCode %in% rm), Precip_Treatment == "NONE", GKR == "GKR", Precinct == "P")
cp <- ddply(cp, .(New.Plot.ID, year, SpeciesCode), summarize, PA = 1)
cp <- pivot_wider(cp, names_from = "year", values_from = "PA")
cp[is.na(cp)] <- 0

#### Create Species List ####

species.list <- list()

for(i in unique(cp$SpeciesCode)) {
  tmp <- filter(cp, SpeciesCode == i)
  tmp <- tmp[,-c(1:2)]
  tmp <- tmp[ , sort(colnames(tmp))]
  if(nrow(tmp) > 8) species.list[[i]] <- tmp 
}

df <- data.frame(names = NA, obs = NA)

for(i in names(species.list)) {
  df <- data.frame(rbind(df, c(names = i, obs = nrow(species.list[[i]]))))
}
df$obs <- as.numeric(as.character(df$obs))

#saveRDS(species.list, "HMMs/Carrizo-Plain/CP_Species-List.RDS")

#### Run HMM ####
# starting values. might not necessarily help much to change these, but we can try
trueParam = numeric(0) #object grouping together the parameters and other quantities which depend on them.
trueParam$c = 0.2      #colonization rate
trueParam$g = 0.5      #germination rate (sigma)
trueParam$r = 1        #reproductive rate (phi) keep this at 1 (if a seed germinates and survives to flower, we assume it produces a seed)
trueParam$s = 0.5      #seed survival
trueParam$p0 = 0.5     #initial state of the seed bank (probability that there were seeds in the soil the year before the first obs of existing flora)

trueParam = makeParametersCalculations(trueParam)

#species.list <- readRDS("CP_Species-List-for-HMM_Site.RDS") # list; each object is a dataframe of PA data for a single species where each row is a site and each column is a year


CP.df <- expand.grid(Species_Name = names(species.list), p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA) # empty df to fill with rates

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
  CP.df[CP.df$Species_Name == j,]$p0 <- EMresult$param$p0
  CP.df[CP.df$Species_Name == j,]$g <- EMresult$param$g
  CP.df[CP.df$Species_Name == j,]$c <- EMresult$param$c
  CP.df[CP.df$Species_Name == j,]$s <- EMresult$param$s
  CP.df[CP.df$Species_Name == j,]$r <- EMresult$param$r
  CP.df[CP.df$Species_Name == j,]$iter <- which.max(EMresult$lllist)
  
}

# remove species that didn't converge in 150 steps
CP.df <- filter(CP.df, iter < 100)

#### Explore Output ####
trait <- read.csv("Data/Long-Term-Datasets/Carrizo-Plain/Carrizo_Traits.csv")
trait$FunGroup <- paste(trait$Status, trait$Habit)
colnames(CP.df)[1] <- "Code"
CP.df <- merge(CP.df, unique(trait[,c(1,2,10,12)]), by.x = "Code", by.y = "SpeciesCode", all.x = T)
CP.df <- filter(CP.df, Annual.Perennial == "Annual")
CP.df <- CP.df[,c(8,1,10,2:7)]

# CHANGE DEPENDING ON WHAT WAS RUN
CP.on.df <- CP.df

# SAVE
saveRDS(CP.all.df, "HMMs/Carrizo-Plain/CP-GKR_HMM.RDS")
saveRDS(CP.on.df, "HMMs/Carrizo-Plain/CP-GKR-Precinct_HMM.RDS")
saveRDS(CP.off.df, "HMMs/Carrizo-Plain/CP-GKR-off-Precinct_HMM.RDS")


CP.off.df <- readRDS("HMMs/Carrizo-Plain/CP-GKR-off-Precinct_HMM.RDS")
CP.on.df <- readRDS("HMMs/Carrizo-Plain/CP-GKR-Precinct_HMM.RDS")
CP.all.df <- readRDS("HMMs/Carrizo-Plain/CP-GKR_HMM.RDS")

CP.off.df$Precinct <- "N"
CP.on.df$Precinct <- "P"
CP.all.df$Precinct <- "both"

CP.df <- rbind(CP.off.df, CP.on.df, CP.all.df)

ggplot(CP.df[CP.df$Precinct != "P",], aes(x = s, y = c, col = Precinct, group = Precinct)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_classic()

ggplot(CP.df, aes(x = s, y = c)) +
  #geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  ylim(0,1) +
  geom_text(aes(label = Code)) +
  facet_wrap(~Precinct) +
  theme_classic()

ggplot(CP.all.df, aes(x = s, y = c)) +
  #geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  geom_text(aes(label = Code)) +
  theme_classic()
