#### Sonoran Desert HMM ####
rm(list=ls())
source("./HMM_Functions.R") #making path relative

#### Load Libraries ####
library(ggplot2)
library(tidyverse)
library(plyr)
#### Prep Data ####
sd <- read.csv("Data/Long-Term-Datasets/Sonoran-Desert/CensusData.csv")
trait <- read.csv("Data/Long-Term-Datasets/Sonoran-Desert/Species_list_2016_2.csv")
trait$Species_Name <- paste(trait$Genus, trait$Species)
trait$FunGroup <- paste(trait$Native.Invasive, trait$Forb.Grass)

remove <- c("CRsp", # removing unidentified species and known perennials
            "DA-SP",
            "eriog",
            "Eriog",
            "ERIOG",
            "LO-AS",
            "LOsp",
            "PHsp",
            "PLsp",
            "plsp",
            "ersp",
            "crsp",
            "losp"
            )

sd <-  filter(sd, !species %in% remove)
unique(sd$species)
sd <- ddply(sd, .(Year, plot.habitat.replicate, species), summarize, PA = 1)
sd <- pivot_wider(sd, names_from = "Year", values_from = "PA")
sd[is.na(sd)] <- 0


#### Create Species List ####
species.list <- list()

for(i in unique(sd$species)) {
  tmp <- filter(sd, species == i)
  tmp <- tmp[,-c(1:2)]
  tmp <- tmp[ , sort(colnames(tmp))]
  if(nrow(tmp) > 8) species.list[[i]] <- tmp # gets rid of species that don't appear in more than 8 plots, can try messing with this
}

df <- data.frame(names = NA, obs = NA)

for(i in names(species.list)) {
  df <- data.frame(rbind(df, c(names = i, obs = nrow(species.list[[i]]))))
}

df$obs <- as.numeric(as.character(df$obs))

#saveRDS(species.list, "SD_Species-List.RDS")

#### Run HMM ####
trueParam = numeric(0) #object grouping together the parameters and other quantities which depend on them.
trueParam$c = 0.2      #colonization rate
trueParam$g = 0.5      #germination rate (sigma)
trueParam$r = 1        #reproductive rate (phi) keep this at 1 (if a seed germinates and survives to flower, we assume it produces a seed)
trueParam$s = 0.5      #seed survival
trueParam$p0 = 0.5     #initial state of the seed bank (probability that there were seeds in the soil the year before the first obs of existing flora)

trueParam = makeParametersCalculations(trueParam)


sd.df <- expand.grid(Species_Name = names(species.list), p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA) # empty df to fill with rates (iter = how many iterations it took the model to converge)

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
    for (k in 1:5) print(i) 
    print(logLikelihood(X, trueParam)) # print the log-likelihood given the "true params"
    EMresult = EMestimation(X, r = 1) # update the model params, stop when new log likelihood is lower than the old log likelihood plus some precision 
    print(logLikelihood(X, trueParam)) 
    # fill in results with new params from the best log likelihood model
    p0Results[i] = EMresult$param$p0 # initial seed bank prob
    gResults[i] = EMresult$param$g # germ
    cResults[i] = EMresult$param$c #col
    sResults[i] = EMresult$param$s #surv
    rResults[i] = EMresult$param$r #prod
    llResults[i] = EMresult$ll # log-likelihood
    lltrueParam[i] = logLikelihood(X,trueParam) # log likelihood of x given "true" parameters
  }

  sd.df[sd.df$Species_Name == j,]$p0 <- EMresult$param$p0
  sd.df[sd.df$Species_Name == j,]$g <- EMresult$param$g
  sd.df[sd.df$Species_Name == j,]$c <- EMresult$param$c
  sd.df[sd.df$Species_Name == j,]$s <- EMresult$param$s
  sd.df[sd.df$Species_Name == j,]$r <- EMresult$param$r
  sd.df[sd.df$Species_Name == j,]$iter <- which.max(EMresult$lllist)
}

# remove species that didn't converge in 100 steps
sd.df2 <- filter(sd.df, iter < 100) # can also change this, look at what species it is throwing away

#### Explore Output ####
colnames(sd.df2)[1] <- "code.old"
sd.df2 <- merge(sd.df2, trait[,c(3,4,13,14)], by.x = "code.old", by.y = "Species.code")

sd.df2 <- sd.df2[,c(9,8,10,2:7)]

#saveRDS(sd.df2, "Sonoran-Desert/SD_HMM.RDS") #tried not to save over Marina's files.  I think the results change a bit each time I run them?  But maybe I changed something...
#saveRDS(sd.df, "Sonoran-Desert/SD_HMM_all.RDS")

#sd.df = readRDS("Sonoran-Desert/SD_HMM.RDS")
ggplot(sd.df2, aes(x = s, y = c)) +
  #geom_point() +
  geom_text(aes(label = Code)) +
  geom_smooth(method = "lm", formula = y ~ x, se = F)

