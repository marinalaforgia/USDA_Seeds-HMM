#### Jasper Ridge HMM ####
rm(list=ls())
source("HMMs/HMM_Functions.R")

#### Load Libraries ####
library(ggplot2)
library(tidyverse)

#### Prep Data ####
jr <- read.csv("Data/Long-Term-Datasets/Jasper-Ridge/JR_Cover_1983-2019.csv")
unique(jr$SpeciesCode)
rm <- c("ROCK", "TRSP", "BR__", "ARDO", "BARE", "POSC", "CHPO","STPU","MECA", "BRTE", "BRTR", "LOUT", "SIJU", "MISP", "THSP", "LASP", "ESCA", "LIAN") # remove perennials and species that don't show up

jr[jr$SpeciesCode == "BRMO",]$SpeciesCode <- "BRHO"
jr[jr$SpeciesCode == "EVSP",]$SpeciesCode <- "HESP"
jr[jr$SpeciesCode == "ORDE",]$SpeciesCode <- "CADE"
jr[jr$SpeciesCode == "EPPA",]$SpeciesCode <- "EPBR"
jr[jr$SpeciesCode == "HELU",]$SpeciesCode <- "HECO"
jr[jr$SpeciesCode == "LIAN",]$SpeciesCode <- "LIPA"
jr[jr$SpeciesCode == "LOSU",]$SpeciesCode <- "LOWR"
jr[jr$SpeciesCode == "TIER",]$SpeciesCode <- "CRCO"
jr[jr$SpeciesCode == "TRTR",]$SpeciesCode <- "TRWI"

jr <-  filter(jr, !(SpeciesCode %in% rm))
jr$PA <- ifelse(jr$PercentCover > 0, 1, 0)

jr.sum <- pivot_wider(jr[,c(1:3,5)], names_from = "SamplingYear", values_from = "PA")
jr.sum[is.na(jr.sum)] <- 0
jr.sum$sumPA <- rowSums(jr.sum[,3:39])
jr.sum <- filter(jr.sum, sumPA > 0) # there are a lot of species that are listed as 0

# TESTing pulling out "independent" quadrats
control.plots <- c("C101", "C103", "C105", "C108", "C110", "C112", "C113", "C115", "C117", "C120", "C122", "C124", "C201", "C203", "C205", "C208", "C210", "C212", "C213", "C215", "C217", "C220", "C222", "C224", "C301", "C303", "C305", "C308", "C310", "C312", "C313", "C315", "C317", "C320", "C322", "C324")

# gopher.plots <- c("G101", "G103", "G105", "G108", "G110", "G112", "G113", "G115", "G117", "G120", "G122", "G124", "G201", "G203", "G205", "G208", "G210", "G212", "G213", "G215", "G217", "G220", "G222", "G224", "G301", "G303", "G305", "G308", "G310", "G312", "G313", "G315", "G317", "G320", "G322", "G324")
# 
# rabbit.plots <- c("R101", "R103", "R105", "R108", "R110", "R112", "R113", "R115", "R117", "R120", "R122", "R124", "R201", "R203", "R205", "R208", "R210", "R212", "R213", "R215", "R217", "R220", "R222", "R224", "R301", "R303", "R305", "R308", "R310", "R312", "R313", "R315", "R317", "R320", "R322", "R324")

jr.sum <- filter(jr.sum, QuadratCode %in% control.plots)
jr.sum <- jr.sum[,-40]

#jr.sum <- filter(jr.sum, QuadratCode %in% gopher.plots)

#### Create Species List ####
species.list <- list()

for(i in unique(jr.sum$SpeciesCode)) {
  tmp <- filter(jr.sum, SpeciesCode == i)
  tmp <- tmp[,-c(1:2)]
  tmp <- tmp[ , sort(colnames(tmp))]
  if(nrow(tmp) > 8) species.list[[i]] <- tmp 
}

#saveRDS(species.list, "HMMs/Jasper-Ridge/JR_Species.RDS")

#### Run HMM ####
# starting values. might not necessarily help much to change these, but we can try
trueParam = numeric(0) #object grouping together the parameters and other quantities which depend on them.
trueParam$c = 0.2      #colonization rate
trueParam$g = 0.5      #germination rate (sigma)
trueParam$r = 1        #reproductive rate (phi) keep this at 1 (if a seed germinates and survives to flower, we assume it produces a seed)
trueParam$s = 0.5      #seed survival
trueParam$p0 = 0.5     #initial state of the seed bank (probability that there were seeds in the soil the year before the first obs of existing flora)

trueParam = makeParametersCalculations(trueParam)

# list; each object is a dataframe of PA data for a single species where each row is a site and each column is a year
jr.df <- expand.grid(Species_Name = names(species.list), p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA) # empty df to fill with rates

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
  
  print(j)
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
  jr.df[jr.df$Species_Name == j,]$p0 <- EMresult$param$p0
  jr.df[jr.df$Species_Name == j,]$g <- EMresult$param$g
  jr.df[jr.df$Species_Name == j,]$c <- EMresult$param$c
  jr.df[jr.df$Species_Name == j,]$s <- EMresult$param$s
  jr.df[jr.df$Species_Name == j,]$r <- EMresult$param$r
  jr.df[jr.df$Species_Name == j,]$iter <- which.max(EMresult$lllist)
}

# remove species that didn't converge in 150 steps
jr.df <- filter(jr.df, iter < 100)

#### Explore Output ####

trait <- read.csv("Data/Long-Term-Datasets/Jasper-Ridge/JR_Traits.csv")
trait$Code.4 <- toupper(trait$Code.4)
trait$FunGroup <- paste(trait$Status, trait$Habit)
jr.df <- merge(jr.df, unique(trait[,c(1,2,3,22)]), by.x = "Species_Name", by.y = "Code.4", all.x = T, all.y = F)
jr.df <- jr.df[,c(8,9,10,2:7)]
colnames(jr.df)[1] <- "Species_Name"
saveRDS(jr.df, "HMMs/Jasper-Ridge/jr_HMM.RDS")

ggplot(jr.df, aes(x = s, y = c, col = FunGroup, group = FunGroup)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  ylim(0,1) +
  theme_classic()
