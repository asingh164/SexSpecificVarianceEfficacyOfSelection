###########################################################################################################################################################
###########################################################################################################################################################
###   Script Details:                                                                                                                                   ###
###   Selection Analysis for Singh, A. A.F., Agrawal. 2020. Sex-Specific Variance in Fitness and the Efficacy of Selection                              ###
###   Script written by Amardeep Singh (amardeep.singh[at]utoronto.ca)                                                                                  ###
###   This script will estimate selection differentials acting on condition in our labratory Drosophila experiment and generate corresponding plots     ###
###########################################################################################################################################################
###########################################################################################################################################################

# This script will generate a plot corresponding to Supplemental Figure S3 (both panels separetly)
# Figure will not look exactly as those presented in the paper. Although the data is identical, I may have switched around the order differently than
#   presented in the script below and I did further processing of fonts, labels etc.
# NOTE: YOU WILL NEED TO SPECIFY THE PATH TO DATA FILES BELOW BY REPLACING "/FILE/PATH/TO/" WITH YOUR ACTUAL FILE PATHS 

###########################################
## 1. Setting up Packages and Data files ##
###########################################
rm(list=ls())

## Loading packages  
require(ggplot2)
require(doBy)
require(RCurl)



# Setting up Functions that I will use for data summary
# Functions for calculating 5% and 95% quantiles for data summary 
# These functions will be called in the summaryBy function below
quantile.function.lower = function(x) {quantile(x, c(0.05))}
quantile.function.upper = function(x) {quantile(x, c(0.95))}

#  Read in fitness data
##  Reading body size data into R
data.github = getURL("https://raw.githubusercontent.com/asingh164/SexSpecificVarianceEfficacyOfSelection/master/FitnessDataNeprojectFinalSept2020.csv")
data=read.csv(text = data.github)
data=subset(data.raw, data.raw$n.total=="32") # Clean data out for any row where we counted less than 32 offspring 

# We could only estimate selection in populations where we had both high- and low-condition flies so remove all data associated with population types 1 and 5
data.hetsub = data[data$population != 1 & data$population != 5,]
high.density.summary.table = summaryBy(n.wt ~ sex + treatment + population.type, data = data.hetsub[data.hetsub$focal.individual.density.treatment == "high",])
colnames(high.density.summary.table) = c("sex", "treatment", "population.type", "n.wt.mean.high.density")
low.density.summary.table = summaryBy(n.wt ~ sex + treatment + population.type, data = data.hetsub[data.hetsub$focal.individual.density.treatment == "low",])
colnames(low.density.summary.table) = c("sex", "treatment", "population.type", "n.wt.mean.low.density")

# Calculating mean fitness for high and low density treatment (i.e., low and high condition respectively)
low.male = mean((data.hetsub[data.hetsub$sex == "m" & data.hetsub$focal.individual.density.treatment == "low",])$n.wt)
high.male = mean((data.hetsub[data.hetsub$sex == "m" & data.hetsub$focal.individual.density.treatment == "high",])$n.wt)

low.female = mean((data.hetsub[data.hetsub$sex == "f" & data.hetsub$focal.individual.density.treatment == "low",])$n.wt)
high.female = mean((data.hetsub[data.hetsub$sex == "f" & data.hetsub$focal.individual.density.treatment == "high",])$n.wt)

# Summarizing Data  and calculating  selection differentials 

# Merge the dataframes for the mean wt fitness and calculate selection differentials for each population type and sex separetly 
hetsub.summary.table = merge(high.density.summary.table, low.density.summary.table, by = c("sex","treatment","population.type"))
hetsub.summary.table$selection.differential = (hetsub.summary.table$n.wt.mean.low.density - hetsub.summary.table$n.wt.mean.high.density) / hetsub.summary.table$n.wt.mean.low.density

# Summarize data and calculate selection differentials for each sex, averaging across populations types 
hetsub.summary.table.by.sex = summaryBy(selection.differential ~ sex, FUN = c(mean), data = hetsub.summary.table)
# Adding sex averaged value 
hetsub.summary.table.by.sex = rbind(hetsub.summary.table.by.sex, c(NA,NA))
hetsub.summary.table.by.sex$sex =as.factor(c("f", "m", "both.sexes"))
hetsub.summary.table.by.sex$selection.differential.mean[3] = as.numeric(mean(hetsub.summary.table$selection.differential)) # adding selection differential for combined 

####################################################################################
## 2. Bootstrapping selection differentials for generate 95% Confidence Intervals ##
####################################################################################

###   A - Bootstrapping each sex*mating treatment*population type combination   ###

## Subset out each category for reasmpling 
# Vials 
male.vial.2.low = (data[data$sex == "m" & data$treatment == "vial" & data$population == 2 & data$focal.individual.density.treatment == "low", ])$n.wt
male.vial.3.low = (data[data$sex == "m" & data$treatment == "vial" & data$population == 3 & data$focal.individual.density.treatment == "low", ])$n.wt
male.vial.4.low = (data[data$sex == "m" & data$treatment == "vial" & data$population == 4 & data$focal.individual.density.treatment == "low", ])$n.wt
male.vial.2.high = (data[data$sex == "m" & data$treatment == "vial" & data$population == 2 & data$focal.individual.density.treatment == "high", ])$n.wt
male.vial.3.high = (data[data$sex == "m" & data$treatment == "vial" & data$population == 3 & data$focal.individual.density.treatment == "high", ])$n.wt
male.vial.4.high = (data[data$sex == "m" & data$treatment == "vial" & data$population == 4 & data$focal.individual.density.treatment == "high", ])$n.wt
female.vial.2.low = (data[data$sex == "f" & data$treatment == "vial" & data$population == 2 & data$focal.individual.density.treatment == "low", ])$n.wt
female.vial.3.low = (data[data$sex == "f" & data$treatment == "vial" & data$population == 3 & data$focal.individual.density.treatment == "low", ])$n.wt
female.vial.4.low = (data[data$sex == "f" & data$treatment == "vial" & data$population == 4 & data$focal.individual.density.treatment == "low", ])$n.wt
female.vial.2.high = (data[data$sex == "f" & data$treatment == "vial" & data$population == 2 & data$focal.individual.density.treatment == "high", ])$n.wt
female.vial.3.high = (data[data$sex == "f" & data$treatment == "vial" & data$population == 3 & data$focal.individual.density.treatment == "high", ])$n.wt
female.vial.4.high = (data[data$sex == "f" & data$treatment == "vial" & data$population == 4 & data$focal.individual.density.treatment == "high", ])$n.wt

# Cages 
male.cage.2.low = (data[data$sex == "m" & data$treatment == "cage" & data$population == 2 & data$focal.individual.density.treatment == "low", ])$n.wt
male.cage.3.low = (data[data$sex == "m" & data$treatment == "cage" & data$population == 3 & data$focal.individual.density.treatment == "low", ])$n.wt
male.cage.4.low = (data[data$sex == "m" & data$treatment == "cage" & data$population == 4 & data$focal.individual.density.treatment == "low", ])$n.wt
male.cage.2.high = (data[data$sex == "m" & data$treatment == "cage" & data$population == 2 & data$focal.individual.density.treatment == "high", ])$n.wt
male.cage.3.high = (data[data$sex == "m" & data$treatment == "cage" & data$population == 3 & data$focal.individual.density.treatment == "high", ])$n.wt
male.cage.4.high = (data[data$sex == "m" & data$treatment == "cage" & data$population == 4 & data$focal.individual.density.treatment == "high", ])$n.wt
female.cage.2.low = (data[data$sex == "f" & data$treatment == "cage" & data$population == 2 & data$focal.individual.density.treatment == "low", ])$n.wt
female.cage.3.low = (data[data$sex == "f" & data$treatment == "cage" & data$population == 3 & data$focal.individual.density.treatment == "low", ])$n.wt
female.cage.4.low = (data[data$sex == "f" & data$treatment == "cage" & data$population == 4 & data$focal.individual.density.treatment == "low", ])$n.wt
female.cage.2.high = (data[data$sex == "f" & data$treatment == "cage" & data$population == 2 & data$focal.individual.density.treatment == "high", ])$n.wt
female.cage.3.high = (data[data$sex == "f" & data$treatment == "cage" & data$population == 3 & data$focal.individual.density.treatment == "high", ])$n.wt
female.cage.4.high = (data[data$sex == "f" & data$treatment == "cage" & data$population == 4 & data$focal.individual.density.treatment == "high", ])$n.wt

# Monogamy 
male.monogamy.3.low = (data[data$sex == "m" & data$treatment == "monogamy" & data$population == 3 & data$focal.individual.density.treatment == "low", ])$n.wt
male.monogamy.3.high = (data[data$sex == "m" & data$treatment == "monogamy" & data$population == 3 & data$focal.individual.density.treatment == "high", ])$n.wt
female.monogamy.3.low = (data[data$sex == "f" & data$treatment == "monogamy" & data$population == 3 & data$focal.individual.density.treatment == "low", ])$n.wt
female.monogamy.3.high = (data[data$sex == "f" & data$treatment == "monogamy" & data$population == 3 & data$focal.individual.density.treatment == "high", ])$n.wt


##  Making empty Dataframes to hold bootstrapped data 
# Vial
#Male 
male.vial.2.low.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(male.vial.2.low.df) = "replicate"
male.vial.2.low.df$mating.regime = "vial"
male.vial.2.low.df$population = 2
male.vial.2.low.df$sex = "m"
male.vial.2.low.df$treatment = "low"
male.vial.2.low.df$n.wt = NA
male.vial.3.low.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(male.vial.3.low.df) = "replicate"
male.vial.3.low.df$mating.regime = "vial"
male.vial.3.low.df$population = 3
male.vial.3.low.df$sex = "m"
male.vial.3.low.df$treatment = "low"
male.vial.3.low.df$n.wt = NA
male.vial.4.low.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(male.vial.4.low.df) = "replicate"
male.vial.4.low.df$mating.regime = "vial"
male.vial.4.low.df$population = 4
male.vial.4.low.df$sex = "m"
male.vial.4.low.df$treatment = "low"
male.vial.4.low.df$n.wt = NA
male.vial.2.high.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(male.vial.2.high.df) = "replicate"
male.vial.2.high.df$mating.regime = "vial"
male.vial.2.high.df$population = 2
male.vial.2.high.df$sex = "m"
male.vial.2.high.df$treatment = "high"
male.vial.2.high.df$n.wt = NA
male.vial.3.high.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(male.vial.3.high.df) = "replicate"
male.vial.3.high.df$mating.regime = "vial"
male.vial.3.high.df$population = 3
male.vial.3.high.df$sex = "m"
male.vial.3.high.df$treatment = "high"
male.vial.3.high.df$n.wt = NA
male.vial.4.high.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(male.vial.4.high.df) = "replicate"
male.vial.4.high.df$mating.regime = "vial"
male.vial.4.high.df$population = 4
male.vial.4.high.df$sex = "m"
male.vial.4.high.df$treatment = "high"
male.vial.4.high.df$n.wt = NA

#female
female.vial.2.low.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(female.vial.2.low.df) = "replicate"
female.vial.2.low.df$mating.regime = "vial"
female.vial.2.low.df$population = 2
female.vial.2.low.df$sex = "f"
female.vial.2.low.df$treatment = "low"
female.vial.2.low.df$n.wt = NA
female.vial.3.low.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(female.vial.3.low.df) = "replicate"
female.vial.3.low.df$mating.regime = "vial"
female.vial.3.low.df$population = 3
female.vial.3.low.df$sex = "f"
female.vial.3.low.df$treatment = "low"
female.vial.3.low.df$n.wt = NA
female.vial.4.low.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(female.vial.4.low.df) = "replicate"
female.vial.4.low.df$mating.regime = "vial"
female.vial.4.low.df$population = 4
female.vial.4.low.df$sex = "f"
female.vial.4.low.df$treatment = "low"
female.vial.4.low.df$n.wt = NA
female.vial.2.high.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(female.vial.2.high.df) = "replicate"
female.vial.2.high.df$mating.regime = "vial"
female.vial.2.high.df$population = 2
female.vial.2.high.df$sex = "f"
female.vial.2.high.df$treatment = "high"
female.vial.2.high.df$n.wt = NA
female.vial.3.high.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(female.vial.3.high.df) = "replicate"
female.vial.3.high.df$mating.regime = "vial"
female.vial.3.high.df$population = 3
female.vial.3.high.df$sex = "f"
female.vial.3.high.df$treatment = "high"
female.vial.3.high.df$n.wt = NA
female.vial.4.high.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(female.vial.4.high.df) = "replicate"
female.vial.4.high.df$mating.regime = "vial"
female.vial.4.high.df$population = 4
female.vial.4.high.df$sex = "f"
female.vial.4.high.df$treatment = "high"
female.vial.4.high.df$n.wt = NA

# Cage
#Male 
male.cage.2.low.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(male.cage.2.low.df) = "replicate"
male.cage.2.low.df$mating.regime = "cage"
male.cage.2.low.df$population = 2
male.cage.2.low.df$sex = "m"
male.cage.2.low.df$treatment = "low"
male.cage.2.low.df$n.wt = NA
male.cage.3.low.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(male.cage.3.low.df) = "replicate"
male.cage.3.low.df$mating.regime = "cage"
male.cage.3.low.df$population = 3
male.cage.3.low.df$sex = "m"
male.cage.3.low.df$treatment = "low"
male.cage.3.low.df$n.wt = NA
male.cage.4.low.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(male.cage.4.low.df) = "replicate"
male.cage.4.low.df$mating.regime = "cage"
male.cage.4.low.df$population = 4
male.cage.4.low.df$sex = "m"
male.cage.4.low.df$treatment = "low"
male.cage.4.low.df$n.wt = NA
male.cage.2.high.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(male.cage.2.high.df) = "replicate"
male.cage.2.high.df$mating.regime = "cage"
male.cage.2.high.df$population = 2
male.cage.2.high.df$sex = "m"
male.cage.2.high.df$treatment = "high"
male.cage.2.high.df$n.wt = NA
male.cage.3.high.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(male.cage.3.high.df) = "replicate"
male.cage.3.high.df$mating.regime = "cage"
male.cage.3.high.df$population = 3
male.cage.3.high.df$sex = "m"
male.cage.3.high.df$treatment = "high"
male.cage.3.high.df$n.wt = NA
male.cage.4.high.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(male.cage.4.high.df) = "replicate"
male.cage.4.high.df$mating.regime = "cage"
male.cage.4.high.df$population = 4
male.cage.4.high.df$sex = "m"
male.cage.4.high.df$treatment = "high"
male.cage.4.high.df$n.wt = NA

# Female 
female.cage.2.low.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(female.cage.2.low.df) = "replicate"
female.cage.2.low.df$mating.regime = "cage"
female.cage.2.low.df$population = 2
female.cage.2.low.df$sex = "f"
female.cage.2.low.df$treatment = "low"
female.cage.2.low.df$n.wt = NA
female.cage.3.low.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(female.cage.3.low.df) = "replicate"
female.cage.3.low.df$mating.regime = "cage"
female.cage.3.low.df$population = 3
female.cage.3.low.df$sex = "f"
female.cage.3.low.df$treatment = "low"
female.cage.3.low.df$n.wt = NA
female.cage.4.low.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(female.cage.4.low.df) = "replicate"
female.cage.4.low.df$mating.regime = "cage"
female.cage.4.low.df$population = 4
female.cage.4.low.df$sex = "f"
female.cage.4.low.df$treatment = "low"
female.cage.4.low.df$n.wt = NA
female.cage.2.high.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(female.cage.2.high.df) = "replicate"
female.cage.2.high.df$mating.regime = "cage"
female.cage.2.high.df$population = 2
female.cage.2.high.df$sex = "f"
female.cage.2.high.df$treatment = "high"
female.cage.2.high.df$n.wt = NA
female.cage.3.high.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(female.cage.3.high.df) = "replicate"
female.cage.3.high.df$mating.regime = "cage"
female.cage.3.high.df$population = 3
female.cage.3.high.df$sex = "f"
female.cage.3.high.df$treatment = "high"
female.cage.3.high.df$n.wt = NA
female.cage.4.high.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(female.cage.4.high.df) = "replicate"
female.cage.4.high.df$mating.regime = "cage"
female.cage.4.high.df$population = 4
female.cage.4.high.df$sex = "f"
female.cage.4.high.df$treatment = "high"
female.cage.4.high.df$n.wt = NA

# Monogamy 
# Male 
male.monogamy.3.low.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(male.monogamy.3.low.df) = "replicate"
male.monogamy.3.low.df$mating.regime = "monogamy"
male.monogamy.3.low.df$population = 3
male.monogamy.3.low.df$sex = "m"
male.monogamy.3.low.df$treatment = "low"
male.monogamy.3.low.df$n.wt = NA
male.monogamy.3.high.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(male.monogamy.3.high.df) = "replicate"
male.monogamy.3.high.df$mating.regime = "monogamy"
male.monogamy.3.high.df$population = 3
male.monogamy.3.high.df$sex = "m"
male.monogamy.3.high.df$treatment = "high"
male.monogamy.3.high.df$n.wt = NA

# Female 
female.monogamy.3.low.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(female.monogamy.3.low.df) = "replicate"
female.monogamy.3.low.df$mating.regime = "monogamy"
female.monogamy.3.low.df$population = 3
female.monogamy.3.low.df$sex = "f"
female.monogamy.3.low.df$treatment = "low"
female.monogamy.3.low.df$n.wt = NA
female.monogamy.3.high.df = data.frame(matrix(seq(from = 1, to = number.bootstraps, by = 1), nrow = number.bootstraps))
names(female.monogamy.3.high.df) = "replicate"
female.monogamy.3.high.df$mating.regime = "monogamy"
female.monogamy.3.high.df$population = 3
female.monogamy.3.high.df$sex = "f"
female.monogamy.3.high.df$treatment = "high"
female.monogamy.3.high.df$n.wt = NA

#Bootstrap Loop
for (i in 1:10000){
  
  # vials
  male.vial.2.low.df[i,6] = mean(sample(male.vial.2.low, length(male.vial.2.low), replace = TRUE))
  male.vial.3.low.df[i,6] = mean(sample(male.vial.3.low, length(male.vial.3.low), replace = TRUE))
  male.vial.4.low.df[i,6] = mean(sample(male.vial.4.low, length(male.vial.4.low), replace = TRUE))
  male.vial.2.high.df[i,6] = mean(sample(male.vial.2.high, length(male.vial.2.high), replace = TRUE))
  male.vial.3.high.df[i,6] = mean(sample(male.vial.3.high, length(male.vial.3.high), replace = TRUE))
  male.vial.4.high.df[i,6] = mean(sample(male.vial.4.high, length(male.vial.4.high), replace = TRUE))
  female.vial.2.low.df[i,6] = mean(sample(female.vial.2.low, length(female.vial.2.low), replace = TRUE))
  female.vial.3.low.df[i,6] = mean(sample(female.vial.3.low, length(female.vial.3.low), replace = TRUE))
  female.vial.4.low.df[i,6] = mean(sample(female.vial.4.low, length(female.vial.4.low), replace = TRUE))
  female.vial.2.high.df[i,6] = mean(sample(female.vial.2.high, length(female.vial.2.high), replace = TRUE))
  female.vial.3.high.df[i,6] = mean(sample(female.vial.3.high, length(female.vial.3.high), replace = TRUE))
  female.vial.4.high.df[i,6] = mean(sample(female.vial.4.high, length(female.vial.4.high), replace = TRUE))
  
  # cage
  male.cage.2.low.df[i,6] = mean(sample(male.cage.2.low, length(male.cage.2.low), replace = TRUE))
  male.cage.3.low.df[i,6] = mean(sample(male.cage.3.low, length(male.cage.3.low), replace = TRUE))
  male.cage.4.low.df[i,6] = mean(sample(male.cage.4.low, length(male.cage.4.low), replace = TRUE))
  male.cage.2.high.df[i,6] = mean(sample(male.cage.2.high, length(male.cage.2.high), replace = TRUE))
  male.cage.3.high.df[i,6] = mean(sample(male.cage.3.high, length(male.cage.3.high), replace = TRUE))
  male.cage.4.high.df[i,6] = mean(sample(male.cage.4.high, length(male.cage.4.high), replace = TRUE))
  female.cage.2.low.df[i,6] = mean(sample(female.cage.2.low, length(female.cage.2.low), replace = TRUE))
  female.cage.3.low.df[i,6] = mean(sample(female.cage.3.low, length(female.cage.3.low), replace = TRUE))
  female.cage.4.low.df[i,6] = mean(sample(female.cage.4.low, length(female.cage.4.low), replace = TRUE))
  female.cage.2.high.df[i,6] = mean(sample(female.cage.2.high, length(female.cage.2.high), replace = TRUE))
  female.cage.3.high.df[i,6] = mean(sample(female.cage.3.high, length(female.cage.3.high), replace = TRUE))
  female.cage.4.high.df[i,6] = mean(sample(female.cage.4.high, length(female.cage.4.high), replace = TRUE))
  
  # #Monogamy
  male.monogamy.3.low.df[i,6] = mean(sample(male.monogamy.3.low, length(male.monogamy.3.low), replace = TRUE))
  male.monogamy.3.high.df[i,6] = mean(sample(male.monogamy.3.high, length(male.monogamy.3.high), replace = TRUE))
  female.monogamy.3.low.df[i,6] = mean(sample(female.monogamy.3.low, length(female.monogamy.3.low), replace = TRUE))
  female.monogamy.3.high.df[i,6] = mean(sample(female.monogamy.3.high, length(female.monogamy.3.high), replace = TRUE))
  
  print(i)
}

# Combine dataframes for the high and low density resampled dataframes 
low.bootstraps = rbind(male.vial.2.low.df,male.vial.3.low.df,male.vial.4.low.df,male.cage.2.low.df,male.cage.3.low.df,male.cage.4.low.df,male.monogamy.3.low.df,
                        female.vial.2.low.df,female.vial.3.low.df,female.vial.4.low.df,female.cage.2.low.df,female.cage.3.low.df,female.cage.4.low.df,female.monogamy.3.low.df)

high.bootstraps = rbind(male.vial.2.high.df,male.vial.3.high.df,male.vial.4.high.df,male.cage.2.high.df,male.cage.3.high.df,male.cage.4.high.df,male.monogamy.3.high.df,
                       female.vial.2.high.df,female.vial.3.high.df,female.vial.4.high.df,female.cage.2.high.df,female.cage.3.high.df,female.cage.4.high.df,female.monogamy.3.high.df)

# Merge dataframes so that selection differentials can be calculated 
bootstrapped.data = merge(high.bootstraps, low.bootstraps, by = c("replicate", "mating.regime", "population", "sex"), sort = FALSE)
bootstrapped.data = bootstrapped.data[,-c(5,7)]
colnames(bootstrapped.data) = c("replicate", "mating.regime", "population", "sex", "mean.wt.high", "mean.wt.low")
bootstrapped.data$sel.diff = (bootstrapped.data$mean.wt.low - bootstrapped.data$mean.wt.high) / bootstrapped.data$mean.wt.low

# Selection differentials for each of the 13 mating regimes and population types 
sel.diff.summary.individual.pops = summaryBy(sel.diff ~  mating.regime + population + sex, FUN = c(mean, quantile.function.lower, quantile.function.upper), data = bootstrapped.data)
sel.diff.summary.by.sex = summaryBy(sel.diff ~ replicate + sex, FUN = c(mean, quantile.function.lower, quantile.function.upper), data = bootstrapped.data)


###   B - Bootstrapping sel differentials avergaing over mating treatment and population type combination   ###
# Averaging across the sexes
male.sel.vector = vector(mode = "numeric", length = number.bootstraps)
female.sel.vector = vector(mode = "numeric", length = number.bootstraps)
both.sexes.sel.vector = vector(mode = "numeric", length = number.bootstraps)

for (i in 1:number.bootstraps){
  bootstrapped.data.subset = bootstrapped.data[bootstrapped.data$replicate == i,]
  male.sel.vector[i] = mean((bootstrapped.data.subset[bootstrapped.data.subset$sex == "m",])$sel.diff)
  female.sel.vector[i] = mean((bootstrapped.data.subset[bootstrapped.data.subset$sex == "f",])$sel.diff)
  both.sexes.sel.vector[i] = mean(bootstrapped.data.subset$sel.diff)
  print(i)
}

# Selection differentials averaged across all population types 
#sel.diff.summary.all.summed.tmp = summaryBy(sel.diff ~ sex+replicate, FUN = c(mean), data = bootstrapped.data)
#sel.diff.summary.all.summed = summaryBy(sel.diff.mean ~ sex, FUN=c(mean, quantile.function.lower, quantile.function.upper), data = sel.diff.summary.all.summed.tmp)
#sel.diff.summary.all.summed

##   Plots

# Plotting all pop and mating regimes combinations 

#######################################################################
## 3. Plotting Selection Differentials with 95% Confidence Intervals ##
#######################################################################

# A - First plot is selection differentials seperately for each sex and population type 
# Merge read means with confidence intervals 
sel.diff.data.merge = merge(hetsub.summary.table, sel.diff.summary.individual.pops, by.x = c("sex", "treatment", "population.type"), by.y = c("sex", "mating.regime", "population"), sort = FALSE)
sel.diff.data.merge$treatment.id = paste(sel.diff.data.merge$sex, sel.diff.data.merge$treatment, sel.diff.data.merge$population, sep = ":")
sel.diff.data.merge = sel.diff.data.merge[,-c(4,5,7)]
colnames(sel.diff.data.merge) = c("sex", "treatment", "population.type", "sel.diff", "lower.ci", "upper.ci", "ID")
selection.plot.all.seperate = ggplot(sel.diff.data.merge, aes(y = sel.diff, x = as.factor(ID))) + 
  geom_point(aes(y = sel.diff, x = as.factor(ID)), size = 5) + 
  geom_errorbar(aes(ymin = sel.diff.data.merge$lower.ci, ymax = sel.diff.data.merge$upper.ci, width = 0.0)) +
  theme_bw() + ylim(-1,1) +
  theme(panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank(),panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank()) +
  geom_hline(yintercept=0)   
selection.plot.all.seperate 

# B - Second plot is selection differentials for each sex and combined sexes averaged over all mating regimes and population types 
# Generate CIs for estimates and add to data frames 
hetsub.summary.table.by.sex$upper.ci = c(as.vector(quantile(female.sel.vector, c(0.95))), as.vector(quantile(male.sel.vector, c(0.95))), as.vector(quantile(both.sexes.sel.vector, c(0.95))))
hetsub.summary.table.by.sex$lower.ci = c(as.vector(quantile(female.sel.vector, c(0.05))), as.vector(quantile(male.sel.vector, c(0.05))), as.vector(quantile(both.sexes.sel.vector, c(0.05))))

selection.plot.all = ggplot(hetsub.summary.table.by.sex, aes(y = as.numeric(selection.differential.mean), x = sex)) + 
  geom_point(aes(y = as.numeric(selection.differential.mean), x = sex),size = 5) + 
  geom_errorbar(aes(ymin = hetsub.summary.table.by.sex$lower.ci, ymax = hetsub.summary.table.by.sex$upper.ci, width = 0.0)) +
  theme_bw() + ylim(-1,1) +
  theme(panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank(),panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank()) +
  geom_hline(yintercept=0)  
selection.plot.all 


# Save both plots
#pdf("/FILE/PATH/TO/Supplemental_Figure3A_SelDiff_Individual_Pops.pdf", width = 10, height = 8)
#selection.plot.all.seperate
#dev.off()
#pdf("/FILE/PATH/TO/Supplemental_Figure3B_SelDiff_Averaged_Across_Pops.pdf", width = 4, height = 10)
#selection.plot.all
#dev.off()


#
