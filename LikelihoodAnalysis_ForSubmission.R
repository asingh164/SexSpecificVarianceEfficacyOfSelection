#
###########################################################################################################################################################
###########################################################################################################################################################
###   Script Details:                                                                                                                                   ###
###   Likelihood Analysis of male:female variance for Singh, A. A.F., Agrawal. 2020. Sex-Specific Variance in Fitness and the Efficacy of Selection     ###
###   Script written by Amardeep Singh (amardeep.singh[at]utoronto.ca)                                                                                  ###
###   This script will run our likelihood procedure to esitimate likelihood functions for both female variance in fitness and the                       ###
###   phi parameter which is the factor by which male fitness variance is greater than male fitness variance                                           ###
###########################################################################################################################################################
###########################################################################################################################################################
rm(list=ls())

# This script will generate plots corresponding to Figure 3, and Supplemental Figures S4-S10. Corresponding figures are denoted in the script below:
#   Figure 3: Line 518 ## Line numbers assume that this line is on line 14, if not add or subtract the amount to get this line onto line 14
#   Figure S4: Line 485
#   Figure S5: Line 434
#   Figure S6: Line 383
#   Figure S7: Line 749
#   Figure S8: Line 698 
#   Figure S9: Line 647
#   Figure S10: Line 782
# Figures may not look exactly as those presented in the paper. Although the data is identical, I may have switched around the order differently than
#   presented in the script below and I did further processing of fonts, labels etc.


###########################################
## 1. Setting up Packages and Data files ##
###########################################

# Loading packages 
require(doBy)
require(Rmpfr)
require(lattice)
require(reshape2)
require(gridExtra)

# Read in fitness data 
data=read.csv("/FILE/PATH/TO/FitnessDataNeprojectFinalSept2020.csv", header = TRUE, sep = ",")
data.summary = summaryBy(n.wt ~ as.factor(treatment) + as.factor(population.type) + as.factor(sex), data = data, FUN = c(mean, var), na.rm = TRUE)


####################################################
## 2. Resampling data and Calculating Likelihoods ##
####################################################

#Sample from fitness data to bring the mean down to 2
v=c()
for (i in 1:nrow(data.summary)){
  if (data.summary$n.wt.mean[i] > 2){
    v[i] = 2/data.summary$n.wt.mean[i]
  } else {
      v[i]=1
  }
}
data.summary$v = v
rownames(data.summary) = seq(from = 1, to = nrow(data.summary), by = 1)

#Binomial Sampling procedure to generate a resampled dataset 
resample.data = subset(data, data$n.wt != "NA") 

rownames(resample.data) = seq(from = 1, to = nrow(resample.data))
resample.data$g = NA

for (i in 1:nrow(resample.data)){
  treatment = as.character(resample.data[i,"treatment"])
  population = as.numeric(as.character(resample.data[i,"population.type"]))
  sex = as.character(resample.data[i,"sex"])
  
  mean.n.wt = (data.summary[data.summary$treatment %in% treatment & data.summary$sex %in% sex & data.summary$population %in% population,])$v
  
  if (as.numeric(as.character(resample.data[i,10])) == 0){
    resample.data$g[i] = 0
  } else {
    resample.data$g[i] = rbinom(n = 1, size = as.numeric(as.character(resample.data$n.wt[i])), prob = mean.n.wt)
  }
}

#Likelihood estimates ## each mating regime seperately  
# Sequence of values for variance in female fitness and omega and phi for our likelihood to be evaluated over
VarF.list = seq(from = 2.01, to = 16, by = 0.02)
omega.list = seq(from = -2, to = 4, by = 0.02)
# phi will be set at 2^omega
phi.list = 2^omega.list

# Setting matricies for saving likelihood values

# First likelihood matrix will be averaging across all population types within a mating regime  
likelihood.matrix.all = matrix(nrow = length(VarF.list), ncol = length(omega.list))
colnames(likelihood.matrix.all) = 2^(omega.list)
rownames(likelihood.matrix.all) = VarF.list
likelihood.matrix.vial = matrix(nrow = length(VarF.list), ncol = length(omega.list))
colnames(likelihood.matrix.vial) = 2^(omega.list)
rownames(likelihood.matrix.vial) = VarF.list
likelihood.matrix.cage = matrix(nrow = length(VarF.list), ncol = length(omega.list))
colnames(likelihood.matrix.cage) = 2^(omega.list)
rownames(likelihood.matrix.cage) = VarF.list
likelihood.matrix.monogamy = matrix(nrow = length(VarF.list), ncol = length(omega.list))
colnames(likelihood.matrix.monogamy) = 2^(omega.list)
rownames(likelihood.matrix.monogamy) = VarF.list

# These likelihoot matricies are breaking down each mating regime into separate population types 
likelihood.matrix.pop1.vial = matrix(nrow = length(VarF.list), ncol = length(omega.list))
colnames(likelihood.matrix.pop1.vial) = 2^(omega.list)
rownames(likelihood.matrix.pop1.vial) = VarF.list
likelihood.matrix.pop2.vial = matrix(nrow = length(VarF.list), ncol = length(omega.list))
colnames(likelihood.matrix.pop2.vial) = 2^(omega.list)
rownames(likelihood.matrix.pop2.vial) = VarF.list
likelihood.matrix.pop3.vial = matrix(nrow = length(VarF.list), ncol = length(omega.list))
colnames(likelihood.matrix.pop3.vial) = 2^(omega.list)
rownames(likelihood.matrix.pop3.vial) = VarF.list
likelihood.matrix.pop4.vial = matrix(nrow = length(VarF.list), ncol = length(omega.list))
colnames(likelihood.matrix.pop4.vial) = 2^(omega.list)
rownames(likelihood.matrix.pop4.vial) = VarF.list
likelihood.matrix.pop5.vial = matrix(nrow = length(VarF.list), ncol = length(omega.list))
colnames(likelihood.matrix.pop5.vial) = 2^(omega.list)
rownames(likelihood.matrix.pop5.vial) = VarF.list

likelihood.matrix.pop1.cage = matrix(nrow = length(VarF.list), ncol = length(omega.list))
colnames(likelihood.matrix.pop1.cage) = 2^(omega.list)
rownames(likelihood.matrix.pop1.cage) = VarF.list
likelihood.matrix.pop2.cage = matrix(nrow = length(VarF.list), ncol = length(omega.list))
colnames(likelihood.matrix.pop2.cage) = 2^(omega.list)
rownames(likelihood.matrix.pop2.cage) = VarF.list
likelihood.matrix.pop3.cage = matrix(nrow = length(VarF.list), ncol = length(omega.list))
colnames(likelihood.matrix.pop3.cage) = 2^(omega.list)
rownames(likelihood.matrix.pop3.cage) = VarF.list
likelihood.matrix.pop4.cage = matrix(nrow = length(VarF.list), ncol = length(omega.list))
colnames(likelihood.matrix.pop4.cage) = 2^(omega.list)
rownames(likelihood.matrix.pop4.cage) = VarF.list
likelihood.matrix.pop5.cage = matrix(nrow = length(VarF.list), ncol = length(omega.list))
colnames(likelihood.matrix.pop5.cage) = 2^(omega.list)
rownames(likelihood.matrix.pop5.cage) = VarF.list

likelihood.matrix.pop1.monogamy = matrix(nrow = length(VarF.list), ncol = length(omega.list))
colnames(likelihood.matrix.pop1.monogamy) = 2^(omega.list)
rownames(likelihood.matrix.pop1.monogamy) = VarF.list
likelihood.matrix.pop3.monogamy = matrix(nrow = length(VarF.list), ncol = length(omega.list))
colnames(likelihood.matrix.pop3.monogamy) = 2^(omega.list)
rownames(likelihood.matrix.pop3.monogamy) = VarF.list
likelihood.matrix.pop5.monogamy = matrix(nrow = length(VarF.list), ncol = length(omega.list))
colnames(likelihood.matrix.pop5.monogamy) = 2^(omega.list)
rownames(likelihood.matrix.pop5.monogamy) = VarF.list

# Likelihood function
# L=dnbinom(x = g[1], size = (4 / (VarF-2)), prob = 2/VarF)

resample.data.vial.pop1.male = resample.data[resample.data$treatment == "vial" & resample.data$sex == "m" & resample.data$population == 1, ] 
resample.data.vial.pop1.female = resample.data[resample.data$treatment == "vial" & resample.data$sex == "f" & resample.data$population == 1, ] 
resample.data.vial.pop2.male = resample.data[resample.data$treatment == "vial" & resample.data$sex == "m" & resample.data$population == 2, ] 
resample.data.vial.pop2.female = resample.data[resample.data$treatment == "vial" & resample.data$sex == "f" & resample.data$population == 2, ] 
resample.data.vial.pop3.male = resample.data[resample.data$treatment == "vial" & resample.data$sex == "m" & resample.data$population == 3, ] 
resample.data.vial.pop3.female = resample.data[resample.data$treatment == "vial" & resample.data$sex == "f" & resample.data$population == 3, ] 
resample.data.vial.pop4.male = resample.data[resample.data$treatment == "vial" & resample.data$sex == "m" & resample.data$population == 4, ] 
resample.data.vial.pop4.female = resample.data[resample.data$treatment == "vial" & resample.data$sex == "f" & resample.data$population == 4, ] 
resample.data.vial.pop5.male = resample.data[resample.data$treatment == "vial" & resample.data$sex == "m" & resample.data$population == 5, ] 
resample.data.vial.pop5.female = resample.data[resample.data$treatment == "vial" & resample.data$sex == "f" & resample.data$population == 5, ] 

resample.data.cage.pop1.male = resample.data[resample.data$treatment == "cage" & resample.data$sex == "m" & resample.data$population == 1, ] 
resample.data.cage.pop1.female = resample.data[resample.data$treatment == "cage" & resample.data$sex == "f" & resample.data$population == 1, ] 
resample.data.cage.pop2.male = resample.data[resample.data$treatment == "cage" & resample.data$sex == "m" & resample.data$population == 2, ] 
resample.data.cage.pop2.female = resample.data[resample.data$treatment == "cage" & resample.data$sex == "f" & resample.data$population == 2, ] 
resample.data.cage.pop3.male = resample.data[resample.data$treatment == "cage" & resample.data$sex == "m" & resample.data$population == 3, ] 
resample.data.cage.pop3.female = resample.data[resample.data$treatment == "cage" & resample.data$sex == "f" & resample.data$population == 3, ] 
resample.data.cage.pop4.male = resample.data[resample.data$treatment == "cage" & resample.data$sex == "m" & resample.data$population == 4, ] 
resample.data.cage.pop4.female = resample.data[resample.data$treatment == "cage" & resample.data$sex == "f" & resample.data$population == 4, ] 
resample.data.cage.pop5.male = resample.data[resample.data$treatment == "cage" & resample.data$sex == "m" & resample.data$population == 5, ] 
resample.data.cage.pop5.female = resample.data[resample.data$treatment == "cage" & resample.data$sex == "f" & resample.data$population == 5, ] 

resample.data.monogamy.pop1.male = resample.data[resample.data$treatment == "monogamy" & resample.data$sex == "m" & resample.data$population == 1, ] 
resample.data.monogamy.pop1.female = resample.data[resample.data$treatment == "monogamy" & resample.data$sex == "f" & resample.data$population == 1, ] 
resample.data.monogamy.pop3.male = resample.data[resample.data$treatment == "monogamy" & resample.data$sex == "m" & resample.data$population == 3, ] 
resample.data.monogamy.pop3.female = resample.data[resample.data$treatment == "monogamy" & resample.data$sex == "f" & resample.data$population == 3, ] 
resample.data.monogamy.pop5.male = resample.data[resample.data$treatment == "monogamy" & resample.data$sex == "m" & resample.data$population == 5, ] 
resample.data.monogamy.pop5.female = resample.data[resample.data$treatment == "monogamy" & resample.data$sex == "f" & resample.data$population == 5, ] 

resample.data.vial.male = resample.data[resample.data$treatment == "vial" & resample.data$sex == "m", ] 
resample.data.cage.male = resample.data[resample.data$treatment == "cage" & resample.data$sex == "m", ] 
resample.data.monogamy.male = resample.data[resample.data$treatment == "monogamy" & resample.data$sex == "m", ] 

resample.data.vial.female = resample.data[resample.data$treatment == "vial" & resample.data$sex == "f", ] 
resample.data.cage.female = resample.data[resample.data$treatment == "cage" & resample.data$sex == "f", ] 
resample.data.monogamy.female = resample.data[resample.data$treatment == "monogamy" & resample.data$sex == "f", ] 

# Looping over each pairwie combination of varF and phi to calculate the likeihood of combinations of both
for (i in 1:length(VarF.list)){
  for (j in 1:length(omega.list)){
    gamma = 2^(omega.list[j])
    
    #Females 
    resample.data.vial.pop1.female$likelihood = apply(resample.data.vial.pop1.female, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i] - 2),prob =  2/VarF.list[i])))
    resample.data.vial.pop2.female$likelihood = apply(resample.data.vial.pop2.female, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i] - 2),prob =  2/VarF.list[i])))
    resample.data.vial.pop3.female$likelihood = apply(resample.data.vial.pop3.female, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i] - 2),prob =  2/VarF.list[i])))
    resample.data.vial.pop4.female$likelihood = apply(resample.data.vial.pop4.female, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i] - 2),prob =  2/VarF.list[i])))
    resample.data.vial.pop5.female$likelihood = apply(resample.data.vial.pop5.female, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i] - 2),prob =  2/VarF.list[i])))
    
    resample.data.cage.pop1.female$likelihood = apply(resample.data.cage.pop1.female, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i] - 2),prob =  2/VarF.list[i])))
    resample.data.cage.pop2.female$likelihood = apply(resample.data.cage.pop2.female, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i] - 2),prob =  2/VarF.list[i])))
    resample.data.cage.pop3.female$likelihood = apply(resample.data.cage.pop3.female, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i] - 2),prob =  2/VarF.list[i])))
    resample.data.cage.pop4.female$likelihood = apply(resample.data.cage.pop4.female, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i] - 2),prob =  2/VarF.list[i])))
    resample.data.cage.pop5.female$likelihood = apply(resample.data.cage.pop5.female, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i] - 2),prob =  2/VarF.list[i])))
    
    resample.data.monogamy.pop1.female$likelihood = apply(resample.data.monogamy.pop1.female, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i] - 2),prob =  2/VarF.list[i])))
    resample.data.monogamy.pop3.female$likelihood = apply(resample.data.monogamy.pop3.female, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i] - 2),prob =  2/VarF.list[i])))
    resample.data.monogamy.pop5.female$likelihood = apply(resample.data.monogamy.pop5.female, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i] - 2),prob =  2/VarF.list[i])))

    if ((VarF.list[i] * gamma) > 2.01){
    
    resample.data.vial.pop1.male$likelihood = apply(resample.data.vial.pop1.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i]*gamma - 2), prob =  2/(VarF.list[i]*gamma))))
    resample.data.vial.pop2.male$likelihood = apply(resample.data.vial.pop2.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i]*gamma - 2), prob =  2/(VarF.list[i]*gamma))))
    resample.data.vial.pop3.male$likelihood = apply(resample.data.vial.pop3.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i]*gamma - 2), prob =  2/(VarF.list[i]*gamma))))
    resample.data.vial.pop4.male$likelihood = apply(resample.data.vial.pop4.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i]*gamma - 2), prob =  2/(VarF.list[i]*gamma))))
    resample.data.vial.pop5.male$likelihood = apply(resample.data.vial.pop5.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i]*gamma - 2), prob =  2/(VarF.list[i]*gamma))))
    
    resample.data.cage.pop1.male$likelihood = apply(resample.data.cage.pop1.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i]*gamma - 2), prob =  2/(VarF.list[i]*gamma))))
    resample.data.cage.pop2.male$likelihood = apply(resample.data.cage.pop2.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i]*gamma - 2), prob =  2/(VarF.list[i]*gamma))))
    resample.data.cage.pop3.male$likelihood = apply(resample.data.cage.pop3.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i]*gamma - 2), prob =  2/(VarF.list[i]*gamma))))
    resample.data.cage.pop4.male$likelihood = apply(resample.data.cage.pop4.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i]*gamma - 2), prob =  2/(VarF.list[i]*gamma))))
    resample.data.cage.pop5.male$likelihood = apply(resample.data.cage.pop5.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i]*gamma - 2), prob =  2/(VarF.list[i]*gamma))))
    
    resample.data.monogamy.pop1.male$likelihood = apply(resample.data.monogamy.pop1.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i]*gamma - 2),prob =  2/(VarF.list[i]*gamma))))
    resample.data.monogamy.pop3.male$likelihood = apply(resample.data.monogamy.pop3.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i]*gamma - 2),prob =  2/(VarF.list[i]*gamma))))
    resample.data.monogamy.pop5.male$likelihood = apply(resample.data.monogamy.pop5.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(VarF.list[i]*gamma - 2),prob =  2/(VarF.list[i]*gamma))))
    
    } else {
    
    resample.data.vial.pop1.male$likelihood = apply(resample.data.vial.pop1.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(2.01 - 2), prob =  2/(2.01))))
    resample.data.vial.pop2.male$likelihood = apply(resample.data.vial.pop2.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(2.01 - 2), prob =  2/(2.01))))
    resample.data.vial.pop3.male$likelihood = apply(resample.data.vial.pop3.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(2.01 - 2), prob =  2/(2.01))))
    resample.data.vial.pop4.male$likelihood = apply(resample.data.vial.pop4.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(2.01 - 2), prob =  2/(2.01))))
    resample.data.vial.pop5.male$likelihood = apply(resample.data.vial.pop5.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(2.01 - 2), prob =  2/(2.01))))
    
    resample.data.cage.pop1.male$likelihood = apply(resample.data.cage.pop1.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(2.01 - 2), prob =  2/(2.01))))
    resample.data.cage.pop2.male$likelihood = apply(resample.data.cage.pop2.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(2.01 - 2), prob =  2/(2.01))))
    resample.data.cage.pop3.male$likelihood = apply(resample.data.cage.pop3.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(2.01 - 2), prob =  2/(2.01))))
    resample.data.cage.pop4.male$likelihood = apply(resample.data.cage.pop4.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(2.01 - 2), prob =  2/(2.01))))
    resample.data.cage.pop5.male$likelihood = apply(resample.data.cage.pop5.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(2.01 - 2), prob =  2/(2.01))))
    
    resample.data.monogamy.pop1.male$likelihood = apply(resample.data.monogamy.pop1.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(2.01 - 2), prob =  2/(2.01))))
    resample.data.monogamy.pop3.male$likelihood = apply(resample.data.monogamy.pop3.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(2.01 - 2), prob =  2/(2.01))))
    resample.data.monogamy.pop5.male$likelihood = apply(resample.data.monogamy.pop5.male, 1, function(x) log(dnbinom(x = as.numeric(x[16]), size = 4/(2.01 - 2), prob =  2/(2.01))))
    
    }
    
    likelihood.matrix.pop1.vial[i,j] = (sum(resample.data.vial.pop1.female$likelihood,resample.data.vial.pop1.male$likelihood))
    likelihood.matrix.pop2.vial[i,j] = (sum(resample.data.vial.pop2.female$likelihood,resample.data.vial.pop2.male$likelihood))
    likelihood.matrix.pop3.vial[i,j] = (sum(resample.data.vial.pop3.female$likelihood,resample.data.vial.pop3.male$likelihood))
    likelihood.matrix.pop4.vial[i,j] = (sum(resample.data.vial.pop4.female$likelihood,resample.data.vial.pop4.male$likelihood))
    likelihood.matrix.pop5.vial[i,j] = (sum(resample.data.vial.pop5.female$likelihood,resample.data.vial.pop5.male$likelihood))
    
    likelihood.matrix.pop1.cage[i,j] = (sum(resample.data.cage.pop1.female$likelihood,resample.data.cage.pop1.male$likelihood))
    likelihood.matrix.pop2.cage[i,j] = (sum(resample.data.cage.pop2.female$likelihood,resample.data.cage.pop2.male$likelihood))
    likelihood.matrix.pop3.cage[i,j] = (sum(resample.data.cage.pop3.female$likelihood,resample.data.cage.pop3.male$likelihood))
    likelihood.matrix.pop4.cage[i,j] = (sum(resample.data.cage.pop4.female$likelihood,resample.data.cage.pop4.male$likelihood))
    likelihood.matrix.pop5.cage[i,j] = (sum(resample.data.cage.pop5.female$likelihood,resample.data.cage.pop5.male$likelihood))
    
    likelihood.matrix.pop1.monogamy[i,j] = (sum(resample.data.monogamy.pop1.female$likelihood,resample.data.monogamy.pop1.male$likelihood))
    likelihood.matrix.pop3.monogamy[i,j] = (sum(resample.data.monogamy.pop3.female$likelihood,resample.data.monogamy.pop3.male$likelihood))
    likelihood.matrix.pop5.monogamy[i,j] = (sum(resample.data.monogamy.pop5.female$likelihood,resample.data.monogamy.pop5.male$likelihood))
    
    
  }
  print(c(i,j))
}
# Save likelhood matricies as RDS objects 
#saveRDS(likelihood.matrix.pop1.vial, "/FILE/PATH/TO/likelihood.matrix.pop1.vial.rds")
#saveRDS(likelihood.matrix.pop2.vial, "/FILE/PATH/TO/likelihood.matrix.pop2.vial.rds")
#saveRDS(likelihood.matrix.pop3.vial, "/FILE/PATH/TO/likelihood.matrix.pop3.vial.rds")
#saveRDS(likelihood.matrix.pop4.vial, "/FILE/PATH/TO/likelihood.matrix.pop4.vial.rds")
#saveRDS(likelihood.matrix.pop5.vial, "/FILE/PATH/TO/likelihood.matrix.pop5.vial.rds")

#saveRDS(likelihood.matrix.pop1.cage, "/FILE/PATH/TO/likelihood.matrix.pop1.cage.rds")
#saveRDS(likelihood.matrix.pop2.cage, "/FILE/PATH/TO/likelihood.matrix.pop2.cage.rds")
#saveRDS(likelihood.matrix.pop3.cage, "/FILE/PATH/TO/likelihood.matrix.pop3.cage.rds")
#saveRDS(likelihood.matrix.pop4.cage, "/FILE/PATH/TO/likelihood.matrix.pop4.cage.rds")
#saveRDS(likelihood.matrix.pop5.cage, "/FILE/PATH/TO/likelihood.matrix.pop5.cage.rds")

#saveRDS(likelihood.matrix.pop1.monogamy, "/FILE/PATH/TO/likelihood.matrix.pop1.monogamy.rds")
#saveRDS(likelihood.matrix.pop3.monogamy, "/FILE/PATH/TO/likelihood.matrix.pop3.monogamy.rds")
#saveRDS(likelihood.matrix.pop5.monogamy, "/FILE/PATH/TO/likelihood.matrix.pop5.monogamy.rds")

#saveRDS(likelihood.matrix.vial, "/FILE/PATH/TO/likelihood.matrix.vial.rds")
#saveRDS(likelihood.matrix.cage, "/FILE/PATH/TO/likelihood.matrix.cage.rds")
#saveRDS(likelihood.matrix.monogamy, "/FILE/PATH/TO/likelihood.matrix.monogamy.rds")


########################################################
## 3. Plotting Likelihood Profiles for VarF and phi  ##
########################################################

####    A. Loading in likelihood matrix RDS files   ####

likelihood.matrix.pop1.vial = readRDS("/FILE/PATH/TO/likelihood.matrix.pop1.vial.rds")
likelihood.matrix.pop2.vial = readRDS("/FILE/PATH/TO/likelihood.matrix.pop2.vial.rds")
likelihood.matrix.pop3.vial = readRDS("/FILE/PATH/TO/likelihood.matrix.pop3.vial.rds")
likelihood.matrix.pop4.vial = readRDS("/FILE/PATH/TO/likelihood.matrix.pop4.vial.rds")
likelihood.matrix.pop5.vial = readRDS("/FILE/PATH/TO/likelihood.matrix.pop5.vial.rds")

likelihood.matrix.pop1.cage = readRDS("/FILE/PATH/TO/likelihood.matrix.pop1.cage.rds")
likelihood.matrix.pop2.cage = readRDS("/FILE/PATH/TO/likelihood.matrix.pop2.cage.rds")
likelihood.matrix.pop3.cage = readRDS("/FILE/PATH/TO/likelihood.matrix.pop3.cage.rds")
likelihood.matrix.pop4.cage = readRDS("/FILE/PATH/TO/likelihood.matrix.pop4.cage.rds")
likelihood.matrix.pop5.cage = readRDS("/FILE/PATH/TO/likelihood.matrix.pop5.cage.rds")

likelihood.matrix.pop1.monogamy = readRDS("/FILE/PATH/TO/likelihood.matrix.pop1.monogamy.rds")
likelihood.matrix.pop3.monogamy = readRDS( "/FILE/PATH/TO/likelihood.matrix.pop3.monogamy.rds")
likelihood.matrix.pop5.monogamy = readRDS( "/FILE/PATH/TO/likelihood.matrix.pop5.monogamy.rds")

# Averaged across population types
likelihood.matrix.vial = likelihood.matrix.pop1.vial+likelihood.matrix.pop2.vial+likelihood.matrix.pop3.vial+likelihood.matrix.pop4.vial+likelihood.matrix.pop5.vial
likelihood.matrix.cage = likelihood.matrix.pop1.cage+likelihood.matrix.pop2.cage+likelihood.matrix.pop3.cage+likelihood.matrix.pop4.cage+likelihood.matrix.pop5.cage
likelihood.matrix.monogamy = likelihood.matrix.pop1.monogamy + likelihood.matrix.pop3.monogamy + likelihood.matrix.pop5.monogamy

likelihood.matrix.polygamy = likelihood.matrix.vial + likelihood.matrix.cage


####    B. Calculating likelihood profiles of phi for each mating regime by population type combination   ####
# Only plotted for a range of values

# Vials 
likelihood.matrix.pop1.vial.phi.profile = as.data.frame((apply(likelihood.matrix.pop1.vial, 2, max))[51:217])
likelihood.matrix.pop1.vial.phi.profile$phi = row.names(likelihood.matrix.pop1.vial.phi.profile)
colnames(likelihood.matrix.pop1.vial.phi.profile)[1] = "log.likelihood"
likelihood.matrix.pop2.vial.phi.profile = as.data.frame((apply(likelihood.matrix.pop2.vial, 2, max))[51:217])
likelihood.matrix.pop2.vial.phi.profile$phi = row.names(likelihood.matrix.pop2.vial.phi.profile)
colnames(likelihood.matrix.pop2.vial.phi.profile)[1] = "log.likelihood"
likelihood.matrix.pop3.vial.phi.profile = as.data.frame((apply(likelihood.matrix.pop3.vial, 2, max))[51:217])
likelihood.matrix.pop3.vial.phi.profile$phi = row.names(likelihood.matrix.pop3.vial.phi.profile)
colnames(likelihood.matrix.pop3.vial.phi.profile)[1] = "log.likelihood"
likelihood.matrix.pop4.vial.phi.profile = as.data.frame((apply(likelihood.matrix.pop4.vial, 2, max))[51:217])
likelihood.matrix.pop4.vial.phi.profile$phi = row.names(likelihood.matrix.pop4.vial.phi.profile)
colnames(likelihood.matrix.pop4.vial.phi.profile)[1] = "log.likelihood"
likelihood.matrix.pop5.vial.phi.profile = as.data.frame((apply(likelihood.matrix.pop5.vial, 2, max))[51:217])
likelihood.matrix.pop5.vial.phi.profile$phi = row.names(likelihood.matrix.pop5.vial.phi.profile)
colnames(likelihood.matrix.pop5.vial.phi.profile)[1] = "log.likelihood"

# Cages 
likelihood.matrix.pop1.cage.phi.profile = as.data.frame((apply(likelihood.matrix.pop1.cage, 2, max))[51:217])
likelihood.matrix.pop1.cage.phi.profile$phi = row.names(likelihood.matrix.pop1.cage.phi.profile)
colnames(likelihood.matrix.pop1.cage.phi.profile)[1] = "log.likelihood"
likelihood.matrix.pop2.cage.phi.profile = as.data.frame((apply(likelihood.matrix.pop2.cage, 2, max))[51:217])
likelihood.matrix.pop2.cage.phi.profile$phi = row.names(likelihood.matrix.pop2.cage.phi.profile)
colnames(likelihood.matrix.pop2.cage.phi.profile)[1] = "log.likelihood"
likelihood.matrix.pop3.cage.phi.profile = as.data.frame((apply(likelihood.matrix.pop3.cage, 2, max))[51:217])
likelihood.matrix.pop3.cage.phi.profile$phi = row.names(likelihood.matrix.pop3.cage.phi.profile)
colnames(likelihood.matrix.pop3.cage.phi.profile)[1] = "log.likelihood"
likelihood.matrix.pop4.cage.phi.profile = as.data.frame((apply(likelihood.matrix.pop4.cage, 2, max))[51:217])
likelihood.matrix.pop4.cage.phi.profile$phi = row.names(likelihood.matrix.pop4.cage.phi.profile)
colnames(likelihood.matrix.pop4.cage.phi.profile)[1] = "log.likelihood"
likelihood.matrix.pop5.cage.phi.profile = as.data.frame((apply(likelihood.matrix.pop5.cage, 2, max))[51:217])
likelihood.matrix.pop5.cage.phi.profile$phi = row.names(likelihood.matrix.pop5.cage.phi.profile)
colnames(likelihood.matrix.pop5.cage.phi.profile)[1] = "log.likelihood"

# Monogamy 
likelihood.matrix.pop1.monogamy.phi.profile = as.data.frame((apply(likelihood.matrix.pop1.monogamy, 2, max))[51:217])
likelihood.matrix.pop1.monogamy.phi.profile$phi = row.names(likelihood.matrix.pop1.monogamy.phi.profile)
colnames(likelihood.matrix.pop1.monogamy.phi.profile)[1] = "log.likelihood"
likelihood.matrix.pop3.monogamy.phi.profile = as.data.frame((apply(likelihood.matrix.pop3.monogamy, 2, max))[51:217])
likelihood.matrix.pop3.monogamy.phi.profile$phi = row.names(likelihood.matrix.pop3.monogamy.phi.profile)
colnames(likelihood.matrix.pop3.monogamy.phi.profile)[1] = "log.likelihood"
likelihood.matrix.pop5.monogamy.phi.profile = as.data.frame((apply(likelihood.matrix.pop5.monogamy, 2, max))[51:217])
likelihood.matrix.pop5.monogamy.phi.profile$phi = row.names(likelihood.matrix.pop5.monogamy.phi.profile)
colnames(likelihood.matrix.pop5.monogamy.phi.profile)[1] = "log.likelihood"

# Each mating regime over all population types
likelihood.matrix.vial.phi.profile = as.data.frame((apply(likelihood.matrix.vial, 2, max))[51:217])
likelihood.matrix.vial.phi.profile$phi = row.names(likelihood.matrix.vial.phi.profile)
colnames(likelihood.matrix.vial.phi.profile)[1] = "log.likelihood"
likelihood.matrix.cage.phi.profile = as.data.frame((apply(likelihood.matrix.cage, 2, max))[51:217])
likelihood.matrix.cage.phi.profile$phi = row.names(likelihood.matrix.cage.phi.profile)
colnames(likelihood.matrix.cage.phi.profile)[1] = "log.likelihood"
likelihood.matrix.monogamy.phi.profile = as.data.frame((apply(likelihood.matrix.monogamy, 2, max))[51:217])
likelihood.matrix.monogamy.phi.profile$phi = row.names(likelihood.matrix.monogamy.phi.profile)
colnames(likelihood.matrix.monogamy.phi.profile)[1] = "log.likelihood"

# Polygamy to constrast with monogamy treatment 
likelihood.matrix.polygamy.phi.profile = as.data.frame((apply(likelihood.matrix.polygamy, 2, max))[51:217])
likelihood.matrix.polygamy.phi.profile$phi = row.names(likelihood.matrix.polygamy.phi.profile)
colnames(likelihood.matrix.polygamy.phi.profile)[1] = "log.likelihood"


####    C. Plotting Likelhiood Profiles for phi    ####

## Plotting for Vials ### FIGURE S6
#pdf("/FILE/PATH/TO/FIGURES6.vials.likelihood.phi..pdf", width = 10, height = 8)  # Run if you want to save the output as a pdf
par(mfrow = c(2,3))
# Vial Pop 1 
plot(y = (likelihood.matrix.pop1.vial.phi.profile$log.likelihood), x = likelihood.matrix.pop1.vial.phi.profile$phi, type = "l", xaxt = "n")
axis(1, at = seq(0.5,5,by = 0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop1.vial.phi.profile.ordered = (likelihood.matrix.pop1.vial.phi.profile[order(likelihood.matrix.pop1.vial.phi.profile$log.likelihood),])[(likelihood.matrix.pop1.vial.phi.profile[order(likelihood.matrix.pop1.vial.phi.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop1.vial.phi.profile[order(likelihood.matrix.pop1.vial.phi.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop1.vial.phi.profile.ordered.range = range(likelihood.matrix.pop1.vial.phi.profile.ordered$phi)
abline(v = (subset(likelihood.matrix.pop1.vial.phi.profile.ordered, likelihood.matrix.pop1.vial.phi.profile.ordered$log.likelihood == max(likelihood.matrix.pop1.vial.phi.profile.ordered$log.likelihood)))$phi)
abline(v = c(likelihood.matrix.pop1.vial.phi.profile.ordered.range[1], likelihood.matrix.pop1.vial.phi.profile.ordered.range[2]), lty = 2)
# Vial Pop 2
plot(y = (likelihood.matrix.pop2.vial.phi.profile$log.likelihood), x = likelihood.matrix.pop2.vial.phi.profile$phi, type = "l", xaxt = "n")
axis(1, at = seq(0.5,5,by = 0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop2.vial.phi.profile.ordered = (likelihood.matrix.pop2.vial.phi.profile[order(likelihood.matrix.pop2.vial.phi.profile$log.likelihood),])[(likelihood.matrix.pop2.vial.phi.profile[order(likelihood.matrix.pop2.vial.phi.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop2.vial.phi.profile[order(likelihood.matrix.pop2.vial.phi.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop2.vial.phi.profile.ordered.range = range(likelihood.matrix.pop2.vial.phi.profile.ordered$phi)
abline(v = (subset(likelihood.matrix.pop2.vial.phi.profile.ordered, likelihood.matrix.pop2.vial.phi.profile.ordered$log.likelihood == max(likelihood.matrix.pop2.vial.phi.profile.ordered$log.likelihood)))$phi)
abline(v = c(likelihood.matrix.pop2.vial.phi.profile.ordered.range[1], likelihood.matrix.pop2.vial.phi.profile.ordered.range[2]), lty = 2)
# Vial Pop 3
plot(y = (likelihood.matrix.pop3.vial.phi.profile$log.likelihood), x = likelihood.matrix.pop3.vial.phi.profile$phi, type = "l", xaxt = "n")
axis(1, at = seq(0.5,5,by = 0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop3.vial.phi.profile.ordered = (likelihood.matrix.pop3.vial.phi.profile[order(likelihood.matrix.pop3.vial.phi.profile$log.likelihood),])[(likelihood.matrix.pop3.vial.phi.profile[order(likelihood.matrix.pop3.vial.phi.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop3.vial.phi.profile[order(likelihood.matrix.pop3.vial.phi.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop3.vial.phi.profile.ordered.range = range(likelihood.matrix.pop3.vial.phi.profile.ordered$phi)
abline(v = (subset(likelihood.matrix.pop3.vial.phi.profile.ordered, likelihood.matrix.pop3.vial.phi.profile.ordered$log.likelihood == max(likelihood.matrix.pop3.vial.phi.profile.ordered$log.likelihood)))$phi)
abline(v = c(likelihood.matrix.pop3.vial.phi.profile.ordered.range[1], likelihood.matrix.pop3.vial.phi.profile.ordered.range[2]), lty = 2)
# Vial Pop 4
plot(y = (likelihood.matrix.pop4.vial.phi.profile$log.likelihood), x = likelihood.matrix.pop4.vial.phi.profile$phi, type = "l", xaxt = "n")
axis(1, at = seq(0.5,5,by = 0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop4.vial.phi.profile.ordered = (likelihood.matrix.pop4.vial.phi.profile[order(likelihood.matrix.pop4.vial.phi.profile$log.likelihood),])[(likelihood.matrix.pop4.vial.phi.profile[order(likelihood.matrix.pop4.vial.phi.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop4.vial.phi.profile[order(likelihood.matrix.pop4.vial.phi.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop4.vial.phi.profile.ordered.range = range(likelihood.matrix.pop4.vial.phi.profile.ordered$phi)
abline(v = (subset(likelihood.matrix.pop4.vial.phi.profile.ordered, likelihood.matrix.pop4.vial.phi.profile.ordered$log.likelihood == max(likelihood.matrix.pop4.vial.phi.profile.ordered$log.likelihood)))$phi)
abline(v = c(likelihood.matrix.pop4.vial.phi.profile.ordered.range[1], likelihood.matrix.pop4.vial.phi.profile.ordered.range[2]), lty = 2)
# Vial Pop 5
plot(y = (likelihood.matrix.pop5.vial.phi.profile$log.likelihood), x = likelihood.matrix.pop5.vial.phi.profile$phi, type = "l", xaxt = "n")
axis(1, at = seq(0.5,5,by = 0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop5.vial.phi.profile.ordered = (likelihood.matrix.pop5.vial.phi.profile[order(likelihood.matrix.pop5.vial.phi.profile$log.likelihood),])[(likelihood.matrix.pop5.vial.phi.profile[order(likelihood.matrix.pop5.vial.phi.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop5.vial.phi.profile[order(likelihood.matrix.pop5.vial.phi.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop5.vial.phi.profile.ordered.range = range(likelihood.matrix.pop5.vial.phi.profile.ordered$phi)
abline(v = (subset(likelihood.matrix.pop5.vial.phi.profile.ordered, likelihood.matrix.pop5.vial.phi.profile.ordered$log.likelihood == max(likelihood.matrix.pop5.vial.phi.profile.ordered$log.likelihood)))$phi)
abline(v = c(likelihood.matrix.pop5.vial.phi.profile.ordered.range[1], likelihood.matrix.pop5.vial.phi.profile.ordered.range[2]), lty = 2)
#dev.off() # Run if you want to save the output as a pdf


## Plotting for Cages ### FIGURE S5
#pdf("/FILE/PATH/TO/FIGURES5.cages.likelihood.phi.pdf", width = 10, height = 8)
par(mfrow = c(2,3))
# cage Pop 1 
plot(y = (likelihood.matrix.pop1.cage.phi.profile$log.likelihood), x = likelihood.matrix.pop1.cage.phi.profile$phi, type = "l", xaxt = "n")
axis(1, at = seq(0.5,5,by = 0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop1.cage.phi.profile.ordered = (likelihood.matrix.pop1.cage.phi.profile[order(likelihood.matrix.pop1.cage.phi.profile$log.likelihood),])[(likelihood.matrix.pop1.cage.phi.profile[order(likelihood.matrix.pop1.cage.phi.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop1.cage.phi.profile[order(likelihood.matrix.pop1.cage.phi.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop1.cage.phi.profile.ordered.range = range(likelihood.matrix.pop1.cage.phi.profile.ordered$phi)
abline(v = (subset(likelihood.matrix.pop1.cage.phi.profile.ordered, likelihood.matrix.pop1.cage.phi.profile.ordered$log.likelihood == max(likelihood.matrix.pop1.cage.phi.profile.ordered$log.likelihood)))$phi)
abline(v = c(likelihood.matrix.pop1.cage.phi.profile.ordered.range[1], likelihood.matrix.pop1.cage.phi.profile.ordered.range[2]), lty = 2)
# cage Pop 2
plot(y = (likelihood.matrix.pop2.cage.phi.profile$log.likelihood), x = likelihood.matrix.pop2.cage.phi.profile$phi, type = "l", xaxt = "n")
axis(1, at = seq(0.5,5,by = 0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop2.cage.phi.profile.ordered = (likelihood.matrix.pop2.cage.phi.profile[order(likelihood.matrix.pop2.cage.phi.profile$log.likelihood),])[(likelihood.matrix.pop2.cage.phi.profile[order(likelihood.matrix.pop2.cage.phi.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop2.cage.phi.profile[order(likelihood.matrix.pop2.cage.phi.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop2.cage.phi.profile.ordered.range = range(likelihood.matrix.pop2.cage.phi.profile.ordered$phi)
abline(v = (subset(likelihood.matrix.pop2.cage.phi.profile.ordered, likelihood.matrix.pop2.cage.phi.profile.ordered$log.likelihood == max(likelihood.matrix.pop2.cage.phi.profile.ordered$log.likelihood)))$phi)
abline(v = c(likelihood.matrix.pop2.cage.phi.profile.ordered.range[1], likelihood.matrix.pop2.cage.phi.profile.ordered.range[2]), lty = 2)
# cage Pop 3
plot(y = (likelihood.matrix.pop3.cage.phi.profile$log.likelihood), x = likelihood.matrix.pop3.cage.phi.profile$phi, type = "l", xaxt = "n")
axis(1, at = seq(0.5,5,by = 0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop3.cage.phi.profile.ordered = (likelihood.matrix.pop3.cage.phi.profile[order(likelihood.matrix.pop3.cage.phi.profile$log.likelihood),])[(likelihood.matrix.pop3.cage.phi.profile[order(likelihood.matrix.pop3.cage.phi.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop3.cage.phi.profile[order(likelihood.matrix.pop3.cage.phi.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop3.cage.phi.profile.ordered.range = range(likelihood.matrix.pop3.cage.phi.profile.ordered$phi)
abline(v = (subset(likelihood.matrix.pop3.cage.phi.profile.ordered, likelihood.matrix.pop3.cage.phi.profile.ordered$log.likelihood == max(likelihood.matrix.pop3.cage.phi.profile.ordered$log.likelihood)))$phi)
abline(v = c(likelihood.matrix.pop3.cage.phi.profile.ordered.range[1], likelihood.matrix.pop3.cage.phi.profile.ordered.range[2]), lty = 2)
# cage Pop 4
plot(y = (likelihood.matrix.pop4.cage.phi.profile$log.likelihood), x = likelihood.matrix.pop4.cage.phi.profile$phi, type = "l", xaxt = "n")
axis(1, at = seq(0.5,5,by = 0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop4.cage.phi.profile.ordered = (likelihood.matrix.pop4.cage.phi.profile[order(likelihood.matrix.pop4.cage.phi.profile$log.likelihood),])[(likelihood.matrix.pop4.cage.phi.profile[order(likelihood.matrix.pop4.cage.phi.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop4.cage.phi.profile[order(likelihood.matrix.pop4.cage.phi.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop4.cage.phi.profile.ordered.range = range(likelihood.matrix.pop4.cage.phi.profile.ordered$phi)
abline(v = (subset(likelihood.matrix.pop4.cage.phi.profile.ordered, likelihood.matrix.pop4.cage.phi.profile.ordered$log.likelihood == max(likelihood.matrix.pop4.cage.phi.profile.ordered$log.likelihood)))$phi)
abline(v = c(likelihood.matrix.pop4.cage.phi.profile.ordered.range[1], likelihood.matrix.pop4.cage.phi.profile.ordered.range[2]), lty = 2)
# cage Pop 5
plot(y = (likelihood.matrix.pop5.cage.phi.profile$log.likelihood), x = likelihood.matrix.pop5.cage.phi.profile$phi, type = "l", xaxt = "n")
axis(1, at = seq(0.5,5,by = 0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop5.cage.phi.profile.ordered = (likelihood.matrix.pop5.cage.phi.profile[order(likelihood.matrix.pop5.cage.phi.profile$log.likelihood),])[(likelihood.matrix.pop5.cage.phi.profile[order(likelihood.matrix.pop5.cage.phi.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop5.cage.phi.profile[order(likelihood.matrix.pop5.cage.phi.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop5.cage.phi.profile.ordered.range = range(likelihood.matrix.pop5.cage.phi.profile.ordered$phi)
abline(v = (subset(likelihood.matrix.pop5.cage.phi.profile.ordered, likelihood.matrix.pop5.cage.phi.profile.ordered$log.likelihood == max(likelihood.matrix.pop5.cage.phi.profile.ordered$log.likelihood)))$phi)
abline(v = c(likelihood.matrix.pop5.cage.phi.profile.ordered.range[1], likelihood.matrix.pop5.cage.phi.profile.ordered.range[2]), lty = 2)
#dev.off()


## Plotting for Monogamy Treatment ###  FIGURE S4
#pdf("/FILE/PATH/TO/FIGURES4.monogamy.likelihood.phi.pdf", width = 10, height = 8)
par(mfrow = c(1,3))
# monogamy Pop 1 
plot(y = (likelihood.matrix.pop1.monogamy.phi.profile$log.likelihood), x = likelihood.matrix.pop1.monogamy.phi.profile$phi, type = "l", xaxt = "n")
axis(1, at = seq(0.5,5,by = 0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop1.monogamy.phi.profile.ordered = (likelihood.matrix.pop1.monogamy.phi.profile[order(likelihood.matrix.pop1.monogamy.phi.profile$log.likelihood),])[(likelihood.matrix.pop1.monogamy.phi.profile[order(likelihood.matrix.pop1.monogamy.phi.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop1.monogamy.phi.profile[order(likelihood.matrix.pop1.monogamy.phi.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop1.monogamy.phi.profile.ordered.range = range(likelihood.matrix.pop1.monogamy.phi.profile.ordered$phi)
abline(v = (subset(likelihood.matrix.pop1.monogamy.phi.profile.ordered, likelihood.matrix.pop1.monogamy.phi.profile.ordered$log.likelihood == max(likelihood.matrix.pop1.monogamy.phi.profile.ordered$log.likelihood)))$phi)
abline(v = c(likelihood.matrix.pop1.monogamy.phi.profile.ordered.range[1], likelihood.matrix.pop1.monogamy.phi.profile.ordered.range[2]), lty = 2)
# monogamy Pop 3
plot(y = (likelihood.matrix.pop3.monogamy.phi.profile$log.likelihood), x = likelihood.matrix.pop3.monogamy.phi.profile$phi, type = "l", xaxt = "n")
axis(1, at = seq(0.5,5,by = 0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop3.monogamy.phi.profile.ordered = (likelihood.matrix.pop3.monogamy.phi.profile[order(likelihood.matrix.pop3.monogamy.phi.profile$log.likelihood),])[(likelihood.matrix.pop3.monogamy.phi.profile[order(likelihood.matrix.pop3.monogamy.phi.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop3.monogamy.phi.profile[order(likelihood.matrix.pop3.monogamy.phi.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop3.monogamy.phi.profile.ordered.range = range(likelihood.matrix.pop3.monogamy.phi.profile.ordered$phi)
abline(v = (subset(likelihood.matrix.pop3.monogamy.phi.profile.ordered, likelihood.matrix.pop3.monogamy.phi.profile.ordered$log.likelihood == max(likelihood.matrix.pop3.monogamy.phi.profile.ordered$log.likelihood)))$phi)
abline(v = c(likelihood.matrix.pop3.monogamy.phi.profile.ordered.range[1], likelihood.matrix.pop3.monogamy.phi.profile.ordered.range[2]), lty = 2)
# monogamy Pop 5
plot(y = (likelihood.matrix.pop5.monogamy.phi.profile$log.likelihood), x = likelihood.matrix.pop5.monogamy.phi.profile$phi, type = "l", xaxt = "n")
axis(1, at = seq(0.5,5,by = 0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop5.monogamy.phi.profile.ordered = (likelihood.matrix.pop5.monogamy.phi.profile[order(likelihood.matrix.pop5.monogamy.phi.profile$log.likelihood),])[(likelihood.matrix.pop5.monogamy.phi.profile[order(likelihood.matrix.pop5.monogamy.phi.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop5.monogamy.phi.profile[order(likelihood.matrix.pop5.monogamy.phi.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop5.monogamy.phi.profile.ordered.range = range(likelihood.matrix.pop5.monogamy.phi.profile.ordered$phi)
abline(v = (subset(likelihood.matrix.pop5.monogamy.phi.profile.ordered, likelihood.matrix.pop5.monogamy.phi.profile.ordered$log.likelihood == max(likelihood.matrix.pop5.monogamy.phi.profile.ordered$log.likelihood)))$phi)
abline(v = c(likelihood.matrix.pop5.monogamy.phi.profile.ordered.range[1], likelihood.matrix.pop5.monogamy.phi.profile.ordered.range[2]), lty = 2)
#dev.off()


# Mating regimes Pops seperate but averaged over all population types   ### FIGURE 3
#pdf("/FILE/PATH/TO/FIGURE3.All.Mating.Regimes.likelihood.phi.pdf", width = 10, height = 8)
par(mfrow = c(1,3))
#par(mai = rep(2,4))
par(mar=c(1,1,1,1))
#par(oma = c(4,4,0.5,0.5))
par("mar")
# Vials -- All Pops
plot(y = (likelihood.matrix.vial.phi.profile$log.likelihood), x = likelihood.matrix.vial.phi.profile$phi, type = "l", xaxt = "n")
axis(1, at = seq(0.5,5,by = 0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.vial.phi.profile.ordered = (likelihood.matrix.vial.phi.profile[order(likelihood.matrix.vial.phi.profile$log.likelihood),])[(likelihood.matrix.vial.phi.profile[order(likelihood.matrix.vial.phi.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.vial.phi.profile[order(likelihood.matrix.vial.phi.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.vial.phi.profile.ordered.range = range(likelihood.matrix.vial.phi.profile.ordered$phi)
abline(v = (subset(likelihood.matrix.vial.phi.profile.ordered, likelihood.matrix.vial.phi.profile.ordered$log.likelihood == max(likelihood.matrix.vial.phi.profile.ordered$log.likelihood)))$phi)
abline(v = c(likelihood.matrix.vial.phi.profile.ordered.range[1], likelihood.matrix.vial.phi.profile.ordered.range[2]), lty = 2)
# Cages -- All Pops
plot(y = (likelihood.matrix.cage.phi.profile$log.likelihood), x = likelihood.matrix.cage.phi.profile$phi, type = "l", xaxt = "n")
axis(1, at = seq(0.5,5,by = 0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.cage.phi.profile.ordered = (likelihood.matrix.cage.phi.profile[order(likelihood.matrix.cage.phi.profile$log.likelihood),])[(likelihood.matrix.cage.phi.profile[order(likelihood.matrix.cage.phi.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.cage.phi.profile[order(likelihood.matrix.cage.phi.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.cage.phi.profile.ordered.range = range(likelihood.matrix.cage.phi.profile.ordered$phi)
abline(v = (subset(likelihood.matrix.cage.phi.profile.ordered, likelihood.matrix.cage.phi.profile.ordered$log.likelihood == max(likelihood.matrix.cage.phi.profile.ordered$log.likelihood)))$phi)
abline(v = c(likelihood.matrix.cage.phi.profile.ordered.range[1], likelihood.matrix.cage.phi.profile.ordered.range[2]), lty = 2)
# Monogamy -- All Pops
plot(y = (likelihood.matrix.monogamy.phi.profile$log.likelihood), x = likelihood.matrix.monogamy.phi.profile$phi, type = "l", xaxt = "n")
axis(1, at = seq(0.5,5,by = 0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.monogamy.phi.profile.ordered = (likelihood.matrix.monogamy.phi.profile[order(likelihood.matrix.monogamy.phi.profile$log.likelihood),])[(likelihood.matrix.monogamy.phi.profile[order(likelihood.matrix.monogamy.phi.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.monogamy.phi.profile[order(likelihood.matrix.monogamy.phi.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.monogamy.phi.profile.ordered.range = range(likelihood.matrix.monogamy.phi.profile.ordered$phi)
abline(v = (subset(likelihood.matrix.monogamy.phi.profile.ordered, likelihood.matrix.monogamy.phi.profile.ordered$log.likelihood == max(likelihood.matrix.monogamy.phi.profile.ordered$log.likelihood)))$phi)
abline(v = c(likelihood.matrix.monogamy.phi.profile.ordered.range[1], likelihood.matrix.monogamy.phi.profile.ordered.range[2]), lty = 2)
#dev.off()

# This figure is not presented in the paper 
# All Polygamy pops
#pdf("/FILE/PATH/TO/polygamy.monogamy.regime.phi.pdf", width = 10, height = 8)
par(mfrow = c(1,2))
#par(mai = rep(2,4))
#par(mar=c(1,1,1,1))
#par(oma = c(4,4,0.5,0.5))
#par("mar")
plot(y = (likelihood.matrix.monogamy.phi.profile$log.likelihood), x = likelihood.matrix.monogamy.phi.profile$phi, type = "l", xaxt = "n")
axis(1, at = seq(0.5,5,by = 0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.monogamy.phi.profile.ordered = (likelihood.matrix.monogamy.phi.profile[order(likelihood.matrix.monogamy.phi.profile$log.likelihood),])[(likelihood.matrix.monogamy.phi.profile[order(likelihood.matrix.monogamy.phi.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.monogamy.phi.profile[order(likelihood.matrix.monogamy.phi.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.monogamy.phi.profile.ordered.range = range(likelihood.matrix.monogamy.phi.profile.ordered$phi)
abline(v = (subset(likelihood.matrix.monogamy.phi.profile.ordered, likelihood.matrix.monogamy.phi.profile.ordered$log.likelihood == max(likelihood.matrix.monogamy.phi.profile.ordered$log.likelihood)))$phi)
abline(v = c(likelihood.matrix.monogamy.phi.profile.ordered.range[1], likelihood.matrix.monogamy.phi.profile.ordered.range[2]), lty = 2)
plot(y = (likelihood.matrix.polygamy.phi.profile$log.likelihood), x = likelihood.matrix.polygamy.phi.profile$phi, type = "l", xaxt = "n")
axis(1, at = seq(0.5,5,by = 0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.polygamy.phi.profile.ordered = (likelihood.matrix.polygamy.phi.profile[order(likelihood.matrix.polygamy.phi.profile$log.likelihood),])[(likelihood.matrix.polygamy.phi.profile[order(likelihood.matrix.polygamy.phi.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.polygamy.phi.profile[order(likelihood.matrix.polygamy.phi.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.polygamy.phi.profile.ordered.range = range(likelihood.matrix.polygamy.phi.profile.ordered$phi)
abline(v = (subset(likelihood.matrix.polygamy.phi.profile.ordered, likelihood.matrix.polygamy.phi.profile.ordered$log.likelihood == max(likelihood.matrix.polygamy.phi.profile.ordered$log.likelihood)))$phi)
abline(v = c(likelihood.matrix.polygamy.phi.profile.ordered.range[1], likelihood.matrix.polygamy.phi.profile.ordered.range[2]), lty = 2)
dev.off()




####    D. Calculating likelihood profiles of VarF for each mating regime by population type combination   ####
# Only plotted for a range of values
# Vials 
likelihood.matrix.pop1.vial.varF.profile = as.data.frame((apply(likelihood.matrix.pop1.vial, 1, max))[1:401])
likelihood.matrix.pop1.vial.varF.profile$varF = row.names(likelihood.matrix.pop1.vial.varF.profile)
colnames(likelihood.matrix.pop1.vial.varF.profile)[1] = "log.likelihood"
likelihood.matrix.pop2.vial.varF.profile = as.data.frame((apply(likelihood.matrix.pop2.vial, 1, max))[1:401])
likelihood.matrix.pop2.vial.varF.profile$varF = row.names(likelihood.matrix.pop2.vial.varF.profile)
colnames(likelihood.matrix.pop2.vial.varF.profile)[1] = "log.likelihood"
likelihood.matrix.pop3.vial.varF.profile = as.data.frame((apply(likelihood.matrix.pop3.vial, 1, max))[1:401])
likelihood.matrix.pop3.vial.varF.profile$varF = row.names(likelihood.matrix.pop3.vial.varF.profile)
colnames(likelihood.matrix.pop3.vial.varF.profile)[1] = "log.likelihood"
likelihood.matrix.pop4.vial.varF.profile = as.data.frame((apply(likelihood.matrix.pop4.vial, 1, max))[1:401])
likelihood.matrix.pop4.vial.varF.profile$varF = row.names(likelihood.matrix.pop4.vial.varF.profile)
colnames(likelihood.matrix.pop4.vial.varF.profile)[1] = "log.likelihood"
likelihood.matrix.pop5.vial.varF.profile = as.data.frame((apply(likelihood.matrix.pop5.vial, 1, max))[1:401])
likelihood.matrix.pop5.vial.varF.profile$varF = row.names(likelihood.matrix.pop5.vial.varF.profile)
colnames(likelihood.matrix.pop5.vial.varF.profile)[1] = "log.likelihood"

# Cages 
likelihood.matrix.pop1.cage.varF.profile = as.data.frame((apply(likelihood.matrix.pop1.cage, 1, max))[1:401])
likelihood.matrix.pop1.cage.varF.profile$varF = row.names(likelihood.matrix.pop1.cage.varF.profile)
colnames(likelihood.matrix.pop1.cage.varF.profile)[1] = "log.likelihood"
likelihood.matrix.pop2.cage.varF.profile = as.data.frame((apply(likelihood.matrix.pop2.cage, 1, max))[1:401])
likelihood.matrix.pop2.cage.varF.profile$varF = row.names(likelihood.matrix.pop2.cage.varF.profile)
colnames(likelihood.matrix.pop2.cage.varF.profile)[1] = "log.likelihood"
likelihood.matrix.pop3.cage.varF.profile = as.data.frame((apply(likelihood.matrix.pop3.cage, 1, max))[1:401])
likelihood.matrix.pop3.cage.varF.profile$varF = row.names(likelihood.matrix.pop3.cage.varF.profile)
colnames(likelihood.matrix.pop3.cage.varF.profile)[1] = "log.likelihood"
likelihood.matrix.pop4.cage.varF.profile = as.data.frame((apply(likelihood.matrix.pop4.cage, 1, max))[1:401])
likelihood.matrix.pop4.cage.varF.profile$varF = row.names(likelihood.matrix.pop4.cage.varF.profile)
colnames(likelihood.matrix.pop4.cage.varF.profile)[1] = "log.likelihood"
likelihood.matrix.pop5.cage.varF.profile = as.data.frame((apply(likelihood.matrix.pop5.cage, 1, max))[1:401])
likelihood.matrix.pop5.cage.varF.profile$varF = row.names(likelihood.matrix.pop5.cage.varF.profile)
colnames(likelihood.matrix.pop5.cage.varF.profile)[1] = "log.likelihood"

# Monogamy 
likelihood.matrix.pop1.monogamy.varF.profile = as.data.frame((apply(likelihood.matrix.pop1.monogamy, 1, max))[1:401])
likelihood.matrix.pop1.monogamy.varF.profile$varF = row.names(likelihood.matrix.pop1.monogamy.varF.profile)
colnames(likelihood.matrix.pop1.monogamy.varF.profile)[1] = "log.likelihood"
likelihood.matrix.pop3.monogamy.varF.profile = as.data.frame((apply(likelihood.matrix.pop3.monogamy, 1, max))[1:401])
likelihood.matrix.pop3.monogamy.varF.profile$varF = row.names(likelihood.matrix.pop3.monogamy.varF.profile)
colnames(likelihood.matrix.pop3.monogamy.varF.profile)[1] = "log.likelihood"
likelihood.matrix.pop5.monogamy.varF.profile = as.data.frame((apply(likelihood.matrix.pop5.monogamy, 1, max))[1:401])
likelihood.matrix.pop5.monogamy.varF.profile$varF = row.names(likelihood.matrix.pop5.monogamy.varF.profile)
colnames(likelihood.matrix.pop5.monogamy.varF.profile)[1] = "log.likelihood"

# Combination of population types in different  
likelihood.matrix.vial.varF.profile = as.data.frame((apply(likelihood.matrix.vial, 1, max))[1:401])
likelihood.matrix.vial.varF.profile$varF = row.names(likelihood.matrix.vial.varF.profile)
colnames(likelihood.matrix.vial.varF.profile)[1] = "log.likelihood"
likelihood.matrix.cage.varF.profile = as.data.frame((apply(likelihood.matrix.cage, 1, max))[1:401])
likelihood.matrix.cage.varF.profile$varF = row.names(likelihood.matrix.cage.varF.profile)
colnames(likelihood.matrix.cage.varF.profile)[1] = "log.likelihood"
likelihood.matrix.monogamy.varF.profile = as.data.frame((apply(likelihood.matrix.monogamy, 1, max))[1:401])
likelihood.matrix.monogamy.varF.profile$varF = row.names(likelihood.matrix.monogamy.varF.profile)
colnames(likelihood.matrix.monogamy.varF.profile)[1] = "log.likelihood"

likelihood.matrix.polygamy.varF.profile = as.data.frame((apply(likelihood.matrix.polygamy, 1, max))[1:401])
likelihood.matrix.polygamy.varF.profile$varF = row.names(likelihood.matrix.polygamy.varF.profile)
colnames(likelihood.matrix.polygamy.varF.profile)[1] = "log.likelihood"

####    E. Plotting Likelhiood Profiles for VarF    ####

## Vial environments for each population type ### Figure S9
#pdf("/FILE/PATH/TO/FigureS9.vials.likelihood.varF.pdf", width = 10, height = 8)
par(mfrow = c(2,3))
# Vial Pop 1 
plot(y = (likelihood.matrix.pop1.vial.varF.profile$log.likelihood), x = likelihood.matrix.pop1.vial.varF.profile$varF, type = "l", xaxt = "n")
axis(1, at = seq(2,10,by=0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop1.vial.varF.profile.ordered = (likelihood.matrix.pop1.vial.varF.profile[order(likelihood.matrix.pop1.vial.varF.profile$log.likelihood),])[(likelihood.matrix.pop1.vial.varF.profile[order(likelihood.matrix.pop1.vial.varF.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop1.vial.varF.profile[order(likelihood.matrix.pop1.vial.varF.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop1.vial.varF.profile.ordered.range = range(likelihood.matrix.pop1.vial.varF.profile.ordered$varF)
abline(v = (subset(likelihood.matrix.pop1.vial.varF.profile.ordered, likelihood.matrix.pop1.vial.varF.profile.ordered$log.likelihood == max(likelihood.matrix.pop1.vial.varF.profile.ordered$log.likelihood)))$varF)
abline(v = c(likelihood.matrix.pop1.vial.varF.profile.ordered.range[1], likelihood.matrix.pop1.vial.varF.profile.ordered.range[2]), lty = 2)
# Vial Pop 2
plot(y = (likelihood.matrix.pop2.vial.varF.profile$log.likelihood), x = likelihood.matrix.pop2.vial.varF.profile$varF, type = "l", xaxt = "n")
axis(1, at = seq(2,10,by=0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop2.vial.varF.profile.ordered = (likelihood.matrix.pop2.vial.varF.profile[order(likelihood.matrix.pop2.vial.varF.profile$log.likelihood),])[(likelihood.matrix.pop2.vial.varF.profile[order(likelihood.matrix.pop2.vial.varF.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop2.vial.varF.profile[order(likelihood.matrix.pop2.vial.varF.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop2.vial.varF.profile.ordered.range = range(likelihood.matrix.pop2.vial.varF.profile.ordered$varF)
abline(v = (subset(likelihood.matrix.pop2.vial.varF.profile.ordered, likelihood.matrix.pop2.vial.varF.profile.ordered$log.likelihood == max(likelihood.matrix.pop2.vial.varF.profile.ordered$log.likelihood)))$varF)
abline(v = c(likelihood.matrix.pop2.vial.varF.profile.ordered.range[1], likelihood.matrix.pop2.vial.varF.profile.ordered.range[2]), lty = 2)
# Vial Pop 3
plot(y = (likelihood.matrix.pop3.vial.varF.profile$log.likelihood), x = likelihood.matrix.pop3.vial.varF.profile$varF, type = "l", xaxt = "n")
axis(1, at = seq(2,10,by=0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop3.vial.varF.profile.ordered = (likelihood.matrix.pop3.vial.varF.profile[order(likelihood.matrix.pop3.vial.varF.profile$log.likelihood),])[(likelihood.matrix.pop3.vial.varF.profile[order(likelihood.matrix.pop3.vial.varF.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop3.vial.varF.profile[order(likelihood.matrix.pop3.vial.varF.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop3.vial.varF.profile.ordered.range = range(likelihood.matrix.pop3.vial.varF.profile.ordered$varF)
abline(v = (subset(likelihood.matrix.pop3.vial.varF.profile.ordered, likelihood.matrix.pop3.vial.varF.profile.ordered$log.likelihood == max(likelihood.matrix.pop3.vial.varF.profile.ordered$log.likelihood)))$varF)
abline(v = c(likelihood.matrix.pop3.vial.varF.profile.ordered.range[1], likelihood.matrix.pop3.vial.varF.profile.ordered.range[2]), lty = 2)
# Vial Pop 4
plot(y = (likelihood.matrix.pop4.vial.varF.profile$log.likelihood), x = likelihood.matrix.pop4.vial.varF.profile$varF, type = "l", xaxt = "n")
axis(1, at = seq(2,10,by=0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop4.vial.varF.profile.ordered = (likelihood.matrix.pop4.vial.varF.profile[order(likelihood.matrix.pop4.vial.varF.profile$log.likelihood),])[(likelihood.matrix.pop4.vial.varF.profile[order(likelihood.matrix.pop4.vial.varF.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop4.vial.varF.profile[order(likelihood.matrix.pop4.vial.varF.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop4.vial.varF.profile.ordered.range = range(likelihood.matrix.pop4.vial.varF.profile.ordered$varF)
abline(v = (subset(likelihood.matrix.pop4.vial.varF.profile.ordered, likelihood.matrix.pop4.vial.varF.profile.ordered$log.likelihood == max(likelihood.matrix.pop4.vial.varF.profile.ordered$log.likelihood)))$varF)
abline(v = c(likelihood.matrix.pop4.vial.varF.profile.ordered.range[1], likelihood.matrix.pop4.vial.varF.profile.ordered.range[2]), lty = 2)
# Vial Pop 5
plot(y = (likelihood.matrix.pop5.vial.varF.profile$log.likelihood), x = likelihood.matrix.pop5.vial.varF.profile$varF, type = "l", xaxt = "n")
axis(1, at = seq(2,10,by=0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop5.vial.varF.profile.ordered = (likelihood.matrix.pop5.vial.varF.profile[order(likelihood.matrix.pop5.vial.varF.profile$log.likelihood),])[(likelihood.matrix.pop5.vial.varF.profile[order(likelihood.matrix.pop5.vial.varF.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop5.vial.varF.profile[order(likelihood.matrix.pop5.vial.varF.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop5.vial.varF.profile.ordered.range = range(likelihood.matrix.pop5.vial.varF.profile.ordered$varF)
abline(v = (subset(likelihood.matrix.pop5.vial.varF.profile.ordered, likelihood.matrix.pop5.vial.varF.profile.ordered$log.likelihood == max(likelihood.matrix.pop5.vial.varF.profile.ordered$log.likelihood)))$varF)
abline(v = c(likelihood.matrix.pop5.vial.varF.profile.ordered.range[1], likelihood.matrix.pop5.vial.varF.profile.ordered.range[2]), lty = 2)
#dev.off()


## Cage environments for each population type   ### FIGURE S8
#pdf("/FILE/PATH/TO/FIGURES8.cages.likelihood.varF.pdf", width = 10, height = 8)
par(mfrow = c(2,3))
# cage Pop 1 
plot(y = (likelihood.matrix.pop1.cage.varF.profile$log.likelihood), x = likelihood.matrix.pop1.cage.varF.profile$varF, type = "l", xaxt = "n")
axis(1, at = seq(2,10,by=0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop1.cage.varF.profile.ordered = (likelihood.matrix.pop1.cage.varF.profile[order(likelihood.matrix.pop1.cage.varF.profile$log.likelihood),])[(likelihood.matrix.pop1.cage.varF.profile[order(likelihood.matrix.pop1.cage.varF.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop1.cage.varF.profile[order(likelihood.matrix.pop1.cage.varF.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop1.cage.varF.profile.ordered.range = range(likelihood.matrix.pop1.cage.varF.profile.ordered$varF)
abline(v = (subset(likelihood.matrix.pop1.cage.varF.profile.ordered, likelihood.matrix.pop1.cage.varF.profile.ordered$log.likelihood == max(likelihood.matrix.pop1.cage.varF.profile.ordered$log.likelihood)))$varF)
abline(v = c(likelihood.matrix.pop1.cage.varF.profile.ordered.range[1], likelihood.matrix.pop1.cage.varF.profile.ordered.range[2]), lty = 2)
# cage Pop 2
plot(y = (likelihood.matrix.pop2.cage.varF.profile$log.likelihood), x = likelihood.matrix.pop2.cage.varF.profile$varF, type = "l", xaxt = "n")
axis(1, at = seq(2,10,by=0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop2.cage.varF.profile.ordered = (likelihood.matrix.pop2.cage.varF.profile[order(likelihood.matrix.pop2.cage.varF.profile$log.likelihood),])[(likelihood.matrix.pop2.cage.varF.profile[order(likelihood.matrix.pop2.cage.varF.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop2.cage.varF.profile[order(likelihood.matrix.pop2.cage.varF.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop2.cage.varF.profile.ordered.range = range(likelihood.matrix.pop2.cage.varF.profile.ordered$varF)
abline(v = (subset(likelihood.matrix.pop2.cage.varF.profile.ordered, likelihood.matrix.pop2.cage.varF.profile.ordered$log.likelihood == max(likelihood.matrix.pop2.cage.varF.profile.ordered$log.likelihood)))$varF)
abline(v = c(likelihood.matrix.pop2.cage.varF.profile.ordered.range[1], likelihood.matrix.pop2.cage.varF.profile.ordered.range[2]), lty = 2)
# cage Pop 3
plot(y = (likelihood.matrix.pop3.cage.varF.profile$log.likelihood), x = likelihood.matrix.pop3.cage.varF.profile$varF, type = "l", xaxt = "n")
axis(1, at = seq(2,10,by=0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop3.cage.varF.profile.ordered = (likelihood.matrix.pop3.cage.varF.profile[order(likelihood.matrix.pop3.cage.varF.profile$log.likelihood),])[(likelihood.matrix.pop3.cage.varF.profile[order(likelihood.matrix.pop3.cage.varF.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop3.cage.varF.profile[order(likelihood.matrix.pop3.cage.varF.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop3.cage.varF.profile.ordered.range = range(likelihood.matrix.pop3.cage.varF.profile.ordered$varF)
abline(v = (subset(likelihood.matrix.pop3.cage.varF.profile.ordered, likelihood.matrix.pop3.cage.varF.profile.ordered$log.likelihood == max(likelihood.matrix.pop3.cage.varF.profile.ordered$log.likelihood)))$varF)
abline(v = c(likelihood.matrix.pop3.cage.varF.profile.ordered.range[1], likelihood.matrix.pop3.cage.varF.profile.ordered.range[2]), lty = 2)
# cage Pop 4
plot(y = (likelihood.matrix.pop4.cage.varF.profile$log.likelihood), x = likelihood.matrix.pop4.cage.varF.profile$varF, type = "l", xaxt = "n")
axis(1, at = seq(2,10,by=0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop4.cage.varF.profile.ordered = (likelihood.matrix.pop4.cage.varF.profile[order(likelihood.matrix.pop4.cage.varF.profile$log.likelihood),])[(likelihood.matrix.pop4.cage.varF.profile[order(likelihood.matrix.pop4.cage.varF.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop4.cage.varF.profile[order(likelihood.matrix.pop4.cage.varF.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop4.cage.varF.profile.ordered.range = range(likelihood.matrix.pop4.cage.varF.profile.ordered$varF)
abline(v = (subset(likelihood.matrix.pop4.cage.varF.profile.ordered, likelihood.matrix.pop4.cage.varF.profile.ordered$log.likelihood == max(likelihood.matrix.pop4.cage.varF.profile.ordered$log.likelihood)))$varF)
abline(v = c(likelihood.matrix.pop4.cage.varF.profile.ordered.range[1], likelihood.matrix.pop4.cage.varF.profile.ordered.range[2]), lty = 2)
# cage Pop 5
plot(y = (likelihood.matrix.pop5.cage.varF.profile$log.likelihood), x = likelihood.matrix.pop5.cage.varF.profile$varF, type = "l", xaxt = "n")
axis(1, at = seq(2,10,by=0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop5.cage.varF.profile.ordered = (likelihood.matrix.pop5.cage.varF.profile[order(likelihood.matrix.pop5.cage.varF.profile$log.likelihood),])[(likelihood.matrix.pop5.cage.varF.profile[order(likelihood.matrix.pop5.cage.varF.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop5.cage.varF.profile[order(likelihood.matrix.pop5.cage.varF.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop5.cage.varF.profile.ordered.range = range(likelihood.matrix.pop5.cage.varF.profile.ordered$varF)
abline(v = (subset(likelihood.matrix.pop5.cage.varF.profile.ordered, likelihood.matrix.pop5.cage.varF.profile.ordered$log.likelihood == max(likelihood.matrix.pop5.cage.varF.profile.ordered$log.likelihood)))$varF)
abline(v = c(likelihood.matrix.pop5.cage.varF.profile.ordered.range[1], likelihood.matrix.pop5.cage.varF.profile.ordered.range[2]), lty = 2)
dev.off()


## Monogamy environments for each population type   ### FIGURE S7
pdf("/FILE/PATH/TO/FIGURES7.monogamy.likelihood.varF.pdf", width = 10, height = 8)
par(mfrow = c(1,3))
# monogamy Pop 1 
plot(y = (likelihood.matrix.pop1.monogamy.varF.profile$log.likelihood), x = likelihood.matrix.pop1.monogamy.varF.profile$varF, type = "l", xaxt = "n")
axis(1, at = seq(2,10,by=0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop1.monogamy.varF.profile.ordered = (likelihood.matrix.pop1.monogamy.varF.profile[order(likelihood.matrix.pop1.monogamy.varF.profile$log.likelihood),])[(likelihood.matrix.pop1.monogamy.varF.profile[order(likelihood.matrix.pop1.monogamy.varF.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop1.monogamy.varF.profile[order(likelihood.matrix.pop1.monogamy.varF.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop1.monogamy.varF.profile.ordered.range = range(likelihood.matrix.pop1.monogamy.varF.profile.ordered$varF)
abline(v = (subset(likelihood.matrix.pop1.monogamy.varF.profile.ordered, likelihood.matrix.pop1.monogamy.varF.profile.ordered$log.likelihood == max(likelihood.matrix.pop1.monogamy.varF.profile.ordered$log.likelihood)))$varF)
abline(v = c(likelihood.matrix.pop1.monogamy.varF.profile.ordered.range[1], likelihood.matrix.pop1.monogamy.varF.profile.ordered.range[2]), lty = 2)
# monogamy Pop 3
plot(y = (likelihood.matrix.pop3.monogamy.varF.profile$log.likelihood), x = likelihood.matrix.pop3.monogamy.varF.profile$varF, type = "l", xaxt = "n")
axis(1, at = seq(2,10,by=0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop3.monogamy.varF.profile.ordered = (likelihood.matrix.pop3.monogamy.varF.profile[order(likelihood.matrix.pop3.monogamy.varF.profile$log.likelihood),])[(likelihood.matrix.pop3.monogamy.varF.profile[order(likelihood.matrix.pop3.monogamy.varF.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop3.monogamy.varF.profile[order(likelihood.matrix.pop3.monogamy.varF.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop3.monogamy.varF.profile.ordered.range = range(likelihood.matrix.pop3.monogamy.varF.profile.ordered$varF)
abline(v = (subset(likelihood.matrix.pop3.monogamy.varF.profile.ordered, likelihood.matrix.pop3.monogamy.varF.profile.ordered$log.likelihood == max(likelihood.matrix.pop3.monogamy.varF.profile.ordered$log.likelihood)))$varF)
abline(v = c(likelihood.matrix.pop3.monogamy.varF.profile.ordered.range[1], likelihood.matrix.pop3.monogamy.varF.profile.ordered.range[2]), lty = 2)
# monogamy Pop 5
plot(y = (likelihood.matrix.pop5.monogamy.varF.profile$log.likelihood), x = likelihood.matrix.pop5.monogamy.varF.profile$varF, type = "l", xaxt = "n")
axis(1, at = seq(2,10,by=0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.pop5.monogamy.varF.profile.ordered = (likelihood.matrix.pop5.monogamy.varF.profile[order(likelihood.matrix.pop5.monogamy.varF.profile$log.likelihood),])[(likelihood.matrix.pop5.monogamy.varF.profile[order(likelihood.matrix.pop5.monogamy.varF.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.pop5.monogamy.varF.profile[order(likelihood.matrix.pop5.monogamy.varF.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.pop5.monogamy.varF.profile.ordered.range = range(likelihood.matrix.pop5.monogamy.varF.profile.ordered$varF)
abline(v = (subset(likelihood.matrix.pop5.monogamy.varF.profile.ordered, likelihood.matrix.pop5.monogamy.varF.profile.ordered$log.likelihood == max(likelihood.matrix.pop5.monogamy.varF.profile.ordered$log.likelihood)))$varF)
abline(v = c(likelihood.matrix.pop5.monogamy.varF.profile.ordered.range[1], likelihood.matrix.pop5.monogamy.varF.profile.ordered.range[2]), lty = 2)
dev.off()


## Polygamy environments for each averaged over all population types ### Figure S10
#pdf("/FILE/PATH/TO/FigureS10.all.mating.regimes.likelihood.varF.pdf", width = 10, height = 8)
par(mfrow = c(1,3))
# monogamy -- All Pops
plot(y = (likelihood.matrix.monogamy.varF.profile$log.likelihood), x = likelihood.matrix.monogamy.varF.profile$varF, type = "l", xaxt = "n")
axis(1, at = seq(2,10,by=0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.monogamy.varF.profile.ordered = (likelihood.matrix.monogamy.varF.profile[order(likelihood.matrix.monogamy.varF.profile$log.likelihood),])[(likelihood.matrix.monogamy.varF.profile[order(likelihood.matrix.monogamy.varF.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.monogamy.varF.profile[order(likelihood.matrix.monogamy.varF.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.monogamy.varF.profile.ordered.range = range(likelihood.matrix.monogamy.varF.profile.ordered$varF)
abline(v = (subset(likelihood.matrix.monogamy.varF.profile.ordered, likelihood.matrix.monogamy.varF.profile.ordered$log.likelihood == max(likelihood.matrix.monogamy.varF.profile.ordered$log.likelihood)))$varF)
abline(v = c(likelihood.matrix.monogamy.varF.profile.ordered.range[1], likelihood.matrix.monogamy.varF.profile.ordered.range[2]), lty = 2)
# Cages -- All Pops
plot(y = (likelihood.matrix.cage.varF.profile$log.likelihood), x = likelihood.matrix.cage.varF.profile$varF, type = "l", xaxt = "n")
axis(1, at = seq(2,10,by=0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.cage.varF.profile.ordered = (likelihood.matrix.cage.varF.profile[order(likelihood.matrix.cage.varF.profile$log.likelihood),])[(likelihood.matrix.cage.varF.profile[order(likelihood.matrix.cage.varF.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.cage.varF.profile[order(likelihood.matrix.cage.varF.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.cage.varF.profile.ordered.range = range(likelihood.matrix.cage.varF.profile.ordered$varF)
abline(v = (subset(likelihood.matrix.cage.varF.profile.ordered, likelihood.matrix.cage.varF.profile.ordered$log.likelihood == max(likelihood.matrix.cage.varF.profile.ordered$log.likelihood)))$varF)
abline(v = c(likelihood.matrix.cage.varF.profile.ordered.range[1], likelihood.matrix.cage.varF.profile.ordered.range[2]), lty = 2)
# Vials -- All Pops
plot(y = (likelihood.matrix.vial.varF.profile$log.likelihood), x = likelihood.matrix.vial.varF.profile$varF, type = "l", xaxt = "n")
axis(1, at = seq(2,10,by=0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.vial.varF.profile.ordered = (likelihood.matrix.vial.varF.profile[order(likelihood.matrix.vial.varF.profile$log.likelihood),])[(likelihood.matrix.vial.varF.profile[order(likelihood.matrix.vial.varF.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.vial.varF.profile[order(likelihood.matrix.vial.varF.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.vial.varF.profile.ordered.range = range(likelihood.matrix.vial.varF.profile.ordered$varF)
abline(v = (subset(likelihood.matrix.vial.varF.profile.ordered, likelihood.matrix.vial.varF.profile.ordered$log.likelihood == max(likelihood.matrix.vial.varF.profile.ordered$log.likelihood)))$varF)
abline(v = c(likelihood.matrix.vial.varF.profile.ordered.range[1], likelihood.matrix.vial.varF.profile.ordered.range[2]), lty = 2)
#dev.off()


# All Polygamy pops
# This plot was not presented in the paper
#pdf("/FILE/PATH/TO/all.likelihood.varF.pdf", width = 10, height = 8)
plot(y = (likelihood.matrix.polygamy.varF.profile$log.likelihood), x = likelihood.matrix.polygamy.varF.profile$varF, type = "l", xaxt = "n")
axis(1, at = seq(2,10,by=0.5))
# Find the x-axis values for the confidence intervals 
# Order the df by likelihood and take values that are within 2 loglikelihood units of the max likelihood 
likelihood.matrix.polygamy.varF.profile.ordered = (likelihood.matrix.polygamy.varF.profile[order(likelihood.matrix.polygamy.varF.profile$log.likelihood),])[(likelihood.matrix.polygamy.varF.profile[order(likelihood.matrix.polygamy.varF.profile$log.likelihood),])$log.likelihood >= (max((likelihood.matrix.polygamy.varF.profile[order(likelihood.matrix.polygamy.varF.profile$log.likelihood),])$log.likelihood) - 2),]
likelihood.matrix.polygamy.varF.profile.ordered.range = range(likelihood.matrix.polygamy.varF.profile.ordered$varF)
abline(v = (subset(likelihood.matrix.polygamy.varF.profile.ordered, likelihood.matrix.polygamy.varF.profile.ordered$log.likelihood == max(likelihood.matrix.polygamy.varF.profile.ordered$log.likelihood)))$varF)
abline(v = c(likelihood.matrix.polygamy.varF.profile.ordered.range[1], likelihood.matrix.polygamy.varF.profile.ordered.range[2]), lty = 2)
#dev.off()



##





