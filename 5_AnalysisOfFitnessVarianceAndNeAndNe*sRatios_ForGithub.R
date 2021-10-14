#################################################################################################################################################################
#################################################################################################################################################################
###   Script Details:                                                                                                                                         ###
###   Analysis of the Variance in Fitness and Ne/Ne*s for Singh, A. and A.F., Agrawal. 2021. Sex-Specific Variance in Fitness and the Efficacy of Selection   ###
###   Script written by Amardeep Singh (amardeep.singh[at]alum.utoronto.ca)                                                                                   ###
#################################################################################################################################################################
#################################################################################################################################################################

# This script will organize and summarize the main variance in sex-specific fitness data and generate Figure 2 and Figure 3 of the main text
# The script below assumes that you have bootstrapped estimates of variance in fitness to obtain 95% CIs for plotting. 
#     Therefore, make sure to run the "" script and save the bootstrapped data  
# Below, I have organized this script into several parts:
#   Part 1  sets up the data, Parts 2 and 3 generates plots corresponding to Figure 2 and Figure 3
#   Part 4 and Part 5 generates the data in Table 2 and Table 3 respectively. 
#   Part 6 will replicate the statistical analyses for comparing means and variances in fitness across the sexes, mating regimes and population types

#   To generate figures there are 'code chunks' that can be run to generate each figure. You need to run through step 1 to set up the data first though
# NOTE: YOU WILL NEED TO SPECIFY THE PATH TO DATA FILES BELOW BY REPLACING "/FILE/PATH/TO/" WITH YOUR ACTUAL FILE PATHS 


###########################################
## 1. Setting up Packages and Data files ##       ## This section needs to be run for everything downstream and these dataframes must not be dumped!
###########################################
rm(list=ls())

##### Initial Data Summary ####

# Loading required packages 
require(doBy)
require(ggplot2)
require(RCurl)

# Reading in fitness data and cleaning it up a bit
fitness.data.github = getURL("https://raw.githubusercontent.com/asingh164/SexSpecificVarianceEfficacyOfSelection/master/FitnessDataNeprojectFinalSept2020.csv")
data=read.csv(text = fitness.data.github)
# Reading in data from bootstrap procedure 
bootstrapped.data.github = getURL("https://raw.githubusercontent.com/asingh164/SexSpecificVarianceEfficacyOfSelection/master/BootstrappedFitnessDataFinalSept2020.csv")
bootstrapped.data=read.csv(text=bootstrapped.data.github)

# Original Data locations 
#data = read.csv("/Users/amardeepsingh/Dropbox/Grad School Stuff/My Research/Environmental Heterogeneity and Mating Regieme Project/Data/Data Files for Dryad/FitnessDataNeprojectFinalSept2020.csv", header = TRUE, sep = ",")
#bootstrapped.data=read.csv("/Users/amardeepsingh/Dropbox/Grad School Stuff/My Research/Environmental Heterogeneity and Mating Regieme Project/Data/Data Files for Dryad/BootstrappedFitnessDataFinalJuly2021.csv")

##NOTE in these datasets population/population.type 1 and 2 refer to the low heterogenity treatments (i.e. 100% high/low condition individuals respectively)
## "high" and "low" refer to larval densities in the focal.individual column

data=as.data.frame(data)
data$n.wt=as.numeric(as.character(data$n.wt))
data=data[data$n.total == "32",] # Remove any observation where we have less than 32 offspring counted ## This caused a drop of 112 replicates across the whole experiment (~5.4% of all replicates)
data=data[!(is.na(data$n.wt)),] # Remove any rows with NA in fitness column 

data.summary=summaryBy(n.wt ~ treatment + population.type + sex, data=data, FUN = c(mean, var))
# Apply the variance adjustment and calculate the adjusted variance (see Methods and Results for details)
#data.summary$adjusted.variance = (2/data.summary$n.wt.mean)* (1 - (2/data.summary$n.wt.mean))*data.summary$n.wt.mean + (2/data.summary$n.wt.mean)^2 * (data.summary$n.wt.var)
data.summary$adjusted.variance = (2 * data.summary$n.wt.var) / data.summary$n.wt.mean  
# For convenience below, convert sex from m and f to male and female
data.summary$sex = gsub('m', 'male', data.summary$sex)
data.summary$sex = gsub('f', 'female', data.summary$sex)

#############################################################
## 2. Figure 2: Adjusted Variance in Sex-Specific Fitness  ##
#############################################################

# Assign 95% CIs for each estimated variance and adjusted variance Values
# Run this loop to get a colum of the lower and upper 95% CI for each variance and adjusted variance value
treatment.levels = levels(as.factor(bootstrapped.data$treatment))
sex.levels = levels(as.factor(bootstrapped.data$sex))
population.type.levels = levels(as.factor((bootstrapped.data$population.type)))
quantile.df = as.data.frame(matrix(NA, nrow = 25, ncol = 9))
names(quantile.df) = c("treatment", "sex", "population", "mean.lower.ci", "mean.upper.ci", "unadj.lower.ci", "unadj.upper.ci", "adj.lower.ci", "adj.upper.ci")
row = 0
for (i in treatment.levels){
  for (j in sex.levels){
    for (k in population.type.levels){
      bootstrapped.data.subset = bootstrapped.data[bootstrapped.data$treatment == i & bootstrapped.data$sex == j & bootstrapped.data$population.type == k ,]
      row = row + 1
      quantile.df[row,"treatment"] = i
      quantile.df[row,"sex"] = j
      quantile.df[row,"population"] = as.integer(k)
      
      
      #calculating upper and lower CIs for mean, var and adj. var
      quantile.df[row,"mean.lower.ci"] = as.numeric(quantile(bootstrapped.data.subset$mean.wt,c(0.025)))
      quantile.df[row,"mean.upper.ci"] = as.numeric(quantile(bootstrapped.data.subset$mean.wt,c(0.975)))
      
      quantile.df[row,"unadj.lower.ci"] = as.numeric(quantile(bootstrapped.data.subset$var.wt,c(0.025)))
      quantile.df[row,"unadj.upper.ci"] = as.numeric(quantile(bootstrapped.data.subset$var.wt,c(0.975)))
      
      quantile.df[row,"adj.lower.ci"] = as.numeric(quantile(bootstrapped.data.subset$adjusted.variance,c(0.025)))
      quantile.df[row,"adj.upper.ci"] = as.numeric(quantile(bootstrapped.data.subset$adjusted.variance,c(0.975)))
      
    }
  }
}
# Combine fitness data summary with CI estimates 
data.summary.figure2 = merge(data.summary, quantile.df, by.x = c("sex", "treatment", "population.type"), by.y = c("sex", "treatment","population"), sort = FALSE)

# Plot of adjusted variances 
# Change labels of data for each population type 
data.summary.figure2$population.type = data.summary$population.type
data.summary.figure2$population.type = gsub(as.numeric(1), "A", data.summary.figure2$population.type)
data.summary.figure2$population.type = gsub(as.numeric(2), "B", data.summary.figure2$population.type)
data.summary.figure2$population.type = gsub(as.numeric(3), "C", data.summary.figure2$population.type)
data.summary.figure2$population.type = gsub(as.numeric(4), "D", data.summary.figure2$population.type)
data.summary.figure2$population.type = gsub(as.numeric(5), "E", data.summary.figure2$population.type)

data.summary.figure2$population.type.label = paste("Population Type ", data.summary.figure2$population.type, sep = "")

# Plotting Adjusted Variances ## Figure 2
var.plot=ggplot(data.summary.figure2, aes(y = adjusted.variance))+
  geom_point(aes(x=sex, shape=treatment, colour=sex), size = 6.5,  position=position_dodge(width=0.5))+ylim(0,12) +
  geom_errorbar(aes(ymin=adj.lower.ci, ymax=adj.upper.ci,x=sex, shape=treatment, colour=sex),size = 1.5, width=0.0,position=position_dodge(width=0.5)) + 
  theme_bw() + 
  theme(panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank(),panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank(),
        text = element_text(size=40, family = "Times"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
        axis.title.y = element_blank(), axis.text.y=element_text(colour = "black"), axis.ticks.y=element_line(colour = "black"),
        strip.background = element_rect(colour = "white", fill = "white"), strip.text = element_text(size = 40), legend.position = "none") +
  facet_wrap(~population.type, nrow = 1, strip.position="bottom", scale = "free_x") + 
  scale_fill_manual(values=c("#F4A582", "#92C5DE")) 
var.plot


###################################################################
## 3. Figure 3: Plotting the ratio r of Ne_mf*s_mf / Ne_ff*s_ff  ##
###################################################################
# Plot of the the ratio of Ne*s in a population of males and females (i.e., Ne_mf*s_mf) relative to a population where males are "female-like" with 
#   respect to the variance in fitness and strength of selection (i.e., Ne_ff*s_ff)

# Loading summary data from Sharp and Agrawal 2013

# Original location of data 
sharp.agrawal.data = read.csv("/Users/amardeepsingh/Dropbox/Grad School Stuff/My Research/Environmental Heterogeneity and Mating Regieme Project/Data/DataFromNathaneil_Sharp_Agrawal_2013.csv")

# Bootstrapping estimates of s_male / s_female from Sharp and Agrawal 2013 data 
selection.ratio.bootstrapped = vector(mode = "numeric", length = 10000)
for (i in 1:10000){
  s.male.sample = sample(sharp.agrawal.data$s.male, replace = TRUE)
  s.female.sample = sample(sharp.agrawal.data$s.female, replace = TRUE)
  selection.ratio.bootstrapped[i] = mean(s.male.sample) / mean(s.female.sample)
}

# Setting up dataframe for plotting 
Ne.df = as.data.frame(matrix(NA, nrow = 13, ncol = 6))
names(Ne.df) = c("treatment", "population","Ne.mf", "Ne.ff", "ratio.of.Ne.mf.to.Ne.ff", "ratio.of.Nes")

treatment.levels = levels(as.factor(data.summary$treatment))
population.type.levels = levels(as.factor((data.summary$population.type)))
# Calculating the ratio r for each population type
row = 0
for (i in treatment.levels){
  for (k in population.type.levels){
    data.summary.subet.male = data.summary[data.summary$treatment == i & data.summary$sex == "male" & data.summary$population.type == k,]
    data.summary.subet.female = data.summary[data.summary$treatment == i & data.summary$sex == "female" & data.summary$population.type == k,]
    
    if (nrow(data.summary.subet.male) > 0 & nrow(data.summary.subet.female) > 0) {
      row = row + 1
      Ne.df[row,"treatment"] = i
      Ne.df[row,"population"] = as.integer(k)
      
      #calculating upper and lower CIs for mean, var and adj. var
      Ne.df[row,"Ne.mf"] = as.numeric((8 *32) / (4 +  data.summary.subet.male$adjusted.variance + data.summary.subet.female$adjusted.variance))
      Ne.df[row,"Ne.ff"] = as.numeric((8 *32) / (4 +  data.summary.subet.female$adjusted.variance + data.summary.subet.female$adjusted.variance))
    }
  }
}
Ne.df$ratio.of.Ne.mf.to.Ne.ff= Ne.df$Ne.mf / Ne.df$Ne.ff 

### Estimate the average strength of selection in males and females in two different ways
# Method 1: Using the average value from Sharp and Agrawal 2013 and Mallet et al 2011.
#Ne.df$ratio.of.Nes.mean = Ne.df$Ne.mf / Ne.df$Ne.ff * 1.25   

# Method 2: Adding error in S via estimates from Sharp and Agrawal 2013
Ne.df$ratio.of.Nes.mean = Ne.df$Ne.mf / Ne.df$Ne.ff * (0.5 * (1 + (mean(sharp.agrawal.data$s.male) / mean(sharp.agrawal.data$s.female)))) 


names(Ne.df)[1:2] = c("treatment", "population")

### Bootstrapping Ne*s ratio values to generate 95% CIs
treatment = c("vial", "cage", "monogamy")
population = c(1,2,3,4,5)

Ne.df$Ne.mf.upper = NA
Ne.df$Ne.mf.lower = NA
Ne.df$Ne.ff.upper = NA
Ne.df$Ne.ff.lower = NA

Ne.df$Ne.s.ratio.upper = NA
Ne.df$Ne.s.ratio.lower = NA
#Ne.df$Ne.s.ff.upper = NA
#Ne.df$Ne.s.ff.lower = NA

for (i in treatment){
  for (j in population){
    tmp.female = subset(bootstrapped.data, bootstrapped.data$treatment == i & bootstrapped.data$sex == "female" & bootstrapped.data$population.type == j)$adjusted.variance
    tmp.male = subset(bootstrapped.data, bootstrapped.data$treatment == i & bootstrapped.data$sex == "male" & bootstrapped.data$population.type == j)$adjusted.variance
    
    Ne.mf.tmp=c(NA, length = 10000)
    Ne.ff.tmp=c(NA, length = 10000)
    #Ne.mf.ratio = c(NA, length = 10000)
    Nes.ratio = c(NA, length = 10000)
    
    for (k in 1:10000){
      Ne.mf.tmp[k]=(8*32)/(tmp.male[k]+tmp.female[k]+4)
      Ne.ff.tmp[k]=(8*32)/(tmp.female[k]+tmp.female[k]+4)
     
      ### Adding error (or not) on the average strength of selection in males and females in two different ways
      # Method 1: Using the average value from Sharp and Agrawal 2013 and Mallet et al 2011.
      #Nes.ratio[k] = ( ((8*32)/(tmp.male[k]+tmp.female[k]+4)) / ((8*32)/(tmp.female[k]+tmp.female[k]+4)) ) * 1.25 
      
      # Method 2: Adding error in S via estimates from Sharp and Agrawal 2013
      Nes.ratio[k] = ( ((8*32)/(tmp.male[k]+tmp.female[k]+4)) / ((8*32)/(tmp.female[k]+tmp.female[k]+4)) ) * (0.5*(1+selection.ratio.bootstrapped[k]))
    }
    
    if (i == "monogamy" & j == 2 | i == "monogamy" & j == 4) {
      
    } else {
      Ne.df[Ne.df$treatment == i & Ne.df$population == j,]$Ne.mf.upper = quantile(Ne.mf.tmp, 0.975)
      Ne.df[Ne.df$treatment == i & Ne.df$population == j,]$Ne.mf.lower = quantile(Ne.mf.tmp, 0.025)
      Ne.df[Ne.df$treatment == i & Ne.df$population == j,]$Ne.ff.upper = quantile(Ne.ff.tmp, 0.975)
      Ne.df[Ne.df$treatment == i & Ne.df$population == j,]$Ne.ff.lower = quantile(Ne.ff.tmp, 0.025)
      
      Ne.df[Ne.df$treatment == i & Ne.df$population == j,]$Ne.s.ratio.upper = quantile(Nes.ratio, 0.975)
      Ne.df[Ne.df$treatment == i & Ne.df$population == j,]$Ne.s.ratio.lower = quantile(Nes.ratio, 0.025)
      #Ne.df[Ne.df$treatment == i & Ne.df$population == j,]$Ne.s.ff.upper = quantile(Nes.ratio, 0.975)
      #Ne.df[Ne.df$treatment == i & Ne.df$population == j,]$Ne.s.ff.lower = quantile(Nes.ratio, 0.025)
    }
  }
}

# Change labels of data for each population type for plotting convenience 
Ne.df$population.type = Ne.df$population
Ne.df$population.type = gsub(as.numeric(1), "A", Ne.df$population.type)
Ne.df$population.type = gsub(as.numeric(2), "B", Ne.df$population.type)
Ne.df$population.type = gsub(as.numeric(3), "C", Ne.df$population.type)
Ne.df$population.type = gsub(as.numeric(4), "D", Ne.df$population.type)
Ne.df$population.type = gsub(as.numeric(5), "E", Ne.df$population.type)
Ne.df$ratio.of.Nes = as.numeric(Ne.df$ratio.of.Nes.mean)
# shapes to be consistent with figure 2:
# Monogamy = 17; Complex environments = 16; Simple Environments = 15
Nes.plot=ggplot(Ne.df, aes(y = as.numeric(ratio.of.Nes))) +
  geom_point(aes(x=treatment, shape = as.factor(treatment)), size = 7.5,  position=position_dodge(width=0.5)) + scale_shape_manual(values=c(16, 17, 15)) +
  geom_point(aes(size=2, x=treatment),  position=position_dodge(0.5))+
  geom_errorbar(aes(ymin = Ne.s.ratio.lower, ymax = Ne.s.ratio.upper,x = as.factor(treatment)),size = 2, width=0.0,position=position_dodge(0.5)) + 
  theme_bw() +
  ylim(0.25,2.5) +
  theme(panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank(), panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank(),
        text = element_text(size=40, family = "Times"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
        axis.title.y = element_blank(), axis.text.y=element_text(colour = "black"), axis.ticks.y=element_line(colour = "black"),
        strip.background = element_rect(colour = "black", fill = "light grey"), strip.text = element_text(size = 40), legend.position = "none") +
  geom_hline(yintercept=1) +
  scale_fill_manual(values = "#FFA500") +
  facet_wrap(~as.factor(population.type), nrow=1, strip.position = "bottom", scale = "free_x") 
Nes.plot


######################################################
## 4.  Table 2: Summary Table for Adjusted variance ##
######################################################
# NOTE: Point estimates reported in Table 2 of the main text are the ACTUAL estimates of variances. Estimates reported in the summary table below are those means from the bootstrapped data.
#       Confidence intervals of estimates reported in Table 2 correspond to the 95% quantiles in the summary table below.

## Function to calculate the adjusted variance in fitness ## NOTE: We didn't end up using this, but and instead found that this adjustment could be used even when the 
#   mean was > 2 (see Supplemental Materials of the paper). Left code here just for legacy to remind me of what we did and why and for posterity.
# If the mean of the resample is < or = 2 then use the obsered variance 
#adjusted.variance.function  = function(mean, variance){
#  ifelse (mean > 2, ((2 / mean)*(1 - (2/mean))*(mean) + ((2 / mean)^2 * (variance))), variance)
#}
## Replace female monogamy with var because the mean was greater than 2
#for (i in 1:nrow(bootstrapped.data)){
#  if (bootstrapped.data[i,"treatment"] == "monogamy" & bootstrapped.data[i,"sex"] == "female") {
#    bootstrapped.data[i,9] = bootstrapped.data[i,7] 
#  } else {
#  }
#} 

# This is the actual function we originally applied but was found later to be inaccurate. I leave it here for posterity.
#adjusted.variance.function  = function(mean, variance){
#  ((2 / mean)*(1 - (2/mean))*(mean)) + ((2 / mean)^2 * (variance))
#}

# Final function -- This is what is reported in the paper
adjusted.variance.function  = function(mean, variance){
  (2*variance)/mean
}

quantile.function.lower = function(x){
  as.numeric(quantile(x, c(0.025))[1])
}
quantile.function.upper = function(x){
  as.numeric(quantile(x, c(0.975))[1])
}


## Summary table for adjusted variance for each sex * treatment combination ## This summary correspoinds to Table 2 in the main text (except the point estimates are not identical! - See above Note)
adj.variance.summary.table = summaryBy(adjusted.variance ~ treatment+population.type+sex, FUN = c(mean, quantile.function.lower,quantile.function.upper), data = bootstrapped.data)


##########################################################
## 5.  Table 3:  Summary table for Ne and Ratios of Ne  ##
##########################################################
# NOTE: Point estimates reported in Table 3 of the main text are the ACTUAL estimates of variances. Estimates reported in the summary table below are those means from the bootstrapped data.
#       Confidence intervals of estimates reported in Table 3 correspond to the 95% quantiles in the summary table below.
# NOTE: There are two summary dataframes generated below, one that corresponds to to columns 1 and 2 and Table 3 and one that corresponds to column 3 of Table 3 - This info is denoted next to commands
# Merge dataframe so that male and female variance for each population type and mating regime are on a single line 
male.tmp = bootstrapped.data[bootstrapped.data$sex == "male", ]
female.tmp = bootstrapped.data[bootstrapped.data$sex == "female",]

Ne.summary.table.tmp = merge(male.tmp, female.tmp, by.x=c("treatment", "population.type", "replicate"), by.y = c("treatment", "population.type", "replicate"), sort = FALSE)
Ne.summary.table.tmp = Ne.summary.table.tmp[,c(1:3,8,13)]
colnames(Ne.summary.table.tmp) = c("treatment", "population.type", "replicate", "male.adj.var", "female.adj.var")

Ne.summary.table.tmp$Ne_mf_N= 8 / (Ne.summary.table.tmp$male.adj.var + Ne.summary.table.tmp$female.adj.var + 4) 
Ne.summary.table.tmp$Ne_ff_N = 8 / (Ne.summary.table.tmp$female.adj.var + Ne.summary.table.tmp$female.adj.var + 4) 
Ne.summary.table.tmp$Ne_mf_Ne_ff = Ne.summary.table.tmp$Ne_mf / Ne.summary.table.tmp$Ne_ff

# Summary of Ne Data ## Columns 1, 2 and 3 in Table 3 in the main text
Ne_mf.summary.table = summaryBy(Ne_mf_N ~ treatment + population.type, FUN = c(mean, quantile.function.lower, quantile.function.upper), data = Ne.summary.table.tmp)
Ne_ff.summary.table = summaryBy(Ne_ff_N ~ treatment + population.type, FUN = c(mean, quantile.function.lower, quantile.function.upper), data = Ne.summary.table.tmp)
Ne_mf_Ne_ff.summary.table = summaryBy(Ne_mf_Ne_ff ~ treatment + population.type, FUN = c(mean, quantile.function.lower, quantile.function.upper), data = Ne.summary.table.tmp)

# Summary for ratio of Ne_mf / Ne_ff ## Column 4 in Table 3 of the main text
Ne.summary.table.tmp$Ne_mf= (8*32) / (Ne.summary.table.tmp$male.adj.var + Ne.summary.table.tmp$female.adj.var + 4) 
Ne.summary.table.tmp$Ne_ff = (8*32) / (Ne.summary.table.tmp$female.adj.var + Ne.summary.table.tmp$female.adj.var + 4) 
Ne.ratio.summary = Ne.summary.table.tmp[,c(1:3,9:10)]
Ne.ratio.summary$ratio = Ne.summary.table.tmp$Ne_mf / Ne.summary.table.tmp$Ne_ff
Ne.ratio.summary$Nes = (Ne.summary.table.tmp$Ne_mf / Ne.summary.table.tmp$Ne_ff) * 1.25

Ne.ratio.summary.table = summaryBy(ratio ~ treatment + population.type, FUN = c(mean, quantile.function.lower, quantile.function.upper), data = Ne.ratio.summary)
## "Ne.Nes.summary.table" corresponds to Column 3 in Table 3 of the main text
Ne.Nes.summary.table = summaryBy(Nes ~ treatment + population.type, FUN = c(mean, quantile.function.lower, quantile.function.upper), data = Ne.ratio.summary)

Nes.vial=c()
Nes.cage=c()
Nes.monogamy=c()

# Comparing averages across mating regimes 

vial.pop1 = Ne.ratio.summary$Nes[Ne.ratio.summary$treatment == "vial" & Ne.ratio.summary$population.type==1]
vial.pop2 = Ne.ratio.summary$Nes[Ne.ratio.summary$treatment == "vial" & Ne.ratio.summary$population.type==2]
vial.pop3 = Ne.ratio.summary$Nes[Ne.ratio.summary$treatment == "vial" & Ne.ratio.summary$population.type==3]
vial.pop4 = Ne.ratio.summary$Nes[Ne.ratio.summary$treatment == "vial" & Ne.ratio.summary$population.type==4]
vial.pop5 = Ne.ratio.summary$Nes[Ne.ratio.summary$treatment == "vial" & Ne.ratio.summary$population.type==5]

cage.pop1 = Ne.ratio.summary$Nes[Ne.ratio.summary$treatment == "cage" & Ne.ratio.summary$population.type==1]
cage.pop2 = Ne.ratio.summary$Nes[Ne.ratio.summary$treatment == "cage" & Ne.ratio.summary$population.type==2]
cage.pop3 = Ne.ratio.summary$Nes[Ne.ratio.summary$treatment == "cage" & Ne.ratio.summary$population.type==3]
cage.pop4 = Ne.ratio.summary$Nes[Ne.ratio.summary$treatment == "cage" & Ne.ratio.summary$population.type==4]
cage.pop5 = Ne.ratio.summary$Nes[Ne.ratio.summary$treatment == "cage" & Ne.ratio.summary$population.type==5]

monogamy.pop1 = Ne.ratio.summary$Nes[Ne.ratio.summary$treatment == "monogamy" & Ne.ratio.summary$population.type==1]
monogamy.pop3 = Ne.ratio.summary$Nes[Ne.ratio.summary$treatment == "monogamy" & Ne.ratio.summary$population.type==3]
monogamy.pop5 = Ne.ratio.summary$Nes[Ne.ratio.summary$treatment == "monogamy" & Ne.ratio.summary$population.type==5]

for (i in 1:10000){
  Nes.vial[i] = mean(c(vial.pop1[i],vial.pop2[i],vial.pop3[i],vial.pop4[i],vial.pop5[i]))
  Nes.cage[i] = mean(c(cage.pop1[i],cage.pop2[i],cage.pop3[i],cage.pop4[i],cage.pop5[i]))
  Nes.monogamy[i] = mean(c(monogamy.pop1[i],monogamy.pop3[i],monogamy.pop5[i]))
}

#

##########################################################################
##  6. Statistical Analysis to Compare Mean and Variance in Fitness     ##
##########################################################################
# The script below will replicate the results presented in the paper for comparisions of means and variances in fitness in the Results text.
# Statistical significance was assessed on the basis of non-overlapping 95% CIs

################ Obtaining true mean differences ################ 

# Subset out data frame into individual sex*pop*mating regime combinations 
# Vials 
male.data.variance.vial.type1 = subset(data.summary[data.summary$sex == "male" & data.summary$treatment == "vial" & data.summary$population.type==1, ])$adjusted.variance
male.data.variance.vial.type2 = subset(data.summary[data.summary$sex == "male" & data.summary$treatment == "vial" & data.summary$population.type==2, ])$adjusted.variance
male.data.variance.vial.type3 = subset(data.summary[data.summary$sex == "male" & data.summary$treatment == "vial" & data.summary$population.type==3, ])$adjusted.variance
male.data.variance.vial.type4 = subset(data.summary[data.summary$sex == "male" & data.summary$treatment == "vial" & data.summary$population.type==4, ])$adjusted.variance
male.data.variance.vial.type5 = subset(data.summary[data.summary$sex == "male" & data.summary$treatment == "vial" & data.summary$population.type==5, ])$adjusted.variance
female.data.variance.vial.type1 = subset(data.summary[data.summary$sex == "female" & data.summary$treatment == "vial" & data.summary$population.type==1, ])$adjusted.variance
female.data.variance.vial.type2 = subset(data.summary[data.summary$sex == "female" & data.summary$treatment == "vial" & data.summary$population.type==2, ])$adjusted.variance
female.data.variance.vial.type3 = subset(data.summary[data.summary$sex == "female" & data.summary$treatment == "vial" & data.summary$population.type==3, ])$adjusted.variance
female.data.variance.vial.type4 = subset(data.summary[data.summary$sex == "female" & data.summary$treatment == "vial" & data.summary$population.type==4, ])$adjusted.variance
female.data.variance.vial.type5 = subset(data.summary[data.summary$sex == "female" & data.summary$treatment == "vial" & data.summary$population.type==5, ])$adjusted.variance
#Cages 
male.data.variance.cage.type1 = subset(data.summary[data.summary$sex == "male" & data.summary$treatment == "cage" & data.summary$population.type==1, ])$adjusted.variance
male.data.variance.cage.type2 = subset(data.summary[data.summary$sex == "male" & data.summary$treatment == "cage" & data.summary$population.type==2, ])$adjusted.variance
male.data.variance.cage.type3 = subset(data.summary[data.summary$sex == "male" & data.summary$treatment == "cage" & data.summary$population.type==3, ])$adjusted.variance
male.data.variance.cage.type4 = subset(data.summary[data.summary$sex == "male" & data.summary$treatment == "cage" & data.summary$population.type==4, ])$adjusted.variance
male.data.variance.cage.type5 = subset(data.summary[data.summary$sex == "male" & data.summary$treatment == "cage" & data.summary$population.type==5, ])$adjusted.variance
female.data.variance.cage.type1 = subset(data.summary[data.summary$sex == "female" & data.summary$treatment == "cage" & data.summary$population.type==1, ])$adjusted.variance
female.data.variance.cage.type2 = subset(data.summary[data.summary$sex == "female" & data.summary$treatment == "cage" & data.summary$population.type==2, ])$adjusted.variance
female.data.variance.cage.type3 = subset(data.summary[data.summary$sex == "female" & data.summary$treatment == "cage" & data.summary$population.type==3, ])$adjusted.variance
female.data.variance.cage.type4 = subset(data.summary[data.summary$sex == "female" & data.summary$treatment == "cage" & data.summary$population.type==4, ])$adjusted.variance
female.data.variance.cage.type5 = subset(data.summary[data.summary$sex == "female" & data.summary$treatment == "cage" & data.summary$population.type==5, ])$adjusted.variance
#Monogamy
male.data.variance.monogamy.type1 = subset(data.summary[data.summary$sex == "male" & data.summary$treatment == "monogamy" & data.summary$population.type==1, ])$adjusted.variance
male.data.variance.monogamy.type3 = subset(data.summary[data.summary$sex == "male" & data.summary$treatment == "monogamy" & data.summary$population.type==3, ])$adjusted.variance
male.data.variance.monogamy.type5 = subset(data.summary[data.summary$sex == "male" & data.summary$treatment == "monogamy" & data.summary$population.type==5, ])$adjusted.variance
female.data.variance.monogamy.type1 = subset(data.summary[data.summary$sex == "female" & data.summary$treatment == "monogamy" & data.summary$population.type==1, ])$adjusted.variance
female.data.variance.monogamy.type3 = subset(data.summary[data.summary$sex == "female" & data.summary$treatment == "monogamy" & data.summary$population.type==3, ])$adjusted.variance
female.data.variance.monogamy.type5 = subset(data.summary[data.summary$sex == "female" & data.summary$treatment == "monogamy" & data.summary$population.type==5, ])$adjusted.variance


#####   Comparison of variances between the sexes    #####
# Vials 
print( mean(c(male.data.variance.vial.type1,male.data.variance.vial.type2,male.data.variance.vial.type3,male.data.variance.vial.type4,male.data.variance.vial.type5)) - mean(c(female.data.variance.vial.type1,female.data.variance.vial.type2,female.data.variance.vial.type3,female.data.variance.vial.type4,female.data.variance.vial.type5)))
# Cage
print( mean(c(male.data.variance.cage.type1,male.data.variance.cage.type2,male.data.variance.cage.type3,male.data.variance.cage.type4,male.data.variance.cage.type5)) - mean(c(female.data.variance.cage.type1,female.data.variance.cage.type2,female.data.variance.cage.type3,female.data.variance.cage.type4,female.data.variance.cage.type5)))
# Monogamy
print( mean(c(male.data.variance.monogamy.type1,male.data.variance.monogamy.type3,male.data.variance.monogamy.type5)) - mean(c(female.data.variance.monogamy.type1,female.data.variance.monogamy.type3,female.data.variance.monogamy.type5)))


####       Comparison of variances across environments (i.e., heterogeniety comparisons)        #####
## Vials
male.vial.low = mean(c(male.data.variance.vial.type1,male.data.variance.vial.type5))
male.vial.high = mean(c(male.data.variance.vial.type3))
female.vial.low = mean(c(female.data.variance.vial.type1,female.data.variance.vial.type5))
female.vial.high = mean(c(female.data.variance.vial.type3))
mf.vial.low = mean(c(male.data.variance.vial.type1,male.data.variance.vial.type5,female.data.variance.vial.type1,female.data.variance.vial.type5))
mf.vial.high = mean(c(male.data.variance.vial.type3,female.data.variance.vial.type3))
print(male.vial.high - male.vial.low)
print(female.vial.high - female.vial.low)
print(mf.vial.high - mf.vial.low)
## Cages
male.cage.low = mean(c(male.data.variance.cage.type1,male.data.variance.cage.type5))
male.cage.high = mean(c(male.data.variance.cage.type3))
female.cage.low = mean(c(female.data.variance.cage.type1,female.data.variance.cage.type5))
female.cage.high = mean(c(female.data.variance.cage.type3))
mf.cage.low = mean(c(male.data.variance.cage.type1,male.data.variance.cage.type5,female.data.variance.cage.type1,female.data.variance.cage.type5))
mf.cage.high = mean(c(male.data.variance.cage.type3,female.data.variance.cage.type3))
print(male.cage.high - male.cage.low)
print(female.cage.high - female.cage.low)
print(mf.cage.high - mf.cage.low)
## Monogamy
male.monogamy.low = mean(c(male.data.variance.monogamy.type1,male.data.variance.monogamy.type5))
male.monogamy.high = mean(c(male.data.variance.monogamy.type3))
female.monogamy.low = mean(c(female.data.variance.monogamy.type1,female.data.variance.monogamy.type5))
female.monogamy.high = mean(c(female.data.variance.monogamy.type3))
mf.monogamy.low = mean(c(male.data.variance.monogamy.type1,male.data.variance.monogamy.type5,female.data.variance.monogamy.type1,female.data.variance.monogamy.type5))
mf.monogamy.high = mean(c(male.data.variance.monogamy.type3,female.data.variance.monogamy.type3))
print(male.monogamy.high - male.monogamy.low)
print(female.monogamy.high - female.monogamy.low)
print(mf.monogamy.high - mf.monogamy.low)

####       Comparison of variances across mating regimes        #####
## Polygamy in Vials 
male.vial=mean(c(male.data.variance.vial.type1,male.data.variance.vial.type2,male.data.variance.vial.type3,male.data.variance.vial.type4,male.data.variance.vial.type5))
female.vial=mean(c(female.data.variance.vial.type1,female.data.variance.vial.type2,female.data.variance.vial.type3,female.data.variance.vial.type4,female.data.variance.vial.type5))
male.cage=mean(c(male.data.variance.cage.type1,male.data.variance.cage.type2,male.data.variance.cage.type3,male.data.variance.cage.type4,male.data.variance.cage.type5))
female.cage=mean(c(female.data.variance.cage.type1,female.data.variance.cage.type2,female.data.variance.cage.type3,female.data.variance.cage.type4,female.data.variance.cage.type5))
male.monogamy = mean(c(male.data.variance.monogamy.type1,male.data.variance.monogamy.type3,male.data.variance.monogamy.type5))
female.monogamy = mean(c(female.data.variance.monogamy.type1,female.data.variance.monogamy.type3,female.data.variance.monogamy.type5))
mf.vial = mean(c(male.data.variance.vial.type1,male.data.variance.vial.type2,male.data.variance.vial.type3,male.data.variance.vial.type4,male.data.variance.vial.type5,
                 female.data.variance.vial.type1,female.data.variance.vial.type2,female.data.variance.vial.type3,female.data.variance.vial.type4,female.data.variance.vial.type5))
mf.cage = mean(c(male.data.variance.cage.type1,male.data.variance.cage.type2,male.data.variance.cage.type3,male.data.variance.cage.type4,male.data.variance.cage.type5,
                 female.data.variance.cage.type1,female.data.variance.cage.type2,female.data.variance.cage.type3,female.data.variance.cage.type4,female.data.variance.cage.type5))
mf.monogamy = mean(c(male.data.variance.monogamy.type1,male.data.variance.monogamy.type2,male.data.variance.monogamy.type3,male.data.variance.monogamy.type4,male.data.variance.monogamy.type5,
                 female.data.variance.monogamy.type1,female.data.variance.monogamy.type2,female.data.variance.monogamy.type3,female.data.variance.monogamy.type4,female.data.variance.monogamy.type5))

print(mean(male.vial - male.cage))
print(mean(female.vial - female.cage))
print(mean(male.vial - male.monogamy))
print(mean(female.vial - female.monogamy))
print(mean(male.cage - male.monogamy))
print(mean(female.cage - female.monogamy))
print(mean(mf.vial - mf.cage))
print(mean(mf.vial - mf.monogamy))
print(mean(mf.cage - mf.monogamy))



################ Obtaining bootstrapped confidence intervals ################ 

# subset out each combination of sex*envronment*mating regime 
male.bootstrapped.variance.vial.type1 = subset(bootstrapped.data[bootstrapped.data$sex == "male" & bootstrapped.data$treatment == "vial" & bootstrapped.data$population.type==1, ])$adjusted.variance
male.bootstrapped.variance.vial.type2 = subset(bootstrapped.data[bootstrapped.data$sex == "male" & bootstrapped.data$treatment == "vial" & bootstrapped.data$population.type==2, ])$adjusted.variance
male.bootstrapped.variance.vial.type3 = subset(bootstrapped.data[bootstrapped.data$sex == "male" & bootstrapped.data$treatment == "vial" & bootstrapped.data$population.type==3, ])$adjusted.variance
male.bootstrapped.variance.vial.type4 = subset(bootstrapped.data[bootstrapped.data$sex == "male" & bootstrapped.data$treatment == "vial" & bootstrapped.data$population.type==4, ])$adjusted.variance
male.bootstrapped.variance.vial.type5 = subset(bootstrapped.data[bootstrapped.data$sex == "male" & bootstrapped.data$treatment == "vial" & bootstrapped.data$population.type==5, ])$adjusted.variance
female.bootstrapped.variance.vial.type1 = subset(bootstrapped.data[bootstrapped.data$sex == "female" & bootstrapped.data$treatment == "vial" & bootstrapped.data$population.type==1, ])$adjusted.variance
female.bootstrapped.variance.vial.type2 = subset(bootstrapped.data[bootstrapped.data$sex == "female" & bootstrapped.data$treatment == "vial" & bootstrapped.data$population.type==2, ])$adjusted.variance
female.bootstrapped.variance.vial.type3 = subset(bootstrapped.data[bootstrapped.data$sex == "female" & bootstrapped.data$treatment == "vial" & bootstrapped.data$population.type==3, ])$adjusted.variance
female.bootstrapped.variance.vial.type4 = subset(bootstrapped.data[bootstrapped.data$sex == "female" & bootstrapped.data$treatment == "vial" & bootstrapped.data$population.type==4, ])$adjusted.variance
female.bootstrapped.variance.vial.type5 = subset(bootstrapped.data[bootstrapped.data$sex == "female" & bootstrapped.data$treatment == "vial" & bootstrapped.data$population.type==5, ])$adjusted.variance
#Cages 
male.bootstrapped.variance.cage.type1 = subset(bootstrapped.data[bootstrapped.data$sex == "male" & bootstrapped.data$treatment == "cage" & bootstrapped.data$population.type==1, ])$adjusted.variance
male.bootstrapped.variance.cage.type2 = subset(bootstrapped.data[bootstrapped.data$sex == "male" & bootstrapped.data$treatment == "cage" & bootstrapped.data$population.type==2, ])$adjusted.variance
male.bootstrapped.variance.cage.type3 = subset(bootstrapped.data[bootstrapped.data$sex == "male" & bootstrapped.data$treatment == "cage" & bootstrapped.data$population.type==3, ])$adjusted.variance
male.bootstrapped.variance.cage.type4 = subset(bootstrapped.data[bootstrapped.data$sex == "male" & bootstrapped.data$treatment == "cage" & bootstrapped.data$population.type==4, ])$adjusted.variance
male.bootstrapped.variance.cage.type5 = subset(bootstrapped.data[bootstrapped.data$sex == "male" & bootstrapped.data$treatment == "cage" & bootstrapped.data$population.type==5, ])$adjusted.variance
female.bootstrapped.variance.cage.type1 = subset(bootstrapped.data[bootstrapped.data$sex == "female" & bootstrapped.data$treatment == "cage" & bootstrapped.data$population.type==1, ])$adjusted.variance
female.bootstrapped.variance.cage.type2 = subset(bootstrapped.data[bootstrapped.data$sex == "female" & bootstrapped.data$treatment == "cage" & bootstrapped.data$population.type==2, ])$adjusted.variance
female.bootstrapped.variance.cage.type3 = subset(bootstrapped.data[bootstrapped.data$sex == "female" & bootstrapped.data$treatment == "cage" & bootstrapped.data$population.type==3, ])$adjusted.variance
female.bootstrapped.variance.cage.type4 = subset(bootstrapped.data[bootstrapped.data$sex == "female" & bootstrapped.data$treatment == "cage" & bootstrapped.data$population.type==4, ])$adjusted.variance
female.bootstrapped.variance.cage.type5 = subset(bootstrapped.data[bootstrapped.data$sex == "female" & bootstrapped.data$treatment == "cage" & bootstrapped.data$population.type==5, ])$adjusted.variance
#Monogamy
male.bootstrapped.variance.monogamy.type1 = subset(bootstrapped.data[bootstrapped.data$sex == "male" & bootstrapped.data$treatment == "monogamy" & bootstrapped.data$population.type==1, ])$adjusted.variance
male.bootstrapped.variance.monogamy.type3 = subset(bootstrapped.data[bootstrapped.data$sex == "male" & bootstrapped.data$treatment == "monogamy" & bootstrapped.data$population.type==3, ])$adjusted.variance
male.bootstrapped.variance.monogamy.type5 = subset(bootstrapped.data[bootstrapped.data$sex == "male" & bootstrapped.data$treatment == "monogamy" & bootstrapped.data$population.type==5, ])$adjusted.variance
female.bootstrapped.variance.monogamy.type1 = subset(bootstrapped.data[bootstrapped.data$sex == "female" & bootstrapped.data$treatment == "monogamy" & bootstrapped.data$population.type==1, ])$adjusted.variance
female.bootstrapped.variance.monogamy.type3 = subset(bootstrapped.data[bootstrapped.data$sex == "female" & bootstrapped.data$treatment == "monogamy" & bootstrapped.data$population.type==3, ])$adjusted.variance
female.bootstrapped.variance.monogamy.type5 = subset(bootstrapped.data[bootstrapped.data$sex == "female" & bootstrapped.data$treatment == "monogamy" & bootstrapped.data$population.type==5, ])$adjusted.variance


#####   Comparison of means and variances between the sexes    #####
male.bootdata.variance.vial=c()
female.bootdata.variance.vial=c()
male.bootdata.variance.cage=c()
female.bootdata.variance.cage=c()
male.bootdata.variance.monogamy=c()
female.bootdata.variance.monogamy=c()

for (i in 1:10000){
  male.bootdata.variance.vial[i] = mean(c(male.bootstrapped.variance.vial.type1[i], male.bootstrapped.variance.vial.type2[i], male.bootstrapped.variance.vial.type3[i], male.bootstrapped.variance.vial.type4[i], male.bootstrapped.variance.vial.type5[i]))
  female.bootdata.variance.vial[i] = mean(c(female.bootstrapped.variance.vial.type1[i], female.bootstrapped.variance.vial.type2[i], female.bootstrapped.variance.vial.type3[i], female.bootstrapped.variance.vial.type4[i], female.bootstrapped.variance.vial.type5[i]))
  male.bootdata.variance.cage[i] = mean(c(male.bootstrapped.variance.cage.type1[i], male.bootstrapped.variance.cage.type2[i], male.bootstrapped.variance.cage.type3[i], male.bootstrapped.variance.cage.type4[i], male.bootstrapped.variance.cage.type5[i]))
  female.bootdata.variance.cage[i] = mean(c(female.bootstrapped.variance.cage.type1[i], female.bootstrapped.variance.cage.type2[i], female.bootstrapped.variance.cage.type3[i], female.bootstrapped.variance.cage.type4[i], female.bootstrapped.variance.cage.type5[i]))
  male.bootdata.variance.monogamy[i] = mean(c(male.bootstrapped.variance.monogamy.type1[i], male.bootstrapped.variance.monogamy.type3[i], male.bootstrapped.variance.monogamy.type5[i]))
  female.bootdata.variance.monogamy[i] = mean(c(female.bootstrapped.variance.monogamy.type1[i], female.bootstrapped.variance.monogamy.type3[i], female.bootstrapped.variance.monogamy.type5[i]))
}
vial.male.female.diff = (male.bootdata.variance.vial-female.bootdata.variance.vial)
cage.male.female.diff = (male.bootdata.variance.cage-female.bootdata.variance.cage)
monogamy.male.female.diff = (male.bootdata.variance.monogamy-female.bootdata.variance.monogamy)
print(c(mean(vial.male.female.diff), quantile(vial.male.female.diff,c(0.025,0.975))))
print(c(mean(cage.male.female.diff), quantile(cage.male.female.diff,c(0.025,0.975))))
print(c(mean(monogamy.male.female.diff), quantile(monogamy.male.female.diff,c(0.025,0.975))))

#####   Comparison of variances across environmental types (i.e., heterogenity treatments    #####

male.bootdata.variance.vial.high=c()
male.bootdata.variance.vial.low=c()
female.bootdata.variance.vial.high=c()
female.bootdata.variance.vial.low=c()
male.bootdata.variance.cage.high=c()
male.bootdata.variance.cage.low=c()
female.bootdata.variance.cage.high=c()
female.bootdata.variance.cage.low=c()
male.bootdata.variance.monogamy.high=c()
male.bootdata.variance.monogamy.low=c()
female.bootdata.variance.monogamy.high=c()
female.bootdata.variance.monogamy.low=c()
mf.bootdata.variance.vial.high=c()
mf.bootdata.variance.vial.low=c()
mf.bootdata.variance.cage.high=c()
mf.bootdata.variance.cage.low=c()
mf.bootdata.variance.monogamy.high=c()
mf.bootdata.variance.monogamy.low=c()
for (i in 1:10000){
  male.bootdata.variance.vial.low[i] = mean(c(male.bootstrapped.variance.vial.type1[i],male.bootstrapped.variance.vial.type5[i]))
  male.bootdata.variance.vial.high[i] = mean(c(male.bootstrapped.variance.vial.type3[i]))
  female.bootdata.variance.vial.low[i] = mean(c(female.bootstrapped.variance.vial.type1[i],female.bootstrapped.variance.vial.type5[i]))
  female.bootdata.variance.vial.high[i] = mean(c(female.bootstrapped.variance.vial.type3[i]))
  male.bootdata.variance.cage.low[i] = mean(c(male.bootstrapped.variance.cage.type1[i],male.bootstrapped.variance.cage.type5[i]))
  male.bootdata.variance.cage.high[i] = mean(c(male.bootstrapped.variance.cage.type3[i]))
  female.bootdata.variance.cage.low[i] = mean(c(female.bootstrapped.variance.cage.type1[i],female.bootstrapped.variance.cage.type5[i]))
  female.bootdata.variance.cage.high[i] = mean(c(female.bootstrapped.variance.cage.type3[i]))
  male.bootdata.variance.monogamy.low[i] = mean(c(male.bootstrapped.variance.monogamy.type1[i],male.bootstrapped.variance.monogamy.type5[i]))
  male.bootdata.variance.monogamy.high[i] = mean(c(male.bootstrapped.variance.monogamy.type3[i]))
  female.bootdata.variance.monogamy.low[i] = mean(c(female.bootstrapped.variance.monogamy.type1[i],female.bootstrapped.variance.monogamy.type5[i]))
  female.bootdata.variance.monogamy.high[i] = mean(c(female.bootstrapped.variance.monogamy.type3[i]))
  
  mf.bootdata.variance.vial.high[i] = mean(c(male.bootstrapped.variance.vial.type3[i],female.bootstrapped.variance.vial.type3[i]))
  mf.bootdata.variance.vial.low[i] = mean(c(male.bootstrapped.variance.vial.type1[i],female.bootstrapped.variance.vial.type1[i],male.bootstrapped.variance.vial.type5[i],female.bootstrapped.variance.vial.type5[i]))
  mf.bootdata.variance.cage.high[i] = mean(c(male.bootstrapped.variance.cage.type3[i],female.bootstrapped.variance.cage.type3[i]))
  mf.bootdata.variance.cage.low[i] = mean(c(male.bootstrapped.variance.cage.type1[i],female.bootstrapped.variance.cage.type1[i],male.bootstrapped.variance.cage.type5[i],female.bootstrapped.variance.cage.type5[i]))
  mf.bootdata.variance.monogamy.high[i] = mean(c(male.bootstrapped.variance.monogamy.type3[i],female.bootstrapped.variance.monogamy.type3[i]))
  mf.bootdata.variance.monogamy.low[i] = mean(c(male.bootstrapped.variance.monogamy.type1[i],female.bootstrapped.variance.monogamy.type1[i],male.bootstrapped.variance.monogamy.type5[i],female.bootstrapped.variance.monogamy.type5[i]))
}
male.vial.high.low = male.bootdata.variance.vial.high - male.bootdata.variance.vial.low
female.vial.high.low = female.bootdata.variance.vial.high - female.bootdata.variance.vial.low
male.cage.high.low = male.bootdata.variance.cage.high - male.bootdata.variance.cage.low
female.cage.high.low = female.bootdata.variance.cage.high - female.bootdata.variance.cage.low
male.monogamy.high.low = male.bootdata.variance.monogamy.high - male.bootdata.variance.monogamy.low
female.monogamy.high.low = female.bootdata.variance.monogamy.high - female.bootdata.variance.monogamy.low
mf.vial.high.low = mf.bootdata.variance.vial.high - mf.bootdata.variance.vial.low
mf.cage.high.low = mf.bootdata.variance.cage.high - mf.bootdata.variance.cage.low
mf.monogamy.high.low = mf.bootdata.variance.monogamy.high - mf.bootdata.variance.monogamy.low

print(c(mean(male.vial.high.low), quantile(male.vial.high.low,c(0.025, 0.975))))
print(c(mean(female.vial.high.low), quantile(female.vial.high.low,c(0.025, 0.975))))
print(c(mean(male.cage.high.low), quantile(male.cage.high.low,c(0.025, 0.975))))
print(c(mean(female.cage.high.low), quantile(female.cage.high.low,c(0.025, 0.975))))
print(c(mean(male.monogamy.high.low), quantile(male.monogamy.high.low,c(0.025, 0.975))))
print(c(mean(female.monogamy.high.low), quantile(female.monogamy.high.low,c(0.025, 0.975))))
print(c(mean(mf.vial.high.low), quantile(mf.vial.high.low, c(0.025,0.975))))
print(c(mean(mf.cage.high.low), quantile(mf.cage.high.low, c(0.025,0.975))))
print(c(mean(mf.monogamy.high.low), quantile(mf.monogamy.high.low, c(0.025,0.975))))



#####   Comparison of variances across mating regimes   #####

male.bootdata.variance.vial=c()
male.bootdata.variance.cage=c()
male.bootdata.variance.monogamy=c()
female.bootdata.variance.vial=c()
female.bootdata.variance.cage=c()
female.bootdata.variance.monogamy=c()
mf.bootdata.variance.vial=c()
mf.bootdata.variance.cage=c()
mf.bootdata.variance.monogamy=c()

for (i in 1:10000){
  male.bootdata.variance.vial[i] = mean(c(male.bootstrapped.variance.vial.type1[i],male.bootstrapped.variance.vial.type2[i],male.bootstrapped.variance.vial.type3[i],male.bootstrapped.variance.vial.type4[i],male.bootstrapped.variance.vial.type5[i]))
  male.bootdata.variance.cage[i] = mean(c(male.bootstrapped.variance.cage.type1[i],male.bootstrapped.variance.cage.type2[i],male.bootstrapped.variance.cage.type3[i],male.bootstrapped.variance.cage.type4[i],male.bootstrapped.variance.cage.type5[i]))
  male.bootdata.variance.monogamy[i] = mean(c(male.bootstrapped.variance.monogamy.type1[i],male.bootstrapped.variance.monogamy.type3[i],male.bootstrapped.variance.monogamy.type5[i]))
  
  female.bootdata.variance.vial[i] = mean(c(female.bootstrapped.variance.vial.type1[i],female.bootstrapped.variance.vial.type2[i],female.bootstrapped.variance.vial.type3[i],female.bootstrapped.variance.vial.type4[i],female.bootstrapped.variance.vial.type5[i]))
  female.bootdata.variance.cage[i] = mean(c(female.bootstrapped.variance.cage.type1[i],female.bootstrapped.variance.cage.type2[i],female.bootstrapped.variance.cage.type3[i],female.bootstrapped.variance.cage.type4[i],female.bootstrapped.variance.cage.type5[i]))
  female.bootdata.variance.monogamy[i] = mean(c(female.bootstrapped.variance.monogamy.type1[i],female.bootstrapped.variance.monogamy.type3[i],female.bootstrapped.variance.monogamy.type5[i]))
 
  mf.bootdata.variance.vial[i] = mean(c(male.bootstrapped.variance.vial.type1[i],male.bootstrapped.variance.vial.type2[i],male.bootstrapped.variance.vial.type3[i],male.bootstrapped.variance.vial.type4[i],male.bootstrapped.variance.vial.type5[i],
                                        female.bootstrapped.variance.vial.type1[i],female.bootstrapped.variance.vial.type2[i],female.bootstrapped.variance.vial.type3[i],female.bootstrapped.variance.vial.type4[i],female.bootstrapped.variance.vial.type5[i]))
  mf.bootdata.variance.cage[i] = mean(c(male.bootstrapped.variance.cage.type1[i],male.bootstrapped.variance.cage.type2[i],male.bootstrapped.variance.cage.type3[i],male.bootstrapped.variance.cage.type4[i],male.bootstrapped.variance.cage.type5[i],
                                        female.bootstrapped.variance.cage.type1[i],female.bootstrapped.variance.cage.type2[i],female.bootstrapped.variance.cage.type3[i],female.bootstrapped.variance.cage.type4[i],female.bootstrapped.variance.cage.type5[i]))
  
  mf.bootdata.variance.monogamy[i] = mean(c(male.bootstrapped.variance.monogamy.type1[i],male.bootstrapped.variance.monogamy.type3[i],male.bootstrapped.variance.monogamy.type5[i],
                                            female.bootstrapped.variance.monogamy.type1[i],female.bootstrapped.variance.monogamy.type3[i],female.bootstrapped.variance.monogamy.type5[i]))

}

# Make vectors of comparisons 
vial.cage.male = male.bootdata.variance.vial - male.bootdata.variance.cage
vial.cage.female = female.bootdata.variance.vial - female.bootdata.variance.cage
vial.monogamy.male = male.bootdata.variance.vial - male.bootdata.variance.monogamy
vial.monogamy.female = female.bootdata.variance.vial - female.bootdata.variance.monogamy
cage.monogamy.male = male.bootdata.variance.cage - male.bootdata.variance.monogamy
cage.monogamy.female = female.bootdata.variance.cage - female.bootdata.variance.monogamy
vial.cage.mf = mf.bootdata.variance.vial- mf.bootdata.variance.cage
vial.monogamy.mf = mf.bootdata.variance.vial- mf.bootdata.variance.monogamy
cage.monogamy.mf = mf.bootdata.variance.cage- mf.bootdata.variance.monogamy

print(c(mean(vial.cage.male), quantile(vial.cage.male, c(0.025,0.975))))
print(c(mean(vial.cage.female), quantile(vial.cage.female, c(0.025,0.975))))
print(c(mean(vial.monogamy.male), quantile(vial.monogamy.male, c(0.025,0.975))))
print(c(mean(vial.monogamy.female), quantile(vial.monogamy.female, c(0.025,0.975))))
print(c(mean(cage.monogamy.male), quantile(cage.monogamy.male, c(0.025,0.975))))
print(c(mean(cage.monogamy.female), quantile(cage.monogamy.female, c(0.025,0.975))))
print(c(mean(vial.cage.mf), quantile(vial.cage.mf, c(0.025,0.975))))
print(c(mean(vial.monogamy.mf), quantile(vial.monogamy.mf, c(0.025,0.975))))
print(c(mean(cage.monogamy.mf), quantile(cage.monogamy.mf, c(0.025,0.975))))



##

