#
###########################################################################################################################################
###########################################################################################################################################
###   Script Details:                                                                                                                   ###
###   Bodysize Analysis and Plots for Singh, A. A.F., Agrawal. 2020. Sex-Specific Variance in Fitness and the Efficacy of Selection     ###
###   Script written by Amardeep Singh (amardeep.singh[at]utoronto.ca)                                                                  ###
###   This script will summarize sex-specific bodymass data from experiment either seperate for each stock, or avageraing across stocks ###
###########################################################################################################################################
###########################################################################################################################################

# NOTE: YOU WILL NEED TO SPECIFY THE PATH TO DATA FILES BELOW BY REPLACING "/FILE/PATH/TO/" WITH YOUR ACTUAL FILE PATHS 
# NOTE: Plots will not look exactly like those in the paper (the data is identical, but I processed images to change fonts, labels etc.)

# Plots generated with this script: Figure S1 and Figure S2

###########################################
## 1. Setting up Packages and Data files ##
###########################################
rm(list=ls())

##  Loading packages
require(doBy)
require(ggplot2)
require(wesanderson) # Don't think I used this, but maybe I'll go back and make the files look nice

##  Reading body size data into R
bodysize.data=read.csv("/FILE/PATH/TO/DATA/BodyMassDataNeprojectFinalSept2020.csv")


###############################################
## 2. Summarizing Data for Visual Inspection ##
###############################################

# First summary is summarizing mean mass per fly by block, stock, sex, treatment and replicate 
bodysize.data.summary = summaryBy(mean.mass.per.fly ~ block+stock+sex+treatment+replicate, data=bodysize.data, FUN=c(length,mean,sd))
bodysize.data.summary$se = bodysize.data.summary$mean.mass.per.fly.mean.sd / sqrt(bodysize.data.summary$mean.mass.per.fly.mean.length)

# Second summary is summarizing mean mass per fly by block, stock, sex, treatment
bodysize.data.summary.averaged.across.blocks = summaryBy(mean.mass.per.fly.mean ~ stock + sex + treatment, data = bodysize.data.summary, FUN=c(length,mean,sd))
bodysize.data.summary.averaged.across.blocks$se = bodysize.data.summary.averaged.across.blocks$mean.mass.per.fly.mean.sd / sqrt(bodysize.data.summary.averaged.across.blocks$mean.mass.per.fly.mean.length)

# Third summary is summarizing mean mass per fly by block, stock and sex 
bodysize.data.summary.averaged.across.blocks.and.treatment = summaryBy(mean.mass.per.fly.mean ~ sex + treatment, data = bodysize.data.summary, FUN=c(length,mean,sd))
bodysize.data.summary.averaged.across.blocks.and.treatment$se = bodysize.data.summary.averaged.across.blocks.and.treatment$mean.mass.per.fly.mean.sd / sqrt(bodysize.data.summary.averaged.across.blocks.and.treatment$mean.mass.per.fly.mean.length)

##############################################
## 3. Data Analysis with Linear Regression  ##
##############################################

# The first model  will analyze at each stock seperately 
model.bodysize.separate.stocks = lm(mean.mass.per.fly.mean ~ sex + stock + treatment+ sex*treatment, data = bodysize.data.summary)
# The second model will ignore stocks and look only at the effect of sex and larval treatment and their interaction 
model.bodysize.all = lm(mean.mass.per.fly.mean ~ sex + treatment+ sex*treatment, data = bodysize.data.summary)

################################
## 4. Plotting Body Size Data ##
################################
# We plotting the LS means from the linear models above 

# Extract LS means from models for plotting 
model.bodysize.lsmeans.separate.stocks = lsmeans(model.bodysize.separate.stocks, ~sex + stock + treatment + sex*treatment)
model.bodysize.lsmeans.all = lsmeans(model.bodysize.all, ~sex + treatment  + sex*treatment)

# Plotting LS mean mass averaged across all stocks (i.e., summing up across different stocks and looking only at differnce b/w males and females in high and low larval density treatments) 

# LS means -- Each stock individually
# This corresponds to Figure S1
lsmeans.plot.by.stock = ggplot(as.data.frame(model.bodysize.lsmeans.separate.stocks), aes(y = lsmean, x = sex, fill = treatment)) + 
  geom_bar(stat="identity", color = "black", position=position_dodge()) +
  geom_errorbar(aes(ymin = (as.data.frame(model.bodysize.lsmeans.separate.stocks))$lower.CL, ymax = (as.data.frame(model.bodysize.lsmeans.separate.stocks))$upper.CL), colour = "black",  position = position_dodge(0.9), width = 0.0, size = 1.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(size = 25), axis.title = element_text(size = 40, colour = "black"),
        axis.title.x = element_text(size = 40), axis.text.x=element_text(colour = "black", family = "Times"), axis.ticks.x=element_line(colour = "black"),
        axis.title.y = element_text(size = 40), axis.text.y=element_text(colour = "black", family = "Times"), axis.ticks.y=element_line(colour = "black"), 
        strip.background = element_rect(colour = "black", fill = "light grey"), strip.text = element_text(size = 40, family = "Times"), 
        legend.text = element_text(color = "black", size = 25, family = "Times"),  legend.key.size = unit(1, "cm")) + 
  scale_y_continuous(limits = c(0,0.4), expand = c(0,0)) +
  scale_fill_manual(name = "", labels = c("High Density", "Low Density"), values = c("#E0E0E0","#808080")) + 
  labs(x = "", y = "") +
  facet_wrap(~stock)
lsmeans.plot.by.stock

# LS means -- Averaged across all stocks 
# This corresponds to Figure S2
lsmeans.plot.all = ggplot(as.data.frame(model.bodysize.lsmeans.all), aes(y = lsmean, x = sex, fill = treatment)) + 
  geom_bar(stat="identity", color = "black", position=position_dodge()) +
  geom_errorbar(aes(ymin = (as.data.frame(model.bodysize.lsmeans.all))$lower.CL, ymax = (as.data.frame(model.bodysize.lsmeans.all))$upper.CL), colour = "black",  position = position_dodge(0.9), width = 0.0, size = 1.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(size = 25), axis.title = element_text(size = 40, colour = "black"),
        axis.title.x = element_text(size = 40), axis.text.x=element_text(colour = "black", family = "Times"), axis.ticks.x=element_line(colour = "black"),
        axis.title.y = element_text(size = 40), axis.text.y=element_text(colour = "black", family = "Times"), axis.ticks.y=element_line(colour = "black"), 
        strip.background = element_rect(colour = "black", fill = "light grey"), strip.text = element_text(size = 40, family = "Times"), 
        legend.text = element_text(color = "black", size = 25, family = "Times"),  legend.key.size = unit(1, "cm")) + 
  scale_y_continuous(limits = c(0,0.35), expand = c(0,0)) +
  scale_fill_manual(name = "", labels = c("High Density", "Low Density"), values = c("#E0E0E0","#808080")) + 
  labs(x = "", y = "")
lsmeans.plot.all

# Saving body size plots 
pdf("FILE/PATH/TO/DATA/Figure S1 - Mean_Mass_per_stock.pdf", width = 16*3, height = 20)
lsmeans.plot.by.stock
dev.off()

pdf("/FILE/PATH/TO/DATA/Figure S2 - Mean_Mass_averaged_across_stocks.pdf", width = 16, height = 20)
lsmeans.plot.all
dev.off()

