#########################################################################################################################################################################
#########################################################################################################################################################################
###   Script Details:                                                                                                                                                 ###
###   Bootstrapping procedure for mean and variance in fitness for Singh, A. and A.F., Agrawal. 2021. Sex-Specific Variance in Fitness and the Efficacy of Selection  ###
###   Script written by Amardeep Singh (amardeep.singh[at]alum.utoronto.ca)                                                                                           ###
#########################################################################################################################################################################
#########################################################################################################################################################################

# NOTE: This code could be written in far fewer lines but I wrote it many many moons ago when I was relatively new to coding/writing in R ¯\_(ツ)_/¯
# This script will take the raw fitness data and bootstrap the mean and variance in sex-specific fitness for each population type and mating regime combination 10,000 times
#   and then save the output to be used downstream.

## Loading in R packages 
require(Rcurl)

# Loading fitness data into R and cleaning it up a bit
data.github <- getURL("https://raw.githubusercontent.com/asingh164/SexSpecificVarianceEfficacyOfSelection/master/FitnessDataNeprojectFinalSept2020.csv")
data=read.csv(text = data.github)
data=as.data.frame(data)
data$n.wt=as.numeric(as.character(data$n.wt))
data=data[data$n.total == "32",] # Remove any observation where we have less than 32 offspring counted
data=data[!(is.na(data$n.wt)),] # Remove any rows with NA in fitness column 

# Subsetting out data
vialsubsetpop1male=subset(data, data$treatment=="vial" & data$population=="1" & data$sex=="m")
vialsubsetpop5male=subset(data, data$treatment=="vial" & data$population=="5" & data$sex=="m")
vialsubsetpop3malelow=subset(data, data$treatment=="vial" & data$population=="3" & data$sex=="m" & data$focal.individual=="low")
vialsubsetpop3malehigh=subset(data, data$treatment=="vial" & data$population=="3" & data$sex=="m" & data$focal.individual=="high")
vialsubsetpop2malelow=subset(data, data$treatment=="vial" & data$population=="2" & data$sex=="m" & data$focal.individual=="low")
vialsubsetpop2malehigh=subset(data, data$treatment=="vial" & data$population=="2" & data$sex=="m" & data$focal.individual=="high")
vialsubsetpop4malelow=subset(data, data$treatment=="vial" & data$population=="4" & data$sex=="m" & data$focal.individual=="low")
vialsubsetpop4malehigh=subset(data, data$treatment=="vial" & data$population=="4" & data$sex=="m" & data$focal.individual=="high")

vialsubsetpop1female=subset(data, data$treatment=="vial" & data$population=="1" & data$sex=="f")
vialsubsetpop5female=subset(data, data$treatment=="vial" & data$population=="5" & data$sex=="f")
vialsubsetpop3femalelow=subset(data, data$treatment=="vial" & data$population=="3" & data$sex=="f" & data$focal.individual=="low")
vialsubsetpop3femalehigh=subset(data, data$treatment=="vial" & data$population=="3" & data$sex=="f" & data$focal.individual=="high")
vialsubsetpop2femalelow=subset(data, data$treatment=="vial" & data$population=="2" & data$sex=="f" & data$focal.individual=="low")
vialsubsetpop2femalehigh=subset(data, data$treatment=="vial" & data$population=="2" & data$sex=="f" & data$focal.individual=="high")
vialsubsetpop4femalelow=subset(data, data$treatment=="vial" & data$population=="4" & data$sex=="f" & data$focal.individual=="low")
vialsubsetpop4femalehigh=subset(data, data$treatment=="vial" & data$population=="4" & data$sex=="f" & data$focal.individual=="high")

cagesubsetpop1male=subset(data, data$treatment=="cage" & data$population=="1" & data$sex=="m")
cagesubsetpop5male=subset(data, data$treatment=="cage" & data$population=="5" & data$sex=="m")
cagesubsetpop3malelow=subset(data, data$treatment=="cage" & data$population=="3" & data$sex=="m" & data$focal.individual=="low")
cagesubsetpop3malehigh=subset(data, data$treatment=="cage" & data$population=="3" & data$sex=="m" & data$focal.individual=="high")
cagesubsetpop2malelow=subset(data, data$treatment=="cage" & data$population=="2" & data$sex=="m" & data$focal.individual=="low")
cagesubsetpop2malehigh=subset(data, data$treatment=="cage" & data$population=="2" & data$sex=="m" & data$focal.individual=="high")
cagesubsetpop4malelow=subset(data, data$treatment=="cage" & data$population=="4" & data$sex=="m" & data$focal.individual=="low")
cagesubsetpop4malehigh=subset(data, data$treatment=="cage" & data$population=="4" & data$sex=="m" & data$focal.individual=="high")

cagesubsetpop1female=subset(data, data$treatment=="cage" & data$population=="1" & data$sex=="f")
cagesubsetpop5female=subset(data, data$treatment=="cage" & data$population=="5" & data$sex=="f")
cagesubsetpop3femalelow=subset(data, data$treatment=="cage" & data$population=="3" & data$sex=="f" & data$focal.individual=="low")
cagesubsetpop3femalehigh=subset(data, data$treatment=="cage" & data$population=="3" & data$sex=="f" & data$focal.individual=="high")
cagesubsetpop2femalelow=subset(data, data$treatment=="cage" & data$population=="2" & data$sex=="f" & data$focal.individual=="low")
cagesubsetpop2femalehigh=subset(data, data$treatment=="cage" & data$population=="2" & data$sex=="f" & data$focal.individual=="high")
cagesubsetpop4femalelow=subset(data, data$treatment=="cage" & data$population=="4" & data$sex=="f" & data$focal.individual=="low")
cagesubsetpop4femalehigh=subset(data, data$treatment=="cage" & data$population=="4" & data$sex=="f" & data$focal.individual=="high")

monogamysubsetpop1male=subset(data, data$treatment=="monogamy" & data$population=="1" & data$sex=="m")
monogamysubsetpop5male=subset(data, data$treatment=="monogamy" & data$population=="5" & data$sex=="m")
monogamysubsetpop3malelow=subset(data, data$treatment=="monogamy" & data$population=="3" & data$sex=="m" & data$focal.individual=="low")
monogamysubsetpop3malehigh=subset(data, data$treatment=="monogamy" & data$population=="3" & data$sex=="m" & data$focal.individual=="high")

monogamysubsetpop1female=subset(data, data$treatment=="monogamy" & data$population=="1" & data$sex=="f")
monogamysubsetpop5female=subset(data, data$treatment=="monogamy" & data$population=="5" & data$sex=="f")
monogamysubsetpop3femalelow=subset(data, data$treatment=="monogamy" & data$population=="3" & data$sex=="f" & data$focal.individual=="low")
monogamysubsetpop3femalehigh=subset(data, data$treatment=="monogamy" & data$population=="3" & data$sex=="f" & data$focal.individual=="high")

#Making Dataframes to store outputs from loop
datavialpop1male = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datavialpop2male = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datavialpop3male = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datavialpop4male = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datavialpop5male = data.frame(matrix(NA, nrow = 10000, ncol = 3))

datavialpop1female = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datavialpop2female = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datavialpop3female = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datavialpop4female = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datavialpop5female = data.frame(matrix(NA, nrow = 10000, ncol = 3))

#Cage
datacagepop1male = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datacagepop2male = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datacagepop3male = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datacagepop4male = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datacagepop5male = data.frame(matrix(NA, nrow = 10000, ncol = 3))

datacagepop1female = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datacagepop2female = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datacagepop3female = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datacagepop4female = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datacagepop5female = data.frame(matrix(NA, nrow = 10000, ncol = 3))

#monogamy
datamonogamypop1male  = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datamonogamypop5male = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datamonogamypop3male = data.frame(matrix(NA, nrow = 10000, ncol = 3))

datamonogamypop1female = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datamonogamypop5female = data.frame(matrix(NA, nrow = 10000, ncol = 3))
datamonogamypop3female = data.frame(matrix(NA, nrow = 10000, ncol = 3))


#Boostrap looping procedure
for (i in 1:10000){
  
  #vials
  vialpop1malesample=vialsubsetpop1male[sample(nrow(vialsubsetpop1male), replace = TRUE), ]
  vialpop5malesample=vialsubsetpop5male[sample(nrow(vialsubsetpop5male), replace = TRUE), ]
  vialpop3malesamplehigh=vialsubsetpop3malehigh[sample(nrow(vialsubsetpop3malehigh), replace = TRUE), ]
  vialpop3malesamplelow=vialsubsetpop3malelow[sample(nrow(vialsubsetpop3malelow), replace = TRUE), ]
  vialpop3malesample=rbind(vialpop3malesamplehigh,vialpop3malesamplelow)
  vialpop2malesamplehigh=vialsubsetpop2malehigh[sample(nrow(vialsubsetpop2malehigh), replace = TRUE), ]
  vialpop2malesamplelow=vialsubsetpop2malelow[sample(nrow(vialsubsetpop2malelow), replace = TRUE), ]
  vialpop2malesample=rbind(vialpop2malesamplehigh,vialpop2malesamplelow)
  vialpop4malesamplehigh=vialsubsetpop4malehigh[sample(nrow(vialsubsetpop4malehigh), replace = TRUE), ]
  vialpop4malesamplelow=vialsubsetpop4malelow[sample(nrow(vialsubsetpop4malelow), replace = TRUE), ]
  vialpop4malesample=rbind(vialpop4malesamplehigh,vialpop4malesamplelow)
  
  
  vialpop1femalesample=vialsubsetpop1female[sample(nrow(vialsubsetpop1female), replace = TRUE), ]
  vialpop5femalesample=vialsubsetpop5female[sample(nrow(vialsubsetpop5female), replace = TRUE), ]
  vialpop3femalesamplehigh=vialsubsetpop3femalehigh[sample(nrow(vialsubsetpop3femalehigh), replace = TRUE), ]
  vialpop3femalesamplelow=vialsubsetpop3femalelow[sample(nrow(vialsubsetpop3femalelow), replace = TRUE), ]
  vialpop3femalesample=rbind(vialpop3femalesamplehigh,vialpop3femalesamplelow)
  vialpop2femalesamplehigh=vialsubsetpop2femalehigh[sample(nrow(vialsubsetpop2femalehigh), replace = TRUE), ]
  vialpop2femalesamplelow=vialsubsetpop2femalelow[sample(nrow(vialsubsetpop2femalelow), replace = TRUE), ]
  vialpop2femalesample=rbind(vialpop2femalesamplehigh,vialpop2femalesamplelow)
  vialpop4femalesamplehigh=vialsubsetpop4femalehigh[sample(nrow(vialsubsetpop4femalehigh), replace = TRUE), ]
  vialpop4femalesamplelow=vialsubsetpop4femalelow[sample(nrow(vialsubsetpop4femalelow), replace = TRUE), ]
  vialpop4femalesample=rbind(vialpop4femalesamplehigh,vialpop4femalesamplelow)
  
  #cage
  
  cagepop1malesample=cagesubsetpop1male[sample(nrow(cagesubsetpop1male), replace = TRUE), ]
  cagepop5malesample=cagesubsetpop5male[sample(nrow(cagesubsetpop5male), replace = TRUE), ]
  cagepop3malesamplehigh=cagesubsetpop3malehigh[sample(nrow(cagesubsetpop3malehigh), replace = TRUE), ]
  cagepop3malesamplelow=cagesubsetpop3malelow[sample(nrow(cagesubsetpop3malelow), replace = TRUE), ]
  cagepop3malesample=rbind(cagepop3malesamplehigh,cagepop3malesamplelow)
  cagepop2malesamplehigh=cagesubsetpop2malehigh[sample(nrow(cagesubsetpop2malehigh), replace = TRUE), ]
  cagepop2malesamplelow=cagesubsetpop2malelow[sample(nrow(cagesubsetpop2malelow), replace = TRUE), ]
  cagepop2malesample=rbind(cagepop2malesamplehigh,cagepop2malesamplelow)
  cagepop4malesamplehigh=cagesubsetpop4malehigh[sample(nrow(cagesubsetpop4malehigh), replace = TRUE), ]
  cagepop4malesamplelow=cagesubsetpop4malelow[sample(nrow(cagesubsetpop4malelow), replace = TRUE), ]
  cagepop4malesample=rbind(cagepop4malesamplehigh,cagepop4malesamplelow)
  
  cagepop1femalesample=cagesubsetpop1female[sample(nrow(cagesubsetpop1female), replace = TRUE), ]
  cagepop5femalesample=cagesubsetpop5female[sample(nrow(cagesubsetpop5female), replace = TRUE), ]
  cagepop3femalesamplehigh=cagesubsetpop3femalehigh[sample(nrow(cagesubsetpop3femalehigh), replace = TRUE), ]
  cagepop3femalesamplelow=cagesubsetpop3femalelow[sample(nrow(cagesubsetpop3femalelow), replace = TRUE), ]
  cagepop3femalesample=rbind(cagepop3femalesamplehigh,cagepop3femalesamplelow)
  cagepop2femalesamplehigh=cagesubsetpop2femalehigh[sample(nrow(cagesubsetpop2femalehigh), replace = TRUE), ]
  cagepop2femalesamplelow=cagesubsetpop2femalelow[sample(nrow(cagesubsetpop2femalelow), replace = TRUE), ]
  cagepop2femalesample=rbind(cagepop2femalesamplehigh,cagepop2femalesamplelow)
  cagepop4femalesamplehigh=cagesubsetpop4femalehigh[sample(nrow(cagesubsetpop4femalehigh), replace = TRUE), ]
  cagepop4femalesamplelow=cagesubsetpop4femalelow[sample(nrow(cagesubsetpop4femalelow), replace = TRUE), ]
  cagepop4femalesample=rbind(cagepop4femalesamplehigh,cagepop4femalesamplelow)
  
  #Monogamy
  
  monogamypop1malesample=monogamysubsetpop1male[sample(nrow(monogamysubsetpop1male), replace = TRUE), ]
  monogamypop1femalesample=monogamysubsetpop1female[sample(nrow(monogamysubsetpop1female), replace = TRUE), ]
  monogamypop5malesample=monogamysubsetpop5male[sample(nrow(monogamysubsetpop5male), replace = TRUE), ]
  monogamypop5femalesample=monogamysubsetpop5female[sample(nrow(monogamysubsetpop5female), replace = TRUE), ]
  monogamypop3malesamplelow=monogamysubsetpop3malelow[sample(nrow(monogamysubsetpop3malelow), replace = TRUE), ]
  monogamypop3malesamplehigh=monogamysubsetpop3malehigh[sample(nrow(monogamysubsetpop3malehigh), replace = TRUE), ]
  monogamypop3malesample=rbind(monogamypop3malesamplelow,monogamypop3malesamplehigh)
  monogamypop3femalesamplelow=monogamysubsetpop3femalelow[sample(nrow(monogamysubsetpop3femalelow), replace = TRUE), ]
  monogamypop3femalesamplehigh=monogamysubsetpop3femalehigh[sample(nrow(monogamysubsetpop3femalehigh), replace = TRUE), ]
  monogamypop3femalesample=rbind(monogamypop3femalesamplelow,monogamypop3femalesamplehigh)
  #   
  
  ##DATA FRAMES
  #vial
  datavialpop1male [i,] = c(i,mean(as.numeric(as.character(vialpop1malesample$n.wt))), var(as.numeric(as.character(vialpop1malesample$n.wt))))
  datavialpop2male [i,] = c(i,mean(as.numeric(as.character(vialpop2malesample$n.wt))), var(as.numeric(as.character(vialpop2malesample$n.wt))))
  datavialpop3male [i,] = c(i,mean(as.numeric(as.character(vialpop3malesample$n.wt))), var(as.numeric(as.character(vialpop3malesample$n.wt))))
  datavialpop4male [i,] = c(i,mean(as.numeric(as.character(vialpop4malesample$n.wt))), var(as.numeric(as.character(vialpop4malesample$n.wt))))
  datavialpop5male [i,] = c(i,mean(as.numeric(as.character(vialpop5malesample$n.wt))), var(as.numeric(as.character(vialpop5malesample$n.wt))))
  
  datavialpop1female [i,] = c(i,mean(as.numeric(as.character(vialpop1femalesample$n.wt))), var(as.numeric(as.character(vialpop1femalesample$n.wt))))
  datavialpop2female [i,] = c(i,mean(as.numeric(as.character(vialpop2femalesample$n.wt))), var(as.numeric(as.character(vialpop2femalesample$n.wt))))
  datavialpop3female [i,] = c(i,mean(as.numeric(as.character(vialpop3femalesample$n.wt))), var(as.numeric(as.character(vialpop3femalesample$n.wt))))
  datavialpop4female [i,] = c(i,mean(as.numeric(as.character(vialpop4femalesample$n.wt))), var(as.numeric(as.character(vialpop4femalesample$n.wt))))
  datavialpop5female [i,] = c(i,mean(as.numeric(as.character(vialpop5femalesample$n.wt))), var(as.numeric(as.character(vialpop5femalesample$n.wt))))
  
  #Cage
  datacagepop1male [i,] = c(i,mean(as.numeric(as.character(cagepop1malesample$n.wt))), var(as.numeric(as.character(cagepop1malesample$n.wt))))
  datacagepop2male [i,] = c(i,mean(as.numeric(as.character(cagepop2malesample$n.wt))), var(as.numeric(as.character(cagepop2malesample$n.wt))))
  datacagepop3male [i,] = c(i,mean(as.numeric(as.character(cagepop3malesample$n.wt))), var(as.numeric(as.character(cagepop3malesample$n.wt))))
  datacagepop4male [i,] = c(i,mean(as.numeric(as.character(cagepop4malesample$n.wt))), var(as.numeric(as.character(cagepop4malesample$n.wt))))
  datacagepop5male [i,] = c(i,mean(as.numeric(as.character(cagepop5malesample$n.wt))), var(as.numeric(as.character(cagepop5malesample$n.wt))))
  
  datacagepop1female [i,] = c(i,mean(as.numeric(as.character(cagepop1femalesample$n.wt))), var(as.numeric(as.character(cagepop1femalesample$n.wt))))
  datacagepop2female [i,] = c(i,mean(as.numeric(as.character(cagepop2femalesample$n.wt))), var(as.numeric(as.character(cagepop2femalesample$n.wt))))
  datacagepop3female [i,] = c(i,mean(as.numeric(as.character(cagepop3femalesample$n.wt))), var(as.numeric(as.character(cagepop3femalesample$n.wt))))
  datacagepop4female [i,] = c(i,mean(as.numeric(as.character(cagepop4femalesample$n.wt))), var(as.numeric(as.character(cagepop4femalesample$n.wt))))
  datacagepop5female [i,] = c(i,mean(as.numeric(as.character(cagepop5femalesample$n.wt))), var(as.numeric(as.character(cagepop5femalesample$n.wt))))
   
  #monogamy
  datamonogamypop1male [i,] = c(i,mean(as.numeric(as.character(monogamypop1malesample$n.wt))), var(as.numeric(as.character(monogamypop1malesample$n.wt))))
  datamonogamypop3male [i,] = c(i,mean(as.numeric(as.character(monogamypop3malesample$n.wt))), var(as.numeric(as.character(monogamypop3malesample$n.wt))))
  datamonogamypop5male [i,] = c(i,mean(as.numeric(as.character(monogamypop5malesample$n.wt))), var(as.numeric(as.character(monogamypop5malesample$n.wt))))
  
  datamonogamypop1female [i,] = c(i,mean(as.numeric(as.character(monogamypop1femalesample$n.wt))), var(as.numeric(as.character(monogamypop1femalesample$n.wt))))
  datamonogamypop3female [i,] = c(i,mean(as.numeric(as.character(monogamypop3femalesample$n.wt))), var(as.numeric(as.character(monogamypop3femalesample$n.wt))))
  datamonogamypop5female [i,] = c(i,mean(as.numeric(as.character(monogamypop5femalesample$n.wt))), var(as.numeric(as.character(monogamypop5femalesample$n.wt))))
  # 
  print(i)
}

# Make a dataframe to store bootstrapped data
vial=rep("vial", 10000)
cage=rep("cage", 10000)
monogamy=rep("monogamy", 10000)

pop1=rep(1, 10000)
pop2=rep(2, 10000)
pop3=rep(3, 10000)
pop4=rep(4, 10000)
pop5=rep(5, 10000)

male=rep("male", 10000)
female=rep("female", 10000)

colnames = c("treatment","population.type","sex","replicate","mean.wt","var.wt")
vialpop1male=cbind(vial,pop1,male, datavialpop1male)
vialpop2male=cbind(vial, pop2, male, datavialpop2male)
vialpop3male=cbind(vial, pop3, male,datavialpop3male)
vialpop4male=cbind(vial, pop4, male, datavialpop4male)
vialpop5male=cbind(vial, pop5, male, datavialpop5male)

vialpop1female=cbind(vial,pop1,female,datavialpop1female)
vialpop2female=cbind(vial, pop2,female,datavialpop2female)
vialpop3female=cbind(vial,pop3,female,datavialpop3female)
vialpop4female=cbind(vial, pop4,female,datavialpop4female)
vialpop5female=cbind(vial, pop5,female,datavialpop5female)

cagepop1male=cbind(cage,pop1,male, datacagepop1male)
cagepop2male=cbind(cage, pop2, male, datacagepop2male)
cagepop3male=cbind(cage,pop3, male,datacagepop3male)
cagepop4male=cbind(cage, pop4, male, datacagepop4male)
cagepop5male=cbind(cage, pop5, male, datacagepop5male)

cagepop1female=cbind(cage,pop1,female,datacagepop1female)
cagepop2female=cbind(cage, pop2,female,datacagepop2female)
cagepop3female=cbind(cage,pop3,female,datacagepop3female)
cagepop4female=cbind(cage, pop4,female,datacagepop4female)
cagepop5female=cbind(cage, pop5,female,datacagepop5female)

monogamypop1male=cbind(monogamy,pop1,male, datamonogamypop1male)
monogamypop3male=cbind(monogamy,pop3, male,datamonogamypop3male)
monogamypop5male=cbind(monogamy, pop5, male, datamonogamypop5male)

monogamypop1female=cbind(monogamy,pop1,female,datamonogamypop1female)
monogamypop3female=cbind(monogamy,pop3,female,datamonogamypop3female)
monogamypop5female=cbind(monogamy, pop5,female,datamonogamypop5female)

colnames(vialpop1male)=colnames
colnames(vialpop2male)=colnames
colnames(vialpop3male)=colnames
colnames(vialpop4male)=colnames
colnames(vialpop5male)=colnames
colnames(vialpop1female)=colnames
colnames(vialpop2female)=colnames
colnames(vialpop3female)=colnames
colnames(vialpop4female)=colnames
colnames(vialpop5female)=colnames
colnames(cagepop1male)=colnames
colnames(cagepop2male)=colnames
colnames(cagepop3male)=colnames
colnames(cagepop4male)=colnames
colnames(cagepop5male)=colnames
colnames(cagepop1female)=colnames
colnames(cagepop2female)=colnames
colnames(cagepop3female)=colnames
colnames(cagepop4female)=colnames
colnames(cagepop5female)=colnames
colnames(monogamypop1male)=colnames
colnames(monogamypop3male)=colnames
colnames(monogamypop5male)=colnames
colnames(monogamypop1female)=colnames
colnames(monogamypop3female)=colnames
colnames(monogamypop5female)=colnames

# Join all dataframes into a single bootstrapped data frame 
bootstrapped.data=rbind(vialpop1male,vialpop2male,vialpop3male,vialpop4male,vialpop5male, 
               vialpop1female,vialpop2female,vialpop3female,vialpop4female,vialpop5female,
               cagepop1male,cagepop2male,cagepop3male,cagepop4male,cagepop5male, 
               cagepop1female,cagepop2female,cagepop3female,cagepop4female,cagepop5female,
               monogamypop1male,monogamypop3male,monogamypop5male, 
               monogamypop1female,monogamypop3female,monogamypop5female)

# Calculate the adjusted variance for each variance values 
bootstrapped.data$adjusted.variance = (2 * bootstrapped.data$var.wt) / bootstrapped.data$mean.wt  

# Write the bootstrapped data file out to disk
write.csv(bootstrapped.data, file="/FILE/PATH/TO/BootstrappedFitnessDataFinalJuly2021.csv")


#
