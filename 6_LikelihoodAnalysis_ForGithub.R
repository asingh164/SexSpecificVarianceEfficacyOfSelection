###########################################################################################################################################################
###########################################################################################################################################################
###   Script Details:                                                                                                                                   ###
###   Selection Analysis for Singh, A. and A.F., Agrawal. 2021. Sex-Specific Variance in Fitness and the Efficacy of Selection                          ###
###   Script written by Aneil Agrawal (a.agrawal[at]utoronto.ca)                                                                                        ###
###   This script will estimate selection differentials acting on condition in our labratory Drosophila experiment and generate corresponding plots     ###
###########################################################################################################################################################
###########################################################################################################################################################


# Load R packages 
library(extraDistr)
library(RCurl)

# Load in data from Github
data.github = getURL("https://raw.githubusercontent.com/asingh164/SexSpecificVarianceEfficacyOfSelection/master/FitnessDataNeprojectFinalSept2020.csv")
data.all=read.csv(text = data.github)

data.all$treatment = as.factor(data.all$treatment)
data.all$population.type = as.factor(data.all$population.type)
data.all$sex = as.factor(data.all$sex)

data.all<-data.all[data.all$n.total == 32,] # Perform the same filtering as the real data (i.e., remove any replicate with fewer than 32 offspring)

data.Vial.PopA.f<- as.numeric(subset(data.all, treatment == "vial" & population.type == "1" & sex == "f" & !is.na(n.total), select= n.wt)[,1])
data.Vial.PopB.f<- as.numeric(subset(data.all, treatment == "vial" & population.type == "2" & sex == "f" & !is.na(n.total), select= n.wt)[,1])
data.Vial.PopC.f<- as.numeric(subset(data.all, treatment == "vial" & population.type == "3" & sex == "f" & !is.na(n.total), select= n.wt)[,1])
data.Vial.PopD.f<- as.numeric(subset(data.all, treatment == "vial" & population.type == "4" & sex == "f" & !is.na(n.total), select= n.wt)[,1])
data.Vial.PopE.f<- as.numeric(subset(data.all, treatment == "vial" & population.type == "5" & sex == "f" & !is.na(n.total), select= n.wt)[,1])
data.Vial.PopA.m<- as.numeric(subset(data.all, treatment == "vial" & population.type == "1" & sex == "m" & !is.na(n.total), select= n.wt)[,1])
data.Vial.PopB.m<- as.numeric(subset(data.all, treatment == "vial" & population.type == "2" & sex == "m" & !is.na(n.total), select= n.wt)[,1])
data.Vial.PopC.m<- as.numeric(subset(data.all, treatment == "vial" & population.type == "3" & sex == "m" & !is.na(n.total), select= n.wt)[,1])
data.Vial.PopD.m<- as.numeric(subset(data.all, treatment == "vial" & population.type == "4" & sex == "m" & !is.na(n.total), select= n.wt)[,1])
data.Vial.PopE.m<- as.numeric(subset(data.all, treatment == "vial" & population.type == "5" & sex == "m" & !is.na(n.total), select= n.wt)[,1])

data.Cage.PopA.f<- as.numeric(subset(data.all, treatment == "cage" & population.type == "1" & sex == "f" & !is.na(n.total), select= n.wt)[,1])
data.Cage.PopB.f<- as.numeric(subset(data.all, treatment == "cage" & population.type == "2" & sex == "f" & !is.na(n.total), select= n.wt)[,1])
data.Cage.PopC.f<- as.numeric(subset(data.all, treatment == "cage" & population.type == "3" & sex == "f" & !is.na(n.total), select= n.wt)[,1])
data.Cage.PopD.f<- as.numeric(subset(data.all, treatment == "cage" & population.type == "4" & sex == "f" & !is.na(n.total), select= n.wt)[,1])
data.Cage.PopE.f<- as.numeric(subset(data.all, treatment == "cage" & population.type == "5" & sex == "f" & !is.na(n.total), select= n.wt)[,1])
data.Cage.PopA.m<- as.numeric(subset(data.all, treatment == "cage" & population.type == "1" & sex == "m" & !is.na(n.total), select= n.wt)[,1])
data.Cage.PopB.m<- as.numeric(subset(data.all, treatment == "cage" & population.type == "2" & sex == "m" & !is.na(n.total), select= n.wt)[,1])
data.Cage.PopC.m<- as.numeric(subset(data.all, treatment == "cage" & population.type == "3" & sex == "m" & !is.na(n.total), select= n.wt)[,1])
data.Cage.PopD.m<- as.numeric(subset(data.all, treatment == "cage" & population.type == "4" & sex == "m" & !is.na(n.total), select= n.wt)[,1])
data.Cage.PopE.m<- as.numeric(subset(data.all, treatment == "cage" & population.type == "5" & sex == "m" & !is.na(n.total), select= n.wt)[,1])

data.Monog.PopA.f<- as.numeric(subset(data.all, treatment == "monogamy" & population.type == "1" & sex == "f" & !is.na(n.total), select= n.wt)[,1])
data.Monog.PopC.f<- as.numeric(subset(data.all, treatment == "monogamy" & population.type == "3" & sex == "f" & !is.na(n.total), select= n.wt)[,1])
data.Monog.PopE.f<- as.numeric(subset(data.all, treatment == "monogamy" & population.type == "5" & sex == "f" & !is.na(n.total), select= n.wt)[,1])
data.Monog.PopA.m<- as.numeric(subset(data.all, treatment == "monogamy" & population.type == "1" & sex == "m" & !is.na(n.total), select= n.wt)[,1])
data.Monog.PopC.m<- as.numeric(subset(data.all, treatment == "monogamy" & population.type == "3" & sex == "m" & !is.na(n.total), select= n.wt)[,1])
data.Monog.PopE.m<- as.numeric(subset(data.all, treatment == "monogamy" & population.type == "5" & sex == "m" & !is.na(n.total), select= n.wt)[,1])


### Just to look at each of the distributions
getTheBasics<-function(data) {hist(data, main = paste("mu = ", round(100*mean(data))/100, ";  var = ", round(100*var(data))/100)); cat(c(mean(data), var(data)))}
getTheBasics(data.Vial.PopA.f)
getTheBasics(data.Vial.PopB.f)
getTheBasics(data.Vial.PopC.f)
getTheBasics(data.Vial.PopD.f)
getTheBasics(data.Vial.PopE.f)

getTheBasics(data.Vial.PopA.m)
getTheBasics(data.Vial.PopB.m)
getTheBasics(data.Vial.PopC.m)
getTheBasics(data.Vial.PopD.m)
getTheBasics(data.Vial.PopE.m)


List.Vial.f = list(data.Vial.PopA.f,data.Vial.PopB.f, data.Vial.PopC.f, data.Vial.PopD.f, data.Vial.PopE.f)
List.Vial.m = list(data.Vial.PopA.m,data.Vial.PopB.m, data.Vial.PopC.m, data.Vial.PopD.m, data.Vial.PopE.m)
List.Cage.f = list(data.Cage.PopA.f,data.Cage.PopB.f, data.Cage.PopC.f, data.Cage.PopD.f, data.Cage.PopE.f)
List.Cage.m = list(data.Cage.PopA.m,data.Cage.PopB.m, data.Cage.PopC.m, data.Cage.PopD.m, data.Cage.PopE.m)
List.Monog.f = list(data.Monog.PopA.f,data.Monog.PopC.f, data.Monog.PopE.f)
List.Monog.m = list(data.Monog.PopA.m, data.Monog.PopC.m, data.Monog.PopE.m)


### The function below is used by the function after it and output is explained there
getTheBasicsForAllPopTypes<-function(listOfDataVecs){
  m = matrix(nrow = 3, ncol = length(listOfDataVecs))
  for(ii in 1:length(listOfDataVecs)){
    m[1,ii] = mean(sapply(listOfDataVecs[ii], as.numeric))
    m[2,ii] = var(sapply(listOfDataVecs[ii], as.numeric))
    m[3,ii] = m[2,ii]*(2/m[1,ii])
  }
  return(m)
}



### The function below is used to give information for each population type of a given mating regime
### The colunmns are the different population types
### The rows are
### 1-3) female mean, obs var, adj var
### 4-6) male mean, obs var, adj var
### 7) Ne(mf)/N  (i.e., using male- and female-specific variances)
### 8) Ne(ff)/N (i.e., setting the male variance to be the same as the female variance)
### 9) Ne(mf)/Ne(ff)
### 10) r calculated assuming alpha = 1.5

getTheBasicsForAllPopTypes.BothSexes<-function(listOfDataVecs.f, listOfDataVecs.m){
  m.f<- getTheBasicsForAllPopTypes(listOfDataVecs.f)
  m.m<- getTheBasicsForAllPopTypes(listOfDataVecs.m)
  m<-rbind(m.f, m.m, matrix(nrow = 4, ncol = length(listOfDataVecs.f)))
  for(ii in 1:length(listOfDataVecs.f)){
    NeOverN.mf<- 8/(4 + m.f[3,ii] + m.m[3,ii])
    NeOverN.ff<- 8/(4 + 2*m.f[3,ii])
    Ne.mf.OVER.Ne.ff<- NeOverN.mf/NeOverN.ff
    rAssumingAlpha1.5 <- (2+m.f[3,ii])*(1 + 1.5)/(4 + m.f[3,ii] + m.m[3,ii])
    m[7,ii] = NeOverN.mf
    m[8,ii] = NeOverN.ff
    m[9,ii] = Ne.mf.OVER.Ne.ff
    m[10,ii] = rAssumingAlpha1.5
  }
  return(m)
}

getTheBasicsForAllPopTypes.BothSexes(List.Vial.f, List.Vial.m)
getTheBasicsForAllPopTypes.BothSexes(List.Cage.f, List.Cage.m)
getTheBasicsForAllPopTypes.BothSexes(List.Monog.f, List.Monog.m)



###########################################################
############# Bootstrapping section #################
###########################################################
getOneBootstrapOneSexOnePop<-function(dt){
  resampled.dt<-sample(dt, size=length(dt), replace = TRUE)
  mu <- mean(resampled.dt)
  v<- var(resampled.dt)
  adjV<- v*(2/mu)
  return(c(mu, v, adjV))
}

getOneBootstrapBothSexesOnePop<-function(dt.f, dt.m){
  summary.f<-getOneBootstrapOneSexOnePop(dt.f)
  summary.m<-getOneBootstrapOneSexOnePop(dt.m)
  SexDiffAdjVar<- summary.m[3] - summary.f[3]
  NeOverN.mf<- 8/(4 + summary.f[3] + summary.m[3])
  NeOverN.ff<- 8/(4 + 2*summary.f[3])
  Ne.mf.OVER.Ne.ff<- NeOverN.mf/NeOverN.ff
  rAssumingAlpha1.5 <- (2+summary.f[3])*(1 + 1.5)/(4 + summary.f[3] + summary.m[3])
  return(c(summary.f, summary.m, SexDiffAdjVar, NeOverN.mf, NeOverN.ff, Ne.mf.OVER.Ne.ff, rAssumingAlpha1.5))
}

getBootstrapDistOneSexOnePop<-function(dt, nBoot = 10000){
  results<-matrix(nrow = nBoot, ncol = length(getOneBootstrapOneSexOnePop(dt)))
  for(i in 1:nBoot) results[i,] <-getOneBootstrapOneSexOnePop(dt)
  return(results)
}

getBootstrapDistBothSexesOnePop<-function(dt.f, dt.m, nBoot = 10000){
  results<-matrix(nrow = nBoot, ncol = length(getOneBootstrapBothSexesOnePop(dt.f, dt.m)))
  for(i in 1:nBoot) results[i,] <-getOneBootstrapBothSexesOnePop(dt.f, dt.m)
  return(results)
}


#### Here I generate bootstrap distributions for each population type for both sexes
#### Distributions are given for 11 variables represented by the columns
### 1-3) female mean, obs var, adj var
### 4-6) male mean, obs var, adj var
### 7) sex diff in (adj) var
### 8) Ne(mf)/N
### 9) Ne(ff)/N
### 10) Ne(mf)/Ne(ff)
### 11) r assuming alpha = 1.5
BootstrapDist.Vial.PopA<-getBootstrapDistBothSexesOnePop(data.Vial.PopA.f, data.Vial.PopA.m)
BootstrapDist.Vial.PopB<-getBootstrapDistBothSexesOnePop(data.Vial.PopB.f, data.Vial.PopB.m)
BootstrapDist.Vial.PopC<-getBootstrapDistBothSexesOnePop(data.Vial.PopC.f, data.Vial.PopC.m)
BootstrapDist.Vial.PopD<-getBootstrapDistBothSexesOnePop(data.Vial.PopD.f, data.Vial.PopD.m)
BootstrapDist.Vial.PopE<-getBootstrapDistBothSexesOnePop(data.Vial.PopE.f, data.Vial.PopE.m)

BootstrapDist.Cage.PopA<-getBootstrapDistBothSexesOnePop(data.Cage.PopA.f, data.Cage.PopA.m)
BootstrapDist.Cage.PopB<-getBootstrapDistBothSexesOnePop(data.Cage.PopB.f, data.Cage.PopB.m)
BootstrapDist.Cage.PopC<-getBootstrapDistBothSexesOnePop(data.Cage.PopC.f, data.Cage.PopC.m)
BootstrapDist.Cage.PopD<-getBootstrapDistBothSexesOnePop(data.Cage.PopD.f, data.Cage.PopD.m)
BootstrapDist.Cage.PopE<-getBootstrapDistBothSexesOnePop(data.Cage.PopE.f, data.Cage.PopE.m)

BootstrapDist.Monog.PopA<-getBootstrapDistBothSexesOnePop(data.Monog.PopA.f, data.Monog.PopA.m)
BootstrapDist.Monog.PopC<-getBootstrapDistBothSexesOnePop(data.Monog.PopC.f, data.Monog.PopC.m)
BootstrapDist.Monog.PopE<-getBootstrapDistBothSexesOnePop(data.Monog.PopE.f, data.Monog.PopE.m)


getBootstrap95CI<-function(BSdist.dt){
  n<-nrow(BSdist.dt)
  nc<-ncol(BSdist.dt)
  bs.results<-matrix(ncol = 3, nrow = nc)
  for(i in 1:nc){
    sortedVec<- sort(BSdist.dt[,i])
    bs.results[i,]<-c(mean(sortedVec), sortedVec[round(0.025*n)], sortedVec[round(0.975*n)])
  }
  return(bs.results)
}


#### Here I get 95% confidence intervals from the bootstrap distributions
#### for 11 variables represented by the columns
### 1-3) female mean, obs var, adj var
### 4-6) male mean, obs var, adj var
### 7) sex diff in (adj) var
### 8) Ne(mf)/N
### 9) Ne(ff)/N
### 10) Ne(mf)/Ne(ff)
### 11) r assuming alpha = 1.5
getBootstrap95CI(BootstrapDist.Vial.PopA)
getBootstrap95CI(BootstrapDist.Vial.PopB)
getBootstrap95CI(BootstrapDist.Vial.PopC)
getBootstrap95CI(BootstrapDist.Vial.PopD)
getBootstrap95CI(BootstrapDist.Vial.PopE)

getBootstrap95CI(BootstrapDist.Cage.PopA)
getBootstrap95CI(BootstrapDist.Cage.PopB)
getBootstrap95CI(BootstrapDist.Cage.PopC)
getBootstrap95CI(BootstrapDist.Cage.PopD)
getBootstrap95CI(BootstrapDist.Cage.PopE)

getBootstrap95CI(BootstrapDist.Monog.PopA)
getBootstrap95CI(BootstrapDist.Monog.PopC)
getBootstrap95CI(BootstrapDist.Monog.PopE)

#####################################################
######### Bootstrap test of heterogeneity effect
### Using only C vs (A & E)
## done by sex and then using sex-averages
#####################################################

## Here are the point estimates
point.estimate.deltaAdjVar.Het1.Vial.f<-getTheBasicsForAllPopTypes(List.Vial.f)[3,3] - (getTheBasicsForAllPopTypes(List.Vial.f)[3,1]+getTheBasicsForAllPopTypes(List.Vial.f)[3,5])/2
point.estimate.deltaAdjVar.Het1.Vial.m<-getTheBasicsForAllPopTypes(List.Vial.m)[3,3] - (getTheBasicsForAllPopTypes(List.Vial.m)[3,1]+getTheBasicsForAllPopTypes(List.Vial.m)[3,5])/2
point.estimate.deltaAdjVar.Het1.Vial.sexAvg<-(point.estimate.deltaAdjVar.Het1.Vial.f + point.estimate.deltaAdjVar.Het1.Vial.m)/2
point.estimate.deltaAdjVar.Het1.Vial.mfDIFF <- point.estimate.deltaAdjVar.Het1.Vial.m - point.estimate.deltaAdjVar.Het1.Vial.f

point.estimate.deltaAdjVar.Het1.Cage.f<-getTheBasicsForAllPopTypes(List.Cage.f)[3,3] - (getTheBasicsForAllPopTypes(List.Cage.f)[3,1]+getTheBasicsForAllPopTypes(List.Cage.f)[3,5])/2
point.estimate.deltaAdjVar.Het1.Cage.m<-getTheBasicsForAllPopTypes(List.Cage.m)[3,3] - (getTheBasicsForAllPopTypes(List.Cage.m)[3,1]+getTheBasicsForAllPopTypes(List.Cage.m)[3,5])/2
point.estimate.deltaAdjVar.Het1.Cage.sexAvg<-(point.estimate.deltaAdjVar.Het1.Cage.f + point.estimate.deltaAdjVar.Het1.Cage.m)/2
point.estimate.deltaAdjVar.Het1.Cage.mfDIFF <- point.estimate.deltaAdjVar.Het1.Cage.m - point.estimate.deltaAdjVar.Het1.Cage.f

point.estimate.deltaAdjVar.Het1.Monog.f<-getTheBasicsForAllPopTypes(List.Monog.f)[3,2] - (getTheBasicsForAllPopTypes(List.Monog.f)[3,1]+getTheBasicsForAllPopTypes(List.Monog.f)[3,3])/2
point.estimate.deltaAdjVar.Het1.Monog.m<-getTheBasicsForAllPopTypes(List.Monog.m)[3,2] - (getTheBasicsForAllPopTypes(List.Monog.m)[3,1]+getTheBasicsForAllPopTypes(List.Monog.m)[3,3])/2
point.estimate.deltaAdjVar.Het1.Monog.sexAvg<-(point.estimate.deltaAdjVar.Het1.Monog.f + point.estimate.deltaAdjVar.Het1.Monog.m)/2
point.estimate.deltaAdjVar.Het1.Monog.mfDIFF <- point.estimate.deltaAdjVar.Het1.Monog.m - point.estimate.deltaAdjVar.Het1.Monog.f


## get bootstrap distributions
deltaAdjVar.Het1.Vial.f<-BootstrapDist.Vial.PopC[,3] - (BootstrapDist.Vial.PopA[,3]+BootstrapDist.Vial.PopE[,3])/2
deltaAdjVar.Het1.Vial.m<-BootstrapDist.Vial.PopC[,6] - (BootstrapDist.Vial.PopA[,6]+BootstrapDist.Vial.PopE[,6])/2
deltaAdjVar.Het1.Cage.f<-BootstrapDist.Cage.PopC[,3] - (BootstrapDist.Cage.PopA[,3]+BootstrapDist.Cage.PopE[,3])/2
deltaAdjVar.Het1.Cage.m<-BootstrapDist.Cage.PopC[,6] - (BootstrapDist.Cage.PopA[,6]+BootstrapDist.Cage.PopE[,6])/2 
deltaAdjVar.Het1.Monog.f<-BootstrapDist.Monog.PopC[,3] - (BootstrapDist.Monog.PopA[,3]+BootstrapDist.Monog.PopE[,3])/2
deltaAdjVar.Het1.Monog.m<-BootstrapDist.Monog.PopC[,6] - (BootstrapDist.Monog.PopA[,6]+BootstrapDist.Monog.PopE[,6])/2

deltaAdjVar.Het1.Vial.sexavg<-(deltaAdjVar.Het1.Vial.f + deltaAdjVar.Het1.Vial.m)/2
deltaAdjVar.Het1.Cage.sexavg<-(deltaAdjVar.Het1.Cage.f + deltaAdjVar.Het1.Cage.m)/2
deltaAdjVar.Het1.Monog.sexavg<-(deltaAdjVar.Het1.Monog.f + deltaAdjVar.Het1.Monog.m)/2

deltaAdjVar.Het1.Vial.mfDiff<-deltaAdjVar.Het1.Vial.m - deltaAdjVar.Het1.Vial.f
deltaAdjVar.Het1.Cage.mfDiff<-deltaAdjVar.Het1.Cage.m - deltaAdjVar.Het1.Cage.f
deltaAdjVar.Het1.Monog.mfDiff<-deltaAdjVar.Het1.Monog.m - deltaAdjVar.Het1.Monog.f


## 95% CIs (note display is in different orientation relative to point estimate outputs)
c(point.estimate.deltaAdjVar.Het1.Vial.f, point.estimate.deltaAdjVar.Het1.Cage.f, point.estimate.deltaAdjVar.Het1.Monog.f)
getBootstrap95CI(cbind(deltaAdjVar.Het1.Vial.f, deltaAdjVar.Het1.Cage.f, deltaAdjVar.Het1.Monog.f))

c(point.estimate.deltaAdjVar.Het1.Vial.m, point.estimate.deltaAdjVar.Het1.Cage.m, point.estimate.deltaAdjVar.Het1.Monog.m)
getBootstrap95CI(cbind(deltaAdjVar.Het1.Vial.m, deltaAdjVar.Het1.Cage.m, deltaAdjVar.Het1.Monog.m))

c(point.estimate.deltaAdjVar.Het1.Vial.sexAvg, point.estimate.deltaAdjVar.Het1.Cage.sexAvg, point.estimate.deltaAdjVar.Het1.Monog.sexAvg)
getBootstrap95CI(cbind(deltaAdjVar.Het1.Vial.sexavg, deltaAdjVar.Het1.Cage.sexavg, deltaAdjVar.Het1.Monog.sexavg))

c(point.estimate.deltaAdjVar.Het1.Vial.mfDIFF, point.estimate.deltaAdjVar.Het1.Cage.mfDIFF, point.estimate.deltaAdjVar.Het1.Monog.mfDIFF)
getBootstrap95CI(cbind(deltaAdjVar.Het1.Vial.mfDiff, deltaAdjVar.Het1.Cage.mfDiff, deltaAdjVar.Het1.Monog.mfDiff))

#####################################################
######### Bootstrap test of mating regime effect
### Use averages across all population types
## done by sex and using sex-averages
#####################################################

## Here are the point estimates
point.estimate.AdjVar.Vial.f<-mean(getTheBasicsForAllPopTypes(List.Vial.f)[3,] )
point.estimate.AdjVar.Vial.m<-mean(getTheBasicsForAllPopTypes(List.Vial.m)[3,] )
point.estimate.AdjVar.Vial.sexavg<-(point.estimate.AdjVar.Vial.f + point.estimate.AdjVar.Vial.m)/2

point.estimate.AdjVar.Cage.f<-mean(getTheBasicsForAllPopTypes(List.Cage.f)[3,] )
point.estimate.AdjVar.Cage.m<-mean(getTheBasicsForAllPopTypes(List.Cage.m)[3,] )
point.estimate.AdjVar.Cage.sexavg<-(point.estimate.AdjVar.Cage.f + point.estimate.AdjVar.Cage.m)/2

point.estimate.AdjVar.Monog.f<-mean(getTheBasicsForAllPopTypes(List.Monog.f)[3,] )
point.estimate.AdjVar.Monog.m<-mean(getTheBasicsForAllPopTypes(List.Monog.m)[3,] )
point.estimate.AdjVar.Monog.sexavg<-(point.estimate.AdjVar.Monog.f + point.estimate.AdjVar.Monog.m)/2

point.estimate.DiffAdjVar.VialCage.f<-point.estimate.AdjVar.Vial.f - point.estimate.AdjVar.Cage.f
point.estimate.DiffAdjVar.VialCage.m<-point.estimate.AdjVar.Vial.m - point.estimate.AdjVar.Cage.m
point.estimate.DiffAdjVar.VialCage.sexavg<-(point.estimate.DiffAdjVar.VialCage.f + point.estimate.DiffAdjVar.VialCage.m)/2
c(point.estimate.DiffAdjVar.VialCage.f, point.estimate.DiffAdjVar.VialCage.m, point.estimate.DiffAdjVar.VialCage.sexavg)

point.estimate.DiffAdjVar.VialMonog.f<-point.estimate.AdjVar.Vial.f - point.estimate.AdjVar.Monog.f
point.estimate.DiffAdjVar.VialMonog.m<-point.estimate.AdjVar.Vial.m - point.estimate.AdjVar.Monog.m
point.estimate.DiffAdjVar.VialMonog.sexavg<-(point.estimate.DiffAdjVar.VialMonog.f + point.estimate.DiffAdjVar.VialMonog.m)/2
c(point.estimate.DiffAdjVar.VialMonog.f, point.estimate.DiffAdjVar.VialMonog.m, point.estimate.DiffAdjVar.VialMonog.sexavg)

point.estimate.DiffAdjVar.CageMonog.f<-point.estimate.AdjVar.Cage.f - point.estimate.AdjVar.Monog.f
point.estimate.DiffAdjVar.CageMonog.m<-point.estimate.AdjVar.Cage.m - point.estimate.AdjVar.Monog.m
point.estimate.DiffAdjVar.CageMonog.sexavg<-(point.estimate.DiffAdjVar.CageMonog.f + point.estimate.DiffAdjVar.CageMonog.m)/2
c(point.estimate.DiffAdjVar.CageMonog.f, point.estimate.DiffAdjVar.CageMonog.m, point.estimate.DiffAdjVar.CageMonog.sexavg)

## get bootstrap distributions
AdjVar.Vial.f<-(BootstrapDist.Vial.PopA[,3] + BootstrapDist.Vial.PopB[,3] + BootstrapDist.Vial.PopC[,3] + BootstrapDist.Vial.PopD[,3] + BootstrapDist.Vial.PopE[,3])/5 
AdjVar.Vial.m<-(BootstrapDist.Vial.PopA[,6] + BootstrapDist.Vial.PopB[,6] + BootstrapDist.Vial.PopC[,6] + BootstrapDist.Vial.PopD[,6] + BootstrapDist.Vial.PopE[,6])/5 
AdjVar.Vial.sexavg<-(AdjVar.Vial.f + AdjVar.Vial.m)/2

AdjVar.Cage.f<-(BootstrapDist.Cage.PopA[,3] + BootstrapDist.Cage.PopB[,3] + BootstrapDist.Cage.PopC[,3] + BootstrapDist.Cage.PopD[,3] + BootstrapDist.Cage.PopE[,3])/5 
AdjVar.Cage.m<-(BootstrapDist.Cage.PopA[,6] + BootstrapDist.Cage.PopB[,6] + BootstrapDist.Cage.PopC[,6] + BootstrapDist.Cage.PopD[,6] + BootstrapDist.Cage.PopE[,6])/5 
AdjVar.Cage.sexavg<-(AdjVar.Cage.f + AdjVar.Cage.m)/2

AdjVar.Monog.f<-(BootstrapDist.Monog.PopA[,3] + BootstrapDist.Monog.PopC[,3] + BootstrapDist.Monog.PopE[,3])/3 
AdjVar.Monog.m<-(BootstrapDist.Monog.PopA[,6] + BootstrapDist.Monog.PopC[,6] + BootstrapDist.Monog.PopE[,6])/3 
AdjVar.Monog.sexavg<-(AdjVar.Monog.f + AdjVar.Monog.m)/2

Diff.AdjVar.VialCage.f<-AdjVar.Vial.f - AdjVar.Cage.f
Diff.AdjVar.VialCage.m<-AdjVar.Vial.m - AdjVar.Cage.m
Diff.AdjVar.VialCage.sexavg<-AdjVar.Vial.sexavg - AdjVar.Cage.sexavg

Diff.AdjVar.VialMonog.f<-AdjVar.Vial.f - AdjVar.Monog.f
Diff.AdjVar.VialMonog.m<-AdjVar.Vial.m - AdjVar.Monog.m
Diff.AdjVar.VialMonog.sexavg<-AdjVar.Vial.sexavg - AdjVar.Monog.sexavg


Diff.AdjVar.CageMonog.f<-AdjVar.Cage.f - AdjVar.Monog.f
Diff.AdjVar.CageMonog.m<-AdjVar.Cage.m - AdjVar.Monog.m
Diff.AdjVar.CageMonog.sexavg<-AdjVar.Cage.sexavg - AdjVar.Monog.sexavg

c(point.estimate.DiffAdjVar.VialCage.f, point.estimate.DiffAdjVar.VialCage.m, point.estimate.DiffAdjVar.VialCage.sexavg)
getBootstrap95CI(cbind(Diff.AdjVar.VialCage.f, Diff.AdjVar.VialCage.m, Diff.AdjVar.VialCage.sexavg))

c(point.estimate.DiffAdjVar.VialMonog.f, point.estimate.DiffAdjVar.VialMonog.m, point.estimate.DiffAdjVar.VialMonog.sexavg)
getBootstrap95CI(cbind(Diff.AdjVar.VialMonog.f, Diff.AdjVar.VialMonog.m, Diff.AdjVar.VialMonog.sexavg))

c(point.estimate.DiffAdjVar.CageMonog.f, point.estimate.DiffAdjVar.CageMonog.m, point.estimate.DiffAdjVar.CageMonog.sexavg)
getBootstrap95CI(cbind(Diff.AdjVar.CageMonog.f, Diff.AdjVar.CageMonog.m, Diff.AdjVar.CageMonog.sexavg))


###########
##### sex differences

point.estimate.DiffAdjVar.mf.Vial<-point.estimate.AdjVar.Vial.m - point.estimate.AdjVar.Vial.f
point.estimate.DiffAdjVar.mf.Cage<-point.estimate.AdjVar.Cage.m - point.estimate.AdjVar.Cage.f
point.estimate.DiffAdjVar.mf.Monog<-point.estimate.AdjVar.Monog.m - point.estimate.AdjVar.Monog.f

c(point.estimate.DiffAdjVar.mf.Vial,point.estimate.DiffAdjVar.mf.Cage,point.estimate.DiffAdjVar.mf.Monog)

## get bootstrap distributions
Diff.AdjVar.mf.Vial<- AdjVar.Vial.m - AdjVar.Vial.f
Diff.AdjVar.mf.Cage<- AdjVar.Cage.m - AdjVar.Cage.f
Diff.AdjVar.mf.Monog<- AdjVar.Monog.m - AdjVar.Monog.f

getBootstrap95CI(cbind(Diff.AdjVar.mf.Vial, Diff.AdjVar.mf.Cage, Diff.AdjVar.mf.Monog))



###########
##### r, the ratio for the efficacy of selection (r = gamma.MF / gamma.FF)

point.estimate.r.Vial<- mean(getTheBasicsForAllPopTypes.BothSexes(List.Vial.f, List.Vial.m)[10,])
point.estimate.r.Cage<-mean(getTheBasicsForAllPopTypes.BothSexes(List.Cage.f, List.Cage.m)[10,])
point.estimate.r.Monog<-mean(getTheBasicsForAllPopTypes.BothSexes(List.Monog.f, List.Monog.m)[10,])

## get bootstrap distributions
r.Vial<-(BootstrapDist.Vial.PopA[,11] + BootstrapDist.Vial.PopB[,11] + BootstrapDist.Vial.PopC[,11] + BootstrapDist.Vial.PopD[,11] + BootstrapDist.Vial.PopE[,11])/5 
r.Cage<-(BootstrapDist.Cage.PopA[,11] + BootstrapDist.Cage.PopB[,11] + BootstrapDist.Cage.PopC[,11] + BootstrapDist.Cage.PopD[,11] + BootstrapDist.Cage.PopE[,11])/5 
r.Monog<-(BootstrapDist.Monog.PopA[,11] + BootstrapDist.Monog.PopC[,11] + BootstrapDist.Monog.PopE[,11])/3 

c(point.estimate.r.Vial, point.estimate.r.Cage, point.estimate.r.Monog)
getBootstrap95CI(cbind(r.Vial, r.Cage, r.Monog))



##########################################################
############ Likelihood analysis section #############################
##########################################################

## The Near-Zero Inflated Possion model
## A mixture of two Poisson distributions 
## With probability pNZ it comes from a Poisson distribution with a mean near zero, lambdaNZ = 0 (default embedded in likelihood function)
## With probability (1-pNZ) it comes from a Poisson distribution with mean lambdaFree
## This has mean mu = lambdaNZ + (1-pNZ)*(lambadFree - lambdaNZ)
## and variance v = lambdaNZ + (1-pNZ)*deltaLambda + pNZ*(1-pNZ)*deltaLambda^2 where deltaLambda = lambdaFree - lambdaNX
## so
## pNZ = (v - mu)/(v + mu*(mu - 1) + lambdaNZ*(lambdaNZ - 2 mu))
## lambdaFree = (v + mu*(mu - lambdaNZ - 1))/(mu - lambdaNZ)
## Note: we only use a simplified version with lambdaNZ = 0, i.e, the "zero-inflated poisson"
dnzip<-function(x, pNZ, lambdaFree, lambdaNZ, useLog = TRUE) pNZ*dpois(x, lambdaNZ) + (1-pNZ)*dpois(x, lambdaFree)


# This analysis is based on a 7 parameter model:
# VarF.Low is the variance of females in a low heterogeneity envirnoment IF competitors had same mean fitness
# psi is the ratio of male to female variance 
# gamma is the ratio of female variance in high vs low heterogeneity environments; 
# note this "gamma" is what we call "lambda" in the manuscript
# other 4 parameters are means of each sex in A and in E

NLLv4.OneSex<-function(distID, listOfDataVecs, parmVec){
  NegLog.NegBinomPDF<-function(data, mu, v) -sum(dnbinom(data, size = (mu^2)/(v - mu), prob = mu/v, log = TRUE))
  NegLog.BetaBinomPDF<-function(data, mu, v) -sum(dbbinom(data, size = 32, alpha = -mu*(mu^2 - 32*mu +v)/(mu^2 - 32*mu + 32*v), beta = (mu -32)*(mu^2 - 32*mu +v)/(mu^2 -32*mu +32*v), log = TRUE))
  NegLog.NearZeroInflatedPoisPDF<-function(data, mu, v, lambdaNZ = 0) -sum(log(dnzip(data, pNZ = (v - mu)/(v + mu*(mu - 1) + lambdaNZ*(lambdaNZ - 2*mu)), 
                                                                                       lambdaFree = (v + mu*(mu - lambdaNZ - 1))/(mu - lambdaNZ), lambdaNZ = lambdaNZ)))

  if(distID == 1) This.Dist =  NegLog.NegBinomPDF ## negative binomial
  if(distID == 2) This.Dist =  NegLog.BetaBinomPDF ## beta binomail
  if(distID == 3) This.Dist =  NegLog.NearZeroInflatedPoisPDF ## "near-zero" inflated poisson

  if(length(listOfDataVecs) == 3) {
    meanVec = c(parmVec[1], sum(parmVec[1:2])/2, parmVec[2])
    varVec = c(parmVec[3]*meanVec[1]/2, parmVec[4]*parmVec[3]*meanVec[2]/2, parmVec[3]*meanVec[3]/2)} 
  if(length(listOfDataVecs) == 5) {
    meanVec = c(parmVec[1], 0.75*parmVec[1] + 0.25*parmVec[2], sum(parmVec[1:2])/2, 0.25*parmVec[1] + 0.75*parmVec[2], parmVec[2])
    varVec = c(parmVec[3]*meanVec[1]/2, 0.5*(1+parmVec[4])*parmVec[3]*meanVec[2]/2, parmVec[4]*parmVec[3]*meanVec[1]/2, 0.5*(1+parmVec[4])*parmVec[3]*meanVec[4]/2,parmVec[3]*meanVec[3]/2)} 
  
  AllVarGreaterThanMeans <- !any(meanVec >=  varVec)
  if(!AllVarGreaterThanMeans)  ans <- 911911911911911
  
  if(AllVarGreaterThanMeans){
    ## function to calculate negative log likelihood based specified distribution
    ans = 0
    for(ii in 1:length(listOfDataVecs)) {
     x <- 9999999999
     try(x <- This.Dist(sapply(listOfDataVecs[ii], as.numeric), meanVec[ii], varVec[ii]) )
     ans = ans + x
    }
  }
  return(ans)
}

## This was a version I used just to explore likelihoods for a single sex (not used for results presented in ms)
## These models suggested that a betabinomial distribution tended to be best
getMLEv4.OneSex<-function(distID, listOfDataVecs){
  n<-50 ## number of starting places for optimization of likelihood
  randparm<- function() {
    needOne<- TRUE; countTries <-0; nMaxTries <-25;
    while(needOne & (countTries < nMaxTries)){
      z.meanVec<-runif(2, min = 1.5, max = 8)
      z.var<-runif(1, min = 2, max = 8)
      z.gamma<-runif(1, min = 0.9, max = 4)
      zParmVec = c(z.meanVec, z.var, z.gamma)
      
      if(length(listOfDataVecs) == 3) {
        z.meanVec = c(zParmVec[1], sum(zParmVec[1:2])/2, zParmVec[2])
        z.varVec = c(zParmVec[3]*z.meanVec[1]/2, zParmVec[4]*zParmVec[3]*z.meanVec[2]/2, zParmVec[3]*z.meanVec[3]/2)} 
      if(length(listOfDataVecs) == 5) {
        z.meanVec = c(zParmVec[1], 0.75*zParmVec[1] + 0.25*zParmVec[2], sum(zParmVec[1:2])/2, 0.25*zParmVec[1] + 0.75*zParmVec[2], zParmVec[2])
        z.varVec = c(zParmVec[3]*z.meanVec[1]/2, 0.5*(1+zParmVec[4])*zParmVec[3]*z.meanVec[2]/2, zParmVec[4]*zParmVec[3]*z.meanVec[1]/2, 0.5*(1+zParmVec[4])*zParmVec[3]*z.meanVec[4]/2,zParmVec[3]*z.meanVec[3]/2)} 
      
      z.AllVarGreaterThanMeans <- !any(z.meanVec >= z.varVec)
      if(z.AllVarGreaterThanMeans) needOne <- FALSE
      countTries<-countTries+1
      if(countTries == nMaxTries) zParmVec<-c(runif(2, min = 2, max = 2.5),  runif(1, min = 2.6, max = 3.1), runif(1, min = 0.9, max = 1.1))
    }
    return(zParmVec)
  }
  TheNLLfn<-function(parmVec) NLLv4.OneSex(distID, listOfDataVecs, parmVec)
  results<-matrix(nrow = n, ncol = length(randparm()) + 1)
  for(i in 1:n){
    oneMLE = results[i,]; j = 0
    while(is.na(oneMLE[1]) & j < 20) {
      randStart<-randparm()
      jj = 0
      while(is.infinite(TheNLLfn(randStart)) & jj < 20) {randStart<-randparm(); jj = jj + 1}
      try(oneMLE<- optim(randStart,TheNLLfn, method = "Nelder-Mead", control=list(maxit = 1000)))
      j = j+1}
    if(j < 20) results[i,]= c(oneMLE$par, oneMLE$value)
  }
  return(results[order(results[,dim(results)[2]]),])
}

NLL.results.Vial.f.Dist1 = getMLEv4.OneSex(1, List.Vial.f); NLL.results.Vial.f.Dist1[1:3,]
NLL.results.Vial.f.Dist2 = getMLEv4.OneSex(2, List.Vial.f); NLL.results.Vial.f.Dist2[1:3,] # 2nd
NLL.results.Vial.f.Dist3 = getMLEv4.OneSex(3, List.Vial.f); NLL.results.Vial.f.Dist3[1:3,] # best

NLL.results.Vial.m.Dist1 = getMLEv4.OneSex(1, List.Vial.m); NLL.results.Vial.m.Dist1[1:3,] # 2nd
NLL.results.Vial.m.Dist2 = getMLEv4.OneSex(2, List.Vial.m); NLL.results.Vial.m.Dist2[1:3,] # best
NLL.results.Vial.m.Dist3 = getMLEv4.OneSex(3, List.Vial.m); NLL.results.Vial.m.Dist3[1:3,]


NLL.results.Cage.f.Dist1 = getMLEv4.OneSex(1, List.Cage.f); NLL.results.Cage.f.Dist1[1:3,] # 2nd
NLL.results.Cage.f.Dist2 = getMLEv4.OneSex(2, List.Cage.f); NLL.results.Cage.f.Dist2[1:3,] # best
NLL.results.Cage.f.Dist3 = getMLEv4.OneSex(3, List.Cage.f); NLL.results.Cage.f.Dist3[1:3,]

NLL.results.Cage.m.Dist1 = getMLEv4.OneSex(1, List.Cage.m); NLL.results.Cage.m.Dist1[1:3,] # 2nd
NLL.results.Cage.m.Dist2 = getMLEv4.OneSex(2, List.Cage.m); NLL.results.Cage.m.Dist2[1:3,] # best
NLL.results.Cage.m.Dist3 = getMLEv4.OneSex(3, List.Cage.m); NLL.results.Cage.m.Dist3[1:3,]


NLL.results.Monog.f.Dist1 = getMLEv4.OneSex(1, List.Monog.f); NLL.results.Monog.f.Dist1[1:3,] # 2nd
NLL.results.Monog.f.Dist2 = getMLEv4.OneSex(2, List.Monog.f); NLL.results.Monog.f.Dist2[1:3,] # best
NLL.results.Monog.f.Dist3 = getMLEv4.OneSex(3, List.Monog.f); NLL.results.Monog.f.Dist3[1:3,]

NLL.results.Monog.m.Dist1 = getMLEv4.OneSex(1, List.Monog.m); NLL.results.Monog.m.Dist1[1:3,] # 2nd
NLL.results.Monog.m.Dist2 = getMLEv4.OneSex(2, List.Monog.m); NLL.results.Monog.m.Dist2[1:3,] # best
NLL.results.Monog.m.Dist3 = getMLEv4.OneSex(3, List.Monog.m); NLL.results.Monog.m.Dist3[1:3,]

## Betabinomial tends to be the best distribution

## Thisis a little function to convert the likelihood parameters to observed means and variances
convertBack<-function(parmVec){
  cat("NLL = ",parmVec[length(parmVec)], "\n")
  meanVec = c(parmVec[1], 0.75*parmVec[1] + 0.25*parmVec[2], sum(parmVec[1:2])/2, 0.25*parmVec[1] + 0.75*parmVec[2], parmVec[2])
  varVec = c(parmVec[3]*meanVec[1]/2, 0.5*(1+parmVec[4])*parmVec[3]*meanVec[2]/2, parmVec[4]*parmVec[3]*meanVec[1]/2, 0.5*(1+parmVec[4])*parmVec[3]*meanVec[4]/2,parmVec[3]*meanVec[3]/2)
  m = matrix(nrow =2, ncol = 5)
  m[1,] = meanVec; m[2,] = varVec;
  return(m)}

convertBack(NLL.results.Vial.f.Dist1[1,])
convertBack(NLL.results.Vial.f.Dist2[1,]) 
convertBack(NLL.results.Vial.f.Dist3[1,]) 
getTheBasicsForAllPopTypes(List.Vial.f)

convertBack(NLL.results.Vial.m.Dist1[1,])
convertBack(NLL.results.Vial.m.Dist2[1,]) ##best
convertBack(NLL.results.Vial.m.Dist3[1,])
getTheBasicsForAllPopTypes(List.Vial.m)


convertBack(NLL.results.Cage.f.Dist1[1,]) ## reasonably close to best
convertBack(NLL.results.Cage.f.Dist2[1,]) ## best
convertBack(NLL.results.Cage.f.Dist3[1,]) 
getTheBasicsForAllPopTypes(List.Cage.f)

convertBack(NLL.results.Cage.m.Dist1[1,])
convertBack(NLL.results.Cage.m.Dist2[1,]) ## best
convertBack(NLL.results.Cage.m.Dist3[1,])
getTheBasicsForAllPopTypes(List.Cage.m)

convertBack(NLL.results.Monog.f.Dist1[1,])
convertBack(NLL.results.Monog.f.Dist2[1,]) ## very close to best
convertBack(NLL.results.Monog.f.Dist3[1,]) ## best
getTheBasicsForAllPopTypes(List.Monog.f)

convertBack(NLL.results.Monog.m.Dist1[1,])
convertBack(NLL.results.Monog.m.Dist2[1,]) ## best
convertBack(NLL.results.Monog.m.Dist3[1,])
getTheBasicsForAllPopTypes(List.Monog.m)





#####################################################
######## *** Likelihood functions to run on data ****
######################################################


NLLv4.BothSexes<-function(distID.f, distID.m, listOfDataVecs.f, listOfDataVecs.m, parmVec){
  
  NegLog.NegBinomPDF<-function(data, mu, v) -sum(dnbinom(data, size = (mu^2)/(v - mu), prob = mu/v, log = TRUE))
  NegLog.BetaBinomPDF<-function(data, mu, v) -sum(dbbinom(data, size = 32, alpha = -mu*(mu^2 - 32*mu +v)/(mu^2 - 32*mu + 32*v), beta = (mu -32)*(mu^2 - 32*mu +v)/(mu^2 -32*mu +32*v), log = TRUE))
  NegLog.NearZeroInflatedPoisPDF<-function(data, mu, v, lambdaNZ = 0) -sum(log(dnzip(data, pNZ = (v - mu)/(v + mu*(mu - 1) + lambdaNZ*(lambdaNZ - 2*mu)), 
                                                                                       lambdaFree = (v + mu*(mu - lambdaNZ - 1))/(mu - lambdaNZ), lambdaNZ = lambdaNZ)))
  
  if(distID.f == 1) This.Dist.f =  NegLog.NegBinomPDF ## negative binomial
  if(distID.f == 2) This.Dist.f =  NegLog.BetaBinomPDF ## beta binomail
  if(distID.f == 3) This.Dist.f =  NegLog.NearZeroInflatedPoisPDF ## "near-zero" inflated poisson
  if(distID.m == 1) This.Dist.m =  NegLog.NegBinomPDF ## negative binomial
  if(distID.m == 2) This.Dist.m =  NegLog.BetaBinomPDF ## beta binomail
  if(distID.m == 3) This.Dist.m =  NegLog.NearZeroInflatedPoisPDF ## "near-zero" inflated poisson
  
  z.means = parmVec[4:7]
  z.VarF.LowHet = parmVec[1]
  z.psi = parmVec[2]
  z.gamma = parmVec[3]
  z.meanVec.f<-c(z.means[1], 0.75*z.means[1] + 0.25*z.means[2], 0.5*sum(z.means[1:2]), 0.25*z.means[1] + 0.75*z.means[2], z.means[2])
  z.meanVec.m<-c(z.means[3], 0.75*z.means[3] + 0.25*z.means[4], 0.5*sum(z.means[3:4]), 0.25*z.means[3] + 0.75*z.means[4], z.means[4])
  z.varVec.f.IfMean2<-c(z.VarF.LowHet, 0.5*(1+z.gamma)*z.VarF.LowHet, z.gamma*z.VarF.LowHet, 0.5*(1+z.gamma)*z.VarF.LowHet, z.VarF.LowHet)
  z.varVec.f<- z.varVec.f.IfMean2*z.meanVec.f/2
  z.varVec.m<-z.psi*z.varVec.f.IfMean2*z.meanVec.m/2
  if(length(listOfDataVecs.f) == 3){
    z.meanVec.f<-z.meanVec.f[c(1,3,5)]; z.meanVec.m<-z.meanVec.m[c(1,3,5)]
    z.varVec.f<-z.varVec.f[c(1,3,5)]; z.varVec.m<-z.varVec.m[c(1,3,5)]
  }

  AllVarGreaterThanMeans <- !any(c(z.meanVec.f, z.meanVec.m) >= c(z.varVec.f, z.varVec.m))
  if(!AllVarGreaterThanMeans)  ans <- 911911911911911
  
  if(AllVarGreaterThanMeans){
    ## function to calculate negative log likelihood based specified distribution
    ans = 0
    for(ii in 1:length(listOfDataVecs.f)) {
      x.f <- 10^10; x.m <- 10^10
      try(x.f <- This.Dist.f(sapply(listOfDataVecs.f[ii], as.numeric), z.meanVec.f[ii], z.varVec.f[ii]) )
      try(x.m <- This.Dist.m(sapply(listOfDataVecs.m[ii], as.numeric), z.meanVec.m[ii], z.varVec.m[ii]) )
      if(is.na(x.f)) x.f <- 10^10
      if(is.na(x.m)) x.m <- 10^10
      ans = ans + x.f + x.m
    }
  }
  return(ans)
}



# The function below is used to get random parameter combinations for searching parameter space for best likelihood
# default values y.varF = -1, y.psi = -1 and y.gamma = -1 indicate these parameters are unconstrained
# otherwise they are constrained to the specified values
randparm.BothSexes<- function(nPopTypes = 5, y.varF = -1, y.psi = -1, y.gamma = -1) {
  needOne<- TRUE; countTries <-0; nMaxTries <-25;
  while(needOne & (countTries < nMaxTries)){
    z.means = runif(4, min = 1.9, max = 8)
    ifelse(y.varF == -1, z.VarF.LowHet <- runif(1, min = 1.9, max = 8), z.VarF.LowHet <- y.varF)
    ifelse(y.psi == -1, z.psi <- runif(1, min = 2.01/z.VarF.LowHet, max = 20.01/z.VarF.LowHet), z.psi <- y.psi)
    ifelse(y.gamma == -1, z.gamma <- runif(1, min = 0.5, max = 2.5), z.gamma <- y.gamma)
    zParmVec<- c(z.VarF.LowHet, z.psi, z.gamma, z.means)
    z.meanVec.f<-c(z.means[1], 0.75*z.means[1] + 0.25*z.means[2], 0.5*sum(z.means[1:2]), 0.25*z.means[1] + 0.75*z.means[2], z.means[2])
    z.meanVec.m<-c(z.means[3], 0.75*z.means[3] + 0.25*z.means[4], 0.5*sum(z.means[3:4]), 0.25*z.means[3] + 0.75*z.means[4], z.means[4])
    z.varVec.f.IfMean2<-c(z.VarF.LowHet, 0.5*(1+z.gamma)*z.VarF.LowHet, z.gamma*z.VarF.LowHet, 0.5*(1+z.gamma)*z.VarF.LowHet, z.VarF.LowHet)
    z.varVec.f<- z.varVec.f.IfMean2*z.meanVec.f/2
    z.varVec.m<-z.psi*z.varVec.f.IfMean2*z.meanVec.m/2
    if(nPopTypes == 3){
      z.meanVec.f<-z.meanVec.f[c(1,3,5)]; z.meanVec.m<-z.meanVec.m[c(1,3,5)]
      z.varVec.f<-z.varVec.f[c(1,3,5)]; z.varVec.m<-z.varVec.m[c(1,3,5)]
    }
    
    z.AllVarGreaterThanMeans <- !any(c(z.meanVec.f, z.meanVec.m) >= c(z.varVec.f, z.varVec.m))
    if(z.AllVarGreaterThanMeans) needOne <- FALSE
    countTries<-countTries+1
    if(countTries == nMaxTries) zParmVec<-c(runif(1, min = 1.9, max = 2.5), runif(2, min = 0.9, max = 1.25),  runif(4, min = 1.9, max = 5))
  }
  return(zParmVec)
}


## This will find maximum likelihood estimates and their negative log likelihoods
## for each of 4 models
## Model 1 - only VarF.Low  and means unconstrained (i.e., 5 parameter)
## Model 2 - 6 free parameters: VarF.Low and psi and means
## Model 3 - 6 free parameters: VarF.Low and gamma and means
## Model 4 - The "full" model: 7 free parameters
## Output is table with results from each model
## Columns are: (1-7) ML estimates of VarF.Low, psi, gamma, meanF.Low, meanF.High, meanM.Low, meanM.High; (8) number of free parameters; (9) NLL
getMLEv4.4Models<-function(distID.f, distID.m, listOfDataVecs.f, listOfDataVecs.m, nStartingPlaces = 50){
  n<-nStartingPlaces ## number of starting places for optimization of likelihood
  nPopTypes = length(listOfDataVecs.f)
  
  nParm = length(randparm.BothSexes())
  
  #Model 1; psi = gamma = 1; (5 parameters: VarFLowHet and 4 means)
  fn.Mod1<-function(parmVec) NLLv4.BothSexes(distID.f, distID.m, listOfDataVecs.f, listOfDataVecs.m,  
                                           c(parmVec[1], 1, 1, parmVec[2], parmVec[3], parmVec[4],parmVec[5]))
  results<-matrix(nrow = n, ncol = 2 + nParm + 2)
  for(i in 1:n){
    oneMLE = results[i,]; j = 0
    while(is.na(oneMLE[1]) & j < 20) {
      randStart<-randparm.BothSexes(nPopTypes = nPopTypes, y.psi = 1, y.gamma = 1)[c(1,4:7)]
      jj = 0
      while(is.infinite(fn.Mod1(randStart)) & jj < 20) {randparm.BothSexes(nPopTypes = nPopTypes, y.psi = 1, y.gamma = 1)[c(1,4:7)]; jj = jj + 1}
      try(oneMLE<- optim(randStart,fn.Mod1, method = "Nelder-Mead", control=list(maxit = 1000)))
      j = j+1}
    if(j < 20) results[i,]= c(distID.f, distID.m, oneMLE$par[1], 1,1,oneMLE$par[2:5], 5, oneMLE$value)
  }
  Model1<-results[order(results[,dim(results)[2]]),]
  cat("best model 1: \n", Model1[1,], "\n")
  
  # Model 2 - 6 free parameters: VarF.Low and psi and means; gamma = 1
  fn.Mod2<-function(parmVec) NLLv4.BothSexes(distID.f, distID.m, listOfDataVecs.f, listOfDataVecs.m,
                                           c(parmVec[1], parmVec[2], 1, parmVec[3], parmVec[4],parmVec[5], parmVec[6]))
  results<-matrix(nrow = n, ncol = 2 + nParm + 2)
  for(i in 1:n){
    oneMLE = results[i,]; j = 0
    while(is.na(oneMLE[1]) & j < 20) {
      randStart<-randparm.BothSexes(nPopTypes = nPopTypes, y.gamma = 1)[c(1:2, 4:7)]
      jj = 0
      while(is.infinite(fn.Mod2(randStart)) & jj < 20) {randStart<-randparm.BothSexes(nPopTypes = nPopTypes, y.gamma = 1)[c(1:2, 4:7)]; jj = jj + 1}
      try(oneMLE<- optim(randStart,fn.Mod2, method = "Nelder-Mead", control=list(maxit = 1000)))
      j = j+1}
    if(j < 20) results[i,]= c(distID.f, distID.m, oneMLE$par[1:2], 1, oneMLE$par[3:6], 6, oneMLE$value)
  }
  Model2<-results[order(results[,dim(results)[2]]),]
  cat("best model 2: \n", Model2[1,], "\n")

  # Model 3 - 6 free parameters: VarF.Low and gamma and means; psi = 1
  fn.Mod3<-function(parmVec) NLLv4.BothSexes(distID.f, distID.m, listOfDataVecs.f, listOfDataVecs.m,
                                                c(parmVec[1], 1, parmVec[2], parmVec[3], parmVec[4], parmVec[5], parmVec[6]))
  results<-matrix(nrow = n, ncol = 2 + nParm + 2)
  for(i in 1:n){
    oneMLE = results[i,]; j = 0
    while(is.na(oneMLE[1]) & j < 20) {
      randStart<-randparm.BothSexes(nPopTypes = nPopTypes, y.psi = 1)[c(1,3:7)]
      jj = 0
      while(is.infinite(fn.Mod3(randStart)) & jj < 20) {randStart<-randparm.BothSexes(nPopTypes = nPopTypes, y.psi = 1)[c(1,3:7)]; jj = jj + 1}
      try(oneMLE<- optim(randStart,fn.Mod3, method = "Nelder-Mead", control=list(maxit = 1000)))
      j = j+1}
    if(j < 20) results[i,]= c(distID.f, distID.m, oneMLE$par[1], 1, oneMLE$par[2:6], 6, oneMLE$value)
  }
  Model3<-results[order(results[,dim(results)[2]]),]
  cat("best model 3: \n", Model3[1,], "\n")
  
  # Model4 - 7 free parameters: VarF.Low, psi and gamma and means
  
  fn.Mod4<-function(parmVec) NLLv4.BothSexes(distID.f, distID.m, listOfDataVecs.f, listOfDataVecs.m,
                                                c(parmVec[1], parmVec[2], parmVec[3], parmVec[4], parmVec[5], parmVec[6], parmVec[7]))
  results<-matrix(nrow = n, ncol = 2+ nParm + 2)
  for(i in 1:n){
    oneMLE = results[i,]; j = 0
    while(is.na(oneMLE[1]) & j < 20) {
      randStart<-randparm.BothSexes(nPopTypes = nPopTypes)
      jj = 0
      while(is.infinite(fn.Mod4(randStart)) & jj < 20) {randStart<-randparm.BothSexes(nPopTypes = nPopTypes); jj = jj + 1}
      try(oneMLE<- optim(randStart,fn.Mod4, method = "Nelder-Mead", control=list(maxit = 1000)))
      j = j+1}
    if(j < 20) results[i,]= c(distID.f, distID.m, oneMLE$par[1:7], 7, oneMLE$value)
  }
  Model4<-results[order(results[,dim(results)[2]]),]
  cat("best model 4: \n", Model4[1,], "\n")

  nbest<-5 #This many of the top results from each model
  return(rbind(Model1[1:nbest,], Model2[1:nbest,], Model3[1:nbest,], Model4[1:nbest,]))
}

FourModels.Vial.F2M2<-getMLEv4.4Models(2, 2, List.Vial.f, List.Vial.m); FourModels.Vial.F2M2
FourModels.Cage.F2M2<-getMLEv4.4Models(2, 2, List.Cage.f, List.Cage.m); FourModels.Cage.F2M2
FourModels.Monog.F2M2<-getMLEv4.4Models(2, 2, List.Monog.f, List.Monog.m); FourModels.Monog.F2M2
#

### Here are the best MLEs from the 4 models
FourModels.Vial.F2M2[c(1, 6, 11, 16),]
FourModels.Cage.F2M2[c(1, 6, 11, 16),]  
FourModels.Monog.F2M2[c(1, 6, 11, 16),]
########################

#####################################################
######## *** Make Likelihood Profile for VarF ****
######################################################
## Use this to generate a likelihood profile for varF after specifying a list of varF values to evaluate
getVarFLikelihoodProfilev4<-function(distID.f, distID.m, listOfDataVecs.f, listOfDataVecs.m, listOfVarFValuesToEvaluate){
  varF.results<-matrix(nrow = 2, ncol = length(listOfVarFValuesToEvaluate))
  n<-20 ## number of starting places for optimization of likelihood
  nPopTypes = length(listOfDataVecs.f)
  nParm = length(randparm.BothSexes())
  
  for(gg in 1:length(listOfVarFValuesToEvaluate)){
    this.varF = listOfVarFValuesToEvaluate[gg]
    fn.Mod.varFspecified<-function(parmVec) NLLv4.BothSexes(distID.f, distID.m, listOfDataVecs.f, listOfDataVecs.m,
                                                            c(this.varF, parmVec[1], parmVec[2], parmVec[3], parmVec[4],parmVec[5], parmVec[6]))
    results<-matrix(nrow = n, ncol = 2 + nParm + 2)
    for(i in 1:n){
      oneMLE = results[i,]; j = 0
      while(is.na(oneMLE[1]) & j < 20) {
        randStart<-randparm.BothSexes(nPopTypes = nPopTypes, y.varF = this.varF)[c(2:7)]
        jj = 0
        while(is.infinite(fn.Mod.varFspecified(randStart)) & jj < 20) {randStart<-randparm.BothSexes(nPopTypes = nPopTypes, y.varF = this.varF)[c(2:7)]; jj = jj + 1}
        try(oneMLE<- optim(randStart,fn.Mod.varFspecified, method = "Nelder-Mead", control=list(maxit = 1000)))
        j = j+1}
      if(j < 20) results[i,]= c(distID.f, distID.m, oneMLE$par[1:2], this.varF, oneMLE$par[3:6], 6, oneMLE$value)
    }
    sorted.results<-results[order(results[,dim(results)[2]]),]
    varF.results[1, gg] = this.varF
    varF.results[2, gg] = sorted.results[1, dim(sorted.results)[2]]
    cat(c(varF.results[1, gg], varF.results[2, gg]), "\n")
  }
  return(varF.results)
}

## Change this as you see fit

UseThisListOfVarFValuesToEvaluate.Vial.F2M2<-c(3.0+ 0.02*(1:100), 2.5, 2.6, 2.7, 2.8, 2.9, 5.1, 5.2, 5.3, 5.4) 
UseThisListOfVarFValuesToEvaluate.Cage.F2M2<-c(2.7+ 0.02*(1:100),2.4, 2.5, 2.6, 4.8, 4.9, 5.1, 5.2)
UseThisListOfVarFValuesToEvaluate.Monog.F2M2<-c(4+ 0.02*(1:150), 3.5, 3.6, 3.7, 3.8, 3.9, 7.1, 7.2, 7.3)

VarFProfile.Vial.F2M2<-getVarFLikelihoodProfilev4(2, 2, List.Vial.f, List.Vial.m, UseThisListOfVarFValuesToEvaluate.Vial.F2M2); VarFProfile.Vial.F2M2
VarFProfile.Cage.F2M2<-getVarFLikelihoodProfilev4(2, 2, List.Cage.f, List.Cage.m, UseThisListOfVarFValuesToEvaluate.Cage.F2M2); VarFProfile.Cage.F2M2
VarFProfile.Monog.F2M2<-getVarFLikelihoodProfilev4(2, 2, List.Monog.f, List.Monog.m, UseThisListOfVarFValuesToEvaluate.Monog.F2M2); VarFProfile.Monog.F2M2

UseThisListOfVarFValuesToEvaluate.Vial.F2M2.more<- c(2.3, 2.8 + 0.02*(1:10))
UseThisListOfVarFValuesToEvaluate.Cage.F2M2.more<- c(4.71 + 0.02*(1:15), 5.4, 5.6, 5.8)
UseThisListOfVarFValuesToEvaluate.Monog.F2M2.more<-c(3.2, 7.4)

VarFProfile.Vial.F2M2.more<-getVarFLikelihoodProfilev4(2, 2, List.Vial.f, List.Vial.m, UseThisListOfVarFValuesToEvaluate.Vial.F2M2.more); VarFProfile.Vial.F2M2.more
VarFProfile.Cage.F2M2.more<-getVarFLikelihoodProfilev4(2, 2, List.Cage.f, List.Cage.m, UseThisListOfVarFValuesToEvaluate.Cage.F2M2.more); VarFProfile.Cage.F2M2.more
VarFProfile.Monog.F2M2.more<-getVarFLikelihoodProfilev4(2, 2, List.Monog.f, List.Monog.m, UseThisListOfVarFValuesToEvaluate.Monog.F2M2.more); VarFProfile.Monog.F2M2.more

VarFProfile.Vial.F2M2.all<-cbind(VarFProfile.Vial.F2M2, VarFProfile.Vial.F2M2.more)
VarFProfile.Vial.F2M2.all<-VarFProfile.Vial.F2M2.all[,order(VarFProfile.Vial.F2M2.all[1,])]
VarFProfile.Cage.F2M2.all<-cbind(VarFProfile.Cage.F2M2, VarFProfile.Cage.F2M2.more)
VarFProfile.Cage.F2M2.all<-VarFProfile.Cage.F2M2.all[,order(VarFProfile.Cage.F2M2.all[1,])]
VarFProfile.Monog.F2M2.all<-cbind(VarFProfile.Monog.F2M2, VarFProfile.Monog.F2M2.more)
VarFProfile.Monog.F2M2.all<-VarFProfile.Monog.F2M2.all[,order(VarFProfile.Monog.F2M2.all[1,])]



VarFProfile.Vial.F2M2.final<-t(rbind(VarFProfile.Vial.F2M2.all, min(FourModels.Vial.F2M2[,11]) - VarFProfile.Vial.F2M2.all[2,]))
VarFProfile.Cage.F2M2.final<-t(rbind(VarFProfile.Cage.F2M2.all, min(FourModels.Cage.F2M2[,11]) - VarFProfile.Cage.F2M2.all[2,]))
VarFProfile.Monog.F2M2.final<-t(rbind(VarFProfile.Monog.F2M2.all, min(FourModels.Monog.F2M2[,11]) - VarFProfile.Monog.F2M2.all[2,]))



{plot(c(VarFProfile.Vial.F2M2.final[,1], VarFProfile.Cage.F2M2.final[,1], VarFProfile.Monog.F2M2.final[,1]), 
      c(VarFProfile.Vial.F2M2.final[,3], VarFProfile.Cage.F2M2.final[,3], VarFProfile.Monog.F2M2.final[,3]), 
      col = c(rep("purple", length(VarFProfile.Vial.F2M2.final[,1])), rep("blue", length(VarFProfile.Cage.F2M2.final[,1])), rep("orange", length(VarFProfile.Monog.F2M2.final[,1]))), ylim = c(-10, 1)); 
  abline(a = 0, b= 0, col = "grey"); abline(a = -2, b= 0, lty = 2, col = "red"); abline(v=2, col = "grey")}


#####################################################
######## *** Make Likelihood Profile for GAMMA (reminder that "GAMMA" is called "lambda" in the mansucript) ****
######################################################
## Use this to generate a likelihood profile for gamma after specifying a list of gamma values to evaluate
getGammaLikelihoodProfilev4<-function(distID.f, distID.m, listOfDataVecs.f, listOfDataVecs.m, listOfGammaValuesToEvaluate){
  gamma.results<-matrix(nrow = 2, ncol = length(listOfGammaValuesToEvaluate))
  n<-20 ## number of starting places for optimization of likelihood
  nPopTypes = length(listOfDataVecs.f)
  nParm = length(randparm.BothSexes())
  
  for(gg in 1:length(listOfGammaValuesToEvaluate)){
    this.gamma = listOfGammaValuesToEvaluate[gg]
    fn.Mod.gammaspecified<-function(parmVec) NLLv4.BothSexes(distID.f, distID.m, listOfDataVecs.f, listOfDataVecs.m,
                                                             c(parmVec[1], parmVec[2], this.gamma, parmVec[3], parmVec[4],parmVec[5], parmVec[6]))
    results<-matrix(nrow = n, ncol = 2 + nParm + 2)
    for(i in 1:n){
      oneMLE = results[i,]; j = 0
      while(is.na(oneMLE[1]) & j < 20) {
        randStart<-randparm.BothSexes(nPopTypes = nPopTypes, y.gamma = this.gamma)[c(1:2, 4:7)]
        jj = 0
        while(is.infinite(fn.Mod.gammaspecified(randStart)) & jj < 20) {randStart<-randparm.BothSexes(nPopTypes = nPopTypes, y.gamma = this.gamma)[c(1:2, 4:7)]; jj = jj + 1}
        try(oneMLE<- optim(randStart,fn.Mod.gammaspecified, method = "Nelder-Mead", control=list(maxit = 1000)))
        j = j+1}
      if(j < 20) results[i,]= c(distID.f, distID.m, oneMLE$par[1:2], this.gamma, oneMLE$par[3:6], 6, oneMLE$value)
    }
    sorted.results<-results[order(results[,dim(results)[2]]),]
    gamma.results[1, gg] = this.gamma
    gamma.results[2, gg] = sorted.results[1, dim(sorted.results)[2]]
    cat(c(gamma.results[1, gg], gamma.results[2, gg]), "\n")
  }
  return(gamma.results)
}

## Change this as you see fit
UseThisListOfGammaValuesToEvaluate.Vial.F2M2 <-sort(c(1.4 + 0.02*(0:60), 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.35, 2.65, 2.7, 2.8, 3))
UseThisListOfGammaValuesToEvaluate.Cage.F2M2<-sort(c(1.2 + 0.02*(0:60), 0.7, 0.9, 1, 1.1, 2.45, 2.5, 2.6, 2.7, 2.9))
UseThisListOfGammaValuesToEvaluate.Monog.F2M2<-sort(c(0.9 + 0.02*(0:50), 0.6, 0.8, 0.85, 1.95, 2, 2.1, 2.2, 2.4, 2.6))

GammaProfile.Vial.F2M2<-getGammaLikelihoodProfilev4(2, 2, List.Vial.f, List.Vial.m, UseThisListOfGammaValuesToEvaluate.Vial.F2M2); GammaProfile.Vial.F2M2
GammaProfile.Cage.F2M2<-getGammaLikelihoodProfilev4(2, 2, List.Cage.f, List.Cage.m, UseThisListOfGammaValuesToEvaluate.Cage.F2M2); GammaProfile.Cage.F2M2
GammaProfile.Monog.F2M2<-getGammaLikelihoodProfilev4(2, 2, List.Monog.f, List.Monog.m, UseThisListOfGammaValuesToEvaluate.Monog.F2M2); GammaProfile.Monog.F2M2

UseThisListOfGammaValuesToEvaluate.Vial.F2M2.more<-c(0.93 +0.02*(0:20), 0.7, 1.37)
UseThisListOfGammaValuesToEvaluate.Cage.F2M2.more<-c(1.15, 1.05)
UseThisListOfGammaValuesToEvaluate.Monog.F2M2.more<-c(0.7, 0.75, 0.81 + 0.02*(0:4))
GammaProfile.Vial.F2M2.more<-getGammaLikelihoodProfilev4(2, 2, List.Vial.f, List.Vial.m, UseThisListOfGammaValuesToEvaluate.Vial.F2M2.more); GammaProfile.Vial.F2M2.more
GammaProfile.Cage.F2M2.more<-getGammaLikelihoodProfilev4(2, 2, List.Cage.f, List.Cage.m, UseThisListOfGammaValuesToEvaluate.Cage.F2M2.more); GammaProfile.Cage.F2M2.more
GammaProfile.Monog.F2M2.more<-getGammaLikelihoodProfilev4(2, 2, List.Monog.f, List.Monog.m, UseThisListOfGammaValuesToEvaluate.Monog.F2M2.more); GammaProfile.Monog.F2M2.more

GammaProfile.Vial.F2M2.all<-cbind(GammaProfile.Vial.F2M2, GammaProfile.Vial.F2M2.more)
GammaProfile.Vial.F2M2.all<-GammaProfile.Vial.F2M2.all[,order(GammaProfile.Vial.F2M2.all[1,])]
GammaProfile.Cage.F2M2.all<-cbind(GammaProfile.Cage.F2M2, GammaProfile.Cage.F2M2.more)
GammaProfile.Cage.F2M2.all<-GammaProfile.Cage.F2M2.all[,order(GammaProfile.Cage.F2M2.all[1,])]
GammaProfile.Monog.F2M2.all<-cbind(GammaProfile.Monog.F2M2, GammaProfile.Monog.F2M2.more)
GammaProfile.Monog.F2M2.all<-GammaProfile.Monog.F2M2.all[,order(GammaProfile.Monog.F2M2.all[1,])]

GammaProfile.Vial.F2M2.final<-t(rbind(GammaProfile.Vial.F2M2.all, min(FourModels.Vial.F2M2[,11]) - GammaProfile.Vial.F2M2.all[2,]))
GammaProfile.Cage.F2M2.final<-t(rbind(GammaProfile.Cage.F2M2.all, min(FourModels.Cage.F2M2[,11]) - GammaProfile.Cage.F2M2.all[2,]))
GammaProfile.Monog.F2M2.final<-t(rbind(GammaProfile.Monog.F2M2.all, min(FourModels.Monog.F2M2[,11]) - GammaProfile.Monog.F2M2.all[2,]))


{plot(c(GammaProfile.Vial.F2M2.final[,1], GammaProfile.Cage.F2M2.final[,1], GammaProfile.Monog.F2M2.final[,1]), 
      c(GammaProfile.Vial.F2M2.final[,3], GammaProfile.Cage.F2M2.final[,3], GammaProfile.Monog.F2M2.final[,3]), 
      col = c(rep("purple", length(GammaProfile.Vial.F2M2.final[,1])), rep("blue", length(GammaProfile.Cage.F2M2.final[,1])), rep("orange", length(GammaProfile.Monog.F2M2.final[,1]))), ylim = c(-10, 1)); 
  abline(a = 0, b= 0, col = "grey"); abline(a = -2, b= 0, lty = 2, col = "red"); abline(v=1, lty = 2, col = "grey")}


#####################################################
######## *** Likelihood profiles for PSI ****
######################################################
## Use this to generate a likelihood profile for psi after specifying a list of psi values to evaluate
getPsiLikelihoodProfilev4<-function(distID.f, distID.m, listOfDataVecs.f, listOfDataVecs.m, listOfPsiValuesToEvaluate){
  psi.results<-matrix(nrow = 2, ncol = length(listOfPsiValuesToEvaluate))
  n<-20 ## number of starting places for optimization of likelihood
  nPopTypes = length(listOfDataVecs.f)
  nParm = length(randparm.BothSexes())
  
  for(gg in 1:length(listOfPsiValuesToEvaluate)){
    this.psi = listOfPsiValuesToEvaluate[gg]
    fn.Mod.psispecified<-function(parmVec) NLLv4.BothSexes(distID.f, distID.m, listOfDataVecs.f, listOfDataVecs.m,
                                                             c(parmVec[1], this.psi, parmVec[2], parmVec[3], parmVec[4],parmVec[5], parmVec[6]))
    results<-matrix(nrow = n, ncol = 2 + nParm + 2)
    for(i in 1:n){
      oneMLE = results[i,]; j = 0
      while(is.na(oneMLE[1]) & j < 20) {
        randStart<-randparm.BothSexes(nPopTypes = nPopTypes, y.psi = this.psi)[c(1, 3:7)]
        jj = 0
        while(is.infinite(fn.Mod.psispecified(randStart)) & jj < 20) {randStart<-randparm.BothSexes(nPopTypes = nPopTypes, y.psi = this.psi)[c(1, 3:7)]; jj = jj + 1}
        try(oneMLE<- optim(randStart,fn.Mod.psispecified, method = "Nelder-Mead", control=list(maxit = 1000)))
        j = j+1}
      if(j < 20) results[i,]= c(distID.f, distID.m, oneMLE$par[1], this.psi, oneMLE$par[2:6], 6, oneMLE$value)
    }
    sorted.results<-results[order(results[,dim(results)[2]]),]
    psi.results[1, gg] = this.psi
    psi.results[2, gg] = sorted.results[1, dim(sorted.results)[2]]
    cat(c(psi.results[1, gg], psi.results[2, gg]), "\n")
  }
  return(psi.results)
}

## Change this as you see fit
UseThisListOfPhiValuesToEvaluate.Vial.F2M2<-sort(c(1.45 + 0.02*(0:50), 0.7, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 2.5, 2.55, 2.6, 2.7, 2.9, 3.1))
UseThisListOfPhiValuesToEvaluate.Cage.F2M2<-sort(c(1.35 + 0.02*(0:50), 0.7, 0.9, 1, 1.1, 1.2, 1.3, 2.4, 2.45, 2.5, 2.6, 2.8, 3))
UseThisListOfPhiValuesToEvaluate.Monog.F2M2<-sort(c(0.9 + 0.02*(0:50), 0.5, 0.6, 0.7, 0.8, 0.85, 1.95, 2, 2.1, 2.2, 2.4))

PsiProfile.Vial.F2M2<-getPsiLikelihoodProfilev4(2, 2, List.Vial.f, List.Vial.m, UseThisListOfPhiValuesToEvaluate.Vial.F2M2); PsiProfile.Vial.F2M2
PsiProfile.Cage.F2M2<-getPsiLikelihoodProfilev4(2, 2, List.Cage.f, List.Cage.m, UseThisListOfPhiValuesToEvaluate.Cage.F2M2); PsiProfile.Cage.F2M2
PsiProfile.Monog.F2M2<-getPsiLikelihoodProfilev4(2, 2, List.Monog.f, List.Monog.m, UseThisListOfPhiValuesToEvaluate.Monog.F2M2); PsiProfile.Monog.F2M2

PsiProfile.Vial.F2M2.final<-t(rbind(PsiProfile.Vial.F2M2, min(FourModels.Vial.F2M2[,11]) - PsiProfile.Vial.F2M2[2,]))
PsiProfile.Cage.F2M2.final<-t(rbind(PsiProfile.Cage.F2M2, min(FourModels.Cage.F2M2[,11]) - PsiProfile.Cage.F2M2[2,]))
PsiProfile.Monog.F2M2.final<-t(rbind(PsiProfile.Monog.F2M2, min(FourModels.Monog.F2M2[,11]) - PsiProfile.Monog.F2M2[2,]))



{plot(c(PsiProfile.Vial.F2M2.final[,1], PsiProfile.Cage.F2M2.final[,1], PsiProfile.Monog.F2M2.final[,1]), 
      c(PsiProfile.Vial.F2M2.final[,3], PsiProfile.Cage.F2M2.final[,3], PsiProfile.Monog.F2M2.final[,3]), 
      col = c(rep("purple", length(PsiProfile.Vial.F2M2.final[,1])), rep("blue", length(PsiProfile.Cage.F2M2.final[,1])), rep("orange", length(PsiProfile.Monog.F2M2.final[,1]))), ylim = c(-10, 1)); 
  abline(a = 0, b= 0, col = "grey"); abline(a = -2, b= 0, lty = 2, col = "red"); abline(v=1, lty = 2, col = "grey")}



#############################
####################################
######### Here we check if it is necessary to add in one additional level 
######### of complexity by allowing the heterogeneity effect to differ between the sexes
######### The male variance is modeled at varM = psi*VarF in low heterogeneity (A and E)
######### The male variance in the other population types (B, C, D)
######### is scaled to varM in low heterogeneity by gamma.M
######### which can differ from gamma.F in the most complex model.
######### If so, that means varM/varF is not the same in B, C, or D as in A and E (which would complicate interpreation of psi parameter)
######### There is no support for the having a separate gamma.M and gamma.F in any 
######### of the mating regimes.
####################################

NLLv5.BothSexes<-function(distID.f, distID.m, listOfDataVecs.f, listOfDataVecs.m, parmVec){
  
  NegLog.NegBinomPDF<-function(data, mu, v) -sum(dnbinom(data, size = (mu^2)/(v - mu), prob = mu/v, log = TRUE))
  NegLog.BetaBinomPDF<-function(data, mu, v) -sum(dbbinom(data, size = 32, alpha = -mu*(mu^2 - 32*mu +v)/(mu^2 - 32*mu + 32*v), beta = (mu -32)*(mu^2 - 32*mu +v)/(mu^2 -32*mu +32*v), log = TRUE))
  NegLog.NearZeroInflatedPoisPDF<-function(data, mu, v, lambdaNZ = 0) -sum(log(dnzip(data, pNZ = (v - mu)/(v + mu*(mu - 1) + lambdaNZ*(lambdaNZ - 2*mu)), 
                                                                                     lambdaFree = (v + mu*(mu - lambdaNZ - 1))/(mu - lambdaNZ), lambdaNZ = lambdaNZ)))
  
  if(distID.f == 1) This.Dist.f =  NegLog.NegBinomPDF ## negative binomial
  if(distID.f == 2) This.Dist.f =  NegLog.BetaBinomPDF ## beta binomail
  if(distID.f == 3) This.Dist.f =  NegLog.NearZeroInflatedPoisPDF ## "near-zero" inflated poisson
  if(distID.m == 1) This.Dist.m =  NegLog.NegBinomPDF ## negative binomial
  if(distID.m == 2) This.Dist.m =  NegLog.BetaBinomPDF ## beta binomail
  if(distID.m == 3) This.Dist.m =  NegLog.NearZeroInflatedPoisPDF ## "near-zero" inflated poisson
  
  z.means = parmVec[5:8]
  z.VarF.LowHet = parmVec[1]
  z.psi.LowHet = parmVec[2]
  z.gamma.F = parmVec[3]
  z.gamma.M = parmVec[4]
  z.meanVec.f<-c(z.means[1], 0.75*z.means[1] + 0.25*z.means[2], 0.5*sum(z.means[1:2]), 0.25*z.means[1] + 0.75*z.means[2], z.means[2])
  z.meanVec.m<-c(z.means[3], 0.75*z.means[3] + 0.25*z.means[4], 0.5*sum(z.means[3:4]), 0.25*z.means[3] + 0.75*z.means[4], z.means[4])
  z.varVec.f.IfMean2<-c(z.VarF.LowHet, 0.5*(1+z.gamma.F)*z.VarF.LowHet, z.gamma.F*z.VarF.LowHet, 0.5*(1+z.gamma.F)*z.VarF.LowHet, z.VarF.LowHet)
  z.varVec.m.IfMean2<-z.psi.LowHet*z.VarF.LowHet*c(1, 0.5*(1+z.gamma.M), z.gamma.M, 0.5*(1+z.gamma.M), 1)
  z.varVec.f<- z.varVec.f.IfMean2*z.meanVec.f/2
  z.varVec.m<-z.varVec.m.IfMean2*z.meanVec.m/2
  if(length(listOfDataVecs.f) == 3){
    z.meanVec.f<-z.meanVec.f[c(1,3,5)]; z.meanVec.m<-z.meanVec.m[c(1,3,5)]
    z.varVec.f<-z.varVec.f[c(1,3,5)]; z.varVec.m<-z.varVec.m[c(1,3,5)]
  }
  
  if(any(is.na(c(z.meanVec.f, z.meanVec.m,z.varVec.f, z.varVec.m)))) c("problem: parmVec = ", 
                                                                       parmVec, "\n z.meanVec.f = ", z.meanVec.f, "\n z.meanVec.m = ",z.meanVec.m,"\n z.varVec.f = ", z.varVec.f, "\n z.varVec.m = ", z.varVec.m, "\n")
  
  AllVarGreaterThanMeans <- !any(c(z.meanVec.f, z.meanVec.m) >= c(z.varVec.f, z.varVec.m))
  if(!AllVarGreaterThanMeans)  ans <- 911911911911911
  
  #cat(z.meanVec.f, "\n", z.varVec.f,"\n", z.meanVec.m, "\n", z.varVec.m,"\n")
  
  if(AllVarGreaterThanMeans){
    ## function to calculate negative log likelihood based specified distribution
    ans = 0
    for(ii in 1:length(listOfDataVecs.f)) {
      x.f <- 10^10; x.m <- 10^10
      try(x.f <- This.Dist.f(sapply(listOfDataVecs.f[ii], as.numeric), z.meanVec.f[ii], z.varVec.f[ii]) )
      try(x.m <- This.Dist.m(sapply(listOfDataVecs.m[ii], as.numeric), z.meanVec.m[ii], z.varVec.m[ii]) )
      if(is.na(x.f)) x.f <- 10^10
      if(is.na(x.m)) x.m <- 10^10
      ans = ans + x.f + x.m
    }
  }
  return(ans)
}


# default values y.varF = -1, y.psi = -1 and y.gamma = -1 indicate these parameters are unconstrained
# y.gamma.M = -999 means that gamma.M is constrained to be equal to gamma.F
# otherwise they are constrained to the specified values
randparmV5.BothSexes<- function(nPopTypes = 5, y.varF = -1, y.psi = -1, y.gamma.F = -1, y.gamma.M = -1) {
  needOne<- TRUE; countTries <-0; nMaxTries <-25;
  while(needOne & (countTries < nMaxTries)){
    z.means = runif(4, min = 1.9, max = 8)
    ifelse(y.varF == -1, z.VarF.LowHet <- runif(1, min = 1.9, max = 8), z.VarF.LowHet <- y.varF)
    ifelse(y.psi == -1, z.psi <- runif(1, min = 2.01/z.VarF.LowHet, max = 20.01/z.VarF.LowHet), z.psi <- y.psi)
    ifelse(y.gamma.F == -1, z.gamma.F <- runif(1, min = 0.5, max = 2.5), z.gamma.F <- y.gamma.F)
    ifelse(y.gamma.M == -1, z.gamma.M <- runif(1, min = 0.5, max = 2.5), z.gamma.M <- y.gamma.M)
    if(y.gamma.M == -999) z.gamma.M<-z.gamma.F
    zParmVec<- c(z.VarF.LowHet, z.psi, z.gamma.F, z.gamma.M, z.means)
    z.meanVec.f<-c(z.means[1], 0.75*z.means[1] + 0.25*z.means[2], 0.5*sum(z.means[1:2]), 0.25*z.means[1] + 0.75*z.means[2], z.means[2])
    z.meanVec.m<-c(z.means[3], 0.75*z.means[3] + 0.25*z.means[4], 0.5*sum(z.means[3:4]), 0.25*z.means[3] + 0.75*z.means[4], z.means[4])
    z.varVec.f.IfMean2<-c(z.VarF.LowHet, 0.5*(1+z.gamma.F)*z.VarF.LowHet, z.gamma.F*z.VarF.LowHet, 0.5*(1+z.gamma.F)*z.VarF.LowHet, z.VarF.LowHet)
    z.varVec.m.IfMean2<-z.psi*z.VarF.LowHet*c(1, 0.5*(1+z.gamma.M), z.gamma.M, 0.5*(1+z.gamma.M), 1)
    z.varVec.f<- z.varVec.f.IfMean2*z.meanVec.f/2
    z.varVec.m<- z.varVec.m.IfMean2*z.meanVec.m/2
    if(nPopTypes == 3){
      z.meanVec.f<-z.meanVec.f[c(1,3,5)]; z.meanVec.m<-z.meanVec.m[c(1,3,5)]
      z.varVec.f<-z.varVec.f[c(1,3,5)]; z.varVec.m<-z.varVec.m[c(1,3,5)]
    }
    
    z.AllVarGreaterThanMeans <- !any(c(z.meanVec.f, z.meanVec.m) >= c(z.varVec.f, z.varVec.m))
    if(z.AllVarGreaterThanMeans) needOne <- FALSE
    countTries<-countTries+1
    if(countTries == nMaxTries) zParmVec<-c(runif(1, min = 1.9, max = 2.5), runif(2, min = 0.9, max = 1.25),  runif(4, min = 1.9, max = 5))
  }
  return(zParmVec)
}


getMLEv5.5Models<-function(distID.f, distID.m, listOfDataVecs.f, listOfDataVecs.m, nStartingPlaces = 50){
  n<-nStartingPlaces ## number of starting places for optimization of likelihood
  nPopTypes = length(listOfDataVecs.f)
  
  nParm = length(randparmV5.BothSexes())
  
  #Model 1; psi = gamma = 1; (5 parameters: VarFLowHet and 4 means)
  fn.Mod1<-function(parmVec) NLLv5.BothSexes(distID.f, distID.m, listOfDataVecs.f, listOfDataVecs.m,  
                                             c(parmVec[1], 1, 1, 1, parmVec[2], parmVec[3], parmVec[4],parmVec[5]))
  results<-matrix(nrow = n, ncol = 2 + nParm + 2)
  for(i in 1:n){
    oneMLE = results[i,]; j = 0
    while(is.na(oneMLE[1]) & j < 20) {
      randStart<-randparmV5.BothSexes(nPopTypes = nPopTypes, y.psi = 1, y.gamma.F = 1, y.gamma.M = 1)[c(1,5:8)]
      jj = 0
      while(is.infinite(fn.Mod1(randStart)) & jj < 20) {randparmV5.BothSexes(nPopTypes = nPopTypes, y.psi = 1, y.gamma.F = 1, y.gamma.M = 1)[c(1,5:8)]; jj = jj + 1}
      try(oneMLE<- optim(randStart,fn.Mod1, method = "Nelder-Mead", control=list(maxit = 1000)))
      j = j+1}
    if(j < 20) results[i,]= c(distID.f, distID.m, oneMLE$par[1], 1,1,1, oneMLE$par[2:5], 5, oneMLE$value)
  }
  Model1<-results[order(results[,dim(results)[2]]),]
  cat("best model 1: \n", Model1[1,], "\n")
  
  # Model 2 - 6 free parameters: VarF.Low and psi and means; gamma = 1
  fn.Mod2<-function(parmVec) NLLv5.BothSexes(distID.f, distID.m, listOfDataVecs.f, listOfDataVecs.m,
                                             c(parmVec[1], parmVec[2], 1, 1, parmVec[3], parmVec[4],parmVec[5], parmVec[6]))
  results<-matrix(nrow = n, ncol = 2 + nParm + 2)
  for(i in 1:n){
    oneMLE = results[i,]; j = 0
    while(is.na(oneMLE[1]) & j < 20) {
      randStart<-randparmV5.BothSexes(nPopTypes = nPopTypes, y.gamma.F = 1, y.gamma.M = 1)[c(1:2, 5:8)]
      jj = 0
      while(is.infinite(fn.Mod2(randStart)) & jj < 20) {randStart<-randparmV5.BothSexes(nPopTypes = nPopTypes, y.gamma.F = 1, y.gamma.M = 1)[c(1:2, 5:8)]; jj = jj + 1}
      try(oneMLE<- optim(randStart,fn.Mod2, method = "Nelder-Mead", control=list(maxit = 1000)))
      j = j+1}
    if(j < 20) results[i,]= c(distID.f, distID.m, oneMLE$par[1:2], 1, 1, oneMLE$par[3:6], 6, oneMLE$value)
  }
  Model2<-results[order(results[,dim(results)[2]]),]
  cat("best model 2: \n", Model2[1,], "\n")
  
  # Model 3 - 6 free parameters: VarF.Low and gamma.F and means; psi = 1
  # gamma.M constrained to be equal to gamma.F
  fn.Mod3<-function(parmVec) NLLv5.BothSexes(distID.f, distID.m, listOfDataVecs.f, listOfDataVecs.m,
                                             c(parmVec[1], 1, parmVec[2], parmVec[2], parmVec[3], parmVec[4], parmVec[5], parmVec[6]))
  results<-matrix(nrow = n, ncol = 2 + nParm + 2)
  for(i in 1:n){
    oneMLE = results[i,]; j = 0
    while(is.na(oneMLE[1]) & j < 20) {
      randStart<-randparmV5.BothSexes(nPopTypes = nPopTypes, y.psi = 1, y.gamma.M = -999)[c(1,3, 5:8)]
      jj = 0
      while(is.infinite(fn.Mod3(randStart)) & jj < 20) {randStart<-randparmV5.BothSexes(nPopTypes = nPopTypes, y.psi = 1, y.gamma.M = -999)[c(1,3, 5:8)]; jj = jj + 1}
      try(oneMLE<- optim(randStart,fn.Mod3, method = "Nelder-Mead", control=list(maxit = 1000)))
      j = j+1}
    if(j < 20) results[i,]= c(distID.f, distID.m, oneMLE$par[1], 1, oneMLE$par[2], oneMLE$par[2], oneMLE$par[3:6], 6, oneMLE$value)
  }
  Model3<-results[order(results[,dim(results)[2]]),]
  cat("best model 3: \n", Model3[1,], "\n")
  
  # Model4 - 7 free parameters: VarF.Low, psi and gamma.F and means
  # gamma.M constrained to be equal to gamma.F
  fn.Mod4<-function(parmVec) NLLv5.BothSexes(distID.f, distID.m, listOfDataVecs.f, listOfDataVecs.m,
                                             c(parmVec[1], parmVec[2], parmVec[3], parmVec[3], parmVec[4], parmVec[5], parmVec[6], parmVec[7]))
  results<-matrix(nrow = n, ncol = 2+ nParm + 2)
  for(i in 1:n){
    oneMLE = results[i,]; j = 0
    while(is.na(oneMLE[1]) & j < 20) {
      randStart<-randparmV5.BothSexes(nPopTypes = nPopTypes,  y.gamma.M = -999)[c(1:3, 5:8)]
      jj = 0
      while(is.infinite(fn.Mod4(randStart)) & jj < 20) {randStart<-randparmV5.BothSexes(nPopTypes = nPopTypes,  y.gamma.M = -999)[c(1:3, 5:8)]; jj = jj + 1}
      try(oneMLE<- optim(randStart,fn.Mod4, method = "Nelder-Mead", control=list(maxit = 1000)))
      j = j+1}
    if(j < 20) results[i,]= c(distID.f, distID.m, oneMLE$par[1:3], oneMLE$par[3], oneMLE$par[4:7], 7, oneMLE$value)
  }
  Model4<-results[order(results[,dim(results)[2]]),]
  cat("best model 4: \n", Model4[1,], "\n")
  
  
  # Model5 - 8 free parameters: VarF.Low, psi and gamma and means
  
  fn.Mod5<-function(parmVec) NLLv5.BothSexes(distID.f, distID.m, listOfDataVecs.f, listOfDataVecs.m,
                                             c(parmVec[1], parmVec[2], parmVec[3], parmVec[4], parmVec[5], parmVec[6], parmVec[7], parmVec[8]))
  results<-matrix(nrow = n, ncol = 2+ nParm + 2)
  for(i in 1:n){
    oneMLE = results[i,]; j = 0
    while(is.na(oneMLE[1]) & j < 20) {
      randStart<-randparmV5.BothSexes(nPopTypes = nPopTypes)
      jj = 0
      while(is.infinite(fn.Mod4(randStart)) & jj < 20) {randStart<-randparmV5.BothSexes(nPopTypes = nPopTypes); jj = jj + 1}
      try(oneMLE<- optim(randStart,fn.Mod5, method = "Nelder-Mead", control=list(maxit = 1000)))
      j = j+1}
    if(j < 20) results[i,]= c(distID.f, distID.m, oneMLE$par[1:8], 8, oneMLE$value)
  }
  Model5<-results[order(results[,dim(results)[2]]),]
  cat("best model 5: \n", Model5[1,], "\n")
  
  nbest<-5 #This many of the top results from each model
  return(rbind(Model1[1:nbest,], Model2[1:nbest,], Model3[1:nbest,], Model4[1:nbest,], Model5[1:nbest,]))
}

FiveModels.Vial.F2M2<-getMLEv5.5Models(2, 2, List.Vial.f, List.Vial.m); FiveModels.Vial.F2M2
FiveModels.Cage.F2M2<-getMLEv5.5Models(2, 2, List.Cage.f, List.Cage.m); FiveModels.Cage.F2M2
FiveModels.Monog.F2M2<-getMLEv5.5Models(2, 2, List.Monog.f, List.Monog.m); FiveModels.Monog.F2M2
# 
# save(FiveModels.Vial.F2M2, FiveModels.Cage.F2M2, FiveModels.Monog.F2M2, file = "/Users/aneil/Google\ Drive/AgrawalLab/Amardeep/SexSpecificVariances/ResultsOf5LikihoodModels.Rdata")
load(file = "/Users/aneil/Google\ Drive/AgrawalLab/Amardeep/SexSpecificVariances/ResultsOf4LikihoodModels.Rdata")


####################################
#############################



#############
########################
################################
######################################
######
## Here we simulate two different ways that focal males could "compete" for fertilizations 
## against competitors and investigate the consequences of inferior competitors on the 
## observed variance among focal males.  We find that when competitors are inferior
## the variance among focal males is exaggerated compared to what it would (should) be if 
## competitors are on average as good as focal males.
######
######################################
################################
########################
#############


### Both models assume a single focal male competes against 16 competitor males for 16 females
### A competitor has an expected fitness of (1 - sComp) of the focal individuals.
### 32 offspring are sampled from each vial

### In the first model, we assume strong last male precedence. 
### Males compete to be the last male to mate with each female.
### This model will tend to have high variance among males because there are only 16 potential
### mating opportunities.
### The argument "nfocal" is the number of focal individuals being measured
### The argument "pv" controls the variance in "quality" among focal males (choose values 5-100).
### i.e., individuals of higher/lower quality will tended to be more/less successful in producing offspring
### Returns: 1) selective disadvantage of competitor; 2) average offspring from focal individual
### 3) variance in offspring among focal individual; 4) variance adjusted by mean using 

## I use this mixed poisson distribution to create overdispersed female fecundity
randMixedPois<-function(n, pHigh, lambdaH, delta){
  lambdaL = lambdaH - delta
  #x<- ifelse(runif(1) < pHigh, rpois(1, lambdaH), rpois(1, lambdaL))
  zz<-sapply(1:n, function(x) ifelse(runif(1) < pHigh, rpois(1, lambdaH), rpois(1, lambdaL)) )
  return(zz)
}

bella<-randMixedPois(5000, 0.5, 100, 30); c(mean(bella), var(bella))

StrongLastMalePrecedence<-function(sComp, pv, nfocal = 5000){
  focalMaleIntrinsic<- rpois(nfocal, pv)
  compMaleIntrinsic<- rpois(5000, (1-sComp)*pv)
  nFemalesPerFocalMale<-sapply(focalMaleIntrinsic, function(x) rbinom(1, 16, prob = x/(x+ sum(sample(compMaleIntrinsic, 15)) + 10^-6 )))
  femaleFec<- randMixedPois(5000, 0.5, 100, 30)
  nOffspringPerFocalMale<- sapply(nFemalesPerFocalMale, function(x) sum(sample(femaleFec, x)))
  nOffspringCompMales<- sapply(nFemalesPerFocalMale, function(x) sum(sample(femaleFec, 16 -x)))
  functToSampleOffspring<-function(focalID){
    if(nOffspringPerFocalMale[focalID] + nOffspringCompMales[focalID] < 32) zz2<-NA;
    if(nOffspringPerFocalMale[focalID]==0 & nOffspringCompMales[focalID] > 0) zz2<-0;
    if(nOffspringPerFocalMale[focalID] > 0 & nOffspringCompMales[focalID] == 0) zz2<-32;
    if(nOffspringPerFocalMale[focalID] > 0 & nOffspringCompMales[focalID] > 0) zz2<-sum(sample( c(rep(1,nOffspringPerFocalMale[focalID]), rep(0,nOffspringCompMales[focalID])), 32))
    return(zz2)
  }
  nFocalOffspring<-sapply(1:nfocal, functToSampleOffspring)
  #nFocalOffspring<-sapply(nFemalesPerFocalMale, function(x) rbinom(1, 32, prob = x/16))
  # varRelIntrinsic<-var(focalMaleIntrinsic/(mean(focalMaleIntrinsic)))
  zbar<-mean(nFocalOffspring[1:nfocal])
  zVar<-var(nFocalOffspring[1:nfocal])
  return(c(sComp, zbar, zVar, zVar*(2/zbar)))
}


StrongLastMalePrecedence(0, 30, 5000)
StrongLastMalePrecedence(0.4, 30, 5000)



### In the second model, we assume high mating rates and no sperm precedence. 
### Essentially the model assumes there is a separate competition among males  
### to fertilize each individual egg.
### This model will tend to have lower variance. 
EveryEggAsSingleTrial<-function(sComp, pv, nfocal = 5000){
  focalMaleIntrinsic<- rpois(max(nfocal, 5000), pv)
  nZygotesPerFocalMale<-sapply(focalMaleIntrinsic, function(x) rbinom(1, 32, prob = x/(x+ (1-sComp)*sum(sample(focalMaleIntrinsic, 15)) )))
  zbar<-mean(nZygotesPerFocalMale[1:nfocal])
  zVar<-var(nZygotesPerFocalMale[1:nfocal])
  # varRelIntrinsic<-var(focalMaleIntrinsic/(mean(focalMaleIntrinsic)))
  return(c(sComp, zbar, zVar, zVar*(2/zbar)))
}

### In the third model, we assume that each female mates with k males.
### Males compete for matings with females.
### The k mates of a female compete to fertilize her zygotes
FinitePolyandry<-function(sComp, pv, nMatesPerFemale = 3, nfocal = 5000){
  nMatesPerFemale.dummyProof = min(c(nMatesPerFemale, 16))
  pvMating = pv; pvPostCopuatory = pv
  sCompMating = sComp; sCompPostCopulatory = 0;
  #sCompMating = 0; sCompPostCopulatory = sComp;
  #sCompMating = 1 - sqrt(1 - sComp); sCompPostCopulatory = 1 - sqrt(1 - sComp);
  nZygotesPerFocalMale<-rep(0, nfocal)
  nZygotesFromCompetitorsOfThisFocalMale<-rep(0, nfocal)
  for(focalID in 1:nfocal){
    MatingAptitudes<- c(rpois(1, pvMating), rpois(15, (1-sCompMating)*pvMating)) + 0.0001 ## adding a small amount to prevent situations where there are no males left with any mating aptitutde
    PostCopulatoryAptitudes<- c(rpois(1, pvPostCopuatory), rpois(15, (1-sCompPostCopulatory)*pvPostCopuatory))
    focalzygotes = 0
    for(femaleID in 1:16){
      winners<-rep(0,16)
      for(mateNumber in 1:nMatesPerFemale.dummyProof){
        TheseMatingAptitudes = MatingAptitudes * (1 - winners) 
        TheseMatingProbs = TheseMatingAptitudes/sum(TheseMatingAptitudes)
        winners<- winners + as.vector(rmultinom(1, 1, prob = TheseMatingProbs))
      }
      femaleFec<- randMixedPois(1, 0.5, 100, 30)
      ThesePostCopAptitudes = PostCopulatoryAptitudes * winners
      focalzygotesFromCurrentFemale = rbinom(1, femaleFec, prob = ThesePostCopAptitudes[1]/sum(ThesePostCopAptitudes))
      compzygotesFromCurrentFemale = femaleFec - focalzygotesFromCurrentFemale
      nZygotesPerFocalMale[focalID] =  nZygotesPerFocalMale[focalID] + focalzygotesFromCurrentFemale
      nZygotesFromCompetitorsOfThisFocalMale[focalID] = nZygotesFromCompetitorsOfThisFocalMale[focalID] + compzygotesFromCurrentFemale
    }
  }
  functToSampleOffspring<-function(focalIDx){
    if(nZygotesPerFocalMale[focalIDx] + nZygotesFromCompetitorsOfThisFocalMale[focalIDx] < 32) zz2<-NA;
    if(nZygotesPerFocalMale[focalIDx]==0 & nZygotesFromCompetitorsOfThisFocalMale[focalIDx] > 0) zz2<-0;
    if(nZygotesPerFocalMale[focalIDx] > 0 & nZygotesFromCompetitorsOfThisFocalMale[focalIDx] == 0) zz2<-32;
    if(nZygotesPerFocalMale[focalIDx] > 0 & nZygotesFromCompetitorsOfThisFocalMale[focalIDx] > 0) zz2<-sum(sample( c(rep(1,nZygotesPerFocalMale[focalIDx]), rep(0,nZygotesFromCompetitorsOfThisFocalMale[focalIDx])), 32))
    return(zz2)
  }
  nFocalOffspring<-sapply(1:nfocal, functToSampleOffspring)
  zbar<-mean(nFocalOffspring[1:nfocal])
  zVar<-var(nFocalOffspring[1:nfocal])
  return(c(sComp, zbar, zVar, zVar*(2/zbar)))
}

FinitePolyandry(0, 10, nMatesPerFemale = 3); FinitePolyandry(0.4, 10, nMatesPerFemale = 3)
FinitePolyandry(0, 10, nMatesPerFemale = 8); FinitePolyandry(0.4, 10, nMatesPerFemale = 8)
FinitePolyandry(0, 10, nMatesPerFemale = 16); FinitePolyandry(0.4, 10, nMatesPerFemale = 16)


### This runs a set of simulations with different values of sComp
GetASimulationSet<-function(MateCompModel, pv, nfocal = 5000){
  z<-t(sapply(0.1*(0:7), MateCompModel, pv = pv, nfocal = nfocal))
  z<-as.data.frame(z)
  names(z)<-c("sComp", "mean", "var", "varAdjByMean")
  z$ObsVarRelVarWithEqualComps<-z[,3]/z[1,3]
  z$AdjVarRelVarWithEqualComps<-z[,4]/z[1,3]
  z$colorByMean<-rep("orange", dim(z)[1])
  for(i in 1:(dim(z)[1])){
    if(z[i,2] < 5) z$colorByMean[i]<-"red"
    if(z[i,2] < 4) z$colorByMean[i]<-"blue"
    if(z[i,2] < 3) z$colorByMean[i]<-"purple"
    if(z[i,2] < 2.2) z$colorByMean[i]<-"black"
  }
  return(z)
}

adjustTable<-function(data){
  z<-as.data.frame(data)
  z$ObsToBaseLine<-z[,3]/z[1,3]
  z$AdjToBaseLine<-z[,4]/z[1,3]
  z$colorByMean<-rep("orange", dim(data)[1])
  for(i in 1:(dim(data)[1])){
    if(z[i,2] < 5) z$colorByMean[i]<-"red"
    if(z[i,2] < 4) z$colorByMean[i]<-"blue"
    if(z[i,2] < 3) z$colorByMean[i]<-"purple"
    if(z[i,2] < 2.2) z$colorByMean[i]<-"black"
  }
  return(z)
}



a1v2<-GetASimulationSet(StrongLastMalePrecedence, 5, nfocal = 20000)
a1v2[,1:6]
b2v2<-GetASimulationSet(EveryEggAsSingleTrial, 10, nfocal = 20000)
b2v2[,1:6]
b2[,1:6]

a1<-GetASimulationSet(StrongLastMalePrecedence, 5)
a2<-GetASimulationSet(StrongLastMalePrecedence, 10)
a3<-GetASimulationSet(StrongLastMalePrecedence, 50)
b1<-GetASimulationSet(EveryEggAsSingleTrial, 5)
b2<-GetASimulationSet(EveryEggAsSingleTrial, 10)
b3<-GetASimulationSet(EveryEggAsSingleTrial, 50)
c1<-GetASimulationSet(FinitePolyandry, 5)
c2<-GetASimulationSet(FinitePolyandry, 10)
c3<-GetASimulationSet(FinitePolyandry, 50)

allResults<-rbind(a1, a2, a3, b1, b2, b3, c1, c2, c3)
#hist(allResults$AdjVarRelVarWithEqualComps)

### This plots the variance relative to the variance under that would be observed
### under identical circumstances but with competitors that have the same average fitness
### Open symbols are using the "observed" variance
### Closed symbols are using observed variance adjusted by the observed mean
### Triangles are from the "every egg as a single trial" model
### Circles are from the "strong last male precedence" model
### Squares are from the "finite polyandry" model
### Points are colured by the mean number offpsring: 
### zbar < 2.2 -- black
### 2.2 < zbar < 3 -- purple
### 4 < zbar < 4 -- blue
### 5 < zbar < 6 -- red
### 6 < zbar -- orange
### This plot shows that for both models of male competition:
### (1) the observed variance of focal males increases 
### with the average fitness of focal individuals relative to competitors
### i.e., higher variance with lighter colours
### (2) adjusting the observed variance by the observed mean
### works well for estimating what the variance would have been with equivalent competitors
### i.e., the ratio is close to 1 for the closed symbols (but not the open symbols)
plot(rep(allResults[,3], 2), c(allResults[,5],allResults[,6]), col = rep(allResults$colorByMean, 2), 
     pch = c(rep(1, length(allResults[,5])/3), rep(2, length(allResults[,5])/3), rep(0, length(allResults[,5])/3), rep(19, length(allResults[,6])/2),rep(17, length(allResults[,6])/2), rep(15, length(allResults[,6])/2)), 
     ylim = c(0, 0.1+max(allResults[,5])),
     xlab = "Observed Variance", ylab = "Ratio",
     main = "Ratio of Var relative to Var with Equal Competitors");abline(a= 1, b=0, lty = "dashed")



###### Simulate data for analysis with likelihood framework



Females.MakeDataToTestWithLikelihoodAnalysis<-function(sComp, pH = 0.5, lambdaH = 100, delta = 30, nfocal = 5000){
  focalFemaleFecundity<- randMixedPois(nfocal, pH, lambdaH, delta)
  competitorFecundity<-randMixedPois(max(nfocal, 5000), pH, (1-sComp)*lambdaH, (1-sComp)*delta)
  nFocalOffspring<-sapply(focalFemaleFecundity, function(x) rbinom(1, 32, prob = x/(x+ sum(sample(competitorFecundity, 15)) )))
  return(nFocalOffspring)
}



StrongLastMalePrecedence.MakeDataToTestWithLikelihoodAnalysis<-function(sComp, pv, pH = 0.5, lambdaH = 100, delta = 30, nfocal = 5000){
  focalMaleIntrinsic<- rpois(nfocal, pv)
  compMaleIntrinsic<- rpois(5000, (1-sComp)*pv)
  nFemalesPerFocalMale<-sapply(focalMaleIntrinsic, function(x) rbinom(1, 16, prob = x/(x+ sum(sample(compMaleIntrinsic, 15)) + 10^-6 )))
  femaleFec<- randMixedPois(5000, pH, lambdaH, delta)
  nOffspringPerFocalMale<- sapply(nFemalesPerFocalMale, function(x) sum(sample(femaleFec, x)))
  nOffspringCompMales<- sapply(nFemalesPerFocalMale, function(x) sum(sample(femaleFec, 16 -x)))
  functToSampleOffspring<-function(focalID){
    if(nOffspringPerFocalMale[focalID] + nOffspringCompMales[focalID] < 32) zz2<-NA;
    if(nOffspringPerFocalMale[focalID]==0 & nOffspringCompMales[focalID] > 0) zz2<-0;
    if(nOffspringPerFocalMale[focalID] > 0 & nOffspringCompMales[focalID] == 0) zz2<-32;
    if(nOffspringPerFocalMale[focalID] > 0 & nOffspringCompMales[focalID] > 0) zz2<-sum(sample( c(rep(1,nOffspringPerFocalMale[focalID]), rep(0,nOffspringCompMales[focalID])), 32))
    return(zz2)
  }
  nFocalOffspring<-sapply(1:nfocal, functToSampleOffspring)
  return(nFocalOffspring)
}

## Note that pH, lambdaH, and delta aren't used but they are listed as arguments so "DoATest" will work
EveryEggAsSingleTrial.MakeDataToTestWithLikelihoodAnalysis<-function(sComp, pv, pH = 0.5, lambdaH = 100, delta = 30, nfocal = 5000){
  focalMaleIntrinsic<- rpois(max(nfocal, 5000), pv)
  nZygotesPerFocalMale<-sapply(focalMaleIntrinsic, function(x) rbinom(1, 32, prob = x/(x+ (1-sComp)*sum(sample(focalMaleIntrinsic, 15)) )))
  return(nZygotesPerFocalMale)
}


FinitePolyandry.MakeDataToTestWithLikelihoodAnalysis<-function(sComp, pv, nMatesPerFemale = 2, pH = 0.5, lambdaH = 100, delta = 30, nfocal = 5000){
  nMatesPerFemale.dummyProof = min(c(nMatesPerFemale, 16))
  pvMating = pv; pvPostCopuatory = pv
  sCompMating = sComp; sCompPostCopulatory = 0;
  #sCompMating = 0; sCompPostCopulatory = sComp;
  #sCompMating = 1 - sqrt(1 - sComp); sCompPostCopulatory = 1 - sqrt(1 - sComp);
  nZygotesPerFocalMale<-rep(0, nfocal)
  nZygotesFromCompetitorsOfThisFocalMale<-rep(0, nfocal)
  for(focalID in 1:nfocal){
    MatingAptitudes<- c(rpois(1, pvMating), rpois(15, (1-sCompMating)*pvMating)) + 0.0001 ## adding a small amount to prevent situations where there are no males left with any mating aptitutde
    PostCopulatoryAptitudes<- c(rpois(1, pvPostCopuatory), rpois(15, (1-sCompPostCopulatory)*pvPostCopuatory))
    focalzygotes = 0
    for(femaleID in 1:16){
      winners<-rep(0,16)
      for(mateNumber in 1:nMatesPerFemale.dummyProof){
        TheseMatingAptitudes = MatingAptitudes * (1 - winners) 
        TheseMatingProbs = TheseMatingAptitudes/sum(TheseMatingAptitudes)
        winners<- winners + as.vector(rmultinom(1, 1, prob = TheseMatingProbs))
      }
      femaleFec<- randMixedPois(1, pH, lambdaH, delta)
      ThesePostCopAptitudes = PostCopulatoryAptitudes * winners
      focalzygotesFromCurrentFemale = rbinom(1, femaleFec, prob = ThesePostCopAptitudes[1]/sum(ThesePostCopAptitudes))
      compzygotesFromCurrentFemale = femaleFec - focalzygotesFromCurrentFemale
      nZygotesPerFocalMale[focalID] =  nZygotesPerFocalMale[focalID] + focalzygotesFromCurrentFemale
      nZygotesFromCompetitorsOfThisFocalMale[focalID] = nZygotesFromCompetitorsOfThisFocalMale[focalID] + compzygotesFromCurrentFemale
    }
  }
  functToSampleOffspring<-function(focalIDx){
    if(nZygotesPerFocalMale[focalIDx] + nZygotesFromCompetitorsOfThisFocalMale[focalIDx] < 32) zz2<-NA;
    if(nZygotesPerFocalMale[focalIDx]==0 & nZygotesFromCompetitorsOfThisFocalMale[focalIDx] > 0) zz2<-0;
    if(nZygotesPerFocalMale[focalIDx] > 0 & nZygotesFromCompetitorsOfThisFocalMale[focalIDx] == 0) zz2<-32;
    if(nZygotesPerFocalMale[focalIDx] > 0 & nZygotesFromCompetitorsOfThisFocalMale[focalIDx] > 0) zz2<-sum(sample( c(rep(1,nZygotesPerFocalMale[focalIDx]), rep(0,nZygotesFromCompetitorsOfThisFocalMale[focalIDx])), 32))
    return(zz2)
  }
  nFocalOffspring<-sapply(1:nfocal, functToSampleOffspring)
  return(nFocalOffspring)
}


getTheBasics(StrongLastMalePrecedence.MakeDataToTestWithLikelihoodAnalysis(0, 50, nfocal = 2000))
getTheBasics(FinitePolyandry.MakeDataToTestWithLikelihoodAnalysis(0, 30, nMatesPerFemale = 2, nfocal = 5000))
getTheBasics(FinitePolyandry.MakeDataToTestWithLikelihoodAnalysis(0.1, 30, nMatesPerFemale = 1, nfocal = 5000))
getTheBasics(FinitePolyandry.MakeDataToTestWithLikelihoodAnalysis(0.3, 30, nMatesPerFemale = 1, nfocal = 5000))
getTheBasics(FinitePolyandry.MakeDataToTestWithLikelihoodAnalysis(0.1, 30, nMatesPerFemale = 2, nfocal = 5000))
getTheBasics(FinitePolyandry.MakeDataToTestWithLikelihoodAnalysis(0.3, 30, nMatesPerFemale = 2, nfocal = 5000))
getTheBasics(FinitePolyandry.MakeDataToTestWithLikelihoodAnalysis(0.1, 30, nMatesPerFemale = 3, nfocal = 5000))
getTheBasics(FinitePolyandry.MakeDataToTestWithLikelihoodAnalysis(0.3, 30, nMatesPerFemale = 3, nfocal = 5000))
getTheBasics(FinitePolyandry.MakeDataToTestWithLikelihoodAnalysis(0.1, 30, nMatesPerFemale = 5, nfocal = 5000))
getTheBasics(FinitePolyandry.MakeDataToTestWithLikelihoodAnalysis(0.3, 30, nMatesPerFemale = 5, nfocal = 5000))

getTheBasics(FinitePolyandry.MakeDataToTestWithLikelihoodAnalysis(0, 30, nMatesPerFemale = 3, nfocal = 5000))
getTheBasics(FinitePolyandry.MakeDataToTestWithLikelihoodAnalysis(0.3, 5, nMatesPerFemale = 2, nfocal = 5000))



DoATest<-function(DistID, MaleCompModID, sCompM, IntraVarM.X, pH.X, lambdaH.X, delta.X, nrep = 1000){
  if(MaleCompModID == 1) MaleCompMod = StrongLastMalePrecedence.MakeDataToTestWithLikelihoodAnalysis
  if(MaleCompModID == 2) MaleCompMod = EveryEggAsSingleTrial.MakeDataToTestWithLikelihoodAnalysis
  if(MaleCompModID == 3) MaleCompMod = FinitePolyandry.MakeDataToTestWithLikelihoodAnalysis
  test.PopA.f <- Females.MakeDataToTestWithLikelihoodAnalysis(0, pH = pH.X, lambdaH = lambdaH.X, delta = delta.X, nfocal = nrep);
  test.PopB.f <- Females.MakeDataToTestWithLikelihoodAnalysis(0, pH = pH.X, lambdaH = lambdaH.X, delta = delta.X, nfocal = nrep);
  test.PopC.f <- Females.MakeDataToTestWithLikelihoodAnalysis(0, pH = pH.X, lambdaH = lambdaH.X, delta = delta.X, nfocal = nrep);
  test.PopD.f <- Females.MakeDataToTestWithLikelihoodAnalysis(0, pH = pH.X, lambdaH = lambdaH.X, delta = delta.X, nfocal = nrep);
  test.PopE.f <- Females.MakeDataToTestWithLikelihoodAnalysis(0, pH = pH.X, lambdaH = lambdaH.X, delta = delta.X, nfocal = nrep);
  test.PopA.m <- MaleCompMod(sCompM, pv = IntraVarM.X, pH = pH.X, lambdaH = lambdaH.X, delta = delta.X, nfocal = nrep);
  test.PopB.m <- MaleCompMod(sCompM, pv = IntraVarM.X, pH = pH.X, lambdaH = lambdaH.X, delta = delta.X, nfocal = nrep);
  test.PopC.m <- MaleCompMod(sCompM, pv = IntraVarM.X, pH = pH.X, lambdaH = lambdaH.X, delta = delta.X, nfocal = nrep);
  test.PopD.m <- MaleCompMod(sCompM, pv = IntraVarM.X, pH = pH.X, lambdaH = lambdaH.X, delta = delta.X, nfocal = nrep);
  test.PopE.m <- MaleCompMod(sCompM, pv = IntraVarM.X, pH = pH.X, lambdaH = lambdaH.X, delta = delta.X, nfocal = nrep);
  
  cat("test.PopA.m: length, mean and var: ", c(length(test.PopA.m), mean(test.PopA.m), var(test.PopA.m)), "\n")
  
  List.test.f = list(test.PopA.f,test.PopB.f, test.PopC.f, test.PopD.f, test.PopE.f)
  List.test.m = list(test.PopA.m,test.PopB.m, test.PopC.m, test.PopD.m, test.PopE.m)
  
  
  summary.f<- getTheBasicsForAllPopTypes(List.test.f)
  summary.m<- getTheBasicsForAllPopTypes(List.test.m)
  xMLEsFromTheFourModels.TestData<- getMLEv4.4Models(DistID, DistID, List.test.f, List.test.m, nStartingPlaces = 25)
  adjustToJustGetBest<-dim(xMLEsFromTheFourModels.TestData)[1]/4
  yMLEsFromTheFourModels.TestData<-xMLEsFromTheFourModels.TestData[1 + adjustToJustGetBest*c(0:3),]
  bestMLEsFromTheFourModels.TestData<-xMLEsFromTheFourModels.TestData[order(xMLEsFromTheFourModels.TestData[,11]),][1,]
  cat(c(MaleCompModID, sCompM, IntraVarM.X, pH.X, lambdaH.X, delta.X, nrep), "\n")
  cat("likelihood \n")
  cat(yMLEsFromTheFourModels.TestData[1,], "\n")
  cat(yMLEsFromTheFourModels.TestData[2,], "\n")
  cat(yMLEsFromTheFourModels.TestData[3,], "\n")
  cat(yMLEsFromTheFourModels.TestData[4,], "\n")
  cat(bestMLEsFromTheFourModels.TestData, "\n")
  cat("females \n")
  cat(summary.f[1,], "\n")
  cat(summary.f[2,], "\n")
  cat(summary.f[3,], "\n")
  cat("males \n")
  cat(summary.m[1,], "\n")
  cat(summary.m[2,], "\n")
  cat(summary.m[3,], "\n")
  return(c(DistID, MaleCompModID, sCompM, IntraVarM.X, pH.X, lambdaH.X, delta.X, nrep, bestMLEsFromTheFourModels.TestData, as.vector(summary.f), as.vector(summary.m)))
  }


testResults.Dist1.MC1.v1<- DoATest(1, 1, 0, 10, 0.5, 50, 15, nrep = 1000)
testResults.Dist2.MC1.v1<- DoATest(2, 1, 0, 10, 0.5, 50, 15, nrep = 1000)
testResults.Dist1.MC1.v2<- DoATest(1, 1, 0, 10, 0.5, 42.5, 0, nrep = 1000)
testResults.Dist2.MC1.v2<- DoATest(2, 1, 0, 10, 0.5, 42.5, 0,nrep = 1000)

testResults.Dist1.MC2.v1<- DoATest(1, 2, 0, 10, 0.5, 50, 15, nrep = 1000)
testResults.Dist2.MC2.v1<- DoATest(2, 2, 0, 10, 0.5, 50, 15, nrep = 1000)
testResults.Dist1.MC2.v2<- DoATest(1, 2, 0, 10, 0.5, 42.5, 0, nrep = 1000)
testResults.Dist2.MC2.v2<- DoATest(2, 2, 0, 10, 0.5, 42.5, 0,nrep = 1000)

testResults.Dist1.MC3.v1<- DoATest(1, 3, 0, 10, 0.5, 50, 15, nrep = 1000)
testResults.Dist2.MC3.v1<- DoATest(2, 3, 0, 10, 0.5, 50, 15, nrep = 1000)
testResults.Dist1.MC3.v2<- DoATest(1, 3, 0, 10, 0.5, 42.5, 0, nrep = 1000)
testResults.Dist2.MC3.v2<- DoATest(2, 3, 0, 10, 0.5, 42.5, 0,nrep = 1000)




sComp.x = 0.35
testResults.Dist1.MC1.v1.poorComp<- DoATest(1, 1, sComp.x, 10, 0.5, 50, 15, nrep = 1000)
testResults.Dist2.MC1.v1.poorComp<- DoATest(2, 1, sComp.x, 10, 0.5, 50, 15, nrep = 1000)
testResults.Dist1.MC1.v2.poorComp<- DoATest(1, 1, sComp.x, 10, 0.5, 42.5, 0, nrep = 1000)
testResults.Dist2.MC1.v2.poorComp<- DoATest(2, 1, sComp.x, 10, 0.5, 42.5, 0,nrep = 1000)

testResults.Dist1.MC2.v1.poorComp<- DoATest(1, 2, sComp.x, 10, 0.5, 50, 15, nrep = 1000)
testResults.Dist2.MC2.v1.poorComp<- DoATest(2, 2, sComp.x, 10, 0.5, 50, 15, nrep = 1000)
testResults.Dist1.MC2.v2.poorComp<- DoATest(1, 2, sComp.x, 10, 0.5, 42.5, 0, nrep = 1000)
testResults.Dist2.MC2.v2.poorComp<- DoATest(2, 2, sComp.x, 10, 0.5, 42.5, 0,nrep = 1000)


testResults.Dist1.MC3.v1.poorComp<- DoATest(1, 3, sComp.x, 10, 0.5, 50, 15, nrep = 1000)
testResults.Dist2.MC3.v1.poorComp<- DoATest(2, 3, sComp.x, 10, 0.5, 50, 15, nrep = 1000)
testResults.Dist1.MC3.v2.poorComp<- DoATest(1, 3, sComp.x, 10, 0.5, 42.5, 0, nrep = 1000)
testResults.Dist2.MC3.v2.poorComp<- DoATest(2, 3, sComp.x, 10, 0.5, 42.5, 0,nrep = 1000)

#nmates = 3
testResults.Dist2.MC1.v3<- DoATest(2, 1, 0, 10, 0.5, 42.5 +40/2, 40, nrep = 1000)
testResults.Dist2.MC2.v3<- DoATest(2, 2, 0, 10, 0.5, 42.5 +40/2, 40, nrep = 1000)
testResults.Dist2.MC3.v3<- DoATest(2, 3, 0, 10, 0.5, 42.5 +40/2, 40, nrep = 1000)
testResults.Dist2.MC1.v3.poorComp<- DoATest(2, 1, sComp.x, 10, 0.5, 42.5 +40/2, 40, nrep = 1000)
testResults.Dist2.MC2.v3.poorComp<- DoATest(2, 2, sComp.x, 10, 0.5, 42.5 +40/2, 40, nrep = 1000)
testResults.Dist2.MC3.v3.poorComp<- DoATest(2, 3, sComp.x, 10, 0.5, 42.5 +40/2, 40, nrep = 1000)

#nmates = 2
testResults.Dist2.MC3.v3b<- DoATest(2, 3, 0, 10, 0.5, 42.5 +40/2, 40, nrep = 1000)
testResults.Dist2.MC3.v3b.poorComp<- DoATest(2, 3, sComp.x, 10, 0.5, 42.5 +40/2, 40, nrep = 1000)



displayTestResults<-function(test){
  cat("DistID = ", test[1], "; MCmod = ", test[2], "; sComp = ", test[3], "\n")
  cat("{IntraVarM.X, pH.X, lambdaH.X, delta.X} = ", test[4:7], "\n")
  cat("female summary: \n", test[20:24], "\n", test[25:29], "\n", test[30:34], "\n")
  cat("male summary: \n", test[35:39], "\n", test[40:44], "\n", test[45:49], "\n")
  cat("MLE: {varF.LowHet, psi, gamma} : ", test[11:13], "\n MLE: means x4: \n", test[14:17], "\n")
  cat("NLL: ", test[19], "\n")
}
displayTestResults(testResults.Dist1.MC1.v1)
displayTestResults(testResults.Dist1.MC1.v1.poorComp)

displayTestResults(testResults.Dist2.MC1.v1)
displayTestResults(testResults.Dist2.MC1.v1.poorComp)

displayTestResults(testResults.Dist1.MC2.v1)
displayTestResults(testResults.Dist1.MC2.v1.poorComp)

displayTestResults(testResults.Dist2.MC2.v1)
displayTestResults(testResults.Dist2.MC2.v1.poorComp)


displayTestResults(testResults.Dist1.MC1.v2)
displayTestResults(testResults.Dist1.MC1.v2.poorComp)

displayTestResults(testResults.Dist2.MC1.v2)
displayTestResults(testResults.Dist2.MC1.v2.poorComp)

displayTestResults(testResults.Dist1.MC2.v2)
displayTestResults(testResults.Dist1.MC2.v2.poorComp)

displayTestResults(testResults.Dist2.MC2.v2)
displayTestResults(testResults.Dist2.MC2.v2.poorComp)

