# Sex Specific Variance and the Efficacy of Selection
This repository contains code and data underlying the figures and analyses associated with Singh, A. and Agrawal, A.F., 2020. Sex Specific Variance and the Efficacy of Selection

## Code 
Code is all written in R. Coding scripts are ordered roughly in the order that the data are presented in the paper. All figures and tables in the paper and supplement are replicable (with the code provided). NOTE: Although the data in the figures is identical to that in the paper + supplement, the figures will not appear exactly the same as I processed figures after generating the basic plot in R (Fonts, labels, panels etc.).

<i>Be sure to read the info at the top of each script. There are several places where you will need to specify a file path either to a data file or for saving outputs</i>

## Data
Raw data is included in csv files. There are three files in total:
  1. [FitnessDataNeprojectFinalSept2020.csv](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/FitnessDataNeprojectFinalSept2020.csv): Raw fitness data from the experiment 
  2. [BodyMassDataNeprojectFinalSept2020.csv](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/BodyMassDataNeprojectFinalSept2020.csv): Body mass data from males and females
  3. [BootstrappedFitnessDataFinalSept2020.csv](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/BootstrappedFitnessDataFinalSept2020.csv): Results from a bootstrapping procedure wherein estimates of the mean and variance in sex-specific fitness was bootstrapped 10,000 times.
  
Metadata describing what each column in each data file mean is described in ["MetaDataForAllDataFiles.txt"](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/MetaDataForAllDataFiles.txt)
