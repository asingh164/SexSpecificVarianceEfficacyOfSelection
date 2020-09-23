# Sex Specific Variance and the Efficacy of Selection
This repository contains code and data underlying the figures and analyses associated with Singh, A. and Agrawal, A.F., 2020. Sex Specific Variance and the Efficacy of Selection

## Code 
Code is all written in R. Coding scripts are ordered roughly in the order that the data are presented in the paper. All figures and tables in the paper and supplement are replicable (with the code provided). NOTE: Although the data in the figures is identical to that in the paper + supplement, the figures will not appear exactly the same as I processed figures after generating the basic plot in R (Fonts, labels, panels etc.).

<i>Be sure to read the info at the top of each script. There are several places where you will need to specify a file path either to a data file or for saving outputs</i>

NOTE: These scripts are some of the first that I had written in R and as such they may not be as pretty, efficient or succinct as they could be.  ¯\_(ツ)_/¯

There are 6 R scripts in total: 
1. [1_TheoreticalModelAnalysis_ForSubmission.R](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/1_TheoreticalModelAnalysis_ForSubmission.R): Script that produces figures to plot results from a simple theoretical model presented in the paper.
2. [2_BodySizeAnalysis_ForSubmission.R](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/2_BodySizeAnalysis_ForSubmission.R): R script analyzing body size data. 
3. [3_SelectionAnalysis_ForSubmission.R](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/3_SelectionAnalysis_ForSubmission.R): Script that estimates of selection differentials acting on male and female body size/condition and plots the results.
4. [4_BootstrappingMeanAndVarianceInFitness_ForSubmission.R](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/4_BootstrappingMeanAndVarianceInFitness_ForSubmission.R): R script that bootstraps estimates of variance and mean fitness from the experiment. This script was used to generate the [BootstrappedFitnessDataFinalSept2020.csv](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/BootstrappedFitnessDataFinalSept2020.csv) file.
5. [5_AnalysisOfFitnessVarianceAndNeAndNe\*sRatios_ForSubmission.R](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/5_AnalysisOfFitnessVarianceAndNeAndNe*sRatios_ForSubmission.R): Script that analyzes sex-specific mean and variance in fitness, plots the data and conducts the statistical analyses presented in the paper
6. [6_LikelihoodAnalysis_ForSubmission.R](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/6_LikelihoodAnalysis_ForSubmission.R): Script that calculated likelihood functions for estimates of female variance in fitness and *phi* which is the factor by which variance in male fitness is inflated relative to female fitness. This script calculates likelihood functions and plots likelihood profiles for the two parameters of interest (variance in female fitness and *phi*.  

## Data
Raw data is included in csv files. Metadata describing what each column in each data file mean is described in ["MetaDataForAllDataFiles.txt"](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/MetaDataForAllDataFiles.txt).There are three files in total:
  1. [FitnessDataNeprojectFinalSept2020.csv](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/FitnessDataNeprojectFinalSept2020.csv): Raw fitness data from the experiment 
  2. [BodyMassDataNeprojectFinalSept2020.csv](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/BodyMassDataNeprojectFinalSept2020.csv): Body mass data from males and females
  3. [BootstrappedFitnessDataFinalSept2020.csv](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/BootstrappedFitnessDataFinalSept2020.csv): Results from a bootstrapping procedure wherein estimates of the mean and variance in sex-specific fitness was bootstrapped 10,000 times.
  

