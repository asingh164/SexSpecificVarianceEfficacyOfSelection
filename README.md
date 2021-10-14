# Sex Specific Variance and the Efficacy of Selection

This repository contains code and data underlying the figures and analyses associated with Singh, A. and Agrawal, A.F., 2021. Sex Specific Variance and the Efficacy of Selection

##### Citation: Singh A., and A.F., Agrawal. 2021. Sex Specific Variance and the Efficacy of Selection. AmNat

## Abstract: 
The variance in fitness is thought to be greater in males than in females in many species. If this is so, there are two potentially contradictory consequences on the efficacy of selection (Nes): greater variance in fitness may allow stronger selection (i.e., increased s) but it will also cause stronger genetic drift (i.e., reduced Ne). We develop a simple model to ask how the stronger condition-dependency of fitness in males than females affects selection and fitness variance in each sex to examine the net effect on the efficacy of selection. We measured the phenotypic variance in fitness for each sex in Drosophila melanogaster in different environmental and mating contexts. The variance in fitness was only ~1.5-2 times higher in males than females; juvenile mortality likely dampens the difference in variation between the sexes. Combining these results with previous studies of sex-specific selection on mutations, we infer that the increased drift due to males counterbalances the stronger selection on males in this species, leaving Nes similar to what would be expected if both sexes where ‘female-like’ with respect to selection and variance in fitness. Reasons why this could differ in other species are discussed![image](https://user-images.githubusercontent.com/43476172/137389618-22e88462-2167-4ec6-8d75-ffb5d966e6e7.png)



## Code 
All code is written in R. Coding scripts are ordered roughly in the order that the data are presented in the paper. All figures and tables in the paper and supplement are replicable (with the code provided). <i>NOTE: Although the data in the figures are identical to those that appear in the paper + supplement, the figures will not be exactly the same as I processed figures after generating the basic plot in R (Fonts, labels, panels etc.).</i>

<i><b>Be sure to read the info a the top of each script. There are several places where you will need to specify a file path either to an input data file or for saving outputs of data or figures</i></b>

<i>NOTE: These scripts are some of the first that I had written in R and as such they may not be as pretty, efficient or succinct as they could be.</i>  ¯\\_(ツ)_/¯

There are 6 R scripts in total: 
1. [1_TheoreticalModelAnalysis_ForGithub.R](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/1_TheoreticalModelAnalysis_ForGithub.R): Script that produces figures to plot results from a simple theoretical model presented in the paper.
2. [2_BodySizeAnalysis_ForGithub.R](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/2_BodySizeAnalysis_ForGithub.R): R script analyzing body size data. 
3. [3_SelectionAnalysis_ForGithub.R](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/3_SelectionAnalysis_ForGithub.R): Script that estimates of selection differentials acting on male and female body size/condition and plots the results.
4. [4_BootstrappingMeanAndVarianceInFitness_ForGithub.R](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/4_BootstrappingVariances_ForGithub.R): R script that bootstraps estimates of variance and mean fitness from the experiment. This script was used to generate the [BootstrappedFitnessDataFinalSept2020.csv](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/BootstrappedFitnessDataFinalJuly2021.csv) file.
5. [5_AnalysisOfFitnessVarianceAndNeAndNe\*sRatios_ForGithub.R](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/5_AnalysisOfFitnessVarianceAndNeAndNe*sRatios_ForGithub.R): Script that analyzes sex-specific mean and variance in fitness, plots the data and conducts the statistical analyses presented in the paper. For this script to work, you need to have either run the "<i>4_BootstrappingMeanAndVarianceInFitness_ForGithub.R</i>" script or use the data in the csv file provided that contains an output of the bootstrapping script.
6. [6_LikelihoodAnalysis_ForGithub.R](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/6_LikelihoodAnalysis_ForGithub.R): Script that calculated likelihood functions for estimates of female variance in fitness and <i>phi</i> which is the factor by which variance in male fitness is inflated relative to female fitness. This script calculates likelihood functions and plots likelihood profiles for the two parameters of interest (variance in female fitness and <i>phi</i>.)

## Data
Raw data is included in csv files. Metadata describing what each column in each data file mean is described in ["MetaDataForAllDataFiles.txt"](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/MetadataForAllDataFiles.txt).There are three files in total:
  1. [FitnessDataNeprojectFinalSept2020.csv](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/FitnessDataNeprojectFinalSept2020.csv): Raw fitness data from the experiment 
  2. [BodyMassDataNeprojectFinalSept2020.csv](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/BodyMassDataNeprojectFinalSept2020.csv): Body mass data from males and females
  3. [BootstrappedFitnessDataFinalJuly2021.csv](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/BootstrappedFitnessDataFinalJuly2021.csv): Results from a bootstrapping procedure wherein estimates of the mean and variance in sex-specific fitness was bootstrapped 10,000 times.
  

