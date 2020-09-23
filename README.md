# Sex Specific Variance and the Efficacy of Selection

This repository contains code and data underlying the figures and analyses associated with Singh, A. and Agrawal, A.F., 2020. Sex Specific Variance and the Efficacy of Selection

##### Citation: Singh A., A.F., Agrawal. 2020. Sex Specific Variance and the Efficacy of Selection.

## Abstract: 
The variance in fitness is thought to be greater in males than in females in many species. If this is so, there are two potentially contradictory consequences on the efficacy of selection (<i>Nes</i>): greater variance in fitness may allow stronger selection (i.e., increased <i>s</i>) but it will also cause stronger genetic drift (i.e., reduced <i>Ne</i>). We develop a simple model to ask how the stronger condition-dependency of fitness in males than females affects selection and fitness variance in each sex in order to examine the net effect on the efficacy of selection. We then present estimates of the phenotypic variance in fitness for each sex in Drosophila melanogaster in different contexts expected to alter sex-specific variances. The variance in fitness was not altered substantially by mating regime or degree of heterogeneity in the rearing environment. Sex differences in the variance in fitness were, at most, modestly higher in males than in females. An explanation for our data is that sexually-independent stochastic processes such as juvenile mortality are a major source of variation that dampen sex differences in fitness variation. In general, strong selection through males can occur without a large increase in the strength of drift in systems where ‘random’ mortality is high.


## Code 
Code is all written in R. Coding scripts are ordered roughly in the order that the data are presented in the paper. All figures and tables in the paper and supplement are replicable (with the code provided). NOTE: Although the data in the figures is identical to that in the paper + supplement, the figures will not appear exactly the same as I processed figures after generating the basic plot in R (Fonts, labels, panels etc.).

<i><b>Be sure to read the info a the top of each script. There are several places where you will need to specify a file path either to an input data file or for saving outputs of data or figures</i></b>

<i>NOTE: These scripts are some of the first that I had written in R and as such they may not be as pretty, efficient or succinct as they could be.</i>  ¯\\_(ツ)_/¯

There are 6 R scripts in total: 
1. [1_TheoreticalModelAnalysis_ForSubmission.R](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/1_TheoreticalModelAnalysis_ForSubmission.R): Script that produces figures to plot results from a simple theoretical model presented in the paper.
2. [2_BodySizeAnalysis_ForSubmission.R](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/2_BodySizeAnalysis_ForSubmission.R): R script analyzing body size data. 
3. [3_SelectionAnalysis_ForSubmission.R](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/3_SelectionAnalysis_ForSubmission.R): Script that estimates of selection differentials acting on male and female body size/condition and plots the results.
4. [4_BootstrappingMeanAndVarianceInFitness_ForSubmission.R](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/4_BootstrappingMeanAndVarianceInFitness_ForSubmission.R): R script that bootstraps estimates of variance and mean fitness from the experiment. This script was used to generate the [BootstrappedFitnessDataFinalSept2020.csv](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/BootstrappedFitnessDataFinalSept2020.csv) file.
5. [5_AnalysisOfFitnessVarianceAndNeAndNe\*sRatios_ForSubmission.R](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/5_AnalysisOfFitnessVarianceAndNeAndNe*sRatios_ForSubmission.R): Script that analyzes sex-specific mean and variance in fitness, plots the data and conducts the statistical analyses presented in the paper
6. [6_LikelihoodAnalysis_ForSubmission.R](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/6_LikelihoodAnalysis_ForSubmission.R): Script that calculated likelihood functions for estimates of female variance in fitness and <i>phi</i> which is the factor by which variance in male fitness is inflated relative to female fitness. This script calculates likelihood functions and plots likelihood profiles for the two parameters of interest (variance in female fitness and <i>phi</i>.)

## Data
Raw data is included in csv files. Metadata describing what each column in each data file mean is described in ["MetaDataForAllDataFiles.txt"](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/MetaDataForAllDataFiles.txt).There are three files in total:
  1. [FitnessDataNeprojectFinalSept2020.csv](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/FitnessDataNeprojectFinalSept2020.csv): Raw fitness data from the experiment 
  2. [BodyMassDataNeprojectFinalSept2020.csv](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/BodyMassDataNeprojectFinalSept2020.csv): Body mass data from males and females
  3. [BootstrappedFitnessDataFinalSept2020.csv](https://github.com/asingh164/SexSpecificVarianceEfficacyOfSelection/blob/master/BootstrappedFitnessDataFinalSept2020.csv): Results from a bootstrapping procedure wherein estimates of the mean and variance in sex-specific fitness was bootstrapped 10,000 times.
  

