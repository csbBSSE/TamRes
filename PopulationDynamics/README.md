# TamRes Population dynamics codes
## 1. Simulation code
The main code for the simulation of the effect of drug on the dynamics of a population of cells is``` Tamres_Population_Dynamics_Simulation.jl ```.

To run this file use the command 

``` julia Tamres_Population_Dynamics_Simulation.jl P_SR P_RS Heterogeneity InitialPopulation DrugInducedPlasticity METInducedPlasticity InitialFracResistant OutputFileName``` 

where 

  - `P_SR` denotes the transition probability from Susceptible to Resistant population
  - `P_RS` denotes the transition probability from Resistant to Susceptible population
  - `Heterogeneity` denotes the variance of the distribution of the TamRes score
  - `InitialPopulation` denotes whether the simulation starts from initial susceptible or resistant population
  - `DrugInducedPlasticity` denotes the probability of swtiching in the presence of drug
  - `METInducedPlasticity` denotes the probability of switching in the presence of MET
  - `InitialFracResistant` denotes the fraction of resistant cells to begin with
  - `OutputFileName` denotes the name of the output file

This can be coupled with the bash code `Multiple_Runs_bash_code.sh` to run the simulation for different values of transition probabilities `P_SR` and `P_RS` and `Heterogeneity`. This data can then be used to plot the different figures including the matrix plot for variation of population size with time. 

This code was written using Julia programming language version 1.5.0.

Requirements for running the codes are the following packages of Julia programming language
  - Random
  - Distributions
  - Statistics
  - StatsBase

## 2. Figure codes
These contains codes for the different plots in the paper. 

The .jl files are written in Julia programming language version 1.5.0. The .py files are written with python version 3.5.8

Requirements for running the codes are the following packages of Julia programming language
  - Random
  - Distributions
  - Statistics
  - StatsBase
  - Plots
  - StasPlots
  - PyPlot (This requires matplotlib to be installed)

Requirements for python codes
  - Numpy
  - Matplotlib
  - Seaborn

For the python codes enter the path to the output files in `path_to_files` variable and the path where figures would be saved in the variable `path_to_savefig`.

