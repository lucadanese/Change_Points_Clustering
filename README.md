# Change_Points_Clustering
This repository contains the code for clustering both time series of survival functions with common change points. It is based on the paper ... and allows to replicate the results in it.

The code is divided in two groups, those starting with TS (time series) and those starting with EPI (epidemiological). 

## Time series data 

- *TS01_functions.R* contains all the functions needed
- *TS02_main_code.R* contains the main algorithm
- *TS03_generate_simdata.R* code to generate simulated data
- *TS04_run_simulations.R* run all the simulations
- *TS05_point_estimate.R* contains the code to get the point estimate of the partition 

## Epidemiological data

- *EPI01_main_code.cpp* contains the main algorithm
- *EPI02_generate_simdata.R* code to generate simulated EPI data
- *EPI03_run_simulations.R* run all the simulations
- *EPI04_point_estimate.R* contains the code to get the point estimate of the partition 
