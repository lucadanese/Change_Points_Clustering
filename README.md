# Change_Points_Clustering
This repository contains the code for replicating the code in the paper *Model-based clustering of time-dependent observations with common structural changes* (2024) by Corradin R., Danese L., KhudaBukhsh W.R. and Ongaro A. . 

The scripts are divided in two groups, those starting with *TS* are for time series data, while those starting with *EPI* are for epidemiological data. 

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
