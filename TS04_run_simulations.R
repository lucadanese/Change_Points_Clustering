##
rm(list = ls())

setwd("")

set.seed(123)

#--LIBRARIES--#

library("Rcpp")
library("BNPmix")
library("parallel")
library("dplyr")
library("doSNOW")
library("purrr")
library("tidyr")
library("matrixStats")
library("LaplacesDemon")
library("MASS")
library("data.table")

#--PARAMETERS--#
alpha_INPUT <- 1
n_iterations_INPUT <- 5000

OU_parameters_INPUT <- list(
  
  a = 1,
  
  b = 1,
  
  c = 0.1,
  
  gamma = 0.1
)

#--IMPORT FUNCTIONS AND DATA--##
source("TS01_functions.R")
source("TS02_main_code.R")

parms_B <- c(1000,10000,100000)
parms_L <- c(1,25,100)
 
matrix_parms <- expand.grid(parms_B,parms_L)

for(j_par in 9:9){
  
  for(rep in 1:50){
    
    source("TS03_generate_simdata.R")
    
	  data_sim <- t(scale(data_sim))
   
    res <- MainFunction(data = data_sim,
	          n_iterations = n_iterations_INPUT,
					  B = matrix_parms[j_par,1], 
						L = matrix_parms[j_par,2], 
						q = 0.5, 
						alpha = alpha_INPUT, 
						OU_parameters = OU_parameters_INPUT)
									
	  path <- c("")
    
    save(data, file = paste0(path,"",matrix_parms[j_par,1],"_" , matrix_parms[j_par,2],"_",rep))
    save(res,  file = paste0(path,"",matrix_parms[j_par,1],"_" , matrix_parms[j_par,2],"_",rep))
    
  }
}

