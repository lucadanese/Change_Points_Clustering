#--LIBRARIES AND FUNCTIONS--##
library(Rcpp)
sourceCpp(file = "EPI01_main_code.cpp")

##--IMPORT DATA--##

# run the following line if data have not been generated yet
#source("EPI02_generate_simdata.R") 

path <- "G:\\My Drive\\Dottorato\\Progetti\\Model Based Clustering of TS\\Simulazioni_EPI\\data\\simEPIdata_24_05_2024" # path with the simulated EPI data  
load(path)

##--SET PARAMETERS--##
maxT <- 200       # maximum obs time 
LowerBound <- 10  # time range lower bound
UpperBound <- 150 # time range upper bound

par_mcrep <- 250 # number MC replications
L <- 1           # split and merge steps
B <- 1000        # dimension of the normalisation constant
sample_prop <- 0.20 # set proportion of sampled infection times

##--RUN SIMULATIONs--##

set.seed(123)

for(r in 1:1){
  
  for(rep in 1:50){

    y <- data_INIT
    
    y_rounded_pre <- lapply(y, function(x) table(sample(ceiling(x), sample_prop * length(x), replace = FALSE)))
    
    y_rounded <- list()
    
    for(i in 1:length(y_rounded_pre)){
      y_rounded[[i]] <- rep(0, length(seq(1,maxT)))
      names(y_rounded[[i]]) <- seq(1,maxT)
      y_rounded[[i]][as.numeric(names(y_rounded_pre[[i]]))] <- as.numeric(y_rounded_pre[[i]])
    }
    
    data_sim <- matrix(NA, 10, maxT)
    
    for(i in 1:length(y_rounded)){
      data_sim[i,1:200] <- y_rounded[[i]]
    }
    
    data_sim <- as.matrix(data_sim[,LowerBound:UpperBound])
    
    est_model <- main_function(data = data_sim, 
                               niter = 5000, nburn = 2000,
                               alpha = 1, q = 0.1, dt = 0.1, 
                               a0 = 4, b0 = 10, c0 = 1, d0 = 1, 
                               gamma = 1/8, MH_var = 0.01, 
                               M = par_mcrep, 
                               B = B,
                               L = L, 
                               S0 = 1, R0 = 0, 
                               p = 0.005, nupd = 250)
	
	  path_to_save <- "G:\\My Drive\\Dottorato\\Progetti\\Model Based Clustering of TS\\Simulazioni_EPI\\11_07_2024\\"
    
    save(est_model, file = paste0(path_to_save,"est_sim_EPI_", "prova", "_", par_mcrep, "_", par_nsplitmerge, "_", par_R,"_",rep ,".RData"))

  }
  
}



