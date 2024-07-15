parms_INIT <- list()
parms_INIT$beta <- list()

parms_INIT$n <- 10 # number of observations 
parms_INIT$partition <- c(rep(1,4),rep(2,3),rep(3,3)) # partition of y_1,....,y_n

parms_INIT$maxT <- 200   # maximum time 
parms_INIT$S0 <- 100000  # number of individuals

# specify time of analysis
parms_INIT$LowerBound <- 10    
parms_INIT$UpperBound <- 150
#


#set parameters

##beta
# k = 1 
parms_INIT$beta[[1]] <- c(0.211, 0.55)
parms_INIT$beta[[2]] <- c(0.221, 0.50)  
parms_INIT$beta[[3]] <- c(0.218, 0.54)   
parms_INIT$beta[[4]] <- c(0.225, 0.51) 

parms_INIT$rho[[1]] <- 115/parms_INIT$S0
parms_INIT$rho[[2]] <- 115/parms_INIT$S0
parms_INIT$rho[[3]] <- 105/parms_INIT$S0
parms_INIT$rho[[4]] <- 100/parms_INIT$S0

#k = 2 
parms_INIT$beta[[5]] <- c(0.213, 0.52) 
parms_INIT$beta[[6]] <- c(0.192, 0.51) 
parms_INIT$beta[[7]] <- c(0.193, 0.57) 

parms_INIT$rho[[5]] <- 120/parms_INIT$S0
parms_INIT$rho[[6]] <- 115/parms_INIT$S0
parms_INIT$rho[[7]] <- 110/parms_INIT$S0

#k = 3
parms_INIT$beta[[8]] <- c(0.195, 0.54) 
parms_INIT$beta[[9]] <- c(0.191, 0.53) 
parms_INIT$beta[[10]] <- c(0.189, 0.51) 

parms_INIT$rho[[8]] <- 105/parms_INIT$S0
parms_INIT$rho[[9]] <- 100/parms_INIT$S0
parms_INIT$rho[[10]] <- 120/parms_INIT$S0
##

##gamma
parms_INIT$gamma <- 1/8 # recovery rate 
##

# position of change points 
cp_INIT <- list()
cp_INIT[[1]] <- c(120)
cp_INIT[[2]] <- c(70)
cp_INIT[[3]] <- c(30)
#

# simulate with Sellke construction 
data_INIT <- sim_epidem_data(parms_INIT, cp_INIT, parallel = TRUE)$infection_times
#

# shrink values bigger than maxT
for(i in 1:length(data_INIT)){
  if(length(which(data_INIT[[i]] > parms_INIT$maxT)) != 0){
    data_INIT[[i]][which(data_INIT[[i]] > parms_INIT$maxT)] =  parms_INIT$maxT 
  }
}
#

# save data
name <- paste0("simEPIdata",today())
path <- ""
save(data_INIT, file = paste0(path,name))
#
