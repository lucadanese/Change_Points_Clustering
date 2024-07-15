MainFunction <- function(data, 
						  n_iterations, 
						  B, 
						  L, 
						  q, 
						  alpha, 
						  OU_parameters, 
						  real_time_diagnostic = TRUE, 
						  step_diagnostic = 250, 
						  coarsening = 1){
  # Perform clustering of time series with common change points. 
  #
  # Args: 
  #	   data: A matrix where the rows are the observations and the columns the realisations of the time series. 
  #    n_iterations: Number of iterations of the MCMC sampler. 
  #    B: Number of orders to sample for each observations in computing the normalisation constant. 
  #    L: Number of split-and-merge steps to perform when proposing a new orders. 
  #    q: Probability of performing a split when a split-and-merge is performed on orders. 
  #    alpha: Concentration parameters for the Dirichlet distributed weights in the distribution of the random orders.
  #    OU_parameters: Parameters for the Ornstein-Uhlenbeck process. 
  #    real_time_diagnostic: Wheter or not print the processing status of the algorithm.
  #    step_diagnostic: How many steps occurs between one diagnostic print and the other. 
  #    coarsening: Coarsening ratio.   

  ## GENERATE STARTING PARTITION AND STARTING ORDERS ## 
  maxT <- ncol(data)
  n <- nrow(data)
  
  # generate initial partition 
  number_starting_groups <- sample(1:n, 
                                   size = 1, 
                                   replace = FALSE)
  
  partition_temp <- sample(c(1:number_starting_groups, 
                  sample(1:number_starting_groups, 
                         size = n - number_starting_groups, 
                         replace = TRUE)))
    
  partition_temp <- SortPartition(partition_temp)
  orders_temp <- lapply(1:number_starting_groups, function(x) RandomPartition(t = maxT, p = 5/maxT))
  
  ## DEFINE CLUSTER FOR PARALLEL COMPUTATION
  cluster <- parallel::makeCluster(parallel::detectCores()-1)
  
  pkg <- c("deSolve", "dplyr", "tidyr", "matrixStats","parallel","parallelRcppArmadillo", "LaplacesDemon")
  
  parallel::clusterExport(cluster, varlist = ls(), envir = environment())
  parallel::clusterExport(cluster, varlist = c("LogLikelihood", "SplitDataPartition", "RandomPartition", "SplitMergeAcceleration", "Indicator1", "Indicator2", "Indicator3"), envir = environment())
  
  invisible(parallel::clusterEvalQ(cluster, lapply(pkg, require, character.only = TRUE)))
  
  registerDoSNOW(cluster)
  ##

  ## COMPUTE STARTING LIKELIHOOD ##  
  likelihood_temp <- as.numeric()
  for(i in 1:n){
	likelihood_temp[i] <- LogLikelihood(x = data[i,],
										partition = orders_temp[[partition_temp[i]]],
										parms = OU_parameters, 
										coarsening = coarsening)
  }

  likelihood_clust <- sapply(1:max(partition_temp), function(x) sum(likelihood_temp[which(partition_temp == x)]))
  
  ## COMPUTE NORMALISATION CONSTANT ##
  
  partition_star_list_nc <- parLapply(1:B,
    cl = cluster,
    fun = function(x) RandomPartition(t = maxT, p = 3/maxT)
  )
  
  print("Normalisation constant: computing...")

  comp_matrix <- matrix(data = NA, nrow = n, ncol = length(partition_star_list_nc))

  grid <- expand.grid(1:nrow(comp_matrix), 1:ncol(comp_matrix))

  parallel::clusterExport(cluster, varlist = c("partition_star_list_nc", "grid", "comp_matrix", "data", "LogLikelihood", "OU_parameters", "SplitDataPartition"), envir = environment())
  
  comp_matrix_2 <- matrix(unlist(parLapply(1:nrow(grid),
    cl = cluster,
    fun = function(x) {
      LogLikelihood(x = data[grid[x, 1],], partition = partition_star_list_nc[[grid[x, 2]]], parms = OU_parameters, coarsening = coarsening) -
  
        dmultinom(
          x = partition_star_list_nc[[grid[x, 2]]],
          prob = rep(
            1 / length(partition_star_list_nc[[grid[x, 2]]]),
            length(partition_star_list_nc[[grid[x, 2]]])
          ),
          log = TRUE
        ) -
        dbinom(length(partition_star_list_nc[[grid[x, 2]]]), size = maxT - 1, prob = 0.5, log = TRUE)
    }
  )), nrow = n)
  
  norm_costs_row <- apply(comp_matrix_2, 1, function(x) logSumExp(logSumExp(x) + log(ncol(comp_matrix_2)) - (maxT - 1) * log(2)))                                      
  
  parallel::clusterExport(cl = cluster,varlist = c("norm_costs_row"), envir = environment())
  
  print("Normalisation constant: computed!") 

  ## MAIN LOOP ##
  
  partitions_matrix <- matrix(NA, nrow = n_iterations, ncol = n)
  entropy_vector <- as.numeric()
  marg_likelihood <- as.numeric()
  
  start_time <- Sys.time()
  
  print(" ----------- MAIN LOOP ----------- " )
  
  for(iter in 1:n_iterations){
    
    # SPLIT-AND-MERGE ON THE PARTITION #

    ij <- sample(1:n, 2, replace = FALSE)
    
    if(partition_temp[ij[1]] == partition_temp[ij[2]]){
      
      S <- which(partition_temp == partition_temp[[ij[1]]] | partition_temp == partition_temp[[ij[2]]])
      
      #randomly assign obs to split_i or split_j
      allocation_index <- as.numeric(rep(NA,length(S)))
      allocation_index[which(S == ij[1])] <- 1 
      allocation_index[which(S == ij[2])] <- 2 
      
      for(i in which(is.na(allocation_index))){
        allocation_index[i] = sample(c(1,2),1)
      }
      
      order <- orders_temp[[partition_temp[ij[1]]]]
      
      ## sample partition 
      
      split_list <- list(S[which(allocation_index == 1)], S[which(allocation_index == 2)])
      
      proposed_list <- list()
        
      for(i in 1:2){
          
          obs_sampled <- sample(S, 1, replace = TRUE)
          
      		proposal_partition_temp <- unlist(orders_temp[partition_temp[obs_sampled]])
  		
      		for(j in 1:L){
      		  
      		  proposal_partition_temp <- SplitMergeAcceleration(y = matrix(data[obs_sampled,], ncol = maxT), 
      		                                                    partition = proposal_partition_temp, 
      		                                                    q = q, parms = OU_parameters, coarsening = coarsening)
  
      		}
      		
        
      		likelihood_proposed_partition <- apply(data, 1, function(x) LogLikelihood(x = x, 
      		                                               partition = proposal_partition_temp, 
      		                                               parms = OU_parameters, 
      		                                               coarsening = coarsening))

      		proposed_list[[i]] <- list(proposal_partition_temp, as.numeric(likelihood_proposed_partition))
          
        }

      lkl_proposal_i <- proposed_list[[1]][[2]]
      lkl_proposal_j <- proposed_list[[2]][[2]]
      
      order_i_split <- proposed_list[[1]][[1]]
      order_j_split <- proposed_list[[2]][[1]]

      lkl_old <- as.numeric(apply(data, 1, function(x) LogLikelihood(x = x, 
                                                           partition = order, 
                                                           parms = OU_parameters, 
                                                           coarsening = coarsening)))
      
      alpha_split <- AlphaSplitPartition(lkl_proposal_i = lkl_proposal_i, 
										 lkl_proposal_j = lkl_proposal_j, 
										 lkl_old = lkl_old, 
                     k = length(table(partition_temp)), 
										 t = maxT,  
										 n = n, 
                     alpha = alpha, 
										 split_i = S[which(allocation_index == 1)], 
										 split_j = S[which(allocation_index == 2)], 
										 norm_costs_row = norm_costs_row)

      if (log(runif(1)) <= alpha_split) {
        
        partition_temp[S[which(allocation_index == 1)]] <- length(table(partition_temp)) + 1
        partition_temp[S[which(allocation_index == 2)]] <- length(table(partition_temp)) + 1
        
        orders_temp[[length(orders_temp)+1]] <- order_i_split
        orders_temp[[length(orders_temp)+1]] <- order_j_split
        
        likelihood_clust[[length(likelihood_clust)+1]] <- sum(lkl_proposal_i[S[which(allocation_index == 1)]])
        likelihood_clust[[length(likelihood_clust)+1]] <- sum(lkl_proposal_j[S[which(allocation_index == 2)]])
        
        orders_temp <- orders_temp[partition_temp[!duplicated(partition_temp)]]
        likelihood_clust <- likelihood_clust[partition_temp[!duplicated(partition_temp)]]
        
		    partition_temp <- SortPartition(partition_temp)
        
		    likelihood_temp[S[which(allocation_index == 1)]] <- lkl_proposal_i[S[which(allocation_index == 1)]]
        likelihood_temp[S[which(allocation_index == 2)]] <- lkl_proposal_j[S[which(allocation_index == 2)]]

      } 
      
    } else {

	  S <- c(which(partition_temp == partition_temp[[ij[1]]]),which(partition_temp == partition_temp[[ij[2]]]))
  
      order_i <- orders_temp[[partition_temp[[ij[1]]]]]
      order_j <- orders_temp[[partition_temp[[ij[2]]]]]
      
      obs_sampled <- sample(S, 1)
      
      proposed_partition <- orders_temp[[partition_temp[obs_sampled]]]
      
      for(j in 1:L){
        
        proposed_partition <- SplitMergeAcceleration(y = matrix(data[obs_sampled,],ncol = maxT), 
                                                 partition = proposed_partition, 
                                                 q = q, parms = OU_parameters, coarsening = coarsening)
      }
      
      lkl_proposal <- apply(data, 1, function(x) LogLikelihood(x = x, partition = proposed_partition, parms = OU_parameters, coarsening = coarsening))
        
      lkl_old_i <- apply(data, 1, function(x) LogLikelihood(x = x, partition = order_i, parms = OU_parameters, coarsening = coarsening))
        
      lkl_old_j <- apply(data, 1, function(x) LogLikelihood(x = x, partition = order_j, parms = OU_parameters, coarsening = coarsening))
      
      alpha_merge <- AlphaMergePartition(lkl_old_i = lkl_old_i, 
										 lkl_old_j = lkl_old_j, 
										 lkl_proposal = lkl_proposal, 
										 k = length(table(partition_temp)), 
										 t = maxT, 
										 n = n, 
										 alpha = alpha, 
										 merge_i = which(partition_temp == partition_temp[[ij[1]]]), 
										 merge_j = which(partition_temp == partition_temp[[ij[2]]]), 
                     norm_costs_row = norm_costs_row)
      
      if (log(runif(1)) <= alpha_merge) {
        
        partition_temp[S] = partition_temp[ij[1]]
        
        orders_temp[[partition_temp[ij[1]]]] <- partition_temp
        orders_temp <- orders_temp[-partition_temp[ij[2]]]
        
        likelihood_clust[[partition_temp[ij[1]]]] <- sum(lkl_proposal[S])
        likelihood_clust <- likelihood_clust[-partition_temp[ij[2]]]
        
        partition_temp <- SortPartition(partition_temp)
        likelihood_temp[S] <- lkl_proposal[S] 
        
      } 
    
    }
    
    # ACCELERATION STEP #
    orders_temp <- lapply(X = 1:length(orders_temp), 
							function(x) SplitMergeAcceleration(y = matrix(data[which(partition_temp == x),], ncol = maxT), partition = orders_temp[[x]], 
                                                                q = q, parms = OU_parameters, coarsening = coarsening))
    
    # SAVE TEMP RESULTS #
    partitions_matrix[iter,] <- partition_temp
    entropy_vector[iter] <- EntropyPartition(partition_temp)
    marg_likelihood[iter] <- sum(likelihood_temp)
	
	# REAL TIME DIAGNOSTIC #
    if(real_time_diagnostic == TRUE){
		if(iter %in% seq(0,n_iterations,step_diagnostic)[-1]){
			print(paste0("Completed: ", iter,"/",n_iterations, " || ", "Elapsed Time: ", round(difftime(Sys.time(), start_time, units='mins'),2), " mins"))
		}
    }
  }
  
  stopCluster(cluster)   
  
  return(list(partitions = partitions_matrix, 
              entropy = entropy_vector,
              likelihod = marg_likelihood))
}
