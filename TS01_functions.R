RandomPartition <- function(t, p = 0.5) {
  # Generate a randm partition of t realisations.
  #  Args:
  #     t: Number of realisations.
  #     p: Parameter of binomial distribution (defaul p = 0.5).
  ngroups <- 1 + rbinom(n = 1, size = t - 1, prob = p)
  par <- as.numeric(rmultinom(1, t, rep(1 / ngroups, ngroups)))
  if (length(which(par == 0) != 0)) {
    par <- par[-which(par == 0)]
  }
  return(par)
}

SplitDataPartition <- function(x, partition) {
  # Split multivariate data in blocks specified by a partition.
  #
  # Args:
  #   x: A matrix where each column is an observation (time series). Rows are the time realisations.
  #
  #   partition: The partition of the data we are considering.
  #
  # Returns:
  #   A list where each element is a block of data.
  blocks <- list()
  x <- as.matrix(x)
  for (i in 1:length(partition)) {
    ifelse(i != 1,
      {
        blocks[[i]] <- x[(sum(partition[1:i - 1]) + 1):(partition[i] + sum(partition[1:i - 1])), ]
      },
      {
        blocks[[i]] <- x[1:partition[i], ]
      }
    )
  }
  return(blocks)
}

LogLikelihood <- function(x, partition, parms, coarsening = 1) {
  # Computes the integrated likelihood for a multivariate set of time series that share the same partition.
  #
  # Args:
  #    x: A matrix where each column is a time series.
  #    partition: The partition of the data we consider.
  #    parms: A list with the parameters of the Ornstein-Uhlenbeck process and
  # the Normal-Inverse-Wishart prior.
  #   coarsening: The coarsening correction. Default is 1 (no coarsening correction).
  output_lkl <- as.numeric()
  
  split_data <- SplitDataPartition(x, partition)
  
  for (j in 1:length(split_data)) {
    
    output_lkl[j] <- marginal_likelihood(y = split_data[[j]], 
                                         a = parms$a, 
                                         b = parms$b, 
                                         c = parms$c, 
                                         gamma = parms$gamma)
    
  }
  
  res <- coarsening * sum(output_lkl)
  return(res)
}


SplitMergeAcceleration <- function(y, partition, q, parms, coarsening = 1) {
  # Performs a split and merge on the partition of observations. 
  #   Args: 
  #      y: Matrix of multivariate data, each column an observation.
  #      partition: The partition of the data we want to update.
  #      q: Probability of performing a split.
  #      parms: A list with the parameters of the Ornstein-Uhlenbeck process and
  # the Normal-Inverse-Wishart prior.
  #      coarsening: The coarsening correction. Default is 1 (no coarsening correction).
  k_r <- length(partition)
  
  t <- sum(partition)
  
  prob_split <- q * Indicator1(k_r, t) + Indicator3(k_r)
  prob_merge <- (1 - q) * Indicator1(k_r, t) + Indicator2(k_r, t)
  
  probs <- c(prob_split, prob_merge)
  probs <- probs / sum(probs)
  
  if (runif(1) <= probs[1]) { # SPLIT
    
    output_list <- Split(partition)
    
    split_partition <- output_list[[1]]
    j_proposal <- output_list[[2]]
    
    u <- log(runif(1))
    
    alpha <- AlphaSplitAcceleration(
      y = as.matrix(y),
      j = j_proposal,
      new_partition = split_partition,
      old_partition = partition,
      parms = parms,
      q = q,
      coarsening = coarsening
    )
    
    if (is.na(alpha) == TRUE) {
      alpha <- u - 1
    }
    
    if (u <= alpha) {
      partition_temp <- split_partition
    } else {
      partition_temp <- partition
    }
    
    k_r <- length(partition_temp)
  } else { # MERGE
    
    output_list <- Merge(partition)
    
    merge_partition <- output_list[[1]]
    j_proposal <- output_list[[2]]
    
    u <- log(runif(1))
    
    alpha <- AlphaMergeAcceleration(
      y = as.matrix(y),
      j = j_proposal,
      new_partition = merge_partition,
      old_partition = partition,
      parms = parms,
      q = q,
      coarsening = coarsening
    )
    
    if (is.na(alpha) == TRUE) {
      alpha <- u - 1
    }
    
    if (u <= alpha) {
      partition_temp <- merge_partition
    } else {
      partition_temp <- partition
    }
    
    k_r <- length(partition_temp)
  }
  
  if (k_r > 1) { # SHUFFLE
    
    shuffle_partition <- Shuffle(partition_temp)
    
    u <- log(runif(1))
    
    alpha <- AlphaShuffleAcceleration(
      y = as.matrix(y),
      new_partition = shuffle_partition,
      old_partition = partition,
      parms = parms,
      coarsening = coarsening
    )
    
    if (is.na(alpha) == TRUE) {
      alpha <- u - 1
    }
    
    if (u <= alpha) {
      partition_temp <- shuffle_partition
    } else {
      partition_temp <- partition
    }
  }
  
  return(partition_temp)
}

AlphaSplitAcceleration <- function(y, j, new_partition, old_partition, parms, q, coarsening) {
  a1 <- log((1 - q) / q) +
    (sum(apply(matrix(y, ncol = sum(new_partition)), 1, function(obs) {
      LogLikelihood(x = obs, partition = new_partition, parms = parms, coarsening = coarsening) -
        LogLikelihood(x = obs, partition = old_partition, parms = parms, coarsening = coarsening)
    }))) +
    log(((length(which(old_partition > 1)) * (old_partition[j] - 1)) / (length(old_partition))))
  
  a2 <- log(1)
  
  return(min(a1, a2))
}


AlphaMergeAcceleration <- function(y, j, new_partition, old_partition, parms, q, coarsening) {
  a1 <- log((q / (1 - q))) + (sum(apply(matrix(y, ncol = sum(new_partition)), 1, function(obs) {
    LogLikelihood(x = obs, partition = new_partition, parms = parms, coarsening = coarsening) -
      LogLikelihood(x = obs, partition = old_partition, parms = parms, coarsening = coarsening)
  }))) +
    log(length(old_partition) - 1) - log((length(which(old_partition > 1)) * (old_partition[j] + old_partition[j+1] - 1)))
  
  a2 <- log(1)
  
  return(min(a1, a2))
}


AlphaShuffleAcceleration <- function(y, new_partition, old_partition, parms, coarsening) {
  a1 <- (sum(apply(matrix(y, ncol = sum(new_partition)), 1, function(obs) {
    LogLikelihood(x = obs, partition = new_partition, parms = parms, coarsening = coarsening) -
      LogLikelihood(x = obs, partition = old_partition, parms = parms, coarsening = coarsening)
  })))
  
  a2 <- log(1)
  
  return(min(a1, a2))
}

AlphaSplitPartition <- function(lkl_proposal_i, lkl_proposal_j, lkl_old, k, t, n, 
								 split_i, split_j, norm_costs_row, alpha) {
  
  n_i_split <- length(split_i)
  n_j_split <- length(split_j)  
  q_split <- -(n_i_split + n_j_split - 2) * log(0.5)
  n_ij_split <- n_i_split + n_j_split

  p_split_1 <- -lgamma(alpha) + lgamma(alpha + n_i_split) + lgamma(alpha + n_j_split) - lgamma(alpha + n_ij_split)


  p_split_2 <- logSumExp(c(log(1), -sum(log(2^
    {
      t - 1
    } - 1:(k + 1)) - log(2^
    {
      t - 1
    } - 0:(k))))) +
    logSumExp(c(log(1), -sum(log(2^
      {
        t - 1
      } - 1:(k + 1)) - log(2^
      {
        t - 1
      } - 0:(k))))) -
    logSumExp(c(log(1), -sum(log(2^
      {
        t - 1
      } - 1:(k)) - log(2^
      {
        t - 1
      } - 0:(k - 1)))))


  p_split <- p_split_1 + p_split_2
  l_split <- sum(lkl_proposal_i[split_i]) + sum(lkl_proposal_j[split_j]) - sum(lkl_old[c(split_i, split_j)])

  f_split_n_1 <- logSumExp(lkl_old - norm_costs_row - log(n))

  f_split_d_1 <- logSumExp(lkl_proposal_i - norm_costs_row - log(n))

  f_split_d_2 <- logSumExp(lkl_proposal_j - norm_costs_row - log(n))

  f_split <- f_split_n_1 - f_split_d_1 - f_split_d_2

  a1 <- q_split + p_split + l_split + f_split

  return(min(a1, log(1)))
}


AlphaMergePartition <- function(lkl_old_i, lkl_old_j, lkl_proposal, k, t, n ,
                                alpha, merge_i, merge_j, norm_costs_row) {
  
  
  n_i_merge <- length(merge_i)
  n_j_merge <- length(merge_j)
  
  q_merge <- (n_i_merge + n_j_merge - 2) * log(0.5)
  n_ij_merge <- n_i_merge + n_j_merge

  p_merge_1 <- lgamma(alpha) + lgamma(alpha + n_ij_merge) - lgamma(alpha + n_i_merge) - lgamma(alpha + n_j_merge)

  p_merge_2 <- logSumExp(c(log(1), -sum(log(2^
    {
      t - 1
    } - 1:(k)) - log(2^
    {
      t - 1
    } - 0:(k - 1))))) -
    logSumExp(c(log(1), -sum(log(2^
      {
        t - 1
      } - 1:(k + 1)) - log(2^
      {
        t - 1
      } - 0:(k))))) +
    logSumExp(c(log(1), -sum(log(2^
      {
        t - 1
      } - 1:(k + 1)) - log(2^
      {
        t - 1
      } - 0:(k)))))

  p_merge <- p_merge_1 + p_merge_2

  l_merge <- sum(lkl_proposal[c(merge_i, merge_j)]) - sum(lkl_old_i[merge_i]) - sum(lkl_old_j[merge_j])

  f_merge_n <- logSumExp(lkl_old_i - norm_costs_row - log(n)) + logSumExp(lkl_old_j - norm_costs_row - log(n))

  f_merge_d <- logSumExp(lkl_proposal - norm_costs_row - log(n))

  f_merge <- f_merge_n - f_merge_d

  a1 <- q_merge + p_merge + l_merge + f_merge

  return(min(a1, log(1)))
  
}

EntropyPartition <- function(partition){
  # Compute the entropy of a partition 
  #   Args: 
  #      partition: A vector with the allocation of each observation 
  vec <- as.numeric()
  
  for(j in as.numeric(names(table(partition)))){

    vec[j] <- (length(which(partition == j)) / length(partition)) * log(length(which(partition == j)) / length(partition))
    
  }
  
  return(- sum(vec))

}

Split <- function(partition) {
  # Perform a random split of the partition 
  #   Args: 
  #      partition: A vector with the allocation of each observation
  #   Output:
  #		 A list of two elements: the split partition and the sampled position where the split is performed
  ifelse(length(which(partition > 1)) != 1,
	{pos <- sample(which(partition > 1), size = 1)},{pos <- which(partition > 1)}
   ) 
   
   l <- sample(1:(partition[pos] - 1), 1)
    
  if (pos != 1) {
    ifelse(pos != length(partition),
      {partition <- c(partition[1:(pos - 1)], l, partition[pos] - l, partition[(pos + 1):length(partition)])},{partition <- c(partition[1:(pos - 1)], l, partition[pos] - l)}
     )  
    }
    
  if (pos == 1) {
    ifelse(pos != length(partition),
             {
               partition <- c(l, partition[pos] - l, partition[(pos + 1):length(partition)])
             },
             {
               partition <- c(l, partition[pos] - l)
             }
      )
    }
  
  return(list(partition, pos))
}

Merge <- function(partition) {
  # Perform a random merge of the partition 
  #   Args: 
  #      partition: A vector with the allocation of each observation
  #   Output:
  #		 A list of two elements: the split partition and the sampled position where the split is performed  
  if (length(partition) != 2) {
    pos <- sample(1:(length(partition) - 1), 1)
    
    u <- runif(1)
    
    ifelse(pos != 1,
           {
             ifelse(pos == (length(partition) - 1),
                    {
                      partition <- c(partition <- c(partition[1:(pos - 1)], partition[pos] + partition[(pos + 1)]))
                    },
                    {
                      partition <- c(partition[1:(pos - 1)], partition[pos] + partition[(pos + 1)], partition[(pos + 2):length(partition)])
                    }
             )
           },
           {
             partition <- c(partition[pos] + partition[(pos + 1)], partition[(pos + 2):length(partition)])
           }
    )
    
  }
  
  if (length(partition) == 2) {
    pos <- 1
    
    partition <- c(partition[pos] + partition[(pos + 1)])
    
    k <- length(partition)
  }

  return(list(partition, pos))
}

Shuffle <- function(partition) {
  # Perform a random shuffle of the observations in two clusters of the partition 
  #   Args: 
  #      partition: A vector with the allocation of each observation
  #   Output:
  #		 The partition after the shuffle  
  pos1 <- sample(1:(length(partition) - 1), 1)
  
  pos2 <- sample(1:(partition[pos1] + partition[pos1 + 1] - 1), 1)
  
  partition[pos1 + 1] <- (partition[pos1] + partition[pos1 + 1] - pos2)
  
  partition[pos1] <- pos2
  
  return(partition)
}

SortPartition <- function(partition){
  # Given a partition (a vector with allocation indexes) this function returns the sorted partition. 
  #
  # Args:
  #   partition: A vector with the allocation indexes of the observations. 
  # Output: 
  #   A vector with the sorted partition.
  partition_fix <- rep(0, length(partition))
  for(i in 1:length(table(partition))){
    partition_fix[which(partition == partition[!duplicated(partition)][i])] <- rep(i, length(which(partition == partition[!duplicated(partition)][i])))
  }
  return(partition_fix)
}


Indicator1 <- function(k, n) ifelse(1 < k & k < n, 1, 0)

Indicator2 <- function(k, n) ifelse(k == n, 1, 0)

Indicator3 <- function(k) ifelse(k == 1, 1, 0)

marginal_likelihood <- function(y, a, b, c, gamma){
  
  if(length(y) != 1){
   Sj <- matrix(0, nrow = length(y), ncol = length(y))
   diag(Sj) <- 1 + gamma^2
   for (i in 2:length(y)) {
     Sj[i, i-1] <- -gamma
   }
   for (i in 1:(length(y)-1)) {
     Sj[i, i+1] <- -gamma
   }
   Sj[1,1] = 1
   Sj[length(y),length(y)] = 1
  } else {
   Sj <- matrix(1,1,1)
  }
  
  r1 <- (a * log(2*b * (1-gamma^2)) + lgamma(length(y)/2 + a)) - ((length(y)/2)*log(pi) + lgamma(a))
  
  #r2 <- 0.5 * (log(c) + log(1+gamma) + log(1-gamma^2) - log(c) - log(length(y)) - log(gamma*(length(y)-c-2)))
  
  r2 <- 0.5 * (log(c * (1-gamma) * (1-gamma^2)) - log(c+length(y)-gamma*(length(y)-c-2)))
  
  if(length(y) != 1){
    r3 <- -(length(y)/2 + a) * log((t(y) %*% Sj %*% y) - (((1-gamma)*(sum(y) - gamma*sum(y[2:(length(y)-1)]))^2) / (c+length(y)-gamma*(length(y)-c-2))) + 2*b*(1-gamma^2))
  } else {
    r3 <- -(length(y)/2 + a) * log((t(y) %*% Sj %*% y) - (((1-0)*(sum(y))^2) / (c+length(y)-0*(length(y)-c-2))) + 2*b*(1-0^2))
  }
  
  return(sum(r1,r2,r3))
  
}
