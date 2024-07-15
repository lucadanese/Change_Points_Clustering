library(BNPmix)
library(coda)
library(salso)

devtools::install_github(“sarawade/mcclust.ext”)
library(mcclust.ext)

ESS_entropy <- as.numeric()
BI_hat <- as.numeric()

# true partition
vec_true <- as.numeric()

vec_true[1:4] <- 1
vec_true[5:7] <- 2
vec_true[8:10] <- 3

name <- ""

for(h in 1:50){
  
  load(paste0(name, "_",h,".RData"))
  
  est_part <- mcclust.ext::minbinder.ext(psm = salso::psm(est_model$clust), cls.draw = est_model$clust)$cl
  
  ESS_entropy[h] <- coda::effectiveSize(apply(est_model$clust, 1, entropy_partition))/nrow(est_model$clust)
  BI_hat[h] <- salso::binder(truth = vec_true, estimate = est_part)
}

mean(ESS_entropy); mean(BI_hat)



