## SIMULATIIVE DATA ##

# CASE 1: Different mean, same variance, time series have similar shapes and 
#         change points

data_sim <- as.data.frame(matrix(NA, nrow = 300, ncol = 10))

gamma <- 0.1

## GROUP 1 ----

j = 1

mu <- 0.5 # 1

data_sim[1,j] <- mu + rnorm(1,0,0.1)

for (i in 2:50){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.1)
}

mu <- 0.85 # 2

for (i in 51:150){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.12)
}

mu <- 0.5 # 3

for (i in 151:195){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.14)
}

mu <- 0.75 # 4

for (i in 196:250){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.13)
}

mu <- 1 # 5

for (i in 251:300){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.15)
}

##

##

j = 2

mu <- 0.15 # 1

data_sim[1,j] <- mu + rnorm(1,0,0.12)

for (i in 2:50){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.12)
}

mu <- 0.75 # 2

for (i in 51:150){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.15)
}

mu <- 0.25 # 3

for (i in 151:195){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.12)
}

mu <- 0 # 4

for (i in 196:250){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.14)
}

mu <- 0.25 # 5

for (i in 251:300){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.13)
}

##

##

j = 3

mu <- 0.25 # 1

data_sim[1,j] <- mu + rnorm(1,0,0.1)

for (i in 2:50){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.1)
}

mu <- 0 # 2

for (i in 51:150){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.12)
}

mu <- 0.15 # 3

for (i in 151:195){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.2)
}

mu <- 0.15 # 4

for (i in 196:250){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.12)
}

mu <- 0.3 # 5

for (i in 251:300){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.14)
}

##

##

j = 4

mu <- 0.75 # 1

data_sim[1,j] <- mu + rnorm(1,0,0.1)

for (i in 2:50){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.1)
}

mu <- 0.4 # 2

for (i in 51:150){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.12)
}

mu <- 0.8 # 3

for (i in 151:195){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.09)
}

mu <- 0.8 # 4

for (i in 196:250){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.24)
}

mu <- 0.4 # 5

for (i in 251:300){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.15)
}


## GROUP 2 ----


j = 5

mu <- 0 # 1

data_sim[1,j] <- mu + rnorm(1,0,0.12)

for (i in 2:40){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.12)
}

mu <- -0.15 # 2

for (i in 41:90){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.13)
}

mu <- 0.15 # 3

for (i in 91:135){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.1)
}

mu <- 0.3  # 4

for (i in 136:180){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.13)
}

mu <- 0.1  # 5

for (i in 181:210){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.14)
}

mu <- 0.3  # 6

for (i in 211:300){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.12)
}

##

j = 6

mu <- 0.5 # 1

data_sim[1,j] <- mu + rnorm(1,0,0.1)

for (i in 2:40){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.1)
}

mu <- 0 # 2

for (i in 41:90){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.24)
}

mu <- -0.5 # 3

for (i in 91:135){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.14)
}

mu <- 0 # 4

for (i in 136:180){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.15)
}

mu <- 0.2  # 5

for (i in 181:210){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.12)
}

mu <- 0  # 6

for (i in 211:300){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.13)
}

##

j = 7

mu <- 0 # 1

data_sim[1,j] <- mu + rnorm(1,0,0.16)

for (i in 2:40){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.16)
}

mu <- 0.2 # 2

for (i in 41:90){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.15)
}

mu <- 0.4 # 3

for (i in 91:135){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.1)
}

mu <- 0.25  # 4

for (i in 136:180){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.13)
}

mu <- -0.1  # 5

for (i in 181:210){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.14)
}

mu <- 0.15  # 5

for (i in 211:300){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.12)
}

## GROUP 3 ----

j = 8

mu <- 0  # 1

data_sim[1,j] <- mu + rnorm(1,0,0.14)

for (i in 2:75){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.14)
}

mu <- -0.25 # 2

for (i in 76:150){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.13)
}

mu <- 0 # 3

for (i in 151:180){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.17)
}

mu <- 0.25 # 4

for (i in 181:200){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.12)
}

mu <- -0.25 # 5

for (i in 201:275){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.14)
}

mu <- 0.1 # 6

for (i in 276:300){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.12)
}

##

j = 9

mu <- 0.25  # 1

data_sim[1,j] <- mu + rnorm(1,0,0.12)

for (i in 2:75){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.12)
}

mu <- 0.25  # 2

for (i in 76:150){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.22)
}

mu <- -0.2 # 3

for (i in 151:180){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.15)
}

mu <- 0.1 # 4

for (i in 181:200){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.14)
}

mu <- 0.3 # 5

for (i in 201:275){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.17)
}

mu <- 0 # 6

for (i in 276:300){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.19)
}

##

j = 10

mu <- 0  # 1

data_sim[1,j] <- mu + rnorm(1,0,0.12)

for (i in 2:75){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.12)
}

mu <- -0.25  # 2

for (i in 76:150){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.13)
}

mu <- 0 # 3

for (i in 151:180){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.15)
}

mu <- 0.25 # 4

for (i in 181:200){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.12)
}

mu <- 0 # 5

for (i in 201:275){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.15)
}

mu <- -0.25 # 6

for (i in 276:300){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.18)
}



