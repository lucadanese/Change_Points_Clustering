data_sim <- as.data.frame(matrix(NA, nrow = 500, ncol = 10))

gamma <- 0.1

## GROUP 1 ----

j = 1

mu <- 1 # 1

data_sim[1,j] <- mu + rnorm(1,0,0.23)

for (i in 2:150){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.23)
}

mu <- 1.25 # 2

for (i in 151:300){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.25)
}

mu <- 0.5 # 3

for (i in 301:500){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.35)
}

##

j = 2

mu <- 1.25 # 1

data_sim[1,j] <- mu + rnorm(1,0,0.21)

for (i in 2:150){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.21)
}

mu <- 1.5 # 2

for (i in 151:300){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.16)
}

mu <- 1 # 3

for (i in 301:500){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.19)
}

##

j = 3

mu <- 0.75 # 1

data_sim[1,j] <- mu + rnorm(1,0,0.25)

for (i in  2:150){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.25)
}

mu <- 0.5 # 2

for (i in 151:300){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.17)
}

mu <- 1 # 3

for (i in 301:500){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.21)
}

##

j = 4

mu <- 0 # 1

data_sim[1,j] <- mu + rnorm(1,0,0.22)

for (i in 2:150){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.22)
}

mu <- 0.5 # 2

for (i in 151:300){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.24)
}

mu <- 0 # 3

for (i in 301:500){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.21)
}

## GROUP 2 ----


j = 5

mu <- 0 # 1

data_sim[1,j] <- mu + rnorm(1,0,0.23)

for (i in 2:100){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.23)
}

mu <- -0.25 # 2

for (i in 101:170){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.19)
}

mu <- 0.25 # 3

for (i in 171:300){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.21)
}

mu <- 0.75  # 4

for (i in 301:400){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.18)
}

mu <- 1.5  # 5

for (i in 401:500){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.17)
}

##

j = 6

mu <- -1 # 1

data_sim[1,j] <- mu + rnorm(1,0,0.24)

for (i in 2:100){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.24)
}

mu <- -0.5 # 2

for (i in 101:170){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.25)
}

mu <- -1 # 3

for (i in 171:300){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.19)
}

mu <- 0  # 4

for (i in 301:400){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.17)
}

mu <- 0.5  # 5

for (i in 401:500){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.14)
}

##

j = 7

mu <- -1 # 1

data_sim[1,j] <- mu + rnorm(1,0,0.20)

for (i in 2:100){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.20)
}

mu <- - 1.25 # 2

for (i in 101:170){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.22)
}

mu <- -1 # 3

for (i in 171:300){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.24)
}

mu <- 0  # 4

for (i in 301:400){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.21)
}

mu <- 1  # 5

for (i in 401:500){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.24)
}



## GROUP 3 ----

j = 8

mu <- 0  # 1

data_sim[1,j] <- mu + rnorm(1,0,0.24)

for (i in 2:180){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.24)
}

mu <- -0.5  # 2

for (i in 181:375){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.20)
}

mu <- 0 # 3

for (i in 376:500){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.19)
}

j = 9

mu <- -1  # 1

data_sim[1,j] <- mu + rnorm(1,0,0.2)

for (i in  2:180){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.2)
}

mu <- -1.5  # 2

for (i in 181:375){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.23)
}

mu <- -0.5 # 3

for (i in 376:500){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.21)
}

j = 10

mu <- 0  # 1

data_sim[1,j] <- mu + rnorm(1,0,0.1)

for (i in 2:180){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.22)
}

mu <- -1  # 2

for (i in 181:375){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.2)
}

mu <- -0.25 # 3

for (i in  376:500){
  data_sim[i,j] <- gamma * data_sim[i-1,j] + (1 - gamma) * mu + rnorm(1,0,(1-gamma)*0.25)
}




