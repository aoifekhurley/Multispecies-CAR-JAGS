library(tidyverse)
library(R2jags)
library(tidybayes)
library(sp)
library(sf)
library(spdep)
library(MASS)
library(raster)
library(prioritizr)
library(mcmcplots)
library(ggpubr)
library(LaplacesDemon)
#### FUNCTIONS ----
  st_rook = function(a, b = a) st_relate(a, b, pattern = "F***1****")
#### SET UP ----
  set.seed(240323)
  max.T <- 5
  R <- 25
  no.c <- 5
  species <- 3
  ##### GRID -----
    x1.length <- sqrt(R)
    x1 <- seq(1:x1.length)
    x2 <- x1
    
    grid <- expand.grid(x = x1, y = x2)
    
    grid.pts <- SpatialPointsDataFrame(coords= grid, 
                                       data=grid, 
                                       proj4string = CRS("+proj=longlat +datum=WGS84")
    )
    gridded(grid.pts) <- TRUE
    
    grid <- as(grid.pts, "SpatialPolygons")
    plot(grid)
    gridspdf <- SpatialPolygonsDataFrame(grid, data=data.frame(ID=row.names(grid),
                                                               row.names=row.names(grid)
    )
    )
    grid.sf <- st_as_sf(gridspdf)
    ggplot(grid.sf) + 
      geom_sf() + 
      geom_sf_text(aes(label = ID))
    
  ##### ADJACENCY MATRIX -----  
    A.list <- st_rook(grid.sf)  
    
    A <- matrix(0, nrow = R, ncol = R)
    for(i in 1:R){
      A[i, A.list[[i]]] <- 1
    }
    
    D <- diag(rowSums(A))    ## Number of neighbours on diagonal
    
#### DATA GENERATION ----    
  alpha.spatial <- 0.6
  #Sigma <- rinvwishart(nu = species+1, S = diag(species))
  tau.spatial <- c(2, 3, 1) #rho12, rho13, rho23
  rho.spat <- c(0.3, -0.2, 0.6)
  
  const.det <- 1 + ((rho.spat[1]^2)*(rho.spat[2]^2)*(rho.spat[3]^2)*(1/tau.spatial[1])*(1/tau.spatial[2])*(1/tau.spatial[3])) - 
                ((rho.spat[1]^2)*(1/tau.spatial[1])*(1/tau.spatial[2])) - ((rho.spat[2]^2)*(1/tau.spatial[1])*(1/tau.spatial[3])) - 
                ((rho.spat[3]^2)*(1/tau.spatial[2])*(1/tau.spatial[3]))
  
  T11 <- tau.spatial[1]*(1 - rho.spat[3]^2*(1/tau.spatial[2])*(1/tau.spatial[3]))
  T12 <- rho.spat[2]*rho.spat[3]*(1/tau.spatial[3]) - rho.spat[1]
  T13 <- rho.spat[1]*rho.spat[3]*(1/tau.spatial[2]) - rho.spat[2]
  T22 <- tau.spatial[2]*(1 - rho.spat[2]^2*(1/tau.spatial[1])*(1/tau.spatial[3]))
  T23 <- rho.spat[1]*rho.spat[2]*(1/tau.spatial[1]) - rho.spat[3]
  T33 <- tau.spatial[3]*(1 - rho.spat[1]^2*(1/tau.spatial[1])*(1/tau.spatial[2]))
  
  Tau.matrix<- matrix(NA, nrow = species, ncol = species)
  for(i in 1:species){
    for(j in 1:species){
      if(i == j){
        Tau.matrix[i,j] <- ifelse(i==1, T11, ifelse(i==2, T22, T33))
      }else if (i < j){
        Tau.matrix[i,j] <- ifelse(i==1 & j==2, T12, ifelse(i==1 & j ==3, T13, T23))
      }else{
        Tau.matrix[i,j] <- Tau.matrix[j,i] 
      }
    }
  }
  Tau.matrix <- Tau.matrix/const.det
  Q <- kronecker(Tau.matrix, D - alpha.spatial*A)
  
  phi <- mvrnorm(n = 1, mu = rep(0, nrow(Q)), Sigma = solve(Q))
  
  phi.1 <- phi[1:R]
  phi.2 <- phi[(1+R):(2*R)]
  phi.3 <- phi[(1+(2*R)):(3*R)]
  
  lambda1 <- exp(phi.1)
  lambda2 <- exp(phi.2)
  lambda3 <- exp(phi.3)
  
  N1 <- rpois(R, lambda1)  
  N2 <- rpois(R, lambda2)
  N3 <- rpois(R, lambda3)
  
  p1 <- 0.3
  p2 <- 0.6
  p3 <- 0.9
  
  y1 <- matrix(rbinom(n = R*max.T, size = N1, prob = p1), nrow = R, ncol = max.T)
  y2 <- matrix(rbinom(n = R*max.T, size = N2, prob = p2), nrow = R, ncol = max.T)
  y3 <- matrix(rbinom(n = R*max.T, size = N3, prob = p3), nrow = R, ncol = max.T)
  
  ##### CULL DATA GENERATION -----
  Z = matrix(0, nrow = R, ncol = no.c) # Matrix, 1 in column for associated county
  for(i in 1:R){
    colh <- sample(c(1:no.c), 1)
    Z[i, colh] <- 1
  }
  cull.p <- matrix(runif(n = no.c*species, min = 0.4, max = 0.6), nrow = no.c, ncol = species)
  countyN1 <- N1%*%Z  
  countyN2 <- N2%*%Z  
  countyN3 <- N3%*%Z  
  NCull <- rbind(countyN1, countyN2, countyN3)
  NCull <- t(NCull)
  cull <- matrix(data = NA, nrow = no.c, ncol = species)
  for(s in 1:species){
    cull[,s] <- rbinom(n = no.c, size = NCull[,s], prob = cull.p[,s])   
  }
  
### DATA FORMATTING FOR JAGS ----  
  R.diag <- diag(R)
  identity.mat <- diag(species)
  y <- list()
  y[[1]] <- y1
  y[[2]] <- y2
  y[[3]] <- y3
  
  BigN <- cbind(N1, N2, N3)
  R_zeros <- rep(0, species*R)
  
  # initialising N 
  max.ys <- matrix(NA, ncol = species, nrow = R)
  for(s in 1:species){
    max.ys[,s] <- apply(y[[s]], MARGIN = 1, FUN = max)+1
  }
#### JAGS ----
  ##### ROYLE -----
    JagsMod.Royle <- '
            model{
                constant.det = 1 + ((rho[1]^2)*(rho[2]^2)*(rho[3]^2)*(sd[1]^2)*(sd[2]^2)*(sd[3]^2)) - 
                               ((rho[1]^2)*(sd[1]^2)*(sd[2])^2) - ((rho[2]^2)*(sd[1]^2)*(sd[3]^2)) - 
                               ((rho[3]^2)*(sd[2]^2)*(sd[3]^2))
                for(i in 1:species){
                  for(j in 1:species){
                    h[i,j] = ifelse(i == j, 
                                    ifelse(i ==1,
                                           tau.s[1]*(1-(rho[3]^2*sd[2]^2*sd[3]^2)),
                                           ifelse(i==2, 
                                                  tau.s[2]*(1 - (rho[2]^2*sd[1]^2*sd[3]^2)),
                                                  tau.s[3]*(1 - (rho[1]^2*sd[1]^2*sd[2]^2))
                                                  )
                                          ),
                                    ifelse(i < j, 
                                           ifelse(i == 1 && j ==2, 
                                                  rho[2]*rho[3]*sd[3]^2 - rho[1], 
                                                  ifelse(i == 1 && j == 3, 
                                                         rho[1]*rho[3]*sd[2]^2 - rho[2],
                                                         rho[1]*rho[2]*sd[1]^2 - rho[3]
                                                         )
                                                  ), 
                                                  h[j,i]
                                            )
                                    )
                  T.mat[i,j] = h[i,j]/constant.det
                  } 
                }
                
                Q0 = D - alpha*A
                
                for(i in 1:(species*R)){
                  ind.i[i] = ifelse(i <= R, 1, ifelse(i>(2*R), 3, 2))
                  row.i[i] = ifelse(i <= R, i, ifelse(i>(2*R), (i - (2*R)),( i - R)))
                  for(j in 1:(species*R)){
                    ind.j[i,j] = ifelse(j <= R, 1, ifelse(j>(2*R), 3, 2))
                    col.j[i,j] = ifelse(j <=R, j, ifelse(j>(2*R), (j - (2*R)), (j - R)))
                    
                    Q[i,j] = Q0[row.i[i], col.j[i,j]]*T.mat[ind.i[i], ind.j[i,j]]
                  }
                }
                
                phi[c(1:(species*R))] ~ dmnorm(R_zeros, inverse(Q)) # Spatial Surface for 3 species
        
                for(s in 1:species){ 
                  S[c(1:R), s] <- ifelse(s==1, 
                                         phi[c(1:R)], 
                                         ifelse(s==2,
                                                phi[c((R+1):(2*R))], 
                                                phi[c(((R*2)+1):(3*R))]
                                                )
                                          )
                }
        
                # Equations
                for(i in 1:R){ # Over Grid Squares
                  N1[i] ~ dpois(lambda[i,1])
                  N2[i] ~ dpois(lambda[i,2])
                  N3[i] ~ dpois(lambda[i,3])
                  for(t in 1:max.T){ # Over Repeated Counts
                    y1[i,t] ~ dpois(p[1] * N1[i])
                    y2[i,t] ~ dpois(p[2] * N2[i])
                    y3[i,t] ~ dpois(p[3] * N3[i])
                  }
                }
                NC[c(1:no.c), 1] = N1[]%*%Z
                NC[c(1:no.c), 2] = N2[]%*%Z
                NC[c(1:no.c), 3] = N3[]%*%Z
                totalN[1] = sum(N1[])
                totalN[2] = sum(N2[])
                totalN[3] = sum(N3[])
                for(s in 1:species){ # Over Species
                  for(c in 1:no.c){
                    culls[c,s] ~ dpois(p.cull[c,s]*NC[c,s])
                    p.cull[c,s] ~ dunif(0.35, 0.65)
                  }
                  p[s] ~ dunif(0.01, 0.99)
                  llambda[c(1:R),s] = S[c(1:R), s] 
                  for(i in 1:R){
                    lambda[i,s] = exp(llambda[i,s])
                  }
                }
                tau ~ dgamma(2, 2)
                alpha ~ dunif(0.2, 0.99)
                
                for(s in 1:species){
                  rho[s] ~ dunif(-1, 1) 
                  sd[s] ~ dunif(0.01, 3)
                }
                
                for(s in 1:species){
                  tau.s[s] = 1/(sd[s]^2)
                }
                
              } 
            '
    jags.data <- list("y1" = y1, "y2" = y2, "y3" = y3,
                      "species" = species, 
                      "D" = D, "A" = A, "R" = R, 
                      "max.T" = max.T,
                      "R_zeros" = R_zeros,
                      "R.diag" = R.diag,
                      "culls" = cull,  "Z" = Z,
                      "no.c" = no.c
    )
    
    jags.inits <- function(){list("N1" = apply(y1, MARGIN = 1, FUN = max)+1,
                                  "N2" = apply(y2, MARGIN = 1, FUN = max)+1,
                                  "N3" = apply(y3, MARGIN = 1, FUN = max)+1,
                                  "rho" = c(-0.5, 0.2, 0.6), # rho12, rho13, rho23
                                  "alpha" = 0.5
    )
    }
    
    jags.pars <- c("p", "lambda", "phi",
                   "N1", "N2", "N3",
                   "alpha",  "tau",
                   "totalN", "p.cull"
    )
    
    mod1 <- jags(data = jags.data, inits = jags.inits,
                 parameters.to.save=jags.pars,
                 model.file = textConnection(JagsMod.Royle),
                 n.iter = 1000,
                 n.burnin = 500,
                 n.chains = 4,
                 n.thin = 1
    )
    
