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
#### FUNCTIONS ----
  st_rook = function(a, b = a) st_relate(a, b, pattern = "F***1****")
#### SET UP ----
  set.seed(20323)
  max.T <- 5
  R <- 25
  no.c <- 5
  species <- 2
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
  X_u <- matrix(rnorm(M*R, 0, 0.5), R, M) # Covariates
  X <- exp(X_u) # Ensuring positive covariates
    
  alpha.spatial <- 0.6
  tau.spatial <- c(3, 2)
  rho <- 0.2
  
  T11 <- tau.spatial[1]/(1 - rho^2)
  T12 <- (-rho*sqrt(tau.spatial[1])*sqrt(tau.spatial[2]))/(1-rho^2)
  T22 <- tau.spatial[2]/(1 - rho^2)
  Tau.matrix<- matrix(NA, nrow = species, ncol = species)
  for(i in 1:species){
    for(j in 1:species){
      if(i == j){
        Tau.matrix[i,j] <- ifelse(i==1, T11, T22)
      }else{
        Tau.matrix[i,j] <- T12
      }
    }
  }
  Q <- kronecker(Tau.matrix, D - alpha.spatial*A)
  
  phi <- mvrnorm(n = 1, mu = rep(0, nrow(Q)), Sigma = solve(Q)) # Spatial Surface for both species

  phi.1 <- phi[1:R] # Spatial Surface for species 1 
  phi.2 <- phi[(1+R):(2*R)] # Spatial Surface for species 2
                    
  lambda1 <- exp(phi.1)
  lambda2 <- exp(phi.2)

  N1 <- rpois(R, lambda1)  # Population species 1
  N2 <- rpois(R, lambda2)  # Population species 2

  p1 <- 0.3
  p2 <- 0.6

  y1 <- matrix(rbinom(n = R*max.T, size = N1, prob = p1), nrow = R, ncol = max.T) # Obs. species 1
  y2 <- matrix(rbinom(n = R*max.T, size = N2, prob = p2), nrow = R, ncol = max.T) # Obs. species 2

  ##### CULL DATA GENERATION -----
    Z = matrix(0, nrow = R, ncol = no.c) # Matrix, 1 in column for associated county
    for(i in 1:R){
      colh <- sample(c(1:no.c), 1)
      Z[i, colh] <- 1
    }
    cull.p <- matrix(runif(n = no.c*species, min = 0.4, max = 0.6), nrow = no.c, ncol = species)
    countyN1 <- N1%*%Z  
    countyN2 <- N2%*%Z  
    NCull <- rbind(countyN1, countyN2)
    NCull <- t(NCull)
    cull <- matrix(data = NA, nrow = no.c, ncol = species)
    for(s in 1:species){
      cull[,s] <- rbinom(n = no.c, size = NCull[,s], prob = cull.p[,s])   
    }

### DATA FORMATTING FOR JAGS ----  
  identity.mat <- diag(species)
  y <- list()
  y[[1]] <- y1
  y[[2]] <- y2

  BigN <- cbind(N1, N2)
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
            constant.det = 1 - rho^2
            for(i in 1:species){
              for(j in 1:species){
                h[i,j] = ifelse(i==j, ifelse(i==1, tau.s[1], tau.s[2]), -1*rho*sqrt(tau.s[1])*sqrt(tau.s[2]))
                T.mat[i,j] = h[i,j]/constant.det
              } 
            }
            
            Q0 = (alpha*(D-A)) + ((1-alpha)*R.diag)
            
            for(i in 1:(species*R)){
              ind.i[i] = ifelse(i <= R, 1, 2)
              row.i[i] = ifelse(i <= R, i, ( i - R))
              for(j in 1:(species*R)){
                ind.j[i,j] = ifelse(j <= R, 1,  2)
                col.j[i,j] = ifelse(j <=R, j, (j - R))
                
                Q[i,j] = Q0[row.i[i], col.j[i,j]]*T.mat[ind.i[i], ind.j[i,j]]
              }
            }
            
            phi[c(1:(species*R))] ~ dmnorm(R_zeros, inverse(Q))
    
            for(s in 1:species){
              S[c(1:R), s] <- ifelse(s==1, phi[c(1:R)], phi[c((R+1):(2*R))])
            }
    
            # Equations
            for(i in 1:R){ # Over Grid Squares
              N1[i] ~ dpois(lambda[i,1])
              N2[i] ~ dpois(lambda[i,2])
              for(t in 1:max.T){ # Over Repeated Counts
                y1[i,t] ~ dpois(p[1] * N1[i])
                y2[i,t] ~ dpois(p[2] * N2[i])
              }
            }
            NC[c(1:no.c), 1] = N1[]%*%Z
            NC[c(1:no.c), 2] = N2[]%*%Z
            totalN[1] = sum(N1[])
            totalN[2] = sum(N2[])
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
            alpha ~ dunif(0.2, 0.99)
            
            rho ~ dunif(-1, 1) 
            for(s in 1:species){
              sd[s] ~ dunif(0.01, 3)
              sd2[s] = sd[s]^2
              tau.s[s] = 1/sd2[s]
            }
            
          } 
        '
    jags.data <- list("y1" = y1, "y2" = y2, 
                      "species" = species, 
                      "D" = D, "A" = A,
                      "R" = R,  "max.T" = max.T,
                      "R_zeros" = R_zeros,
                      "culls" = cull, 
                      "Z" = Z, "no.c" = no.c
    )
    
    jags.inits <- function(){list(N1 = apply(y1, MARGIN = 1, FUN = max)+1,
                                  N2 = apply(y2, MARGIN = 1, FUN = max)+1,
                                  rho = 0.5, 
                                  alpha = 0.5
    )
    }
    
    jags.pars <- c("p",
                   "lambda",
                   "phi",
                   "N1", "N2",
                   "alpha", 
                   "tau.s",
                   "totalN",
                   "p.cull"
    )
    
    mod1 <- jags(data = jags.data, inits = jags.inits,
                 parameters.to.save=jags.pars,
                 model.file = textConnection(JagsMod.Royle),
                 n.iter = 1000,
                 n.burnin = 500,
                 n.chains = 4,
                 n.thin = 1
    )
    
