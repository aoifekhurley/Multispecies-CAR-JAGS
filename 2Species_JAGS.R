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
  M <- 3      # Number of Covariates
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
  # Sigma <- rinvwishart(nu = species+1, S = diag(species))
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
  #Q <- kronecker(alpha.spatial*(D-A) + diag((1 - alpha.spatial), nrow = nrow(A)), solve(Sigma))
  
  phi <- mvrnorm(n = 1, mu = rep(0, nrow(Q)), Sigma = solve(Q))

  phi.1 <- phi[1:R]
  phi.2 <- phi[(1+R):(2*R)] 
                    
  lambda1 <- exp(phi.1)
  lambda2 <- exp(phi.2)

  N1 <- rpois(R, lambda1)  
  N2 <- rpois(R, lambda2)  

  p1 <- 0.3
  p2 <- 0.6

  y1 <- matrix(rbinom(n = R*max.T, size = N1, prob = p1), nrow = R, ncol = max.T)
  y2 <- matrix(rbinom(n = R*max.T, size = N2, prob = p2), nrow = R, ncol = max.T)

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
  R.diag <- diag(R)
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
    jags.data <- list("y1" = y1,
                      "y2" = y2, 
                      "species" = species, 
                      "D" = D,
                      "A" = A,
                      "R" = R, 
                      "max.T" = max.T,
                      "R_zeros" = R_zeros,
                      "R.diag" = R.diag,
                      "culls" = cull, 
                      "Z" = Z,
                      "no.c" = no.c
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
    tictoc::tic()
    mod1 <- jags(data = jags.data, inits = jags.inits,
                 parameters.to.save=jags.pars,
                 model.file = textConnection(JagsMod.Royle),
                 n.iter = 1000,
                 n.burnin = 500,
                 n.chains = 4,
                 n.thin = 1
    )
    tictoc::toc()
    
#### Visualising Priors ----
  m1 <- mod1$BUGSoutput$sims.matrix
  post.samples <- as_tibble(m1)
    
  # plot detection probability parameter
  ggplot(post.samples, aes(x = `p[1]`)) +
    stat_halfeye() +
    geom_vline(xintercept = p1, colour = "blue", linetype = "dashed") +
    theme_classic()
    
  ggplot(post.samples, aes(x = `lambda[1,1]`)) +
    stat_halfeye() +
    geom_vline(xintercept = lambda1[1], colour = "blue", linetype = "dashed") +
    theme_classic()
    
  ggplot(post.samples, aes(x = `lambda[2,1]`)) +
    stat_halfeye() +
    geom_vline(xintercept = lambda2[1], colour = "blue", linetype = "dashed") +
    theme_classic()
    
  ggplot(post.samples, aes(x = `alpha`)) +
    stat_halfeye() +
    geom_vline(xintercept = alpha.spatial, colour = "blue", linetype = "dashed") +
    theme_classic()
    
  ggplot(post.samples, aes(x = `tau.s[1]`)) +
    stat_halfeye() +
    geom_vline(xintercept = tau.spatial[1], colour = "blue", linetype = "dashed") +
    theme_classic()
    
  quantiles.samples <- apply(X = post.samples, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.975, 0.995, 1))
    
  p1.quantiles <- quantiles.samples[ , grepl("^p\\[1]$", colnames(quantiles.samples))]
  p1.quantiles <- c(p1.quantiles, p1)
  p1.quantiles <- as.data.frame(t(p1.quantiles))
  colnames(p1.quantiles)[9] <- "Actual"
  sum(p1.quantiles$Actual > p1.quantiles$`2.5%` & p1.quantiles$Actual < p1.quantiles$`97.5%`)
  p1.quantiles$Variable <- "p1"
  
  ggplot(p1.quantiles, aes(x = Variable)) + 
    geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`)) + 
    geom_point(aes(y=`50%`)) + 
    geom_point(aes(y=Actual), color = "red", shape = 8) + 
    labs(x="", y = "") + 
    theme_classic()
  
  p2.quantiles <- quantiles.samples[ , grepl("^p\\[2]$", colnames(quantiles.samples))]
  p2.quantiles <- c(p2.quantiles, p2)
  p2.quantiles <- as.data.frame(t(p2.quantiles))
  colnames(p2.quantiles)[9] <- "Actual"
  sum(p2.quantiles$Actual > p2.quantiles$`2.5%` & p2.quantiles$Actual < p2.quantiles$`97.5%`)
  p2.quantiles$Variable <- "p2"
  
  ggplot(p2.quantiles, aes(x = Variable)) + 
    geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`)) + 
    geom_point(aes(y=`50%`)) + 
    geom_point(aes(y=Actual), color = "red", shape = 8) + 
    labs(x="", y = "") + 
    theme_classic()
  
  p.cull.q <- quantiles.samples[ , grepl("^p.cull", colnames(quantiles.samples))]
  
  p.cull.1 <- p.cull.q[,1:no.c]
  p.cull.1 <- rbind(p.cull.1, cull.p[,1])
  p.cull.1 <- as.data.frame(t(p.cull.1))
  colnames(p.cull.1)[9] <- "Actual"
  sum(p.cull.1$Actual > p.cull.1$`2.5%` & p.cull.1$Actual < p.cull.1$`97.5%`)
  p.cull.1$Variable <- paste0("Cull Percentages ", seq(1:no.c))
  
  ggplot(p.cull.1, aes(x = Variable)) + 
    geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`)) + 
    geom_point(aes(y=`50%`)) + 
    geom_point(aes(y=Actual), color = "red", shape = 8) + 
    labs(x="", y = "") + 
    theme_classic()
  
  p.cull.2 <- p.cull.q[, (no.c+1):(2*no.c)]
  p.cull.2 <- rbind(p.cull.2, cull.p[,2])
  p.cull.2 <- as.data.frame(t(p.cull.2))
  colnames(p.cull.2)[9] <- "Actual"
  sum(p.cull.2$Actual > p.cull.2$`2.5%` & p.cull.2$Actual < p.cull.2$`97.5%`)
  p.cull.2$Variable <- paste0("Cull Percentages ", seq(1:no.c))
  
  ggplot(p.cull.2, aes(x = Variable)) + 
    geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`)) + 
    geom_point(aes(y=`50%`)) + 
    geom_point(aes(y=Actual), color = "red", shape = 8) + 
    labs(x="", y = "") + 
    theme_classic()
  
  phi.q <- quantiles.samples[ , grepl("^phi\\[", colnames(quantiles.samples))]
  phi.q <- rbind(phi.q, phi)
  phi.q <- as.data.frame(t(phi.q))
  colnames(phi.q)[9] <- "Actual"
  sum(phi.q$Actual > phi.q$`2.5%` & phi.q$Actual < phi.q$`97.5%`)
  sum(phi.q$Actual > phi.q$`2.5%` & phi.q$Actual < phi.q$`97.5%`)/nrow(phi.q)
  mean(phi.q$Actual - phi.q$`50%`)
  
  BigN.q <- quantiles.samples[, grepl("^totalN", colnames(quantiles.samples))]
  BigN.q <- rbind(BigN.q,c(sum(N1), sum(N2)))
  BigN.q <- as.data.frame(t(BigN.q))
  colnames(BigN.q)[9] <- "Actual"
  sum(BigN.q$Actual > BigN.q$`2.5%` & BigN.q$Actual < BigN.q$`97.5%`)
  BigN.q$Variable <- c("N1", "N2")
  
  ggplot(BigN.q, aes(x = Variable)) + 
    geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`)) + 
    geom_point(aes(y=`50%`)) + 
    geom_point(aes(y=Actual), color = "red", shape = 8) + 
    labs(x="", y = "") + 
    theme_classic()
  
  