### Title: Simulating dynamic occupancy data to test model 
### Author: Abbey Feuka
### Date: 19OCT22
### Notes:
####################
dynOccSim <- function(seed, #random seed for simulation reproducabilty
                      mat.dim.x=10,mat.dim.y=5, #dimensions of simulated raster
                      nyrs=5, #nubmer of years 
                      extinct, #extinction probability (1-persistence)
                      p, #detection probability/ies
                      trials.mn, #mean number of site visits/trials (intensity of Poisson)
                      nb=F, #neighborhood site effect
                      inf.dist=F){ #distance to nearest occupied site effect
  require(sf)
  require(boot)
  require(raster)
  require(nimble)
  
  ntyps <- length(p)
  #make simulated raster landscape
  r <- raster(extent(matrix(c(0,0,mat.dim.x,mat.dim.x),nrow=2)), 
              nrow=mat.dim.x, ncol=mat.dim.x)            
  r[] <- 1:ncell(r)
  
  #covert raster to polygons
  r.sp <- rasterToPolygons(r)
  r.sf <- st_as_sf(r.sp)
  
  #calculate adjacent cells
  site.neighbors <- st_intersects(r.sf,remove_self=T)
  adj <- as.carAdjacency(site.neighbors) # remove extra section creaeted by site.neighbors
  num<-adj$num
  adj<-adj$adj
  
  #distance matrix
  site.cent <- st_centroid(r.sf)
  site.dist <- st_distance(r.sf)
  
  #number of sites
  nsites <- nrow(r.sf)
  
  #landscape covariates and intercept
  x.static <- data.frame(site=rep(1:nsites,nyrs),
                         yr=sort(rep(1:nyrs,nsites)),
                         int=rep(1,nyrs*nsites),
                         x1=rnorm(nyrs*nsites,0,1))
  
  #number of trials/visits
  ntrials <- data.frame(site=sort(rep(1:nsites,nyrs*ntyps)),
                        yr=rep(sort(rep(1:nyrs,ntyps)),nsites),
                        typ=rep(1:ntyps,nsites*nyrs),
                        n=rpois(nsites*nyrs,trials.mn))
  #occupancy data structure
  y <- data.frame(site=sort(rep(1:nsites,nyrs*ntyps)),
                  yr=rep(sort(rep(1:nyrs,ntyps)),nsites),
                  typ=rep(1:ntyps,nsites*nyrs),
                  y=rep(NA,ntyps*nsites*nyrs))
  #decay paramter for distance function
  decay <- 0.5
  #calculate number of covariates
  if(nb==F & inf.dist==F){
    pX <- dim(x.static[,-c(1:2)])[2]
  } else if(nb==T & inf.dist==F){
    pX <- dim(x.static[,-c(1:2)])[2]+1
  }else if(nb==F & inf.dist==T){
    pX <- dim(x.static[,-c(1:2)])[2]+1
  }else if(nb==T & inf.dist==T){
    pX <- dim(x.static[,-c(1:2)])[2]+2
  } 
  set.seed(seed)
  
  #covariate relatinoships
  beta <- rnorm(pX,0,0.5)
  
  #first year
  psi <- gam <- z <- matrix(NA,nsites,nyrs)
  
  #initialize occupancy
  eps <- extinct
  psi[,1] <- rbeta(nsites,1,1)
  z[,1] <- rbinom(nsites,1,psi[,1])
  for(i in 1:nrow(y[y$yr==1,])){
    if(length(p)>1){
      y$y[y$yr==1][i] <- rbinom(1,ntrials$n[ntrials$yr==1][i],
                                p[y$typ[y$yr==1][i]]*z[y$site[y$yr==1][i],1])
    }else {
      y$y[y$yr==1][i] <- rbinom(1,ntrials$n[ntrials$yr==1][i],
                              p*z[y$site[y$yr==1][i],1])
    }
  }

  #following years
  for(t in 1:nyrs){
    if(t>1){
      psi[,t] <- (1-z[,t-1])*gam[,t-1] + (1-eps)*z[,t-1]
      z[,t] <- rbinom(nsites,1,psi[,t])
      for(i in 1:nrow(y[y$yr==t,])){
        if(length(p)>1){
          y$y[y$yr==t][i] <- rbinom(1,ntrials$n[ntrials$yr==t][i],
                                    p[y$typ[y$yr==t][i]]*z[y$site[y$yr==t][i],t])
        }else {
          y$y[y$yr==t][i] <- rbinom(1,ntrials$n[ntrials$yr==t][i],
                                    p*z[y$site[y$yr==t][i],t])
        }
      }
    }
    for(i in 1:nsites){
      if(inf.dist==T){ #if using distance effect
        infect.dist <- exp((-decay)*min(site.dist[-i,i][z[-i,t]==1]))
      }
      if(nb==T){ #if using neighborhood effect
        if(i==1){
          mn.nb <- mean(z[adj[1:(sum(num[1:i]))],t])
        } else {
          mn.nb <- mean(z[adj[(sum(num[1:(i-1)])+1):(sum(num[1:i]))],t])
        }
      }
      if(nb==T & inf.dist==T){ #adjusting x for spatial, model-estimated covariates
        x <- do.call("c",(c(x.static[x.static$site==i & x.static$yr==t,-c(1:2)],
                            infect.dist,mn.nb)))
      } else if(nb==T & inf.dist==F){
        x <- do.call("c",(c(x.static[x.static$site==i & x.static$yr==t,-c(1:2)],
                            mn.nb)))
      } else if(nb==F & inf.dist==T){
        x <- do.call("c",(c(x.static[x.static$site==i & x.static$yr==t,-c(1:2)],
                            infect.dist)))
      } else if(nb==F & inf.dist==F){
        x <- as.matrix(x.static[x.static$site==i & x.static$yr==t,-c(1:2)])
      }
      gam[i,t] <- inv.logit(x%*%beta)
    }
  }

  #add occupancy values to original dataframe
  for(t in 1:ncol(psi)){
    r.sf <- cbind(r.sf,psi[,t])
  }
  plot(r.sf["psi...t..2"])
  
  #create matrix of covariates
  x.mat <- matrix(NA,dim(psi)[1],dim(psi)[2])
  for(i in 1:nrow(x.static)){
    x.mat[x.static$site[i],x.static$yr[i]] <- x.static$x1[i]
  }

  list(psi.field=r.sf,y=y,gam=gam,psi=psi,eps=eps,beta=beta,decay=decay,
       p=p,x=x.static,x.mat=x.mat,z=z,ntrials=ntrials,nsites=nsites,
       nyrs=nyrs,ntyps=ntyps,adj=adj,num=num,site.dist=site.dist)
}


