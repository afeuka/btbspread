### Title: Leave one out cross validation using dynamic occupancy MCMC
### Author: Abbey Feuka
### Date: 15DEC22
### Notes: saves posterior predictions automatically
######################
fitDynOccModLoo <- function(loo, #year index to hold out from model fitting
                            seed, #random seed for initial values
                            y, #number of positive tests in site in each year (nxt)
                            ntrials, #number of tests in each site in each year
                            x=NA, #covariates on colonization prob
                            adj, #adjacency matrix, vectorized
                            num, #index vector for adjacency matrix
                            site.dist, #distance matrix between observations
                            
                            nmcmc, #number of mcmc iterations
                            burnProp, #proportion of samples as burn in
                            beta.tune=0.6,#tuning parameter for beta (occupancy)
                            
                            nb.gam=F, #include neighborhood coefficient
                            inf.dist.gam=F,#include infection distance coefficient
                            decay="estimated"#fix decay rate on distance effect or estimate it
                            
){
  require(boot)
  require(tidyverse)
  
  btune <- beta.acc <- NA
  
  #year to leave out of data specified by loo
  y.loo <- y
  y.loo$y[y.loo$yr==loo] <- rep(NA,length(y.loo$y[y.loo$yr==loo]))
  
  #function for covariate calculations
  inv.logit.mat.beta<-function(cov,coeff,na.rm=F){
    #creates nsite x nyr matrix of coeffs
    cov.prod <- sapply(1:length(unique(cov[,2])),function(t){
      sapply(1:length(unique(cov[,1])),function(i){
        inv.logit(as.numeric(cov[cov[,1]==i & cov[,2]==t,-c(1:2)])%*%coeff)
      })})
    if(na.rm==T){
      cov.prod <- c(cov.prod)
      cov.prod <- cov.prod[!is.na(cov.prod)]
    }
    return(cov.prod)
  }
  
  nsites <- length(unique(y$site)) #number of sites
  nyrs <- length(unique(y$yr)) #number of years (full dataset)
  ntyps <- length(unique(y$typ)) #number of surveillance/detection types
  
  nmcmc.eff <- nmcmc - (nmcmc*burnProp) #number of effective mcmc samples
  ntsamp <- 1000 #tuning samples
  
  #count covariates
  if(!is.na(sum(x))){
    x.static <- subset(x,yr!=nyrs)
    if(dim(x)[2]>3){ #covariates
      pX <- dim(x[,-c(1:2)])[2]
    } else {
      pX <- 1
    }
    if(nb.gam==T & inf.dist.gam==F){
      pX<-pX+1
    } else if(nb.gam==F & inf.dist.gam==T){
      pX<-pX+1
    } else if(nb.gam==T & inf.dist.gam==T){
      pX<-pX+2
    } 
    #colonization coefficients
    beta.save <- matrix(NA,nmcmc.eff,pX)
    beta.acc <- numeric(pX)
    b.tune.save <- matrix(NA,nmcmc,pX)
  }
  #detection probabilities
  p <- rbeta(ntyps,1,1)
  
  #sample storage
  ypred.save <- list()
  for(i in 1:nmcmc.eff){
    ypred.save[[i]] <- cbind.data.frame(y[,1:3],param=rep(NA,nrow(y)))
  }
  pocc.save <- ypred.save
  ypred <- y
  colnames(ypred)[colnames(ypred)=="y"] <- "ypred"
  
  z.save <- array(NA,dim=c(nsites,nyrs,nmcmc.eff))
  psi.save <- z.save <- array(NA,dim=c(nsites,nyrs,nmcmc.eff))
  gam.save <- array(NA,dim=c(nsites,nyrs-1,nmcmc.eff))
  eps.save <- rep(NA,nmcmc.eff)
  p.save <- matrix(NA,nmcmc.eff,ntyps)
  loglik.oos.save <- matrix(NA,nmcmc.eff,length(subset(y.loo,yr==loo)$y))
  
  #priors
  beta.mn <- 0
  beta.sd <- 1
  alpha.p <- 1
  beta.p <-1 
  
  #initial values
  z <- matrix(0,nsites,nyrs)
  y.site.yr <- as.data.frame(suppressMessages(y %>% 
                                                group_by(site,yr) %>% 
                                                summarise(y = sum(y,na.rm=T))))
  ntrials.site.yr <- as.data.frame(suppressMessages(ntrials %>% 
                                                      group_by(site,yr) %>% 
                                                      summarise(n = sum(n,na.rm=T))))
  possites <- y.site.yr$site[y.site.yr$y>0]
  posyears <- y.site.yr$yr[y.site.yr$y>0]
  for(i in 1:length(possites)){z[possites[i],posyears[i]] <-1}
  
  if(!is.na(sum(x))){
    if(ncol(x)==3){
      beta <- as.vector(glm(c(z) ~ x[,-(1:2)],family=binomial())$coefficients)[1:pX]
    } else if (ncol(x)==5){
      beta <- as.vector(glm(c(z) ~ x[,4] + x[,5],family=binomial())$coefficients)[1:pX]
    } else if (ncol(x)==6){
      beta <- as.vector(glm(c(z) ~ x[,4] + x[,5] + x[,6],family=binomial())$coefficients)[1:pX]
    } else {
      beta <- as.vector(glm(c(z) ~ x[,-(1:3)],family=binomial())$coefficients)[1:pX]
    }
    beta <- beta[!is.na(beta)]
    if(nb.gam==T & inf.dist.gam==T){
      beta <- c(beta,rnorm(3,0,1))
    } else if(nb.gam==F & inf.dist.gam==T){
      beta <- c(beta,rnorm(2,0,1))
    } else if(nb.gam==T & inf.dist.gam==F){
      beta <- c(beta,rnorm(1,0,1))
    }
    beta <- rnorm(pX,beta,0.05)
    btune <- beta.tune*abs(beta)
    beta.save <- matrix(NA,nmcmc.eff,pX)
  }

  #pocc = occupancy * detection
  pocc <- y
  colnames(pocc)[length(colnames(pocc))] <- "pocc"
  pocc$pocc <- NA
  pocc$pocc <- sapply(1:nrow(pocc),function(i)p[pocc$typ[i]]*z[pocc$site[i],pocc$yr[i]])
  
  set.seed(seed)
  psi <- matrix(rbeta(nsites*nyrs,1,1),nsites,nyrs)
  gam <- matrix(rbeta(nsites*nyrs,1,1),nsites,nyrs)
  eps <- rbeta(1,1,nsites,nyrs)
  
  #initialize spatial covariates
  if(!is.na(sum(x))){
    if(inf.dist.gam==T){
      if(decay=="fixed"){
        infect.dist <- sapply(2:nyrs,function(t){ #verified same as for loop
          sapply(1:nsites,function(i){
            exp((-0.2)*min(site.dist[-i,i][z[-i,t-1]==1]))})})
      } else {
        infect.dist <- sapply(2:nyrs,function(t){ #verified same as for loop
          sapply(1:nsites,function(i){
            exp((-beta[pX])*min(site.dist[-i,i][z[-i,t-1]==1]))})})
      }
    }
    if(nb.gam==T){
      #mean infection status of neighboring sites
      mn.nb <- sapply(2:nyrs,function(t){ #verified same as for loop
        sapply(1:nsites,function(i){
          if(i==1){
            mean(z[adj[1:(sum(num[1:i]))],t-1])
          } else {
            mean(z[adj[(sum(num[1:(i-1)])+1):(sum(num[1:i]))],t-1])
          }})})
    }
    if(nb.gam==T & inf.dist.gam==T){
      x <- cbind(x.static,c(t(infect.dist)),c(t(mn.nb)))
    } else if(nb.gam==T & inf.dist.gam==F){
      x <- cbind(x.static,c(t(mn.nb)))
    } else if(nb.gam==F & inf.dist.gam==T){
      x <- cbind(x.static,c(t(infect.dist)))
    } else {
      x <- x.static
    }
  }
  
  ### MCMC ###
  for(k in 1:nmcmc){
    if(k%%1000==0)cat(k," ")
    
    #psi 1
    alpha.psi1 <- sum(z[,1])+1
    beta.psi1 <- sum(1-z[,1])+1
    psi[,1] <- rbeta(nsites,alpha.psi1,beta.psi1)
    
    #epsilon (for z[,t-1]==1)
    alpha.tmp <-sum(sapply(2:nyrs,function(t)(sum(z[z[,t-1]==1,t]))))
    beta.tmp <-sum(sapply(2:nyrs,function(t)(sum(1-z[z[,t-1]==1,t]))))
    eps <- rbeta(1,alpha.tmp+1,beta.tmp+1)
    
    ##betas 
    #(gamma coeffs) for z[,t-1]==0 (z[,t] estimated using gam[,t-1])
    if(is.na(sum(x))){
      alpha.tmp <- sum(sapply(2:nyrs,function(t)(sum(z[z[,t-1]==0,t]))))
      beta.tmp <- sum(sapply(2:nyrs,function(t)(sum(1-z[z[,t-1]==0,t]))))
      gam <- matrix(rbeta(1,alpha.tmp+1,beta.tmp+1),nsites,nyrs-1)
    } else{
      #beta update (each individually, vectorized)
      beta.star <- rnorm(pX,beta,btune)
      if(inf.dist.gam==T){
        if(decay=="fixed"){
          infect.dist.star <- sapply(2:nyrs,function(t){ #verified same as for loop
            sapply(1:nsites,function(i){
              exp((-0.2)*min(site.dist[-i,i][z[-i,t-1]==1]))})})
        } else {
          infect.dist.star <- sapply(2:nyrs,function(t){ #verified same as for loop
            sapply(1:nsites,function(i){
              exp((-beta.star[pX])*min(site.dist[-i,i][z[-i,t-1]==1]))})})
        }
      }
      if(nb.gam==F & inf.dist.gam==T){
        x.star <- cbind(x.static,c(t(infect.dist.star)))
        x.tmp.star <- lapply(2:nyrs,function(t){
          as.matrix(x.star[x$site%in%which(z[,t-1]==0) & x.star$yr==(t-1),-c(1:2)])})
      }
      if(nb.gam==T & inf.dist.gam==T){
        x.star <- cbind(x.static,c(t(infect.dist.star)),c(t(mn.nb)))
        x.tmp.star <- lapply(2:nyrs,function(t){
          as.matrix(x.star[x$site%in%which(z[,t-1]==0) & x.star$yr==(t-1),-c(1:2)])})
      }
      x.tmp <- lapply(2:nyrs,function(t){
        as.matrix(x[x$site%in%which(z[,t-1]==0) & x$yr==(t-1),-c(1:2)])})
      
      if(inf.dist.gam==T){ #mh update for distance
        if(decay=="fixed"){
          mh1 <- rowSums(sapply(2:nyrs,function(t){
            sum(dbinom(z[z[,t-1]==0,t],1,inv.logit(x.tmp.star[[t-1]]%*%beta.star),log=T))+
              dnorm(beta.star,beta.mn,beta.sd,log=T)}))
          mh2 <- rowSums(sapply(2:nyrs,function(t){
            sum(dbinom(z[z[,t-1]==0,t],1,inv.logit(x.tmp[[t-1]]%*%beta),log=T))+
              dnorm(beta,beta.mn,beta.sd,log=T)}))
        } else {
          mh1 <- rowSums(sapply(2:nyrs,function(t){
            sum(dbinom(z[z[,t-1]==0,t],1,inv.logit(x.tmp.star[[t-1]]%*%beta.star[1:(pX-1)]),log=T))+
              dnorm(beta.star,beta.mn,beta.sd,log=T)}))
          mh2 <- rowSums(sapply(2:nyrs,function(t){
            sum(dbinom(z[z[,t-1]==0,t],1,inv.logit(x.tmp[[t-1]]%*%beta[1:(pX-1)]),log=T))+
              dnorm(beta,beta.mn,beta.sd,log=T)}))
        }
      } else { #mh update for no distance
        if(pX!=1){ #model with int + covariates
          mh1 <- rowSums(sapply(2:nyrs,function(t){
            sum(dbinom(z[z[,t-1]==0,t],1,inv.logit(x.tmp[[t-1]]%*%beta.star),log=T))+
              dnorm(beta.star,beta.mn,beta.sd,log=T)}))
          mh2 <- rowSums(sapply(2:nyrs,function(t){
            sum(dbinom(z[z[,t-1]==0,t],1,inv.logit(x.tmp[[t-1]]%*%beta),log=T))+
              dnorm(beta,beta.mn,beta.sd,log=T)}))
        } else { #model with int only
          mh1 <- sum(sapply(2:nyrs,function(t){
            sum(dbinom(z[z[,t-1]==0,t],1,inv.logit(x.tmp[[t-1]]%*%beta.star),log=T))+
              dnorm(beta.star,beta.mn,beta.sd,log=T)}))
          mh2 <- sum(sapply(2:nyrs,function(t){
            sum(dbinom(z[z[,t-1]==0,t],1,inv.logit(x.tmp[[t-1]]%*%beta),log=T))+
              dnorm(beta,beta.mn,beta.sd,log=T)}))
        }
      }
      mh <- exp(mh1-mh2)
      keep.idx <- mh>runif(pX)
      beta[keep.idx] <- beta.star[keep.idx]
      beta.acc[keep.idx] <- beta.acc[keep.idx]+1
      
      if(inf.dist.gam==T){
        #transmission distance
        if(decay=="fixed"){
          infect.dist <- sapply(2:nyrs,function(t){ #verified same as for loop
            sapply(1:nsites,function(i){
              exp((-0.2)*min(site.dist[-i,i][z[-i,t-1]==1]))})})
        } else {
          infect.dist <- sapply(2:nyrs,function(t){ #verified same as for loop
            sapply(1:nsites,function(i){
              exp((-beta[pX])*min(site.dist[-i,i][z[-i,t-1]==1]))})})
        }
      }
      if(nb.gam==T & inf.dist.gam==T){
        x <- cbind(x.static,c(t(infect.dist)),c(t(mn.nb)))
      } else if(nb.gam==T & inf.dist.gam==F){
        x <- cbind(x.static,c(t(mn.nb)))
      } else if(nb.gam==F & inf.dist.gam==T){
        x <- cbind(x.static,c(t(infect.dist)))
      }
      if(inf.dist.gam==T){
        gam <- inv.logit.mat.beta(as.matrix(x),beta[1:(pX-1)])
      } else {
        gam <- inv.logit.mat.beta(as.matrix(x),beta)
      }
    }#end regression logic (MH update for betas)

    ##p - detection occupied sites only
    occ.idx <- sapply(1:nrow(y),function(i){z[y$site[i],y$yr[i]]==1})
    pocc.tmp <- pocc[occ.idx,] 
    pocc.tmp <- pocc.tmp[pocc.tmp$yr!=loo,]
    y.tmp <- y[occ.idx,]
    y.tmp <- y.tmp[y.tmp$yr!=loo,]
    ntrials.tmp <- ntrials[occ.idx,]
    ntrials.tmp <- ntrials.tmp[ntrials.tmp$yr!=loo,]
    for(i in 1:ntyps){
      alpha.p.tmp <- sum(y.tmp$y[y.tmp$typ==i])+alpha.p
      beta.p.tmp <- sum(ntrials.tmp$n[ntrials.tmp$typ==i]-y.tmp$y[y.tmp$typ==i])+beta.p
      p[i] <- rbeta(1,alpha.p.tmp,beta.p.tmp)
    }
    pocc$pocc <- sapply(1:nrow(pocc),function(i)p[pocc$typ[i]]*z[pocc$site[i],pocc$yr[i]])
    
    #z for sites with no detections
    for(i in 1:nsites){
      for(t in 1:nyrs){
        if(sum(y.site.yr$site==i & y.site.yr$yr==t)>0){ #site-year sampled
          if(y.site.yr$y[y.site.yr$site==i & y.site.yr$yr==t]==0){
            if(t==1){ #first timestep
              if(z[i,t+1]==0){
                psi.star <- (1-eps)/((1-eps)+(1-gam[i,t]))
              } else if(z[i,t+1]==1){
                psi.star <- eps/(eps+gam[i,t])
              }
            } else if(t>1 & t<nyrs){ #middle timesteps
              if(z[i,t-1]==1 & z[i,t+1]==1){
                psi.star <- (eps*eps)/((eps*eps)+(1-eps)*gam[i,t])
              } else if(z[i,t-1]==1 & z[i,t+1]==0){
                psi.star <- (eps*(1-eps))/(eps*(1-eps)+(1-eps)*(1-gam[i,t]))
              } else if(z[i,t-1]==0 & z[i,t+1]==1){
                psi.star <- (gam[i,t-1]*eps)/(gam[i,t-1]*eps+(1-gam[i,t-1])*gam[i,t])
              } else if(z[i,t-1]==0 & z[i,t+1]==0){
                psi.star <- (gam[i,t-1]*(1-eps))/(gam[i,t-1]*(1-eps)+(1-gam[i,t-1])*(1-gam[i,t]))
              }
            } else if(t==nyrs){ #last timestep
              if(z[i,t-1]==0){
                psi.star <- gam[i,t-1]/(gam[i,t-1]+(1-gam[i,t-1]))
              } else if(z[i,t-1]==1){
                psi.star <- eps/(eps+(1-eps))
              }
            }
            z.tmp <- prod(unlist(sapply(1:ntyps,function(j){
              ((1-p[j])^ntrials$n[ntrials$site==i & ntrials$yr==t & ntrials$typ==j])})))*psi.star
            z.p<- z.tmp/(z.tmp+(1-psi.star))
            z[i,t] <- rbinom(1,1,z.p)  #only z[y==0] updated
          }
          if(t>1){
            psi[i,t] <- gam[i,t-1]*(1-z[i,t-1]) + z[i,t-1]*eps
          }
        }# end site-year sampled condition
      }# end year loop
    }# end site loop
    
    if(k > (nmcmc*burnProp)){ #don't save burn in
      #posterior predictions
      ypred$ypred <- rbinom(nrow(ypred),ntrials$n,pocc$pocc)
      loglik.oos <- dbinom(subset(y,yr==loo)$y,subset(ntrials,yr==loo)$n,
                           subset(pocc,yr==loo)$pocc)
      q <- k-(nmcmc*burnProp)
      #save samples
      psi.save[,,q] <- psi
      z.save[,,q] <- z
      pocc.save[[q]]$param <- pocc$pocc
      p.save[q,] <- p
      ypred.save[[q]]$param <- ypred$ypred
      loglik.oos.save[q,] <- loglik.oos
    } #end burn in condition
  }#end mcmc loop
  
  if(!is.na(sum(x))){beta.acc <- beta.acc/nmcmc}
  
  list(psi.save=psi.save,
    loglik.oos.save=loglik.oos.save,
    z.save=z.save,p.save=p.save,
    pocc.save=pocc.save,ypred.save=ypred.save,
    beta.acc=beta.acc,nmcmc=nmcmc.eff)
}
