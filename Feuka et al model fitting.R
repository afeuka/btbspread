### Title: Fitting dynamic occupancy models in parallel
### Author: Abbey Feuka
### Date: 06JAN22
### Notes: uses multiple p specification
### calculates AUC and bayesian p values post hoc
################
library(parallel)
library(pROC)

#spatial data (outputs of spatDatPrep_fun.R)
load("./Data/spatData 2021 clean.RData")

mod.names <- c(
  "null","dens","agri","temp",
  "dens agri", "dens temp", "agri temp",
  "dist","nb","nb dist",
  "dens dist", "agri dist", "temp dist",
  "dens nb", "agri nb", "temp nb","dens agri temp"
)
dens.logic <- agri.logic <- temp.logic <- 
  nb.logic <- dist.logic <- rep(FALSE,length(mod.names))
dens.logic[grep("dens",mod.names)] <- TRUE
agri.logic[grep("agri",mod.names)] <- TRUE
temp.logic[grep("temp",mod.names)] <- TRUE
nb.logic[grep("nb",mod.names)] <- TRUE
dist.logic[grep("dist",mod.names)] <- TRUE

#set up data according to model specification 
source("./Functions/modDataPrep_fun.R")
modDat <- modDataPrep("site","quarter",surv,site.info,samp.info,
                      dens.occ=dens.logic[idx],
                      agri.occ=agri.logic[idx],
                      temp.occ=temp.logic[idx],
                      surv.det=T) 

nsites <- modDat$n.sites #number of sites
nyrs <- modDat$n.yr #number of years
ntyps <- modDat$n.typs #number of surveillance/detection types
ntrials <- modDat$ntrials #number of tests

#data
y <- modDat$y #number of positives observed 
x <- modDat$x #covariates

nb <- nb.logic[idx] #neighborhood effect
dist <- dist.logic[idx] #distance effect
site.dist.sc <- scale(site.dist,center=F) #scaled distances

nChains <- 5 #number of chains
seed <- 1:nChains #random seeds
nmcmc <- 50000 #number of mcmc iterations
burnProp <- 0.1 #proportion of chain as burn in
beta.tune <- 0.1 #tuning parameter on beta coefficients

#fit chains in parallel
source("./Functions/fitDynOccMod_fun.R")
mcmcPar <- function(j){
  samp <- fitDynOccMod(seed,y,ntrials,x=x,adj,num,site.dist.sc, 
                       nmcmc,burnProp, beta.tune=beta.tune,
                       nb.gam=nb,inf.dist.gam=dist,postPred=T,
                       chain.idx=j)
}

#one chain per node
cl <- makeCluster(nChains, "SOCK")
clusterExport(cl, list("seed","y","x","ntrials","adj", "num","beta.tune",
                       "nb","dist","site.dist.sc","fitDynOccMod",
                       "nmcmc","burnProp","beta.tune","nb","dist"))

system.time(
  parSamples<- clusterApply(cl, 1:nChains, mcmcPar)
)
stopCluster(cl)

#save samples
save(parSamples,file=paste("./Model Outputs/pj/modfit mult p",mod.names[idx],"col all psi.RData"))
##############

#post hoc within sample validation
nmcmc <- nmcmc-(nmcmc*burnProp)
psi <- array(NA,dim=c(nsites,nyrs,nmcmc*nChains))
for(j in 1:nChains){
  psi[,,(nmcmc*(j-1)+1):(nmcmc*j)] <- parSamples[[j]]$psi.save
}
nChains <- length(parSamples)
pocc <- lapply(1:nChains,function(i){parSamples[[i]]$pocc.save})
pocc <- do.call("rbind",pocc)
ypred <- lapply(1:nChains,function(i){parSamples[[i]]$ypred.save})
ypred <- do.call("rbind",ypred)

# AUC for prevalance process
pocc.mn <- pocc %>% group_by(site,yr,typ) %>% 
  summarise(mn = mean(param))
pocc.mn<-as.data.frame(pocc.mn)
ybin <- ifelse(modDat$y$y>0,1,0)
auc.prev.pocc <- as.numeric(roc(ybin,pocc.mn$mn)$auc)

# AUC for occupancy process
psi.mn <- apply(psi,c(1,2),mean)
psi.mn.vec<-numeric(nrow(y))
for(i in 1:nrow(y)){
  psi.mn.vec[i] <- psi.mn[y$site[i],y$yr[i]]
}
auc.occ <- suppressMessages(as.numeric(roc(ybin,psi.mn.vec)$auc))

dev.y <- pocc %>% group_by(samp,chain) %>% 
  summarise(ll = sum(dbinom(y$y,ntrials$n,param,log=T)))
dev.y$ll <- (-2)*dev.y$ll
sapply(1:3,function(i){identical(pocc[,i],ypred[,i])})
ypred.pocc <- cbind(pocc,ypred$param)

dev.pred <- ypred.pocc %>% group_by(samp,chain) %>% 
  summarise(ll = sum(dbinom(`ypred$param`,ntrials$n,param,log=T)))
dev.pred$ll <- (-2)*dev.pred$ll

pVal <- sum(dev.y$ll>dev.pred$ll)/(nrow(dev.y))

save(pocc.mn,psi.mn,auc.prev.pocc,auc.occ,dev.y,dev.pred,pVal,
     file=paste0("./Model checks and comparisons/Within Sample/modval mult p ",
                 mod.names[idx]," col.RData"))

