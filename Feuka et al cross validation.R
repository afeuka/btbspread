### Title: Leave one out cross validation for dynamic occupancy models
### Author: Abbey Feuka
### Date: 15DEC22
### Notes: 
#######################
library(parallel)

load("./Data/spatData 2020.RData")

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

source("./Functions/modDataPrep_fun.R")
for(idx in 1:length(mod.names)){
modDat <- modDataPrep("site","quarter",surv,site.info,samp.info,
                      dens.occ=dens.logic[idx],
                      agri.occ=agri.logic[idx],
                      temp.occ=temp.logic[idx],
                      surv.det=T) 

nsites <- modDat$n.sites #number of sites
nyrs <- modDat$n.yr #number of years (full data set)
ntyps <- modDat$n.typs #number of surveillance/detection types
ntrials <- modDat$ntrials #number of tests/site visits per year

y <- modDat$y #number of positive observations
x <- modDat$x #covariates

#used tuning parameters from model fit
btune <- 0.1

nb <- nb.logic[idx] #neighborhood effect
dist <- dist.logic[idx] #distance effect

#fit cross validation models in parallel
#one model per year using multiple nodes
source("./Functions/fitDynOccMod_pj_loo_fun.R")
loo.par <- function(t){
  nmcmc <- 50000
  burnProp <- 0.1
  n.eff <- nmcmc-(nmcmc*burnProp)
  
  modfit.loo <- fitDynOccModLoo(loo=t,seed=1,y=y, x=x, ntrials=ntrials,
                                adj=adj,num=num,site.dist=site.dist.sc,
                                nmcmc=nmcmc, burnProp=0.1,
                                nb.gam=nb,inf.dist.gam=dist,decay="fixed")
}

cl <- makeCluster(5, "SOCK")
clusterExport(cl, list("y","x","ntrials","adj", "num","btune",
                       "nb","dist","site.dist.sc","fitDynOccModLoo"))

system.time(
  parSamplesLoo<- clusterApply(cl, 1:nyrs, loo.par)
)
stopCluster(cl)

##-2xloglikelihood
cv <- sapply(1:nyrs,function(i){
  -2*sum(log(colMeans(parSamplesLoo[[i]]$loglik.oos.save)))})
cv.mn <- mean(cv)
cv.tot <- sum(cv)

##out of sample AUC
z.fit <- lapply(1:nyrs,function(i)parSamplesLoo[[i]]$z.save)

nmcmc <- 45000
pocc <- lapply(1:nyrs,function(t){lapply(1:nmcmc,function(i){
  subset(parSamplesLoo[[t]]$pocc.save[[i]],yr==t)})})
pocc <- lapply(1:nyrs,function(t){do.call("rbind",pocc[[t]])})

# AUC for prevalance process from pocc samples 
pocc.mn <- lapply(1:nyrs,function(t){
  suppressMessages(pocc[[t]] %>% group_by(site,typ) %>% 
                     summarise(mn = mean(param)))})
identical(ntrials[ntrials$yr==1,]$site,pocc.mn[[1]]$site)
identical(ntrials[ntrials$yr==1,]$typ,pocc.mn[[1]]$typ)
ybin <- ifelse(modDat$y$y>0,1,0)
pocc.mn.vec <- do.call("rbind",pocc.mn)
auc.prev.oos.pocc <- as.numeric(roc(ybin,pocc.mn.vec$mn)$auc)

# AUC for occupancy process
psi <- array(NA,dim=c(nsites,nyrs,nmcmc))
for(t in 1:nyrs){ # out of sample psi's
  psi[,t,] <- parSamplesLooPsi[[t]]$psi.save[,t,]
}
psi.mn <- apply(psi,c(1,2),mean)
psi.mn.vec<-numeric(nrow(y))
for(i in 1:nrow(y)){
  psi.mn.vec[i] <- psi.mn[y$site[i],y$yr[i]]
}
ybin <- ifelse(y$y>0,1,0)
auc.occ.oos <- suppressMessages(as.numeric(roc(ybin,psi.mn.vec)$auc))

#oos log like
save(cv,cv.mn,cv.tot,
  file=paste0("./Model checks and comparisons/LOO/All Years/CV Scores/Colonization/LogLik/cv scores loo 10yrs mult p ",mod.names[idx]," col.RData"))

#oos auc
save(pocc.mn,auc.occ.oos,auc.prev.oos.pocc,
     file=paste("./Model checks and comparisons/LOO/All Years/CV Scores/Colonization/AUC/AUC cv scores loo 10yrs mult p",
                mod.names[idx],"col.RData"))
#samples
save(parSamplesLoo,
  file=paste0("./Model checks and comparisons/LOO/All Years/Samples/Colonization/samples loo 10yrs mult p ",mod.names[idx]," col.RData"))
}
