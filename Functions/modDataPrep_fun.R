### Title: bTB transmission model data prep, using spatial data
### Author: Abbey Feuka
### Date: 18NOV22
### Notes:
#########################################
modDataPrep <- function(siteLevel, # scale of occupancy site (section,quarter,or township)
                        surv, #surveillance data, sample level (from spatDataPrep)
                        site.info, #all sites, regardless of surveyed, landscape covs (from spatDataPrep)
                        samp.info, #only surveyed sites, sample-level covs (from spatDataPrep)
                        dens.occ=F, #colonization as a function of host density
                        agri.occ=F, #colonization as a function of proportion of agriculture
                        temp.occ=F, #colonization as a function of days below 0C
                        surv.det=T){ #separate by survey type

  require(tidyverse)
  require(sf)
  
  #data attributes
  years <- sort(unique(site.info$Year))
  site.info$year.idx <- as.numeric(as.factor(site.info$Year))
  sites <- sort(unique(site.info$SITE))
  n.sites <- length(unique(site.info$site.idx))
  n.yrs <- length(unique(site.info$Year))
  n.typs <- length(unique(site.info$survey.idx[!is.na(site.info$survey.idx)]))
  
  #summarised by site, year, and survey type
  site.info.survey <- site.info %>% 
      st_drop_geometry() %>% 
      dplyr::select(site.idx,year.idx,survey.idx,TB,n,dens.sc,agri.sc,
                    hab.sc,temp.sc,propMale,propAdult)
    #sites not sampled
    site.info.survey$TB[is.na(site.info.survey$TB)]<-0 
    site.info.survey$n[is.na(site.info.survey$n)]<-0 
    
    site.level.covs <- suppressMessages(as.data.frame(site.info.survey %>%
                            group_by(year.idx,site.idx) %>%
                            summarise(TB = sum(TB,na.rm=T),
                                      n = sum(n,na.rm=T),
                                      dens.sc = max(dens.sc),
                                      agri.sc = max(agri.sc),
                                      temp.sc = max(temp.sc),
                                      hab.sc = max(hab.sc)) %>%
                            distinct() %>% 
                            arrange(year.idx,site.idx)))
    if(surv.det==F){
      surv.site.level <- suppressMessages(as.data.frame(surv %>%
                                                          group_by(year.idx,site.idx) %>%
                                                          summarise(TB = sum(TB,na.rm=T),
                                                                    n = n(),
                                                                    SITE = SITE[1],
                                                                    DMUname = DMUname[1],
                                                                    Year = Year) %>%
                                                          distinct() %>% 
                                                          arrange(year.idx,site.idx)))
      y <- data.frame(
        site=sort(rep(1:n.sites,n.yrs)),
        yr=rep(sort(rep(1:n.yrs)),n.sites),
        typ=rep(1,n.sites*n.yrs),
        y=rep(NA,n.sites*n.yrs))
      n <- data.frame(
        site=sort(rep(1:n.sites,n.yrs)),
        yr=rep(sort(rep(1:n.yrs,)),n.sites),
        typ=rep(1,n.sites*n.yrs),
        n=rep(NA,n.sites*n.yrs))
      
      for(i in 1:nrow(surv.site.level)){
        site <- surv.site.level$site.idx[i]
        year <- surv.site.level$year.idx[i]
        y$y[y$site==site & y$yr==year] <-surv.site.level$TB[i]
        n$n[n$site==site & n$yr==year] <-surv.site.level$n[i]
      }
    }else if(surv.det==T){ #separate by survey type, sampled sites
      y <- data.frame(
        site=sort(rep(1:n.sites,n.yrs*n.typs)),
        yr=rep(sort(rep(1:n.yrs,n.typs)),n.sites),
        typ=rep(1:n.typs,n.sites*n.yrs),
        y=rep(NA,n.typs*n.sites*n.yrs))
      n <- data.frame(
        site=sort(rep(1:n.sites,n.yrs*n.typs)),
        yr=rep(sort(rep(1:n.yrs,n.typs)),n.sites),
        typ=rep(1:n.typs,n.sites*n.yrs),
        n=rep(NA,n.typs*n.sites*n.yrs))
      y.array <- n.array  <- propMale <- propAdult <-array(NA,dim=c(n.sites,n.yrs,n.typs))
      for(i in 1:nrow(samp.info)){
        site <- samp.info$site.idx[i]
        year <- samp.info$year.idx[i]
        survey <- samp.info$survey.idx[i]
        y$y[y$site==site & y$yr==year & y$typ==survey] <-samp.info$TB[i]
        n$n[n$site==site & n$yr==year & n$typ==survey] <-samp.info$n[i]
        y.array[site,year,survey] <- samp.info$TB[i]
        n.array[site,year,survey] <- samp.info$n[i]
      }
    }
    #observed y only
    y.obs <- y[!is.na(y$y),]
    ntrials.obs <- n[!is.na(n$n),]
    
    #reorder factor** not original factor numbers
    y.obs$siteorig <- y.obs$site
    y.obs$yrorig <- y.obs$yr
    ntrials.obs$siteorig <- ntrials.obs$site
    ntrials.obs$yrorig <- ntrials.obs$yr
    
    y.obs$site <- as.numeric(as.factor(y.obs$site))
    y.obs$yr <- as.numeric(as.factor(y.obs$yr))
    ntrials.obs$site <- as.numeric(as.factor(ntrials.obs$site))
    ntrials.obs$yr <- as.numeric(as.factor(ntrials.obs$yr))
    y.obs<-y.obs[!duplicated(y.obs),]
    ntrials.obs<-ntrials.obs[!duplicated(ntrials.obs),]
  
  #covariates - occupancy (all sites)
  if(dens.occ==F & agri.occ==F & temp.occ==F){
    x <- site.level.covs[,c("site.idx","year.idx")]
    x <- cbind(x$site.idx,x$year.idx,rep(1,nrow(x)))
    colnames(x) <- c("site","yr","int")
  } else if(dens.occ==T & agri.occ==F & temp.occ==F){
    x <- site.level.covs[,c("site.idx","year.idx","dens.sc")]
    x <- cbind(x$site.idx,x$year.idx,rep(1,nrow(x)),x$dens.sc)
    colnames(x) <- c("site","yr","int","dens")
  } else if(dens.occ==F & agri.occ==T & temp.occ==F){
    x <- site.level.covs[,c("site.idx","year.idx","agri.sc")]
    x <- cbind(x$site.idx,x$year.idx,rep(1,nrow(x)),x$agri.sc)
    colnames(x) <- c("site","yr","int","agri")
  }else if(dens.occ==F & agri.occ==F & temp.occ==T){
    x <- site.level.covs[,c("site.idx","year.idx","temp.sc")]
    x <- cbind(x$site.idx,x$year.idx,rep(1,nrow(x)),x$temp.sc)
    colnames(x) <- c("site","yr","int","temp")
  }else if(dens.occ==T & agri.occ==T & temp.occ==F){
    x <- site.level.covs[,c("site.idx","year.idx","dens.sc","agri.sc")]
    x <- cbind(x$site.idx,x$year.idx,rep(1,nrow(x)),x$dens.sc,x$agri.sc)
    colnames(x) <- c("site","yr","int","dens","agri")
  }else if(dens.occ==F & agri.occ==T & temp.occ==T){
    x <- site.level.covs[,c("site.idx","year.idx","agri.sc","temp.sc")]
    x <- cbind(x$site.idx,x$year.idx,rep(1,nrow(x)),x$agri.sc,x$temp.sc)
    colnames(x) <- c("site","yr","int","agri","temp")
  }else if(dens.occ==T & agri.occ==F & temp.occ==T){
    x <- site.level.covs[,c("site.idx","year.idx","dens.sc","temp.sc")]
    x <- cbind(x$site.idx,x$year.idx,rep(1,nrow(x)),x$dens.sc,x$temp.sc)
    colnames(x) <- c("site","yr","int","dens","temp")
  }else if(dens.occ==T & agri.occ==T & temp.occ==T){
    x <- site.level.covs[,c("site.idx","year.idx","dens.sc","agri.sc","temp.sc")]
    x <- cbind(x$site.idx,x$year.idx,rep(1,nrow(x)),x$dens.sc,x$agri.sc,x$temp.sc)
    colnames(x) <- c("site","yr","int","dens","agri","temp")
  }
  x <- as.data.frame(x) #site and year only (occupancy)
  
  #refactor site
  if(surv.det==T){
    #refactor site
    x.samp <- subset(x, site%in%y.obs$site)
    x$siteorig <- x$site
    x$yrorig <- x$yr
    x.samp$site <- as.numeric(as.factor(x.samp$site))
    x.samp$yr <- as.numeric(as.factor(x.samp$yr))
  } else {
    x.samp <-x
  }

  #covariates - detection (smapled sites only)
  w <- site.info.survey[,c("site.idx","year.idx","survey.idx")]
  w <- cbind(w$site.idx,w$year.idx,w$survey.idx)
  colnames(w) <- c("site","yr","survey")
  w <- as.data.frame(w)
  w <- w[!(duplicated(w)),]

  w <- w[-which(is.na(w$survey)),]
  w$siteorig <- w$site
  w$yrorig <- w$yr
  w$site <- as.numeric(as.factor(w$site))
  w$yr <- as.numeric(as.factor(w$yr))
  w$typ <- w$survey
  w <- as.data.frame(pivot_wider(w,names_from="survey",
                                 names_prefix="survey",
                                 values_from = "survey"))
  w$survey1[is.na(w$survey1)] <- 0
  w$survey2[is.na(w$survey2)] <- 0
  w$survey3[is.na(w$survey3)] <- 0
  w$survey2[w$survey2==2] <- 1
  w$survey3[w$survey3==3] <- 1
  w.refact <- w[,c("site","yr","typ",
                     "survey1","survey2","survey3")]
  w.refact[,1] <- as.integer(w.refact[,1])
  w.refact[,2] <- as.integer(w.refact[,2])
  w.refact[,3] <- as.integer(w.refact[,3])
  w.refact <- w.refact %>% arrange(site,yr)

  list(n.yrs=n.yrs,n.sites=n.sites,n.typs=n.typs,
         years=years,sites=sites,
         site.idx = site.info.survey$site.idx,
         year.idx = site.info.survey$year.idx,
         ntrials=ntrials.obs,
         ntrials.allsites=n,
         x.all=x, #occupancy covariates (site and year) for all sites in alpena
         x=x.samp, #occupancy covariates for samples sites only
         w=w.refact, #detection level covariates for sampled sites only
         w.allsites=w,#detection level covariates for all sites in alpena
         y=y.obs,#surveillance for sampled sites only
         y.allstites=y#surveillance for all sites in alpena
)
}
