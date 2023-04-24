### Title: Prep spatial data for dyn occ mod
### Author: Abbey Feuka
### Date: 18NOV22
### Notes: 
###################
spatDatPrep <- function(siteLevel, # scale of occupancy site (section,quarter,or township)
                        end.yr) # cutoff year for data (2020 or 2021)
  { 
  require(nimble)
  require(tidyverse)
  require(sf)
  require(terra)
  
  ###surveillance all 5 counties
  surv.allsites <- read.csv("./Data/bTB_surveillance_dmu_5co_2011_2021.csv")
  colnames(surv.allsites)[colnames(surv.allsites)=="TB."]<-"TB"
  #spatial data
  sec <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/bTB/GIS Data/DMUs/section5.shp",quiet=T)
  dmu <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/bTB/GIS Data/DMUs/dmus.shp",quiet=T)
  
  sec <- st_transform(sec,"epsg:3078")
  dmu <- st_transform(dmu,"epsg:3078")
  alp.dmu.sec <- st_intersection(sec,subset(subset(dmu,Name=="Alpena")))
  alp.count.sec <- st_intersection(subset(sec,CNTYNAME=="Alpena"),
                                   subset(subset(dmu,Name=="bTB Core Area")))
  identical(colnames(alp.count.sec),colnames(alp.dmu.sec))
  alp.sec.all <- rbind(alp.dmu.sec,alp.count.sec)
  alp.sec.all <- subset(alp.sec.all, CNTYNAME=="Alpena")
  alp.sec.all$Area <- st_area(alp.sec.all)
  
  #object id
  dups <- alp.sec.all[duplicated(alp.sec.all$OBJECTID),]
  for(i in 1:nrow(dups)){
    x <- subset(alp.sec.all,OBJECTID==unique(dups$OBJECTID)[i])
    x.min <- min(x$Area)
    x.max <- max(x$Area)
    #remove smaller section row
    alp.sec.all <- alp.sec.all[-(which(alp.sec.all$Area==x.min &
                                         alp.sec.all$OBJECTID==x$OBJECTID[1])),]
    #combine section geometries
    st_geometry(alp.sec.all[which(alp.sec.all$Area==x.max &
                                    alp.sec.all$OBJECTID==x$OBJECTID[1]),]) <- st_union(x)
  }
  
  #sections
  dups <- alp.sec.all[duplicated(alp.sec.all$TWNRNGSEC),]
  for(i in 1:nrow(dups)){
    x <- subset(alp.sec.all,TWNRNGSEC==unique(dups$TWNRNGSEC))
    x.min <- min(x$Area)
    x.max <- max(x$Area)
    #remove smaller section row
    alp.sec.all <- alp.sec.all[-(which(alp.sec.all$Area==x.min)),]
    #combine section geometries
    st_geometry(alp.sec.all[which(alp.sec.all$Area==x.max),]) <- st_union(x)
  }
  
  #4-section quads
  Q1<-c(1,2,3,10,11,12,13,14,15)
  Q2<-c(4,5,6,7,8,9,16,17,18)
  Q3<-c(19,20,21,28,29,30,31,32,33)
  Q4<-c(22,23,24,25,26,27,34,35,36)
  
  alp.sec.all$QUAD <- numeric(nrow(alp.sec.all))
  for(i in 1:nrow(alp.sec.all)){
    if(alp.sec.all$SEC[i] %in% Q1){
      alp.sec.all$QUAD[i] <- "Q1"
    } else if(alp.sec.all$SEC[i] %in% Q2){
      alp.sec.all$QUAD[i] <- "Q2"
    } else if(alp.sec.all$SEC[i] %in% Q3){
      alp.sec.all$QUAD[i] <- "Q3"
    } else if(alp.sec.all$SEC[i] %in% Q4){
      alp.sec.all$QUAD[i] <- "Q4"
    }
  }
  alp.sec.all$TWNRNGQUAD <- paste0(alp.sec.all$TWNRNG, alp.sec.all$QUAD)
  
  if(siteLevel=="section"){
    alp.sec.all$SITE <- alp.sec.all$TWNRNGSEC
  } else if(siteLevel=="quarter"){
    alp.sec.all$SITE <- alp.sec.all$TWNRNGQUAD
  }
  
  ###site info
  #includes ALL SITES IN COUNTY, sampled or not
  site.info <- alp.sec.all[,c("SITE","CNTYNAME","Name")]
  site.info <- site.info %>% 
    rename(DMUname="Name") %>% 
    arrange(SITE) # IN ORDER BY SITE
  site.info$Area <- st_area(site.info)
  site.info$ObjectID <- 1:nrow(site.info)
  
  if(siteLevel=="quarter"){
    dups <- unique(site.info[duplicated(site.info$SITE),]$SITE)
    for(i in 1:length(dups)){
      x <- subset(site.info,SITE==dups[i])
      x.max <- x$ObjectID[which.max(x$Area)]
      x.min <- x$ObjectID[x$ObjectID!=x.max]
      #remove smaller geometries
      site.info <- site.info[!(site.info$ObjectID%in%x.min),]
      #combine section geometries
      st_geometry(site.info[which(site.info$ObjectID==x.max),]) <- st_union(x)
    }
  }
  #remove small sites
  site.info$site.idx.orig <- as.numeric(as.factor(site.info$SITE))
  site.info$Area<-st_area(site.info)
  site.info$area.mi <- st_area(site.info)/2.59e6
  smol <- unique(site.info$SITE[which(as.numeric(site.info$area.mi)<0.5)])
  site.info <- site.info[-which(site.info$SITE%in%smol),]
  
  # density estimates 
  dens <- read.csv("C:/Users/Abigail.Feuka/OneDrive - USDA/bTB/Spatial_Risk_Assessment/deerabundance/nimble/Model Outputs/density est 29SEP22.csv")
  #align density with surveillance years
  dens <- dens[,c(which(colnames(dens)=="Name"),which(colnames(dens)=="year"),grep("dens",colnames(dens)))]
  dens <- subset(dens, Name%in%c("Alpena","bTB Core Area") & year!=2010 & year<=end.yr) #2011-2020
  colnames(dens)[colnames(dens)=="Name"] <- "DMUname"
  
  #join density estimates to section data by DMU name
  site.info <- left_join(site.info,dens,by=c("DMUname"))
  colnames(site.info)[colnames(site.info)=="year"] <- "Year"
  
  # agriculture
  if(siteLevel=="section"){
    covs <- read.csv("./Data/dmu_landscape_covs_section.csv")
    agri <- covs %>% filter(source=="nlcd" & 
                              Ecosystem%in%c("Cultivated Crops")&
                              year==2019) %>% #very little variation among years 
      dplyr::select(TWNRNGSEC,prop) %>% 
      rename(propAg = prop,
             SITE=TWNRNGSEC)
  }else if(siteLevel=="quarter"){
    covs <- read.csv("./Data/dmu_landscape_covs_quarter.csv")
    class(covs)
    agri <- covs %>% 
      filter(source=="nlcd" &
               Ecosystem%in%c("Cultivated Crops")&
               year==2019) %>% #very little variation among years
      dplyr::select(TWNRNGQUAD,prop) %>% 
      rename(propAg = prop,
             SITE=TWNRNGQUAD)
  }
  site.info <- left_join(site.info, agri, by=c("SITE"))
  site.info$propAg[is.na(site.info$propAg)] <- 0 #no ag = 0 ag
  
  #deer habitat 
  hab <- rast("C:/Users/Abigail.Feuka/OneDrive - USDA/bTB/GIS Data/FWFP/fwfp_dmu_epsg3078.tif")
  site.info$habPot <- numeric(nrow(site.info))
  for(i in 1:nrow(site.info)){
    x <- values(crop(hab,site.info[i,]))
    x[is.nan(x)] <- 0
    site.info$habPot[i] <- mean(x)
  }
  
  #days below freezing
  noaa <- read.csv("C:/Users/Abigail.Feuka/OneDrive - USDA/bTB/Spatial_Risk_Assessment/deerabundance/nimble/noaa_weather_all_dmus_pooled.csv")
  dt32 <- noaa %>% 
    filter(type=="DT32" & DATE>2009 & DATE<=end.yr) %>% 
    rename(Year=DATE,dt32=mn) %>% 
    dplyr::select(Year,dt32)
  site.info <- left_join(site.info,dt32,by="Year")
  
  #indices 
  site.info$site.idx <- as.numeric(as.factor(site.info$SITE))
  
  #adjacency matrix (remove multiples per year from site.info df)
  site.idx.df <- site.info %>% dplyr::select(site.idx,geometry) %>% arrange(site.idx)
  site.idx.df <- site.idx.df[!duplicated(site.idx.df$site.idx),]
  site.neighbors <- st_intersects(site.idx.df,remove_self=T)
  adj <- as.carAdjacency(site.neighbors) # remove extra section creaeted by site.neighbors
  identical(as.numeric(max(adj$adj)),as.numeric(max(site.info$site.idx)))
  
  #distance matrix
  site.cent <- suppressWarnings(st_centroid(site.idx.df))
  site.dist <- st_distance(site.cent)/1000 #in km
  #remove units
  attr(site.dist,"units") <- NULL
  class(site.dist) <- setdiff(class(site.dist),"units")
  colnames(site.dist)<-NULL
  
  #scale covariates
  site.info$dens.sc <- scale(site.info$mn.dens.mi) #mean density
  site.info$agri.sc <- scale(site.info$propAg) #proportion agriculture
  site.info$hab.sc <- scale(site.info$habPot) #fall winter food potential
  site.info$temp.sc <- scale(site.info$dt32) #days below 0C
  
  years <- sort(unique(site.info$Year))
  site.info$year.idx <- as.numeric(as.factor(site.info$Year))
  sites <- sort(unique(site.info$SITE))
  n.sites <- length(unique(site.info$site.idx))
  n.yr <- length(unique(site.info$Year))

  ### sample info 
  if(siteLevel=="section"){
    temp <- site.info %>% 
      rename(TWNRNGSEC=SITE) %>% 
      dplyr::select(TWNRNGSEC,site.idx,DMUname,Year) %>% 
      st_drop_geometry()
    surv <- left_join(surv.allsites,temp,by=c("TWNRNGSEC","Year")) %>% 
      dplyr::select(TWNRNGSEC,CNTYNAME,DMUname,Year,Survey_Type,Sex,Age,TB,site.idx) %>% 
      rename(SITE=TWNRNGSEC) %>% 
      filter(!is.na(site.idx))
  } else {
    temp <- site.info %>% 
      rename(TWNRNGQUAD=SITE) %>% 
      dplyr::select(TWNRNGQUAD,site.idx,DMUname,Year) %>% 
      st_drop_geometry()
    surv.temp <- left_join(surv.allsites,alp.sec.all,by="TWNRNGSEC") %>% 
      filter(!is.na(TWNRNGQUAD))
    surv <- left_join(surv.temp,temp,by=c("TWNRNGQUAD","Year")) %>% 
      dplyr::select(TWNRNGQUAD,DMUname,Year,Survey_Type,Sex,Age,TB,site.idx) %>% 
      rename(SITE=TWNRNGQUAD) %>% 
      filter(!is.na(site.idx))
  }
  
  surv$year.idx <- as.numeric(as.factor(surv$Year))
  surv <- surv %>% arrange(year.idx,site.idx)

  surv$TB[surv$TB=="Neg"] <- 0
  surv$TB[surv$TB=="Pos"] <- 1
  surv <- subset(surv, TB!="IS") #removes sections with no samples
  surv$TB <- as.numeric(surv$TB)
  
  surv <- surv %>% arrange(year.idx,site.idx)
  surv$samp.idx <- 1:nrow(surv) #one for each sample

  #surv is one line per data point
  #samp.info is grouped by site and year
  samp.info <- suppressMessages(as.data.frame(surv %>%
                                    group_by(year.idx,site.idx,survey.idx) %>%
                                    summarise(TB = sum(TB),
                                              n = n(),
                                              propMale = sum(Sex,na.rm=T)/n,
                                              propAdult = sum(Age,na.rm=T)/n,
                                              SITE = SITE[1],
                                              DMUname = DMUname[1],
                                              Year = Year) %>%
                                    distinct() %>% 
                                    arrange(year.idx,site.idx)))
  
  n.total <- nrow(surv) #including unsampled sites
  n.samp <- nrow(subset(surv,!is.na(TB))) #only samples
  identical(max(surv$samp.idx),n.total)
  
  #sample-level information
  samp.info.survey <- surv %>% 
    dplyr::select(samp.idx,site.idx,year.idx,survey.idx,Sex,Age,TB) %>% 
    arrange(samp.idx)

  site.info <- left_join(site.info,samp.info,
                           by=c("year.idx","site.idx","SITE","DMUname","Year"))
  all.info <- left_join(samp.info[,c("year.idx","site.idx","survey.idx",
                                     "TB","n","Year")],
                        site.info[,c("year.idx","site.idx","survey.idx",
                                     "SITE","CNTYNAME","DMUname","area.mi",
                                     "site.idx.orig","dens.sc","agri.sc","temp.sc")],
                        by=c("year.idx","site.idx","survey.idx"))

  list(alp.sec.all=alp.sec.all, #all sections in alpena county
       site.info=site.info, #data on all sites in Alpena, w/ or w/o TB surveillance
       surv=surv, #surveillance data by individual
       adj=adj$adj,num=adj$num,site.dist=site.dist, #adjacency matrix
       samp.info=samp.info, #data on surveyed sites only, w/ covariates
       all.info=all.info) #everything in one dataframe
}

