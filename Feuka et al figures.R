### Title: Plots from dyn occ model results
### Author: Abbey Feuka
### Date: 09JAN23
### Notes:
##################
library(tidyverse)
library(HDInterval)
library(boot)
library(abind)

#top model
load("./Model Outputs/pj/modfit mult p dens col all psi.RData")
nmcmc <- length(parSamples[[1]]$eps.save) #number of saved mcmc iterations
n.chains <- length(parSamples) #number of chains

#load original data
load("./Data/spatData sec 2020.RData")
harv.year.max <- 2020
source("./Functions/modDataPrep_fun.R")
modDat <- modDataPrep("site","quarter",surv,site.info,samp.info,
                      dens.occ=T,agri.occ=F,temp.occ=F,
                      surv.det=T) 

nsites <- modDat$n.sites #number of sites
nyrs <- modDat$n.yr #number of years
ntyps <- length(unique(modDat$y$typ)) #number of surveillance/detection types
ntrials <- modDat$ntrials #number of tests per site each year

y <- modDat$y #number of postiive samples 
x <- modDat$x #covariates

#remove sites less than 0.5 mi2 for plotting
site.info$area.mi <- site.info$Area/(2.59e6)
site.info <- subset(site.info,as.numeric(area.mi)>0.5)
site.info.pt <- st_centroid(site.info)

y <- y[y$site%in%unique(site.info$site.idx),]
ntrials <- ntrials[ntrials$site%in%unique(site.info$site.idx),]
x <- x[x$site%in%unique(site.info$site.idx),]

y.yr <- y %>% group_by(yr) %>% 
  summarise(totPos = sum(y))
ntrials.yr <- ntrials %>% group_by(yr) %>% 
  summarise(tot = sum(n))
y.yr$prev <- y.yr$totPos/ntrials.yr$tot

#############
#trace plots
b.names <- c("Colonization - Intercept","Colonization - Density")
for(j in 1:ncol(parSamples[[1]]$beta.save)){
  png(file=paste0("./Model Outputs/Plots/Traceplots/trace beta",j,"dens all psi.png"),
      width=600, height=350)
  {plot(parSamples[[1]]$beta.save[,j],typ="l",ylab=b.names[j])
  for(i in 2:5){points(parSamples[[i]]$beta.save[,j],typ="l",col=i)}}
  dev.off()
}

p.names<- c("Detection - Public Harvest","Detection - Opportunistic",
            "Detection - Agency Removals")
for(j in 1:ncol(parSamples[[1]]$p.save)){
  png(file=paste0("./Model Outputs/Plots/Traceplots/trace p",j,"dens all psi.png"),
      width=600, height=350)
  {plot(parSamples[[1]]$p.save[,j],typ="l",ylab=p.names[j])
  for(i in 2:5){points(parSamples[[i]]$p.save[,j],typ="l",col=i)}}
  dev.off()
}

png(file=paste0("./Model Outputs/Plots/Traceplots/trace eps dens all psi.png"),
    width=600, height=350)
{plot(parSamples[[1]]$eps.save,typ="l",ylab="Extinction")
for(i in 2:5){points(parSamples[[i]]$eps.save,typ="l",col=i)}}
dev.off()

plot(parSamples[[1]]$psi.save[10,1,],typ="l",ylab="psi")
for(i in 2:5){points(parSamples[[i]]$psi.save[10,1,],typ="l",col=i)}

#############
#parameter estimates 
#extinction/peristence probability 
eps <- unlist(lapply(1:n.chains,function(i)parSamples[[i]]$eps.save))
mean(1-eps)
hdi(1-eps)

#colonization
gam <- lapply(1:n.chains,function(i)parSamples[[i]]$gam.save[unique(site.info$site.idx),,])
gam <- do.call("abind",gam)
gam.mn <- apply(gam,c(1,2),mean)
site.info$site.idx.new <- as.numeric(as.factor(site.info$site.idx))
gam.df <- site.info[,c("year.idx","site.idx.new","mn.dens.mi","dens.sc")]
gam.df$gam.mn <- rep(NA,nrow(gam.df))
for(i in 1:length(unique(site.info$site.idx))){
  for(t in 2:nyrs){
    gam.df$gam.mn[gam.df$year.idx==t&gam.df$site.idx.new==i] <- gam.mn[i,t-1]
  }
}
range(gam.df$gam.mn,na.rm=T)

#beta coefficients
beta.orig <- do.call("rbind",lapply(1:n.chains,function(i){parSamples[[i]]$beta.save}))
beta <- cbind(c(beta.orig),sort(rep(1:ncol(beta.orig),nrow(beta.orig))))
colnames(beta)<-c("value","beta")
beta <- as.data.frame(beta)

beta.df <- beta %>% group_by(beta) %>% summarise(mn = mean(value),
                                                 lci=quantile(value,prob=0.025),
                                                 uci=quantile(value,prob=0.975),
                                                 ldi=hdi(value,0.95)[1],
                                                 udi=hdi(value,0.95)[2])
gbeta<- ggplot(beta.df,aes(x=mn,y=as.factor(beta)))+
  geom_point()+
  geom_errorbar(aes(xmin=ldi,xmax=udi))+
  geom_vline(xintercept = 0,col="grey")+
  scale_y_discrete(labels=c("Intercept","Density"))+
  ylab("Colonization probability coefficient") +
  xlab("Value")+
  xlim(c(-2.3,1.3))
ggsave(gbeta,filename = paste("./Model Outputs/Plots/col coeffs dens dist", harv.year.max,"hdpi all psi.tiff"),
       device="tiff",width=unit(5,"in"),height=unit(5,"in"))
ggsave(gbeta,filename = paste("./Model Outputs/Plots/col coeffs dens dist", harv.year.max,"hdpi all psi.png"),
       device="png",width=unit(5,"in"),height=unit(5,"in"))

#density effect
test <- unique(sort(site.info$dens.sc))
test.dens <- unique(sort(site.info$mn.dens.mi))
dens.df <- t(sapply(1:nrow(beta.orig),function(i)inv.logit(beta.orig[i,1]+beta.orig[i,2]*test)))
dens.mn <- colMeans(dens.df)
dens.lci <- apply(dens.df,2,function(i)quantile(i,probs=0.025))
dens.uci <- apply(dens.df,2,function(i)quantile(i,probs=0.975))
dens.ldi <- apply(dens.df,2,function(i)hdi(i)[1])
dens.udi <- apply(dens.df,2,function(i)hdi(i)[2])
dens.df <- as.data.frame(cbind.data.frame(dens = test.dens,
                           mn = dens.mn,
                           lci=dens.lci,
                           uci=dens.uci,
                           ldi=dens.ldi,
                           udi=dens.udi))
gdens<-
  ggplot(dens.df,aes(y=mn,x=dens)) +geom_line(aes(color="Posterior Mean"))+
  geom_ribbon(aes(ymax=lci,ymin=uci,x=dens,fill="95% HPDI"),alpha=0.4)+
  ylab("Colonization probability")+
  xlab("Deer per km2")+
  scale_fill_manual(name="",label="95% HPDI",values="grey12")+
  scale_color_manual(name="",label="Posterior Mean",values="black")
ggsave(gdens,filename = paste("./Model Outputs/Plots/dens density effect", harv.year.max,"hdpi all psi.png"),
       device="png",width=unit(7,"in"),height=unit(5,"in"))

####FIGURE 3
fig3 <- gridExtra::grid.arrange(gbeta,gdens,nrow=1,ncol=2)
ggsave(fig3,filename = paste("./Model Outputs/Plots/Final/Figure_3.tiff"),
       device="tiff",width=unit(8,"in"),height=unit(3,"in"))

#detection and samples sizes to presence
n.chains <- length(parSamples)
p <- lapply(1:n.chains,function(i){parSamples[[i]]$p.save})
p <- do.call("rbind.data.frame",p)
eff <- log(0.05)/log(1-p)

colnames(p) <- colnames(eff)<-c("Public Harvest","Opportunistic Sampling","Agency Removals")

p.df <- p %>% 
  pivot_longer(cols=1:3,names_to="type",values_to="p") %>% 
  group_by(type) %>% 
  summarise(mn=mean(p),
            md=quantile(p,probs=0.5),
            lci=quantile(p,probs=0.025),
            uci=quantile(p,probs=0.975),
            ldi=hdi(p,0.95)[1],
            udi=hdi(p,0.95)[2])
eff.df <- eff %>% pivot_longer(cols=1:3,names_to="type",values_to="eff")%>% 
  group_by(type) %>% 
  summarise(mn=mean(eff),
            md=quantile(eff,probs=0.5),
            lci=quantile(eff,probs=0.025),
            uci=quantile(eff,probs=0.975),
            ldi=hdi(eff,0.95)[1],
            udi=hdi(eff,0.95)[2])

g.p <- ggplot(p.df,aes(y=md,x=reorder(type,md,decreasing=T)))+
  geom_point() +
  geom_errorbar(data=p.df,aes(ymin=ldi,ymax=udi))+
  ylab("Detection Probability") + xlab("Survey Type") +
  theme(axis.text.x=element_text(size=8))
harv.year.max<-2020
ggsave(g.p,filename = paste("./Model Outputs/Plots/detection probabilities",harv.year.max,"hpdi all psi.tiff"),
       device="tiff",width=unit(6,"in"),height=unit(5,"in"))
ggsave(g.p,filename = paste("./Model Outputs/Plots/detection probabilities",harv.year.max,"hpdi all psi.png"),
       device="png",width=unit(6,"in"),height=unit(5,"in"))

g.eff<-ggplot(eff.df,aes(y=md,x=reorder(type,md,decreasing=F)))+
  geom_point() +
  geom_errorbar(data=eff.df,aes(ymin=ldi,ymax=udi))+
  ylab("No. Samples to 95% Detection given Presence") + xlab("Survey Type")+
  theme(axis.text.x=element_text(size=8))

ggsave(g.eff,filename = paste("./Model Outputs/Plots/samples to 95 detection",harv.year.max,"hpdi all psi.tiff"),
       device="tiff",width=unit(6,"in"),height=unit(5,"in"))
ggsave(g.eff,filename = paste("./Model Outputs/Plots/samples to 95 detection",harv.year.max,"hpdi all psi.png"),
       device="png",width=unit(6,"in"),height=unit(5,"in"))

####FIGURE 2
fig2 <- gridExtra::grid.arrange(g.p,g.eff,nrow=1,ncol=2)
ggsave(fig2,filename = paste("./Model Outputs/Plots/Final/Figure_2.tiff"),
       device="tiff",width=unit(8,"in"),height=unit(4,"in"))

#occupancy
psi.list <- list()
nsites <-length(unique(site.info$site.idx))
for(i in 1:n.chains){
  psi.list[[i]] <-list()
  for(k in 1:nmcmc){ #remove burn in
    psi.list[[i]][[k]] <- parSamples[[i]]$psi.save[unique(site.info$site.idx),,k]
    psi.list[[i]][[k]] <- cbind.data.frame(psi.list[[i]][[k]],
                                           chain=rep(i,nsites),
                                           idx=rep(k,nsites),
                                           site=unique(site.info$site.idx))
  }
}
psi.list <- do.call("rbind",psi.list)
psi.df <- do.call("rbind",psi.list)
colnames(psi.df[,1:nyrs]) <- sapply(1:nyrs,function(i)paste0("year",i))

post.psi <- psi.df %>% 
  pivot_longer(cols=colnames(psi.df)[1:nyrs],names_to="year",values_to="psi") %>% 
  group_by(site,year) %>% 
  summarise(mn = mean(psi),
            md = quantile(psi,prob=0.5),
            lci = quantile(psi,prob=0.025),
            uci= quantile(psi,prob=0.975),
            ldi = hdi(psi)[1],
            udi = hdi(psi)[2]) %>% 
  rename(site.idx=site,
         year.idx=year)
post.psi$year.idx <- as.numeric(post.psi$year.idx)
post.psi$range <- post.psi$uci-post.psi$lci

post.psi.sf <- left_join(site.info,post.psi,by=c("site.idx","year.idx"))
psi.mn <- post.psi.sf %>% 
  select(SITE,Year,mn) %>% 
  rename(year=Year,psi.mn=mn,TWNRNGQUAD=SITE)

post.psi.sf %>% 
  group_by(DMUname,Year) %>% 
  summarise(mn = mean(mn)) %>% 
  group_by(DMUname) %>% 
  summarise(grndmn = mean(mn),
            min = min(mn),
            max = max(mn))

psi.area <- post.psi.sf %>% 
  group_by(DMUname,Year) %>% 
  summarise(grndmn = mean(mn),
            ldi=hdi(mn)[1],
            udi=hdi(mn)[2])
psi.area$DMUname[psi.area$DMUname=="Alpena"] <- "Non-Core"

psi.area %>% group_by(DMUname) %>% 
  summarise(mn=mean(grndmn),
            min=min(grndmn),
            max=max(grndmn))

post.psi.sf %>% group_by(DMUname) %>% 
  summarise(mean = mean(mn),
            min = min(mn),
            max = max(mn))

psi.all <- post.psi.sf %>% 
  group_by(Year) %>% 
  summarise(grndmn = mean(mn),
            ldi=hdi(mn)[1],
            udi=hdi(mn)[2])
samples <- ntrials %>% group_by(yr,typ) %>% 
  summarise(tot = sum(n))
samples$Year <- modDat$years[samples$yr]

#average occupancy over time
axis.scale <- 10000
gocc <-
  ggplot()+
    geom_line(data=psi.all,aes(y=grndmn,x=factor(Year),size="x",col="x"),group=1)+
    geom_line(data=psi.area,aes(y=grndmn,x=factor(Year),
                                group=DMUname,col=DMUname),alpha=0.5)+
  geom_bar(data=samples,aes(x=factor(Year),y=tot/axis.scale,fill=factor(typ)),stat="identity") +
  xlab("Year")+
  scale_size_manual(values=c(1))+
  guides(size="none")+
  labs(color="Area")+
    scale_y_continuous(name="Mean bTB occupancy probability",
                       sec.axis = sec_axis(trans=~.*axis.scale,name="Samples"))+
    scale_fill_manual(name="Survey Type",values=c("skyblue","lightblue","navy"),
                      labels=c("Public Hunt","Opportunistic","Agency Removals"))+
  scale_color_manual(labels=c("bTB Core Area","Non-Core","All"),
                     values=c("magenta4","slateblue1","purple4"))

ggsave(paste0("dens_alpena_occ_yr_samples all psi.png"),
       plot=gocc,device="png",
       width=unit(7,"in"), height=unit(5,"in"),
       dpi="print",path="./Model Outputs/Plots")

####FIGURE 5
ggsave(gocc,filename = paste("./Model Outputs/Plots/Final/Figure_5.tiff"),
       device="tiff",width=unit(7,"in"),height=unit(5,"in"))

####FIGURE 4
#occupancy maps 
alpena.psi <- psi.df %>% 
  pivot_longer(cols=1:10,names_to="Year",values_to="psi") %>% 
  filter(site%in%unique(site.info$site.idx)) %>% 
  group_by(Year,idx,chain) %>% 
  summarise(alp.mn = mean(psi)) %>%
  group_by(Year) %>% 
  summarise(mn = mean(alp.mn),
            md = quantile(alp.mn,prob=0.5),
            lci = quantile(alp.mn,prob=0.025),
            uci = quantile(alp.mn,prob=0.975),
            ldi = hdi(alp.mn)[1],
            udi = hdi(alp.mn)[2]) %>% 
  arrange(as.numeric(Year))
colnames(alpena.psi)[colnames(alpena.psi)=="site"] <- "site.idx"

dmus <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/bTB/GIS Data/DMUs/dmus.shp")
counties <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/bTB/GIS Data/DMUs/county5.shp")
counties <- st_transform(counties, st_crs(dmus))
btb.area <- st_intersection(subset(dmus,Name=="bTB Core Area"),subset(counties,CNTYNAME=="Alpena"))

g <- list()
for(yr in 1:10){
  post.psi.sf.yr <- post.psi.sf %>% filter(year.idx==yr)
  site.info.pt.yr <- site.info.pt %>% filter(year.idx==yr&TB>0)
  g[[yr]] <-
    ggplot(post.psi.sf.yr)+geom_sf(aes(fill=mn)) +
      geom_sf(data=btb.area,aes(col="Name"),alpha=0)+
      scale_color_manual(name=" ",labels="bTB Core Area",values="black") +
    theme(panel.background = element_blank(),
          axis.text=element_blank(),
          axis.ticks = element_blank(),
          legend.position="none")+
    scale_fill_continuous(name="Probabilty of TB",type="viridis",limits=c(0,1)) +
    ggtitle(unique(post.psi.sf.yr$Year))+
    geom_sf(data=site.info.pt.yr,aes(size=TB))+
    scale_size(breaks=c(0,1,2,3,4),name="No. TB+ Tests",range=c(1,3))
}

for(yr in 1:nyrs){
  ggsave(paste0("dens_dist_col_tb_occ_postmn_",modDat$years[yr],"blank_core all psi.tiff"),
         plot=g[[yr]],device="tiff",
         width=unit(4,"in"), height=unit(3,"in"),
         dpi="print",path="./Model Outputs/Plots")
}

#raw data
g.samp <- list()
for(yr in 1:nyrs){
  g.samp[[yr]] <-
    site.info %>% filter(year.idx==yr) %>% group_by(site.idx) %>% 
    summarise(tot.TB = sum(TB)) %>% 
    ggplot()+geom_sf(aes(fill=factor(tot.TB))) +
    geom_sf(data=btb.area,aes(col="Name"),alpha=0)+
    theme(panel.background = element_blank(),
          axis.text=element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none")+
          # legend.direction="horizontal")+
    ggtitle(modDat$years[yr])+
    scale_color_manual(values="black")+
    scale_fill_viridis_d(name="No. TB+ Samples",limits=as.character(c(0,1,2,3,4,NA)),
                         labels=as.character(c(0,1,2,3,4,"Not sampled")),na.value = "grey88")
  ggsave(paste0("pos_samps_",modDat$years[yr],".tiff"),
         plot=g.samp[[yr]],device="tiff",
         width=unit(4,"in"), height=unit(3,"in"),
         dpi="print",path="./Model Outputs/Plots")
}

