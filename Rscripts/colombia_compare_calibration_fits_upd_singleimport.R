##
## Plot of RMSE or MAE
##
library(jsonlite)
## Comparison of fit to departmental time series b/w single node and admin2 node calibrations
## 1. cumulative cases (with uncertainty), peak week, length of outbreak
## 2. weekly error over course of outbreak

## Comparison of fit to municipality time series for admin2 calibrations 
## 1. cumulative cases
## 2. timing of first case and timing of peak 

fig.dir="d:/dtk-aedes_scripts/CHIKV/calibration_figures/upd/"

##
## Function to convert daily timeseries to weekly
##
daily.to.weekly=function(daily.ts){
  week.inds=seq(1,length(daily.ts),by=7)
  week.end.inds=c(week.inds[-1]-1,week.inds[length(week.inds)]+6)
  week.ts=rep(NA,length(week.inds))
  for(i in 1:length(week.inds)){
    week.ts[i]=sum(daily.ts[week.inds[i]:week.end.inds[i]])
  }
  return(round(week.ts,0))
}

##
## Single-node
##
##Observations
setwd("d:/dtk-aedes_scripts/CHIKV/inputs")
obs.f=list.files(pattern=".json")
obs.f=obs.f[-grep("example",obs.f)]
obs.list=list()
for(i in 1:length(obs.f)){
  ij=fromJSON(txt=obs.f[i])
  obs.list[[i]]=ij$weekly_reported_cases$`Reported Cases`
}
obs.len=159 #unique(unlist(lapply(obs.list, length)))
obs.names=lapply(strsplit(obs.f,"_"), function(x) x[1])
##weekly time series
# year.v=c(rep(2013,52),rep(2014,52),rep(2015,52),rep(2016,length(obs.ts)-156))
# week.v=c(1:52,1:52,1:52,1:(length(obs.ts)-156))
# week.year=paste(year.v,ifelse(week.v>10,week.v,paste("0",week.v,sep="")),sep="-")

##Simulations
setwd("d:/dtk-aedes_scripts/CHIKV/Single_node_calibration_output_decay")
sim.dir=list.dirs(recursive=F,full.names = F)
#sim.dir=sim.dir[-which(sim.dir=="Colombia")]
sim.names=sim.dir
sim.list=list()
for(i in 1:length(sim.dir)){
  idirs=list.dirs(path=sim.dir[i])
  idirs=idirs[-1]
  sim.infections=matrix(NA,nrow=length(idirs),ncol=obs.len)
  rownames(sim.infections)=unlist(lapply(strsplit(idirs,"_"),function(x)x[2]))
  for(j in 1:length(idirs)){
    jsim=fromJSON(txt=paste(idirs[j],"InsetChart.json",sep="/"))
    #sim.i=floor(as.integer(unlist(strsplit(dept.sim.files[i],"_"))[6])/2)+1
    #sample.i=as.integer(substring(unlist(strsplit(dept.sim.files[i],"_"))[7],1,1))
    rep.infs=jsim$Channels$`New Reported Infections`$Data
    #Pop size
    pop.i=unique(jsim$Channels$`Statistical Population`$Data)
    ## For now multiply cases by 10
    rep.infs=rep.infs*10
    ##Convert daily infections to weekly
    # if(length(rep.infs)/7!=length(ts.j$weekly_reported_cases$`Reported Cases`)){
    #   print(paste("Simulated time series for",dept.sim.files[i],"not the same length as observed time series"))
    #   next()
    # }
    week.infs=daily.to.weekly(rep.infs)
    sim.infections[j,]=week.infs
  }
  sim.list[[i]]=sim.infections
  print(i)
}

##Calculate MAE and RMSE
sims.rmse=list()
sims.mae=list()
sims.mase=sims.mae
sims.rmse.norm=sims.rmse
best.sims.rmse=list()
best.sims.mae=list()
best.sims.mase=sims.mae
best.sims.rmse.norm=sims.rmse
for(xx in 1:length(sim.names)){
  xx.obs=obs.list[[which(obs.names==sim.names[xx])]]
  xx.sims=sim.list[[xx]]
  best.sims=which(rownames(xx.sims) %in% as.numeric(seq(0,39)))
  #RMSE
  sims.rmse[[xx]]=apply(xx.sims,1,function(ii) sum((xx.obs-ii)^2)/length(xx.obs))
  sims.rmse.norm[[xx]]=sims.rmse[[xx]]/mean(xx.obs)
  #MAE
  sims.mae[[xx]]=apply(xx.sims,1,function(ii) sum(abs(xx.obs-ii))/length(xx.obs))
  ##MASE
  m1=apply(xx.sims,1,function(ii) sum(abs(xx.obs-ii)))
  m2=sum(abs(diff(xx.obs)))
  sims.mase[[xx]]=m1/(m2*(length(xx.obs)/(length(xx.obs)-1)))
  ##Best sims only
  best.sims.rmse[[xx]]=sims.rmse[[xx]][best.sims]
  best.sims.rmse.norm[[xx]]=sims.rmse.norm[[xx]][best.sims]
  best.sims.mae[[xx]]=sims.mae[[xx]][best.sims]
  best.sims.mase[[xx]]=sims.mase[[xx]][best.sims]
}

## Level 2 sims
##Simulations
setwd("d:/dtk-aedes_scripts/CHIKV/chikv_admin2_spatial_output")
sim.dir=list.dirs(recursive=F,full.names = F)
#sim.dir=sim.dir[-which(sim.dir=="old")]
sim2.names=sim.dir
sim2.list=list()

for(i in 1:length(sim.dir)){
  if(sim2.names[i] %in% c("Vaupes","Amazonas","Guainia")){
    setwd("d:/dtk-aedes_scripts/CHIKV/chikv_admin2_spatial_output")
    idirs=list.dirs(path=sim.dir[i])
  }else{
    setwd("d:/dtk-aedes_scripts/CHIKV/Dept_admin2_nodes_calibration_output")
    idirs=list.dirs(path=sim.dir[i])
  }

  idirs=idirs[-1]
  sim.infections=matrix(NA,nrow=length(idirs),ncol=obs.len)
  rownames(sim.infections)=unlist(lapply(strsplit(idirs,"_"),function(x)x[2]))
  for(j in 1:length(idirs)){
    jsim=fromJSON(txt=paste(idirs[j],"InsetChart.json",sep="/"))
    #sim.i=floor(as.integer(unlist(strsplit(dept.sim.files[i],"_"))[6])/2)+1
    #sample.i=as.integer(substring(unlist(strsplit(dept.sim.files[i],"_"))[7],1,1))
    rep.infs=jsim$Channels$`New Reported Infections`$Data
    #Pop size
    pop.i=unique(jsim$Channels$`Statistical Population`$Data)
    ## For now multiply cases by 10
    rep.infs=rep.infs
    ##Convert daily infections to weekly
    # if(length(rep.infs)/7!=length(ts.j$weekly_reported_cases$`Reported Cases`)){
    #   print(paste("Simulated time series for",dept.sim.files[i],"not the same length as observed time series"))
    #   next()
    # }
    week.infs=daily.to.weekly(rep.infs)
    sim.infections[j,]=week.infs
  }
  sim2.list[[i]]=sim.infections
  print(i)
}

##Calculate MAE and RMSE

#sim2.names=sim2.names[-which(sim2.names=="old")]
a2.sims.rmse=list()
a2.sims.mae=list()
a2.sims.mase=list()
a2.sims.rmse.norm=list()
a2.best.sims.rmse=list()
a2.best.sims.mae=list()
a2.best.sims.mase=list()
a2.best.sims.rmse.norm=list()
for(xx in 1:length(sim2.names)){
  if(sim2.names[xx]=="old") next()
  xx.obs=obs.list[[which(obs.names==sim2.names[xx])]]
  xx.sims=sim2.list[[xx]]
  best.sims=which(rownames(xx.sims) %in% as.character(seq(0,39)))
  #RMSE
  a2.sims.rmse[[xx]]=apply(xx.sims,1,function(ii) sum((xx.obs-ii)^2)/length(xx.obs))
  a2.sims.rmse.norm[[xx]]=a2.sims.rmse[[xx]]/mean(xx.obs)
  #MAE
  a2.sims.mae[[xx]]=apply(xx.sims,1,function(ii) sum(abs(xx.obs-ii))/length(xx.obs))
  #MASE
  m1=apply(xx.sims,1,function(ii) sum(abs(xx.obs-ii)))
  m2=sum(abs(diff(xx.obs)))
  a2.sims.mase[[xx]]=m1/(m2*(length(xx.obs)/(length(xx.obs)-1)))
  ##Best-fit only
  a2.best.sims.rmse[[xx]]=a2.sims.rmse[[xx]][best.sims]
  a2.best.sims.mae[[xx]]=a2.sims.mae[[xx]][best.sims]
  a2.best.sims.rmse.norm[[xx]]=a2.sims.rmse.norm[[xx]][best.sims]
  a2.best.sims.mase[[xx]]=a2.sims.mase[[xx]][best.sims]
}


##Generate data frames
for(i in 1:length(sim.names)){
  i.name=sim.names[i]
  a2.ind=which(sim2.names==i.name)
  i.df=data.frame(mae=sims.mae[[i]],rmse=sims.rmse[[i]],mase=sims.mase[[i]],rmse.n=sims.rmse.norm[[i]])
  i.df$dept=i.name
  i.df$lev=1
  if(length(a2.ind)>0){
    a2.df=data.frame(mae=a2.sims.mae[[a2.ind]],rmse=a2.sims.rmse[[a2.ind]],mase=a2.sims.mase[[a2.ind]],rmse.n=a2.sims.rmse.norm[[a2.ind]])
    a2.df$dept=sim2.names[a2.ind]
    a2.df$lev=2
    i.df=rbind(i.df,a2.df)
  }
  if(i==1){
    sim.df=i.df
  }else{
    sim.df=rbind(sim.df,i.df)
  }
}

#Best only
for(i in 1:length(sim.names)){
  i.name=sim.names[i]
  a2.ind=which(sim2.names==i.name)
  i.df=data.frame(mae=best.sims.mae[[i]],rmse=best.sims.rmse[[i]],
                  mase=best.sims.mase[[i]],rmse.n=best.sims.rmse.norm[[i]])
  i.df$dept=i.name
  i.df$lev=1
  if(length(a2.ind)>0){
    a2.df=data.frame(mae=a2.best.sims.mae[[a2.ind]],rmse=a2.best.sims.rmse[[a2.ind]],
                     mase=a2.best.sims.mase[[a2.ind]],rmse.n=a2.best.sims.rmse.norm[[a2.ind]])
    a2.df$dept=sim2.names[a2.ind]
    a2.df$lev=2
    i.df=rbind(i.df,a2.df)
  }
  if(i==1){
    bestsim.df=i.df
  }else{
    bestsim.df=rbind(bestsim.df,i.df)
  }
}

library(ggplot2)
library(plyr)

ggplot(sim.df,aes(x=as.factor(lev),y=mae,fill=as.factor(lev)))+geom_boxplot()
ggplot(sim.df,aes(x=dept,y=mae,fill=as.factor(lev)))+geom_boxplot()

##Comparison of level 1 to level 2 results
dept.cnt=ddply(sim.df,.(dept,lev),summarize,
      cnt=length(mae))
dept.cnt=ddply(dept.cnt,.(dept),summarize,
               cnt=min(cnt))
for(i in 1:length(sim.names)){
  id=sim.names[i]
  dcnt=dept.cnt$cnt[dept.cnt$dept==id]
  i1=sim.df[sim.df$dept==id&sim.df$lev==1,]
  i2=sim.df[sim.df$dept==id&sim.df$lev==2,]
  i1=i1[1:dcnt,]
  i2=i2[1:dcnt,]
  if(i==1){
    msim.df=rbind(i1,i2)
  }else{
    msim.df=rbind(msim.df,i1,i2)
  }
}

bestsim.df1=bestsim.df[bestsim.df$lev==1,]
bestsim.df2=bestsim.df[bestsim.df$lev==2,]
for(i in 1:length(sim.names)){
  bdf1=bestsim.df1[bestsim.df1$dept==sim.names[i],]
  bdf2=bestsim.df2[bestsim.df2$dept==sim.names[i],]
  if(nrow(bdf2)==0){
    next()
  }
  imae=bdf1$mae-bdf2$mae
  imase=bdf1$mase-bdf2$mase
  iRmase=bdf2$mase/bdf1$mase
  irmse=bdf1$rmse-bdf2$rmse
  irmsen=bdf1$rmse.n-bdf2$rmse.n
  idf=data.frame(dept=sim.names[i],mae=imae,mase=imase,rel.mase=iRmase,rmse=irmse,rmse.n=irmsen)
  if(i==1){
    bs.comp.df=idf
  }else{
    bs.comp.df=rbind(bs.comp.df,idf)
  }
}
bs.m=ddply(bs.comp.df,.(dept),summarize,mae=median(mae),mase=median(mase),rmse=median(rmse.n),rel.mase=median(rel.mase))
bs.m=bs.m[(order(bs.m$rel.mase)),]
bs.comp.df$dept.rank=NA
bs.comp.df$dept.rank2=NA
for(i in 1:length(unique(bs.m$dept))){
  bs.comp.df$dept.rank[which(bs.comp.df$dept==bs.m$dept[i])]=i
  bs.comp.df$dept.rank2[which(bs.comp.df$dept==bs.m$dept[i])]=i
}
bs.comp.df=bs.comp.df[order(bs.comp.df$dept.rank),]
dept.levels=unique(bs.comp.df$dept)
bs.comp.df$deptf=factor(bs.comp.df$dept,levels=dept.levels)

##PLOT OF DIFFERENCE
#pdf(file="d:/dtk-aedes_scripts/CHIKV/calibration_figures/combo/admin1_vs_admin2_fits_mase.pdf",width=14,height=8)
ggplot(bs.comp.df,aes(x=deptf,y=mase))+geom_boxplot(fill="lightblue")+
  theme(axis.text.x=element_text(angle=60,hjust=1))+ylim(-14,14)+
  labs(x="Department",y="Relative mean absolute scaled error (rel MASE)")+
  geom_abline(slope=0,intercept=0)
dev.off()

##PLOT OF DIFFERENCE
pdf(file="d:/dtk-aedes_scripts/CHIKV/calibration_figures/combo/upd/admin1_vs_admin2_fits_mase_rev2_singleimport.pdf",width=12,height=8,useDingbats = F)
#png(file="d:/dtk-aedes_scripts/CHIKV/calibration_figures/combo/admin1_vs_admin2_fits_mase.png",width=12,height=8,units="in",res=300)
ggplot(bs.comp.df,aes(x=deptf,y=mase))+geom_boxplot(fill="lightblue")+
  theme(axis.text.x=element_text(angle=60,hjust=1),panel.background = element_rect(fill="white",color="black"))+ylim(-12.5,30)+
  labs(x="Department",y="Multi-patch fit vs. single-patch fit (relative MASE)")+
  geom_abline(slope=0,intercept=0)

dev.off()

# ##PLOT OF DIFFERENCE
# pdf(file="d:/dtk-aedes_scripts/CHIKV/calibration_figures/combo/upd/admin1_vs_admin2_fits_mase_rev2_singleimport.pdf",width=12,height=8,useDingbats = F)
# #png(file="d:/dtk-aedes_scripts/CHIKV/calibration_figures/combo/admin1_vs_admin2_fits_mase.png",width=12,height=8,units="in",res=300)
# ggplot(bs.comp.df,aes(x=deptf,y=mase))+geom_boxplot(fill="lightblue")+
#   theme(axis.text.x=element_text(angle=60,hjust=1),panel.background = element_rect(fill="white",color="black"))+
#   labs(x="Department",y="Multi-patch fit vs. single-patch fit (relative MASE)")+
#   geom_abline(slope=0,intercept=0)
# 
# dev.off()

pdf(file="d:/dtk-aedes_scripts/CHIKV/calibration_figures/combo/upd/admin1_vs_admin2_fits_relative_mase_rev2_singleimport.pdf",width=12,height=8,useDingbats = F)
#png(file="d:/dtk-aedes_scripts/CHIKV/calibration_figures/combo/admin1_vs_admin2_fits_mase.png",width=12,height=8,units="in",res=300)
ggplot(bs.comp.df,aes(x=deptf,y=rel.mase))+geom_boxplot(fill="lightblue")+
  theme(axis.text.x=element_text(angle=60,hjust=1),panel.background = element_rect(fill="white",color="black"))+
  labs(x="Department",y="Multi-patch fit vs. single-patch fit (relative MASE)")+
  geom_abline(slope=0,intercept=0)+scale_y_log10()

dev.off()


###
## Stats comparing results
##
quantile(sim.df$mase[sim.df$lev==1],c(0.025,0.5,0.975))
quantile(sim.df$mase[sim.df$lev==2],c(0.025,0.5,0.975))

quantile(bestsim.df1$mase,c(0.025,0.5,0.975))
quantile(bestsim.df2$mase,c(0.025,0.5,0.975))
quantile(bs.comp.df$rel.mase,c(0.025,0.5,0.975))

dept.mase.c=ddply(bs.comp.df,.(dept),summarize,
      mase.m=mean(mase),
      mase.lo=quantile(mase,0.025),
      mase.md=quantile(mase,0.5),
      mase.hi=quantile(mase,0.975))

dept.rel.mase=ddply(bs.comp.df,.(dept),summarize,
                  mase.m=mean(rel.mase),
                  mase.lo=quantile(rel.mase,0.025),
                  mase.md=quantile(rel.mase,0.5),
                  mase.hi=quantile(rel.mase,0.975))

for(i in 1:length(sim2.names)){
  print(sim2.names[i])
  #print(chisq.test(bestsim.df1$mase[bestsim.df1$dept==sim2.names[i]],bestsim.df2$mase[bestsim.df2$dept==sim2.names[i]]))
  print(t.test(log(bestsim.df1$mase[bestsim.df1$dept==sim2.names[i]]),
                   log(bestsim.df2$mase[bestsim.df2$dept==sim2.names[i]]),var.equal = F))
  
}
t.est <- t.test(log(bestsim.df1$mase[bestsim.df1$dept==sim2.names[i]]),
                       log(bestsim.df2$mase[bestsim.df2$dept==sim2.names[i]]),var.equal = F)$stat
# t = 5.6598, df = 255.185, p-value = 4.066e-08


for(i in 1:length(sim2.names)){
  print(sim2.names[i])
  print(quantile(bestsim.df1$mase[bestsim.df1$dept==sim2.names[i]],c(0.025,.5,.975)))
  print(quantile(bestsim.df2$mase[bestsim.df2$dept==sim2.names[i]],c(0.025,.5,.975)))
  #print(wilcox.test(bestsim.df1$mase[bestsim.df1$dept==sim2.names[i]],bestsim.df2$mase[bestsim.df2$dept==sim2.names[i]])$p.value)
  #print(t.test(bestsim.df1$mase[bestsim.df1$dept==sim2.names[i]],bestsim.df2$mase[bestsim.df2$dept==sim2.names[i]])$p.value)
  
}

b <- function(){
  A <- sample(log(bestsim.df1$mase[bestsim.df1$dept==sim2.names[i]]), 40, replace=T)  
  B <- sample( log(bestsim.df2$mase[bestsim.df2$dept==sim2.names[i]]), 40, replace=T) 
  stud_test <- t.test(A, B, var.equal=FALSE)
  stud_test$stat
}
t.stat.vect = vector(length=1000)
t.vect <- replicate(1000, b())

1 - mean(t.est>t.vect)
###
### National-level results
### Use first 100 simulations
### 

##National-level sim
setwd("d:/dtk-aedes_scripts/CHIKV/Single_node_calibration_output")
sim.dir=list.dirs(recursive=F,full.names = F)
sim.dir="Colombia"
sim.names=sim.dir
natsim.list=list()
idirs=list.dirs(path=sim.dir)
idirs=idirs[-1]
sim.infections=matrix(NA,nrow=length(idirs),ncol=obs.len)
rownames(sim.infections)=unlist(lapply(strsplit(idirs,"_"),function(x)x[2]))
for(j in 1:length(idirs)){
    jsim=fromJSON(txt=paste(idirs[j],"InsetChart.json",sep="/"))
    #sim.i=floor(as.integer(unlist(strsplit(dept.sim.files[i],"_"))[6])/2)+1
    #sample.i=as.integer(substring(unlist(strsplit(dept.sim.files[i],"_"))[7],1,1))
    rep.infs=jsim$Channels$`New Reported Infections`$Data
    #Pop size
    pop.i=unique(jsim$Channels$`Statistical Population`$Data)
    ## For now multiply cases by 10
    rep.infs=rep.infs*10
    ##Convert daily infections to weekly
    # if(length(rep.infs)/7!=length(ts.j$weekly_reported_cases$`Reported Cases`)){
    #   print(paste("Simulated time series for",dept.sim.files[i],"not the same length as observed time series"))
    #   next()
    # }
    week.infs=daily.to.weekly(rep.infs)
    sim.infections[j,]=week.infs
  }
natsim.list[[1]]=sim.infections

##National-level simulation results for all 3 levels
top.sims=1:100
for(i in 1:length(sim.list)){
  i1.mat=sim.list[[i]]
  i2.mat=sim2.list[[i]]
  i1.mat=i1.mat[which(as.integer(rownames(i1.mat)) %in% top.sims),]
  i2.mat=i2.mat[which(as.integer(rownames(i1.mat)) %in% top.sims),]
  if(i==1){
    nat1.mat=i1.mat
    nat2.mat=i2.mat
  }else{
    nat1.mat=nat1.mat+i1.mat
    nat2.mat=nat2.mat+i2.mat
  }
}
nat0.mat=natsim.list[[1]][which(rownames(natsim.list[[1]]) %in% top.sims),]

nat0.m=unlist(apply(nat0.mat,2,mean))
nat0.md=unlist(apply(nat0.mat,2,function(x) quantile(x,0.5)))
nat0.l=unlist(apply(nat0.mat,2,function(x) quantile(x,0.025)))
nat0.h=unlist(apply(nat0.mat,2,function(x) quantile(x,0.975)))
nat1.m=unlist(apply(nat1.mat,2,mean))
nat1.md=unlist(apply(nat1.mat,2,function(x) quantile(x,0.5)))
nat1.l=unlist(apply(nat1.mat,2,function(x) quantile(x,0.025)))
nat1.h=unlist(apply(nat1.mat,2,function(x) quantile(x,0.975)))
nat2.m=unlist(apply(nat2.mat,2,mean))
nat2.md=unlist(apply(nat2.mat,2,function(x) quantile(x,0.5)))
nat2.l=unlist(apply(nat2.mat,2,function(x) quantile(x,0.025)))
nat2.h=unlist(apply(nat2.mat,2,function(x) quantile(x,0.975)))
nat1.df=data.frame(week=1:159,mean=nat1.m,med=nat1.md,lo=nat1.l,hi=nat1.h)
nat2.df=data.frame(week=1:159,mean=nat2.m,med=nat2.md,lo=nat2.l,hi=nat2.h)
nat0.df=data.frame(week=1:159,mean=nat0.m,med=nat0.md,lo=nat0.l,hi=nat0.h)
nat1.df$cases=nat1.df$med
nat2.df$cases=nat2.df$med
nat0.df$cases=nat0.df$med

##Get national time series
chik.dat=read.csv(file="d:/Colombia/Chikv_timeseries.csv",header=T,as.is=T)
nat.ts=c(rep(0,73),rowSums(chik.dat[,-1]))
nat.df=data.frame(week=1:159,cases=nat.ts)
plot(nat.ts,typ="p",pch=19,col="red",ylim=c(0,20000))
lines(nat.ts,lwd=0.5,col="red")

year.v=c(rep(2013,52),rep(2014,52),rep(2015,52),rep(2016,nrow(nat.df)-156))
week.v=c(1:52,1:52,1:52,1:(nrow(nat.df)-156))
week.year=paste(year.v,ifelse(week.v>10,week.v,paste("0",week.v,sep="")),sep="-")
x.at=seq(1:nrow(nat.df))[seq(53,nrow(nat.df),13)]
x.lab=week.year[seq(53,nrow(nat.df),13)]

#png(file="national_time_series_fit.png",width=12,height=8,units="in",res=300)
pdf(file=paste(fig.dir,"national_time_series_fit.pdf",sep=""),width=12,height=8,useDingbats = F,pointsize = 14)
ggplot(nat.df[-(1:52),],aes(x=week,y=cases))+geom_point(color="black")+geom_line(color="black")+
 geom_ribbon(aes(ymin=nat0.l[-(1:52)],ymax=nat0.h[-(1:52)]),fill="green",alpha=0.25)+geom_line(data=nat0.df[-(1:52),],aes(x=week,y=med),col="green")+
  geom_ribbon(aes(ymin=nat1.l[-(1:52)],ymax=nat1.h[-(1:52)]),fill="blue",alpha=0.25)+geom_line(data=nat1.df[-(1:52),],aes(x=week,y=med),col="blue")+
  geom_ribbon(aes(ymin=nat2.l[-(1:52)],ymax=nat2.h[-(1:52)]),fill="red",alpha=0.25)+geom_line(data=nat2.df[-(1:52),],aes(x=week,y=med),col="red")+
  
      theme(axis.text.x=element_text(angle=60,hjust=1),panel.background = element_rect(fill="white",color="black"))+
  scale_y_continuous(name="Weekly cases")+scale_x_continuous(name="Week",breaks=x.at,labels=x.lab)+coord_cartesian(ylim = c(0, 50000)) 
##Need to add legend
dev.off()

##
## Split into 2 figures
##
#png(file="national_time_series_fit.png",width=12,height=8,units="in",res=300)
pdf(file=paste(fig.dir,"national_time_series_fit_upd.pdf",sep=""),width=10,height=8,useDingbats = F,pointsize = 18)
ggplot(nat.df[-(1:52),],aes(x=week,y=cases))+geom_point(color="black")+geom_line(color="black")+
  geom_ribbon(aes(ymin=nat0.l[-(1:52)],ymax=nat0.h[-(1:52)]),fill="green",alpha=0.25)+geom_line(data=nat0.df[-(1:52),],aes(x=week,y=med),col="green")+
 # geom_ribbon(aes(ymin=nat1.l[-(1:52)],ymax=nat1.h[-(1:52)]),fill="blue",alpha=0.25)+geom_line(data=nat1.df[-(1:52),],aes(x=week,y=med),col="blue")+
#  geom_ribbon(aes(ymin=nat2.l[-(1:52)],ymax=nat2.h[-(1:52)]),fill="red",alpha=0.25)+geom_line(data=nat2.df[-(1:52),],aes(x=week,y=med),col="red")+
  
  theme(axis.text.x=element_text(angle=60,hjust=1),panel.background = element_rect(fill="white",color="black"))+
  scale_y_continuous(name="Weekly cases")+scale_x_continuous(name="Week",breaks=x.at,labels=x.lab)+coord_cartesian(ylim = c(0, 50000)) 
##Need to add legend
dev.off()

pdf(file=paste(fig.dir,"national_time_series_dept_fit_upd_singleimport.pdf",sep=""),width=10,height=8,useDingbats = F,pointsize = 14)
ggplot(nat.df[-(1:52),],aes(x=week,y=cases))+geom_point(color="black")+geom_line(color="black")+
  #geom_ribbon(aes(ymin=nat0.l[-(1:52)],ymax=nat0.h[-(1:52)]),fill="green",alpha=0.25)+geom_line(data=nat0.df[-(1:52),],aes(x=week,y=med),col="green")+
  geom_ribbon(aes(ymin=nat1.l[-(1:52)],ymax=nat1.h[-(1:52)]),fill="blue",alpha=0.25)+geom_line(data=nat1.df[-(1:52),],aes(x=week,y=med),col="blue")+
  geom_ribbon(aes(ymin=nat2.l[-(1:52)],ymax=nat2.h[-(1:52)]),fill="red",alpha=0.25)+geom_line(data=nat2.df[-(1:52),],aes(x=week,y=med),col="red")+
  
  theme(axis.text.x=element_text(angle=60,hjust=1),panel.background = element_rect(fill="white",color="black"))+
  scale_y_continuous(name="Weekly cases")+scale_x_continuous(name="Week",breaks=x.at,labels=x.lab)+coord_cartesian(ylim = c(0, 50000)) 
##Need to add legend
dev.off()

##
##national summary stats
##
##National-level model
quantile(rowSums(nat0.mat),c(0.025,0.5,0.975))
nat.sum=rowSums(nat0.mat)
nat.sum=nat.sum[which(nat.sum!=0)]
quantile(nat.sum,c(0.025,0.5,0.975))


quantile(rowSums(nat1.mat),c(0.025,0.5,0.975))
nat.sum=rowSums(nat1.mat)
nat.sum=nat.sum[which(nat.sum!=0)]
quantile(nat.sum,c(0.025,0.5,0.975))

nat1.df[which.max(nat.df$cases),2:5]/max(nat.df$cases)


quantile(rowSums(nat2.mat),c(0.025,0.5,0.975))
nat.sum=rowSums(nat2.mat)
nat.sum=nat.sum[which(nat.sum!=0)]
quantile(nat.sum,c(0.025,0.5,0.975))
quantile(nat.sum,c(0.025,0.5,0.975))/481284

##Cases during Observed and estimtaed peaks as pct of observed peak
nat2.df[which.max(nat.df$cases),2:5]/max(nat.df$cases)
nat2.df[which.max(nat2.df$cases),2:5]/max(nat.df$cases)

##summarize
sim.df.sum=ddply(sim.df,.(dept,lev),summarize,
                 mae.m=mean(mae),
                 mae.md=median(mae),
                 mae.l=quantile(mae,c(0.025)),
                 mae.h=quantile(mae,c(0.975)))

sim.df.sum=ddply(sim.df,.(dept),summarize,
                 mae1.m=mean(mae[lev==1]),
                 mae1.md=median(mae[lev==1]),
                 mae1.l=quantile(mae[lev==1],c(0.025)),
                 mae1.h=quantile(mae[lev==1],c(0.975)),
                 mae1.25=quantile(mae[lev==1],c(0.25)),
                 mae1.75=quantile(mae[lev==1],c(0.75)),
                 mae2.25=quantile(mae[lev==2],c(0.25)),
                 mae2.75=quantile(mae[lev==2],c(0.75)),
                 mae2.m=mean(mae[lev==2]),
                 mae2.md=median(mae[lev==2]),
                 mae2.l=quantile(mae[lev==2],c(0.025)),
                 mae2.h=quantile(mae[lev==2],c(0.975)))

bestsim.df.sum=ddply(bestsim.df,.(dept),summarize,
                 mae1.m=mean(mae[lev==1]),
                 mae1.md=median(mae[lev==1]),
                 mae1.l=quantile(mae[lev==1],c(0.025)),
                 mae1.h=quantile(mae[lev==1],c(0.975)),
                 mae1.25=quantile(mae[lev==1],c(0.25)),
                 mae1.75=quantile(mae[lev==1],c(0.75)),
                 mae2.25=quantile(mae[lev==2],c(0.25)),
                 mae2.75=quantile(mae[lev==2],c(0.75)),
                 mae2.m=mean(mae[lev==2]),
                 mae2.md=median(mae[lev==2]),
                 mae2.l=quantile(mae[lev==2],c(0.025)),
                 mae2.h=quantile(mae[lev==2],c(0.975)))

#Remove those without lev2 yet
sim.df.sum=sim.df.sum[!is.na(sim.df.sum$mae2.h),]
bestsim.df.sum=bestsim.df.sum[!is.na(sim.df.sum$mae2.h),]

#Quartile ranges
admin2.imp=sim.df.sum$dept[which(sim.df.sum$mae1.md>sim.df.sum$mae2.75)]
admin1.imp=sim.df.sum$dept[which(sim.df.sum$mae2.md>sim.df.sum$mae1.75)]
ggplot(sim.df[which(sim.df$dept %in% admin2.imp),],aes(x=dept,y=mae,fill=as.factor(lev)))+
  geom_boxplot()+labs(x="Department",y="Mean absolute error",fill="Simulation\n Level")
ggplot(sim.df[which(sim.df$dept %in% admin1.imp),],aes(x=dept,y=mae,fill=as.factor(lev)))+
  geom_boxplot()+labs(x="Department",y="Mean absolute error",fill="Simulation\n Level")

#Quartile ranges
admin2.impb=bestsim.df.sum$dept[which(bestsim.df.sum$mae1.md>bestsim.df.sum$mae2.75)]
admin1.impb=bestsim.df.sum$dept[which(bestsim.df.sum$mae2.md>bestsim.df.sum$mae1.75)]
pdf(file="chikv_colombia_admin2_improvements_mae.pdf",width=12,height=9)
ggplot(bestsim.df[which(bestsim.df$dept %in% admin2.impb),],aes(x=dept,y=sqrt(mae),fill=as.factor(lev)))+
  geom_boxplot()+labs(x="Department",y="Mean absolute error",fill="Simulation\n Level")
dev.off()
pdf(file="chikv_colombia_admin2_no_improvement_mae.pdf",width=12,height=9)
ggplot(bestsim.df[which(bestsim.df$dept %in% admin1.impb),],aes(x=dept,y=sqrt(mae),fill=as.factor(lev)))+
  geom_boxplot()+labs(x="Department",y="Mean absolute error",fill="Simulation\n Level")
dev.off()