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
v.list=list()
b.list=list()
for(i in 1:length(sim.dir)){
  idirs=list.dirs(path=sim.dir[i])
  idirs=idirs[-1]
  jsim=fromJSON(txt=paste(idirs[1],"InsetChart.json",sep="/"))
  sim.length=length(jsim$Channels$`Daily Bites per Human`$Data)
  sim.v=matrix(NA,nrow=length(idirs),ncol=sim.length)
  sim.b=sim.v
  rownames(sim.v)=unlist(lapply(strsplit(idirs,"_"),function(x)x[2]))
  for(j in 1:length(idirs)){
    jsim=fromJSON(txt=paste(idirs[j],"InsetChart.json",sep="/"))
    #sim.i=floor(as.integer(unlist(strsplit(dept.sim.files[i],"_"))[6])/2)+1
    #sample.i=as.integer(substring(unlist(strsplit(dept.sim.files[i],"_"))[7],1,1))
    bites=jsim$Channels$`Daily Bites per Human`$Data
    vectors=jsim$Channels$`Adult Vectors`$Data
    #Pop size
    #pop.i=unique(jsim$Channels$`Statistical Population`$Data)
    
    sim.v[j,]=vectors
    sim.b[j,]=bites
  }
  v.list[[i]]=sim.v
  b.list[[i]]=sim.b
  print(i)
}


## Level 2 sims
##Simulations
setwd("d:/dtk-aedes_scripts/CHIKV/chikv_admin2_spatial_output")
sim.dir=list.dirs(recursive=F,full.names = F)
#sim.dir=sim.dir[-which(sim.dir=="old")]
sim2.names=sim.dir

v2.list=list()
b2.list=list()
for(i in 1:16) { #length(sim.dir)){
  if(sim2.names[i] %in% c("Vaupes","Amazonas","Guainia")){
    setwd("d:/dtk-aedes_scripts/CHIKV/chikv_admin2_spatial_output")
    idirs=list.dirs(path=sim.dir[i])
  }else{
    setwd("d:/dtk-aedes_scripts/CHIKV/Dept_admin2_nodes_calibration_output")
    idirs=list.dirs(path=sim.dir[i])
  }
  
  idirs=idirs[-1]
  jsim=fromJSON(txt=paste(idirs[1],"InsetChart.json",sep="/"))
  sim.length=length(jsim$Channels$`Daily Bites per Human`$Data)
  sim.v=matrix(NA,nrow=length(idirs),ncol=sim.length)
  sim.b=sim.v
  rownames(sim.v)=unlist(lapply(strsplit(idirs,"_"),function(x)x[2]))
  rownames(sim.b)=unlist(lapply(strsplit(idirs,"_"),function(x)x[2]))
  
  
  for(j in 1:length(idirs)){
    
    jsim=fromJSON(txt=paste(idirs[j],"InsetChart.json",sep="/"))
    #sim.i=floor(as.integer(unlist(strsplit(dept.sim.files[i],"_"))[6])/2)+1
    #sample.i=as.integer(substring(unlist(strsplit(dept.sim.files[i],"_"))[7],1,1))
    jsim=fromJSON(txt=paste(idirs[j],"InsetChart.json",sep="/"))
    #sim.i=floor(as.integer(unlist(strsplit(dept.sim.files[i],"_"))[6])/2)+1
    #sample.i=as.integer(substring(unlist(strsplit(dept.sim.files[i],"_"))[7],1,1))
    bites=jsim$Channels$`Daily Bites per Human`$Data
    vectors=jsim$Channels$`Adult Vectors`$Data
    #Pop size
    #pop.i=unique(jsim$Channels$`Statistical Population`$Data)
    
    sim.v[j,]=vectors
    sim.b[j,]=bites
  }
  v2.list[[i]]=sim.v
  b2.list[[i]]=sim.b
  print(i)
 
}

##Top 10
v1.mat=matrix(NA,nrow=length(sim2.names),ncol=ncol(v.list[[1]]))
b1.mat=v2.mat=b2.mat=v1.mat
for (i in 1:length(sim2.names)){
  l1.ind=which(sim.names==sim2.names[i])
  lev1.v=v.list[[l1.ind]]
  lev1.b=b.list[[l1.ind]]
  lev2.v=v2.list[[i]]
  lev2.b=b2.list[[i]]
  
  lev1.sims=as.numeric(rownames(lev1.v))
  lev2.sims=as.numeric(rownames(lev2.v))
  l1=order(lev1.sims)[1:40]
  l2=order(lev2.sims)[1:40]
  v1.mat[i,]=colMeans(lev1.v[l1,])
  b1.mat[i,]=colMeans(lev1.b[l1,])
  v2.mat[i,]=colMeans(lev2.v[l2,])
  b2.mat[i,]=colMeans(lev2.b[l2,])
  
}


##Plots
setwd("C:/users/perkinslab/Dropbox/CHIKV/other_figures/")
library(RColorBrewer)
reds=brewer.pal(9,"Reds")
blues=brewer.pal(9,"Blues")

x.at=seq(1,1100,by=365)
x.labs=2014:2017

dname="NortedeSantander"
dplot=which(sim2.names==dname)

pdf(file=paste(dname,"_bites.pdf",sep=""),width=8,height=6,pointsize=12,useDingbats = F)
plot(b1.mat[dplot,-(1:365)],typ="l",cex.lab=1.2,cex.axis=1.1,col=blues[7],
     xlab="",ylab="Daily bites per human",axes=F,ylim=c(0,30))
axis(1,at=x.at,labels=x.labs)
axis(2)
box()
lines(b2.mat[dplot,-(1:365)],col=reds[7])
legend("topleft",c("Single-patch model","Multi-patch model"),lwd=2,col=c(blues[7],reds[7]),bty="n")

dev.off()