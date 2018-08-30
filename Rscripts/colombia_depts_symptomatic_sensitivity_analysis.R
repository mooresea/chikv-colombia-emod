##
## Compare EMOD results to time series for Colombia departments
##
library(jsonlite)
#args = commandArgs(trailingOnly=TRUE)
#dept.nm=args[1]

setwd("D:/dtk-aedes_scripts/CHIKV")
dept.nm="Antioquia"

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

##Observations
dept.ts.f=paste("inputs/",dept.nm,"_2013-2016_timeseries.json",sep="")
ts.j=fromJSON(txt=dept.ts.f)
obs.ts=ts.j$weekly_reported_cases$`Reported Cases`

##weekly time series
year.v=c(rep(2013,52),rep(2014,52),rep(2015,52),rep(2016,length(obs.ts)-156))
week.v=c(1:52,1:52,1:52,1:(length(obs.ts)-156))
week.year=paste(year.v,ifelse(week.v>10,week.v,paste("0",week.v,sep="")),sep="-")


## Admin1 Results
dept.sim.dir=paste("Single_node_calibration_output_decay/",dept.nm,sep="")
dept.sim.files=paste(dept.sim.dir,list.files(path=dept.sim.dir,pattern="Inset",recursive = T),sep="/")

## Admin1 Results Low
dept.sim.dir.lo=paste("Single_node_calibration_output_decay_symp_low/",dept.nm,sep="")
dept.sim.files.lo=paste(dept.sim.dir.lo,list.files(path=dept.sim.dir.lo,pattern="Inset",recursive = T),sep="/")

## Admin1 Results High
dept.sim.dir.hi=paste("Single_node_calibration_output_decay_symp_hi/",dept.nm,sep="")
dept.sim.files.hi=paste(dept.sim.dir.hi,list.files(path=dept.sim.dir.hi,pattern="Inset",recursive = T),sep="/")

## Admin2 Results 
dept2.sim.dir=paste("chikv_admin2_spatial_output/",dept.nm,sep="")
dept2.sim.files=paste(dept2.sim.dir,list.files(path=dept2.sim.dir,pattern="Inset",recursive = T),sep="/")

## Admin2 Results Low
dept2.sim.dir.lo=paste("Multi_node_calibration_output_decay_symp_low/",dept.nm,sep="")
dept2.sim.files.lo=paste(dept2.sim.dir.lo,list.files(path=dept2.sim.dir.lo,pattern="Inset",recursive = T),sep="/")

## Admin2 Results High
dept2.sim.dir.hi=paste("Multi_node_calibration_output_decay_symp_hi/",dept.nm,sep="")
dept2.sim.files.hi=paste(dept2.sim.dir.hi,list.files(path=dept2.sim.dir.hi,pattern="Inset",recursive = T),sep="/")

##First 10 only
top10.1=paste("simulation_",seq(0,39),"_",sep="")
#top10.2=paste("simulation_",seq(0,39),"/",sep="")
if(any(grepl(top10.1[1],dept.sim.files.lo))){
  dept.sim.files=dept.sim.files[(sapply(top10.1,function(x) grep(x,dept.sim.files)))]
  dept.sim.files.lo=dept.sim.files.lo[(sapply(top10.1,function(x) grep(x,dept.sim.files.lo)))]
  dept.sim.files.hi=dept.sim.files.hi[(sapply(top10.1,function(x) grep(x,dept.sim.files.hi)))]
}else{
  dept.sim.files=dept.sim.files[(sapply(top10.2,function(x) grep(x,dept.sim.files)))]
  dept.sim.files.lo=dept.sim.files.lo[(sapply(top10.2,function(x) grep(x,dept.sim.files.lo)))]
  dept.sim.files.hi=dept.sim.files.hi[(sapply(top10.2,function(x) grep(x,dept.sim.files.hi)))]  
}

dept2.sim.files=dept2.sim.files[(sapply(top10.1,function(x) grep(x,dept2.sim.files)))]
dept2.sim.files.lo=dept2.sim.files.lo[(sapply(top10.1,function(x) grep(x,dept2.sim.files.lo)))]
dept2.sim.files.hi=dept2.sim.files.hi[(sapply(top10.1,function(x) grep(x,dept2.sim.files.hi)))]


## Get sim results
sim.infections=list()
for (i in 1:length(dept.sim.files)){
  json.i=fromJSON(txt=dept.sim.files[i])
  sim.i=floor(as.integer(unlist(strsplit(dept.sim.files[i],"_"))[6]))+1
  #sample.i=as.integer(substring(unlist(strsplit(dept.sim.files[i],"_"))[7],1,1))
  rep.infs=json.i$Channels$`New Reported Infections`$Data
  #Pop size
  pop.i=unique(json.i$Channels$`Statistical Population`$Data)
  ## For now multiply cases by 10
  rep.infs=rep.infs*10
  ##Convert daily infections to weekly
  if(length(rep.infs)/7!=length(ts.j$weekly_reported_cases$`Reported Cases`)){
    print(paste("Simulated time series for",dept.sim.files[i],"not the same length as observed time series"))
    next()
  }
  week.infs=daily.to.weekly(rep.infs)
  sim.infections[[i]]=list("simulation"=sim.i,"reported_infections"=week.infs)
  
}

# ## Get sim results
# sim2.infections=list()
# for (i in 1:length(dept2.sim.files)){
#   json.i=fromJSON(txt=dept2.sim.files[i])
#   i1=unlist(strsplit(dept2.sim.files[i],"_"))[6]
#   sim.i=floor(as.integer(unlist(strsplit(i1,"/"))[1]))+1
#   #sample.i=as.integer(substring(unlist(strsplit(dept.sim.files[i],"_"))[7],1,1))
#   rep.infs=json.i$Channels$`New Reported Infections`$Data
#   #Pop size
#   pop.i=unique(json.i$Channels$`Statistical Population`$Data)
#   ## For now multiply cases by 10
#   rep.infs=rep.infs*1
#   ##Convert daily infections to weekly
#   if(length(rep.infs)/7!=length(ts.j$weekly_reported_cases$`Reported Cases`)){
#     print(paste("Simulated time series for",dept2.sim.files[i],"not the same length as observed time series"))
#     next()
#   }
#   week.infs=daily.to.weekly(rep.infs)
#   sim2.infections[[i]]=list("simulation"=sim.i,"reported_infections"=week.infs)
#   
# }
## Get sim results
sim.infections.lo=list()
for (i in 1:length(dept.sim.files.lo)){
  json.i=fromJSON(txt=dept.sim.files.lo[i])
  sim.i=floor(as.integer(unlist(strsplit(dept.sim.files.lo[i],"_"))[8]))+1
  #sample.i=as.integer(substring(unlist(strsplit(dept.sim.files[i],"_"))[7],1,1))
  rep.infs=json.i$Channels$`New Reported Infections`$Data
  #Pop size
  pop.i=unique(json.i$Channels$`Statistical Population`$Data)
  ## For now multiply cases by 10
  rep.infs=rep.infs*10
  ##Convert daily infections to weekly
  if(length(rep.infs)/7!=length(ts.j$weekly_reported_cases$`Reported Cases`)){
    print(paste("Simulated time series for",dept.sim.files.lo[i],"not the same length as observed time series"))
    next()
  }
  week.infs=daily.to.weekly(rep.infs)
  sim.infections.lo[[i]]=list("simulation"=sim.i,"reported_infections"=week.infs)
  
}

## Get sim results
sim.infections.hi=list()
for (i in 1:length(dept.sim.files.hi)){
  json.i=fromJSON(txt=dept.sim.files.hi[i])
  sim.i=floor(as.integer(unlist(strsplit(dept.sim.files.hi[i],"_"))[8]))+1
  #sample.i=as.integer(substring(unlist(strsplit(dept.sim.files[i],"_"))[7],1,1))
  rep.infs=json.i$Channels$`New Reported Infections`$Data
  #Pop size
  pop.i=unique(json.i$Channels$`Statistical Population`$Data)
  ## For now multiply cases by 10
  rep.infs=rep.infs*10
  ##Convert daily infections to weekly
  if(length(rep.infs)/7!=length(ts.j$weekly_reported_cases$`Reported Cases`)){
    print(paste("Simulated time series for",dept.sim.files.hi[i],"not the same length as observed time series"))
    next()
  }
  week.infs=daily.to.weekly(rep.infs)
  sim.infections.hi[[i]]=list("simulation"=sim.i,"reported_infections"=week.infs)
  
}

## Get sim results
sim2.infections=list()
for (i in 1:length(dept2.sim.files)){
  json.i=fromJSON(txt=dept2.sim.files[i])
  sim.i=unlist(strsplit(dept2.sim.files[i],"_"))[8]
  #sim.i=floor(as.integer(unlist(strsplit(i1,"/"))[1]))+1
  #sample.i=as.integer(substring(unlist(strsplit(dept.sim.files[i],"_"))[7],1,1))
  rep.infs=json.i$Channels$`New Reported Infections`$Data
  #Pop size
  pop.i=unique(json.i$Channels$`Statistical Population`$Data)
  ## For now multiply cases by 10
  rep.infs=rep.infs*1
  ##Convert daily infections to weekly
  if(length(rep.infs)/7!=length(ts.j$weekly_reported_cases$`Reported Cases`)){
    print(paste("Simulated time series for",dept2.sim.files.lo[i],"not the same length as observed time series"))
    next()
  }
  week.infs=daily.to.weekly(rep.infs)
  sim2.infections[[i]]=list("simulation"=sim.i,"reported_infections"=week.infs)
  
}

## Get sim results
sim2.infections.lo=list()
for (i in 1:length(dept2.sim.files.lo)){
  json.i=fromJSON(txt=dept2.sim.files.lo[i])
  sim.i=unlist(strsplit(dept2.sim.files.lo[i],"_"))[8]
  #sim.i=floor(as.integer(unlist(strsplit(i1,"/"))[1]))+1
  #sample.i=as.integer(substring(unlist(strsplit(dept.sim.files[i],"_"))[7],1,1))
  rep.infs=json.i$Channels$`New Reported Infections`$Data
  #Pop size
  pop.i=unique(json.i$Channels$`Statistical Population`$Data)
  ## For now multiply cases by 10
  rep.infs=rep.infs*1
  ##Convert daily infections to weekly
  if(length(rep.infs)/7!=length(ts.j$weekly_reported_cases$`Reported Cases`)){
    print(paste("Simulated time series for",dept2.sim.files.lo[i],"not the same length as observed time series"))
    next()
  }
  week.infs=daily.to.weekly(rep.infs)
  sim2.infections.lo[[i]]=list("simulation"=sim.i,"reported_infections"=week.infs)
  
}

## Get sim results
sim2.infections.hi=list()
for (i in 1:length(dept2.sim.files.hi)){
  json.i=fromJSON(txt=dept2.sim.files.hi[i])
  sim.i=unlist(strsplit(dept2.sim.files.hi[i],"_"))[8]
  #sim.i=floor(as.integer(unlist(strsplit(i1,"/"))[1]))+1
  #sample.i=as.integer(substring(unlist(strsplit(dept.sim.files[i],"_"))[7],1,1))
  rep.infs=json.i$Channels$`New Reported Infections`$Data
  #Pop size
  pop.i=unique(json.i$Channels$`Statistical Population`$Data)
  ## For now multiply cases by 10
  rep.infs=rep.infs*1
  ##Convert daily infections to weekly
  if(length(rep.infs)/7!=length(ts.j$weekly_reported_cases$`Reported Cases`)){
    print(paste("Simulated time series for",dept2.sim.files.hi[i],"not the same length as observed time series"))
    next()
  }
  week.infs=daily.to.weekly(rep.infs)
  sim2.infections.hi[[i]]=list("simulation"=sim.i,"reported_infections"=week.infs)
  
}

##
##Plots
##
##Max infections in simulated or observed time series
#sim.inf.m=rep(NA,length(sim.infections[[1]]$reported_infections))
for(i in 1:length(sim.infections.lo)){
  i1=sim.infections[[i]]$reported_infections
  i1l=sim.infections.lo[[i]]$reported_infections
  i1h=sim.infections.hi[[i]]$reported_infections
  i2=sim2.infections[[i]]$reported_infections
  i2l=sim2.infections.lo[[i]]$reported_infections
  i2h=sim2.infections.hi[[i]]$reported_infections
  if(i==1){
    sim.inf.m=i1
    sim.inf.ml=i1l
    sim.inf.mh=i1h
    sim2.inf.m=i2
    sim2.inf.ml=i2l
    sim2.inf.mh=i2h
  }else{
    sim.inf.m=cbind(sim.inf.m,i1)
    sim.inf.ml=cbind(sim.inf.ml,i1l)
    sim.inf.mh=cbind(sim.inf.mh,i1h)
    sim2.inf.m=cbind(sim2.inf.m,i2)
    sim2.inf.ml=cbind(sim2.inf.ml,i2l)
    sim2.inf.mh=cbind(sim2.inf.mh,i2h)
  }
  
}
sim.inf.m=rowSums(sim.inf.m)/40  ##Why not 40?
sim.inf.ml=rowSums(sim.inf.ml)/40  ##Why not 40?
sim.inf.mh=rowSums(sim.inf.mh)/40  ##Why not 40?
sim2.inf.m=rowSums(sim2.inf.m)/40  ##Why not 40?
sim2.inf.ml=rowSums(sim2.inf.ml)/40  ##Why not 40?
sim2.inf.mh=rowSums(sim2.inf.mh)/40  ##Why not 40?
#sim2.inf.m=rowSums(sim2.inf.m)/40
cases.max=max(obs.ts,
    unlist(lapply(sim.infections, function(x) max(x$reported_infections))),
    unlist(lapply(sim.infections.hi, function(x) max(x$reported_infections))),
    unlist(lapply(sim.infections.lo, function(x) max(x$reported_infections))))
x.at=seq(1:length(obs.ts))[seq(53,length(obs.ts),13)]
x.lab=week.year[seq(53,length(obs.ts),13)]
cases.max2=max(obs.ts,
              unlist(lapply(sim2.infections, function(x) max(x$reported_infections))),
              unlist(lapply(sim2.infections.hi, function(x) max(x$reported_infections))),
              unlist(lapply(sim2.infections.lo, function(x) max(x$reported_infections))))
if(dept.nm=="Atlantico"){
  cases.max=20000
}


pdf(file=paste("calibration_figures/combo/upd2/",dept.nm,"_symptomatic_sensitivity_singlenode.pdf",sep=""),pointsize = 14,width=7,height=7,useDingbats = F)

plot(obs.ts,xlim=c(53,length(obs.ts)),ylim=c(0,cases.max+10),
     pch=19,col="black",xlab="Week",ylab="Reported Cases",axes=F)
axis(2)
axis(1,at=x.at,labels=x.lab)
box()
##plot simulations
for (i in 1:length(sim.infections)){
  lines(sim.infections[[i]]$reported_infections,lwd=1.5,col=rgb(0,0,200/255,alpha=.1))
  lines(sim.infections.lo[[i]]$reported_infections,lwd=1.5,col=rgb(0,200/255,0,alpha=.1))
  lines(sim.infections.hi[[i]]$reported_infections,lwd=1.5,col=rgb(250/255,0,0,alpha=.1))
}
lines(sim.inf.m,lwd=2,col=rgb(0,0,200/255,alpha=.8))
lines(sim.inf.ml,lwd=2,col=rgb(0,250/255,0,alpha=.8))
lines(sim.inf.mh,lwd=2,col=rgb(250/255,0,0,alpha=.8))
#lines(sim2.inf.m,lwd=2,col=rgb(250/255,0,0,alpha=.8))

legend("topleft",c("Observed","% symptomatic = 90%","% symptomatic = 72%","% symptomatic = 54%"),pch=c(19,NA,NA,NA),
       lty=c(NA,1,1,1),col=c("black","red","blue","green"),lwd=c(NA,2,2,2),bty='n')



dev.off()

pdf(file=paste("calibration_figures/combo/upd2/",dept.nm,"_symptomatic_sensitivity_multinode.pdf",sep=""),pointsize = 14,width=7,height=7,useDingbats = F)

plot(obs.ts,xlim=c(53,length(obs.ts)),ylim=c(0,cases.max2+10),
     pch=19,col="black",xlab="Week",ylab="Reported Cases",axes=F)
axis(2)
axis(1,at=x.at,labels=x.lab)
box()
##plot simulations
for (i in 1:length(sim2.infections)){
  lines(sim2.infections[[i]]$reported_infections,lwd=1.5,col=rgb(0,0,200/255,alpha=.1))
  lines(sim2.infections.lo[[i]]$reported_infections,lwd=1.5,col=rgb(0,200/255,0,alpha=.1))
  lines(sim2.infections.hi[[i]]$reported_infections,lwd=1.5,col=rgb(250/255,0,0,alpha=.1))
}
lines(sim2.inf.m,lwd=2,col=rgb(0,0,200/255,alpha=.8))
lines(sim2.inf.ml,lwd=2,col=rgb(0,250/255,0,alpha=.8))
lines(sim2.inf.mh,lwd=2,col=rgb(250/255,0,0,alpha=.8))
#lines(sim2.inf.m,lwd=2,col=rgb(250/255,0,0,alpha=.8))

legend("topleft",c("Observed","% symptomatic = 90%","% symptomatic = 72%","% symptomatic = 54%"),pch=c(19,NA,NA,NA),
       lty=c(NA,1,1,1),col=c("black","red","blue","green"),lwd=c(NA,2,2,2),bty='n')


dev.off()