##
## Compare EMOD results to time series for Colombia departments
##
library(jsonlite)
#args = commandArgs(trailingOnly=TRUE)
#dept.nm=args[1]

setwd("D:/dtk-aedes_scripts/CHIKV")
dept.nm="Atlantico"

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

## Admin2 Results
dept2.sim.dir=paste("Dept_admin2_nodes_calibration_output/",dept.nm,sep="")
dept2.sim.files=paste(dept2.sim.dir,list.files(path=dept2.sim.dir,pattern="Inset",recursive = T),sep="/")

##First 10 only
top10.1=paste("simulation_",seq(0,39),"_",sep="")
top10.2=paste("simulation_",seq(0,39),"/",sep="")
if(any(grepl(top10.1[1],dept.sim.files))){
  dept.sim.files=dept.sim.files[(sapply(top10.1,function(x) grep(x,dept.sim.files)))]
}else{
  dept.sim.files=dept.sim.files[(sapply(top10.2,function(x) grep(x,dept.sim.files)))]  
}

dept2.sim.files=dept2.sim.files[(sapply(top10.2,function(x) grep(x,dept2.sim.files)))]

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

## Get sim results
sim2.infections=list()
for (i in 1:length(dept2.sim.files)){
  json.i=fromJSON(txt=dept2.sim.files[i])
  i1=unlist(strsplit(dept2.sim.files[i],"_"))[6]
  sim.i=floor(as.integer(unlist(strsplit(i1,"/"))[1]))+1
  #sample.i=as.integer(substring(unlist(strsplit(dept.sim.files[i],"_"))[7],1,1))
  rep.infs=json.i$Channels$`New Reported Infections`$Data
  #Pop size
  pop.i=unique(json.i$Channels$`Statistical Population`$Data)
  ## For now multiply cases by 10
  rep.infs=rep.infs*1
  ##Convert daily infections to weekly
  if(length(rep.infs)/7!=length(ts.j$weekly_reported_cases$`Reported Cases`)){
    print(paste("Simulated time series for",dept2.sim.files[i],"not the same length as observed time series"))
    next()
  }
  week.infs=daily.to.weekly(rep.infs)
  sim2.infections[[i]]=list("simulation"=sim.i,"reported_infections"=week.infs)
  
}

##
##Plots
##
##Max infections in simulated or observed time series
#sim.inf.m=rep(NA,length(sim.infections[[1]]$reported_infections))
for(i in 1:length(sim.infections)){
  i1=sim.infections[[i]]$reported_infections
  i2=sim2.infections[[i]]$reported_infections
  if(i==1){
    sim.inf.m=i1
    sim2.inf.m=i2
  }else{
    sim.inf.m=cbind(sim.inf.m,i1)
    sim2.inf.m=cbind(sim2.inf.m,i2)
  }
  
}
sim.inf.m=rowSums(sim.inf.m)/39
sim2.inf.m=rowSums(sim2.inf.m)/39
cases.max=max(obs.ts,
    unlist(lapply(sim.infections, function(x) max(x$reported_infections))),
    unlist(lapply(sim2.infections, function(x) max(x$reported_infections))))
x.at=seq(1:length(obs.ts))[seq(53,length(obs.ts),13)]
x.lab=week.year[seq(53,length(obs.ts),13)]

if(dept.nm=="Atlantico"){
  cases.max=10000
}


pdf(file=paste("calibration_figures/combo/upd/",dept.nm,"_top10_simulations_singleimport.pdf",sep=""),pointsize = 14,width=7,height=7,useDingbats = F)

plot(obs.ts,xlim=c(53,length(obs.ts)),ylim=c(0,cases.max+10),
     pch=19,col="black",xlab="Week",ylab="Reported Cases",axes=F)
axis(2)
axis(1,at=x.at,labels=x.lab)
box()
##plot simulations
for (i in 1:length(sim.infections)){
  lines(sim.infections[[i]]$reported_infections,lwd=1.5,col=rgb(0,0,200/255,alpha=.1))
  lines(sim2.infections[[i]]$reported_infections,lwd=1.5,col=rgb(250/255,0,0,alpha=.1))
}
lines(sim.inf.m,lwd=2,col=rgb(0,0,200/255,alpha=.8))
lines(sim2.inf.m,lwd=2,col=rgb(250/255,0,0,alpha=.8))

legend("topleft",c("","Observed","Single-patch","Multi-patch"),pch=c(NA,19,NA,NA),
       lty=c(NA,1,1,1),col=c(NA,"black","blue","red"),lwd=c(NA,NA,2,2),bty='n')
legend("topleft",c(dept.nm),cex=1.1,bty='n')
#points(obs.ts,pch=19,col="red")
# lines(sim.infections[[1]]$reported_infections,lwd=2,col="blue")
# if(sim.infections[[2]]$sample==1 | is.na(sim.infections[[2]]$sample)){
#   lines(sim.infections[[2]]$reported_infections,lwd=2,col="blue")
#   lines(sim.infections[[3]]$reported_infections,lwd=2,col="blue")
#   lines(sim.infections[[4]]$reported_infections,lwd=2,col="blue")
# } 

dev.off()