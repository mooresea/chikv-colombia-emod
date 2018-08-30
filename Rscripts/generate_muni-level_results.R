##
## Comparison of simulated and observed municipality data
##
setwd("d:/")

##
## Function to convert daily timeseries to weekly
##
daily.to.weekly=function(daily.mat){
  week.inds=seq(1,ncol(daily.mat),by=7)
  week.end.inds=c(week.inds[-1]-1,week.inds[length(week.inds)]+6)
  week.mat=matrix(NA,ncol=length(week.inds),nrow=nrow(daily.mat))
  for(i in 1:nrow(week.mat)){
    for(j in 1:length(week.inds)){
      week.mat[i,j]=sum(daily.mat[i,week.inds[j]:week.end.inds[j]])
    }    
  }
  return(round(week.mat,0))
}

#
# Observed data
#
chik.dat=read.csv(file="Colombia/laptop/Chikv_timeseries.csv",header=T,as.is=T)
mcodes=read.table(file="Colombia/Unique_divipola_codes.txt",sep=":",header=T,as.is = T)

##
## Siumulated data
dept.names=list.dirs(path="d:/dtk-aedes_scripts/CHIKV/chikv_admin2_spatial_output/",recursive=F,full.names = F)
#dept.names=dept.names[-which(dept.names=="old")]
#dept.names=c("Vichada")

for(dd in 1:length(dept.names)){
#Dept
  dept=dept.names[dd]
  print(paste("Start",dept,":",dd,"out of",length(dept.names),"depts"))
  dept.dir=paste("d:/dtk-aedes_scripts/CHIKV/chikv_admin2_spatial_output/",dept,sep="")
  sim.dirs=list.dirs(path=dept.dir,recursive=F)
  sim.mat=read.csv(file=paste(sim.dirs[1],"spatial_infections_upd.csv",sep="/"),header=T,as.is=T)
  sim.ids=as.numeric(unlist(lapply(strsplit(sim.dirs,"_"),function(x) x[6])))
  #Get spatial results for each sim
  ##Array
  sim.a=array(NA,dim=c(length(sim.dirs),nrow(sim.mat),160))
  for(i in 1:length(sim.dirs)){
    sim.mat=read.csv(file=paste(sim.dirs[i],"spatial_infections_upd.csv",sep="/"),header=T,as.is=T)
    sim.a[i,,]=cbind(sim.mat[,1],daily.to.weekly(sim.mat[,-1]))
    if(i%%10==0) print(i)
    
  }
  ##Mean of best-fit sim (4 samples)
  bestsim.m=t(apply(sim.a[which(sim.ids<4),,-1],2,colMeans))
  ##Mean from array
  sim.m=t(apply(sim.a[,,-1],2,colMeans))
  
  sim.nodes=substring(sim.a[1,,1],first=3)
  #Dept id
  if(max(nchar(sim.nodes))==5){
    dept.id=substring(sim.nodes[1],first=1,last=2)
  }else{
    dept.id=substring(sim.nodes[1],first=1,last=1)
  }
  save(sim.a,bestsim.m,sim.m,sim.nodes,file=paste(dept,dept.id,"spatial_sim_results_upd.rda",sep="_"))
}