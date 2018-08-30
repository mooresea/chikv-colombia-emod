##
##  Plot observed vs. simulated for one department
##
library(ggplot2)
library(jsonlite)
setwd("d:/")

#
# Observed data
#
chik.dat=read.csv(file="Colombia/laptop/Chikv_timeseries.csv",header=T,as.is=T)
mcodes=read.table(file="Colombia/Unique_divipola_codes.txt",sep=":",header=T,as.is = T)

##
## Siumulated data
dept.names=list.dirs(path="d:/dtk-aedes_scripts/CHIKV/chikv_admin2_spatial_output/",recursive=F,full.names = F)
#dept.names=dept.names[-which(dept.names=="old")]
dept.files=list.files(path="d:/Colombia/muni_output/upd/",full.names = T)
###
##Calc metrics and plot fits for each department
dept="Bolivar"
  #LOAD Data
  deptf=dept.files[grep(paste("/",dept,"_",sep=""),dept.files)]
  #deptf=paste("d:/Colombia/muni_output/upd/",deptf,sep="")
  load(deptf)
  #deptf2=paste("d:/Colombia/muni_output/",dept,"__spatial_sim_results_incidence.rda",sep="")
  #load(deptf2)
  
  ##Get pop data from demographics file
  deptdf=paste("d:/Eradication/InputDataFiles/Colombia_dept_nodes_admin2/",dept,"/",dept,"_admin2_nodes_demographics.json",sep="")
  dept.j=fromJSON(txt=deptdf)
  
  sim.nodes=substring(sim.a[1,,1],first=3)
  j.nodes=substring(dept.j$Nodes$NodeID,first=3)
  j.pop=dept.j$Nodes$NodeAttributes$InitialPopulation
  j.pop=j.pop[unlist(lapply(sim.nodes,function(x) which(j.nodes==x)))]
  #Dept id
  if(max(nchar(sim.nodes))==5){
    dept.id=substring(sim.nodes[1],first=1,last=2)
    #j.nodes=substring(j.nodes,first=3)
  }else{
    dept.id=substring(sim.nodes[1],first=1,last=1)
    #j.nodes=substring(j.nodes,first=2)
  }
  
  #Select appropriate columns from obs data
  obs.c=which(names(chik.dat) %in% paste("X",sim.nodes,sep=""))
  if(length(obs.c)<length(sim.nodes)){
    print(paste("Only", length(obs.c), "municipalities in observed match the",
                length(sim.nodes),"simulated municipalities",sep=" "))
  }
  
  obs.n=substring(names(chik.dat)[obs.c],first=2)
  #re-order if necessary
  obs.c=obs.c[unlist(lapply(obs.n,function(x) which(sim.nodes==x)))]
  obs.m=cbind(matrix(0,nrow=length(obs.c),ncol=73),t(chik.dat[,obs.c]))
  
  # ##weekly time series
  year.v=c(rep(2013,52),rep(2014,52),rep(2015,52),rep(2016,ncol(obs.m)-156))
  week.v=c(1:52,1:52,1:52,1:(ncol(obs.m)-156))
  week.year=paste(year.v,ifelse(week.v>10,week.v,paste("0",week.v,sep="")),sep="-")
  x.at=seq(1:ncol(obs.m))[seq(53,ncol(obs.m),13)]
  x.lab=week.year[seq(53,ncol(obs.m),13)]
  ##Order by size of observed outbreak?
  obs.m2=obs.m[rev(order(rowSums(obs.m))),]
  sim.m2=sim.m[rev(order(rowSums(obs.m))),]
  bestsim.m2=bestsim.m[rev(order(rowSums(obs.m))),]
  sim.a2=sim.a[,rev(order(rowSums(obs.m))),]
  sim.nodes2=sim.nodes[rev(order(rowSums(obs.m)))]
  
  plot(obs.m2[1,])
  lines(obs.m2[2,])
  lines(obs.m2[3,],lty=2)
  lines(obs.m2[4,],lty=3)
  lines(obs.m2[5,],lty=4)
  par(mfrow=c(1,1))
  # pdf(file=paste("d:/Colombia/muni_output/",dept,"_municipality_timeseries_fit_upd2.pdf",sep=""),useDingbats = F)
  # for(i in 1:nrow(obs.m)){
  #   #maxv=max(max(obs.m[i,]),max(sim.a[,i,-1]))
  #   maxv=max(max(obs.m2[i,]),max(sim.m2[i,]),max(bestsim.m2[i,]))
  #   plot(obs.m2[i,],pch=19,col="red",ylim=c(0,maxv+100),xlim=c(53,ncol(obs.m)),
  #        xlab="Week",ylab="Reported Cases",axes=F)
  #   axis(2)
  #   axis(1,at=x.at,labels=x.lab)
  #   box()
  #   ##plot simulations
  #   for (j in 1:(dim(sim.a)[1])){
  #     lines(sim.a2[j,i,-1],lwd=1.5,col=rgb(208/255,240/255,242/255,alpha=.25))
  #   }
  #   legend("topleft",c("Observed","Simulations"),pch=c(19,NA),
  #          lty=c(NA,1,1),col=c("red","blue"),lwd=c(NA,2),bty='n')
  #   points(obs.m2[i,],pch=19,col="red")
  #   #lines(sim.m2[i,],lwd=2,col=rgb(147/255,249/255,190/255))
  #   lines(sim.m2[i,],lwd=2,col="blue")
  #   legend("topright",mcodes$Mun_NAME[which(mcodes$Municipality==sim.nodes2[i])],bty="n")
  #   if(i%%10==0) print(i)
  # }
  # dev.off()
  
  pdf(file=paste("d:/Colombia/muni_output/figs/",dept,"_municipality_timeseries_fit_upd_ribbon.pdf",sep=""),useDingbats = F)
  
  for(i in 1:nrow(obs.m2)){
    #maxv=max(max(obs.m[i,]),max(sim.a[,i,-1]))
    obs1.df=data.frame(week=1:ncol(obs.m2),cases=obs.m2[i,])
    obs2.df=data.frame(week=1:ncol(obs.m2),cases=obs.m2[2,])
    obs3.df=data.frame(week=1:ncol(obs.m2),cases=obs.m2[3,])
    obs4.df=data.frame(week=1:ncol(obs.m2),cases=obs.m2[4,])
    maxv=max(max(obs.m2[i,]),max(sim.m2[i,]),max(bestsim.m2[i,]))
    sim1.m=unlist(apply(sim.a2[,i,-1],2,mean))
    sim1.md=unlist(apply(sim.a2[,i,-1],2,function(x) quantile(x,0.5)))
    sim1.l=unlist(apply(sim.a2[,i,-1],2,function(x) quantile(x,0.025)))
    sim1.h=unlist(apply(sim.a2[,i,-1],2,function(x) quantile(x,0.975)))
    sim1.df=data.frame(week=1:159,mean=sim1.m,cases=sim1.md,lo=sim1.l,hi=sim1.h)
    # sim2.m=unlist(apply(sim.a2[,2,-1],2,mean))
    # sim2.md=unlist(apply(sim.a2[,2,-1],2,function(x) quantile(x,0.5)))
    # sim2.l=unlist(apply(sim.a2[,2,-1],2,function(x) quantile(x,0.025)))
    # sim2.h=unlist(apply(sim.a2[,2,-1],2,function(x) quantile(x,0.975)))
    # sim2.df=data.frame(week=1:159,mean=sim2.m,cases=sim2.md,lo=sim2.l,hi=sim2.h)
    # sim3.m=unlist(apply(sim.a2[,3,-1],2,mean))
    # sim3.md=unlist(apply(sim.a2[,3,-1],2,function(x) quantile(x,0.5)))
    # sim3.l=unlist(apply(sim.a2[,3,-1],2,function(x) quantile(x,0.025)))
    # sim3.h=unlist(apply(sim.a2[,3,-1],2,function(x) quantile(x,0.975)))
    # sim3.df=data.frame(week=1:159,mean=sim3.m,cases=sim3.md,lo=sim3.l,hi=sim3.h)
    # sim4.m=unlist(apply(sim.a2[,4,-1],2,mean))
    # sim4.md=unlist(apply(sim.a2[,4,-1],2,function(x) quantile(x,0.5)))
    # sim4.l=unlist(apply(sim.a2[,4,-1],2,function(x) quantile(x,0.025)))
    # sim4.h=unlist(apply(sim.a2[,4,-1],2,function(x) quantile(x,0.975)))
    # sim4.df=data.frame(week=1:159,mean=sim4.m,cases=sim4.md,lo=sim4.l,hi=sim4.h)
    maxv=max(maxv,max(sim1.h))

    print(ggplot(obs1.df[-(1:52),],aes(x=week,y=cases))+geom_point(col="red")+#geom_line(color="blue",lty=1,lwd=1)+
      #geom_point(data=obs2.df[-(1:52),],aes(x=week,y=cases),col="red")+
      #geom_point(data=obs3.df[-(1:52),],aes(x=week,y=cases),col="green")+
      #geom_point(data=obs4.df[-(1:52),],aes(x=week,y=cases),col="orange")+
      #geom_ribbon(aes(ymin=sim2.l[-(1:52)],ymax=sim2.h[-(1:52)]),fill="red",alpha=0.25)+
      geom_ribbon(aes(ymin=sim1.l[-(1:52)],ymax=sim1.h[-(1:52)]),fill="blue",alpha=0.25)+
      geom_line(data=sim1.df[-(1:52),],aes(x=week,y=mean),color="blue",lwd=1,lty=1)+
      geom_line(data=sim1.df[-(1:52),],aes(x=week,y=cases),color="purple",lwd=1,lty=1)+
      theme(axis.text.x=element_text(angle=60,hjust=1),panel.background = element_rect(fill="white",color="black"))+
      scale_y_continuous(name="Weekly cases")+scale_x_continuous(name="Week",breaks=x.at,labels=x.lab)+
      coord_cartesian(ylim = c(0, maxv))) 
    #maxv2=max(sim2.h,sim3.h,sim4.h)
    # ggplot(obs2.df[-(1:52),],aes(x=week,y=cases))+geom_point(col="red")+geom_line(color="red",lty=1,lwd=1)+
    #   geom_point(data=obs3.df[-(1:52),],aes(x=week,y=cases),col="purple")+geom_line(data=obs3.df[-(1:52),],aes(x=week,y=cases),col="purple")+
    #   #geom_point(data=obs4.df[-(1:52),],aes(x=week,y=cases),col="green")+      geom_line(data=obs4.df[-(1:52),],aes(x=week,y=cases),col="green")+
    #   #geom_point(data=obs4.df[-(1:52),],aes(x=week,y=cases),col="orange")+
    #   geom_ribbon(aes(ymin=sim2.l[-(1:52)],ymax=sim2.h[-(1:52)]),fill="red",alpha=0.1)+
    #   geom_ribbon(aes(ymin=sim3.l[-(1:52)],ymax=sim3.h[-(1:52)]),fill="purple",alpha=0.1)+
    #   #geom_ribbon(aes(ymin=sim4.l[-(1:52)],ymax=sim4.h[-(1:52)]),fill="green",alpha=0.25)+
    #   #geom_line(data=sim1.df[-(1:52),],aes(x=week,y=cases),color="blue",lwd=1,lty=1)+
    # theme(axis.text.x=element_text(angle=60,hjust=1),panel.background = element_rect(fill="white",color="black"))+
    #   scale_y_continuous(name="Weekly cases")+scale_x_continuous(name="Week",breaks=x.at,labels=x.lab)+coord_cartesian(ylim = c(0, maxv2)) 
    # 

    
  }
  dev.off()    

