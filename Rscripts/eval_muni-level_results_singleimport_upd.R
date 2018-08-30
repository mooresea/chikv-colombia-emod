##
## Comparison of simulated and observed municipality data
##
setwd("d:/")

#
# Observed data
#
chik.dat=read.csv(file="Colombia/laptop/Chikv_timeseries.csv",header=T,as.is=T)
mcodes=read.table(file="Colombia/Unique_divipola_codes.txt",sep=":",header=T,as.is = T)

##
## Siumulated data
dept.names=list.dirs(path="d:/dtk-aedes_scripts/CHIKV/Dept_admin2_nodes_calibration_output/",recursive=F,full.names = F)
dept.names=dept.names[-which(dept.names=="old")]

###
### Function to calcuate various eval metrics
###
eval.calcs=function(obs.dat,sim.dat,calc="mae"){
  #Remove first year since no introductions
  obs.dat=obs.dat[-(1:52)]
  sim.dat=sim.dat[-(1:52)]
  if(calc=="mae"){
    mval=sum(abs(obs.dat-sim.dat))/length(sim.dat)
  }else{
    if(calc=="rmse"){
      mval=sqrt(sum((obs.dat-sim.dat)^2)/length(sim.dat))
    }else{
      if(calc=="mase"){
        m1=sum(abs(obs.dat-sim.dat))
        m2=sum(abs(diff(obs.dat)))
        m2=ifelse(m2>0,m2,1)
        mval=m1/((length(sim.dat)/(length(sim.dat)-1))*m2)
        malt=sum(abs(obs.dat-sim.dat)/(m2/(length(sim.dat)-1)))/length(sim.dat)
      }else{
        if(calc=="rmse.norm"){
          mval=sqrt(sum((obs.dat-sim.dat)^2)/length(sim.dat))
          mval=mval/ifelse(mean(obs.dat)>0,mean(obs.dat),1/length(obs.dat))
        }
        if(calc=="rmse.norm2"){
          mval=sqrt(sum((obs.dat-sim.dat)^2)/length(sim.dat))
          mval=mval/ifelse(max(obs.dat)>0,max(obs.dat),1)
        }
      }
    }
  }
  return(mval)
}

##
## Convert case count ts to incidence ts
## 
calc.inc=function(dmat){
  for(rind in 1:nrow(dmat)){
    dmat[rind,]=dmat[rind,]*1000/j.pop[rind]
  }
  return(dmat)
}

dept.files=list.files(path="d:/Colombia/muni_output/results")
dept.names=unlist(lapply(strsplit(dept.files,"_"),function(x) x[1]))
##Calc metrics and plot fits for each department
sims.sd=list()
for(dd in 1:length(dept.names)){
  dept=dept.names[dd]
  #LOAD Data
  deptf=paste("d:/Colombia/muni_output/results/",dept.files[dd],sep="")
  load(deptf)
  #deptf2=paste("d:/Colombia/muni_output/",dept,"__spatial_sim_results_incidence.rda",sep="")
  #load(deptf2)
  
  ##Get pop data from demographics file
  deptdf=paste("d:/Eradication/InputDataFiles/Colombia_dept_nodes_admin2/",dept,"/",dept,"_admin2_nodes_demographics.json",sep="")
  dept.j=fromJSON(txt=deptdf)
  
  print(paste("Begin eval for",dept,":",dd,"out of",length(dept.names),"departments"))
  
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

  ###Incidence
  obs.i.m=calc.inc(obs.m)
  sim.i.m=calc.inc(sim.m)
  bestsim.i.m=calc.inc(bestsim.m)
  sim.i.a=array(NA,dim=dim(sim.a))
  for(aa in 1:(dim(sim.a)[1])){
    aa.mat=calc.inc(sim.a[aa,,-1])
    sim.i.a[aa,,]=cbind(sim.a[aa,,1],aa.mat)
    if(aa%%50==0) print(aa)
  }
  
  ##Calculate different error measurements by node
  best.mae=rep(NA,nrow(obs.m))
  best.rmse=best.mae
  best.mase=best.mae
  best.rmse.n1=best.mae
  best.rmse.n2=best.mae
  all.mae=best.mae
  all.rmse=best.mae
  all.rmse.n1=best.mae
  all.rmse.n2=best.mae
  all.mase=best.mae
  
  for(i in 1:nrow(obs.m)){
    #MAE
    all.mae[i]=eval.calcs(obs.m[i,],sim.m[i,],calc="mae")
    best.mae[i]=eval.calcs(obs.m[i,],bestsim.m[i,],calc="mae")
    #RMSE
    all.rmse[i]=eval.calcs(obs.m[i,],sim.m[i,],calc="rmse")
    best.rmse[i]=eval.calcs(obs.m[i,],bestsim.m[i,],calc="rmse")
    #Normalized RMSE
    all.rmse.n1[i]=eval.calcs(obs.m[i,],sim.m[i,],calc="rmse.norm")
    best.rmse.n1[i]=eval.calcs(obs.m[i,],bestsim.m[i,],calc="rmse.norm")
    all.rmse.n2[i]=eval.calcs(obs.m[i,],sim.m[i,],calc="rmse.norm2")
    best.rmse.n2[i]=eval.calcs(obs.m[i,],bestsim.m[i,],calc="rmse.norm2")
    #MASE
    all.mase[i]=eval.calcs(obs.m[i,],sim.m[i,],calc="mase")
    best.mase[i]=eval.calcs(obs.m[i,],bestsim.m[i,],calc="mase")
  }
  ##Plots
  
  ##Cor for each simulation
  sim.sums=(apply(sim.a,1,function(x) rowSums(x[,-1])))
  sim.sums.i=sim.sums/j.pop
  sim.cor=apply(sim.sums,2,function(x) cor(x,rowSums(obs.m)))
  sim.cor.i=apply(sim.sums.i,2,function(x) cor(x,rowSums(obs.m)/j.pop))
  sim.cor.df=data.frame(dept=dept,cor=sim.cor)
  sim.cori.df=data.frame(dept=dept,cor=sim.cor.i)
  
  ##Get standard deviation for muni-level incidence for each simulation
  sims.sd[[dd]]=apply(sim.sums.i,2, sd)
  
  if(dd==1){
    all.cor.df=sim.cor.df
    all.cori.df=sim.cori.df
  }else{
    all.cor.df=rbind(all.cor.df,sim.cor.df)
    all.cori.df=rbind(all.cori.df,sim.cori.df)
  }
  
  ##Correlation between obs and simulated muni-level final totals
  corp=cor((rowSums(sim.m)+1),(rowSums(obs.m)+1))
  corpb=cor((rowSums(bestsim.m)+1),(rowSums(obs.m)+1))
  ##Correlation between obs and simulated muni-level final totals
  corp.i=cor((rowSums(sim.i.m)+1),(rowSums(obs.i.m)+1))
  corpb.i=cor((rowSums(bestsim.i.m)+1),(rowSums(obs.i.m)+1))
  cora=rep(NA,dim(sim.a)[1])
  cora.i=cora
  for(yy in 1:(dim(sim.a)[1])){
    cora[yy]=cor((rowSums(sim.a[yy,,-1])+1),(rowSums(obs.m)+1))
    cora.i[yy]=cor((rowSums(sim.i.a[yy,,-1])+1),(rowSums(obs.i.m)+1))
  }
  
  ## Correlations for log data
  ##Correlation between obs and simulated muni-level final totals
  corp=cor(log10(rowSums(sim.m)+1),log10(rowSums(obs.m)+1))
  corpb=cor(log10(rowSums(bestsim.m)+1),log10(rowSums(obs.m)+1))
  ##Correlation between obs and simulated muni-level final totals
  corp.i=cor(log10(rowSums(sim.i.m)+1),log10(rowSums(obs.i.m)+1))
  corpb.i=cor(log10(rowSums(bestsim.i.m)+1),log10(rowSums(obs.i.m)+1))


  dd.df=data.frame(dept=dept,dept.id=dept.id,muni=substring(row.names(obs.m),first=2),
                   mae=all.mae,mase=all.mase,rmse=all.rmse,rmse.n1=all.rmse.n1,rmse.n2=all.rmse.n2,
                   mae.best=best.mae,mase.best=best.mase,rmse.best=best.rmse,rmse.n1.best=best.rmse.n1,rmse.n2.best=best.rmse.n2,
                  cor=corp,cor.best=corpb,cori=corp.i,cor.best.i=corpb.i,
                  obs.cases=rowSums(obs.m),obs.inc=rowSums(obs.i.m),
                  sim.cases=rowSums(sim.m),sim.inc=rowSums(sim.i.m),
                  sim.inc.m=apply((apply(sim.i.a,1,function(x) rowSums(x[,-1]))),1,function(x) quantile(x,0.5)),
                  sim.inc.l=apply((apply(sim.i.a,1,function(x) rowSums(x[,-1]))),1,function(x) quantile(x,0.025)),
                  sim.inc.h=apply((apply(sim.i.a,1,function(x) rowSums(x[,-1]))),1,function(x) quantile(x,0.975)))
  
  if(dd==1){
    eval.df=dd.df
  }else{
    eval.df=rbind(eval.df,dd.df)
  }

}

##Infections
s_r=0.72
r_r=0.08
eval.df$iar=eval.df$sim.inc/(1000*s_r*r_r)

library(ggplot2)

ggplot(eval.df,aes(x=dept,y=mase))+geom_boxplot()+
  theme(axis.text.x=element_text(angle=60,hjust=1))

#pdf(file="d:/Colombia/muni_output/correlation_obs_vs_sim_cases_and_incid.pdf",width=10,height=6)
par(mfrow=c(1,2))
hist(eval.df$cor,breaks=20,ylim=c(0,300),main="Cases")
abline(v=mean(eval.df$cor),lty=2,col="red")
legend("topleft","mean Pearson's r",lty=2,col="red",bty="n")
hist(eval.df$cori,breaks=20,ylim=c(0,300),main="Incidence")
abline(v=mean(eval.df$cori),lty=2,col="red")
legend("topleft","mean Pearson's r",lty=2,col="red",bty="n")

cor.df1=data.frame(dept=eval.df$dept,cor=eval.df$cor)
cor.df2=data.frame(dept=eval.df$dept,cor=eval.df$cori)
cor.df1$meas="cases"
cor.df2$meas="incid"
cor.df=rbind(cor.df1,cor.df2)
ggplot(cor.df,aes(x=cor,fill=meas))+geom_density(alpha=0.3)
ggplot(cor.df,aes(x=cor,fill=meas))+
  geom_histogram(alpha=0.5,position = "identity",binwidth = 0.05)
#dev.off()

##
##
##
library(plyr)

####
#### Generate muni-level results from single-node models 
####
#load("d:/dtk-aedes_scripts/CHIKV/calibration_figures/combo/upd/model_fit_measures_admin1+2.rda")

load("d:/dtk-aedes_scripts/CHIKV/calibration_figures/combo/upd/model_fit_measures_admin1+2.rda")

### Get mean sd for each dept
dept.sims.sd=unlist(lapply(sims.sd, mean))
#
#
#
dept.list=sim.list
dept.list[[which(sim.names=="BogotaDC")]]=NULL
for(i in 1:length(sim2.names)){
  imat=dept.list[[i]]
  itotals=rowSums(imat)
  dept=sim2.names[i]
  
  deptf=paste("d:/Colombia/muni_output/results/",dept.files[i],sep="")
  load(deptf)
  #deptf2=paste("d:/Colombia/muni_output/",dept,"__spatial_sim_results_incidence.rda",sep="")
  #load(deptf2)
  
  ##Get pop data from demographics file
  deptdf=paste("d:/Eradication/InputDataFiles/Colombia_dept_nodes_admin2/",dept,"/",dept,"_admin2_nodes_demographics.json",sep="")
  dept.j=fromJSON(txt=deptdf)
  
  print(paste("Begin eval for",dept,":",dd,"out of",length(dept.names),"departments"))
  
  sim.nodes=substring(sim.a[1,,1],first=3)
  j.nodes=substring(dept.j$Nodes$NodeID,first=3)
  j.pop=dept.j$Nodes$NodeAttributes$InitialPopulation
  j.pop=j.pop[unlist(lapply(sim.nodes,function(x) which(j.nodes==x)))]
  pop.sum=sum(j.pop)
  pop.pct=j.pop/pop.sum
  
  muni.mat=matrix(NA,nrow=length(pop.pct),ncol=length(itotals))
  for(j in 1:nrow(muni.mat)){
    muni.mat[j,]=itotals*pop.pct[j]
  }
  
  ###
  ### Random draw from multinomial
  ### 
  muni.mat.multi=rmultinom(length(itotals),size=itotals,prob=pop.pct)
  
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

  ##Cor for each simulation
  #sim.sums=(apply(sim.a,1,function(x) rowSums(x[,-1])))
  ##Null expectations

  sim.sums=muni.mat.multi
  sim.sums.i=muni.mat.multi/j.pop
  ##Use SD from multi-node
  sim.sums.i2=sim.sums.i
  sim.sums.i2=matrix(rnorm(length(sim.sums.i2),mean=sim.sums.i2,sd=dept.sims.sd[i]),nrow=nrow(sim.sums.i),ncol=ncol(sim.sums.i))
  sim.sums.i2=ifelse(sim.sums.i2<0,0,sim.sums.i2)
  ##Add measurement error
  #sim.sums.i=sim.sums.i+rnorm(length(sim.sums.i),mean=0,sd=mean(sim.sums.i)/20)
  #sim.sums.i=ifelse(sim.sums.i<0,0,sim.sums.i)
  sim.cor=apply(sim.sums,2,function(x) cor(x,rowSums(obs.m)))
  sim.cor.i=apply(sim.sums.i,2,function(x) cor(round(x,8),rowSums(obs.m)/j.pop))
  sim.cor.i2=apply(sim.sums.i2,2,function(x) cor(round(x,8),rowSums(obs.m)/j.pop))
  sim.cor.df=data.frame(dept=dept,cor=sim.cor)
  sim.cori.df=data.frame(dept=dept,cor=sim.cor.i)
  sim.cori.df2=data.frame(dept=dept,cor=sim.cor.i2)
  if(i==1){
    s.cor.df=sim.cor.df
    s.cori.df=sim.cori.df
    s.cori.df2=sim.cori.df2
  }else{
    s.cor.df=rbind(s.cor.df,sim.cor.df)
    s.cori.df=rbind(s.cori.df,sim.cori.df)
    s.cori.df2=rbind(s.cori.df2,sim.cori.df2)
  }
  
}

#sim1=sim.df[sim.df$lev==1,]
sim.df=sim.df[sim.df$lev==2,]
sim.df.sum=ddply(sim.df,.(dept), summarize,
                 mae=mean(mae),
                 rmse=mean(rmse),
                 mase=mean(mase),
                 rmse.n=mean(rmse.n))
mase.tri=quantile(sim.df.sum$mase,c(0.33,0.66))
mase.qtr=quantile(sim.df.sum$mase,c(0.25,0.5,0.75))
rmse.qtr=quantile(sim.df.sum$rmse.n,c(0.25,0.5,0.75))

sim.df.sum$tri=ifelse(sim.df.sum$mase<mase.tri[1],1,ifelse(sim.df.sum$mase>mase.tri[2],3,2))
sim.df.sum$qtr=ifelse(sim.df.sum$mase<mase.qtr[1],1,
                      ifelse(sim.df.sum$mase<mase.qtr[2],2,
                      ifelse(sim.df.sum$mase<mase.qtr[3],3,4)))
sim.df.sum$qtr.rmse=ifelse(sim.df.sum$rmse.n<rmse.qtr[1],1,
                      ifelse(sim.df.sum$rmse.n<rmse.qtr[2],2,
                             ifelse(sim.df.sum$rmse.n<rmse.qtr[3],3,4)))


all(eval.df$dept %in% sim.df.sum$dept)
eval.df$tri=NA
eval.df$qtr=NA
all.cor.df$tri=NA
all.cor.df$qtr=NA
all.cor.df$qtr.rmse=NA
all.cori.df$tri=NA
all.cori.df$qtr=NA
all.cori.df$qtr.rmse=NA
eval.df$dept.mase=NA
eval.df$dept.rmse=NA
for (i in 1:nrow(sim.df.sum)){
  idept=sim.df.sum$dept[i]
  eval.df$dept.mase[which(eval.df$dept==idept)]=sim.df.sum$mase[i]
  eval.df$dept.rmse[which(eval.df$dept==idept)]=sim.df.sum$rmse.n[i]
  eval.df$tri[which(eval.df$dept==idept)]=sim.df.sum$tri[i]
  eval.df$qtr[which(eval.df$dept==idept)]=sim.df.sum$qtr[i]
  all.cor.df$tri[which(all.cor.df$dept==idept)]=sim.df.sum$tri[i]
  all.cor.df$qtr[which(all.cor.df$dept==idept)]=sim.df.sum$qtr[i]
  all.cor.df$qtr.rmse[which(all.cor.df$dept==idept)]=sim.df.sum$qtr.rmse[i]
  all.cori.df$tri[which(all.cori.df$dept==idept)]=sim.df.sum$tri[i]
  all.cori.df$qtr[which(all.cori.df$dept==idept)]=sim.df.sum$qtr[i]
  all.cori.df$qtr.rmse[which(all.cori.df$dept==idept)]=sim.df.sum$qtr.rmse[i]
}

##Histogram
ggplot(cor.df,aes(x=cor,fill=meas))+geom_density(alpha=0.3)
ggplot(cor.df,aes(x=cor,fill=meas))+
  geom_histogram(alpha=0.5,position = "identity",binwidth = 0.1)

all.cor.df$mase=ifelse(all.cor.df$qtr<3,1,0)
all.cori.df$mase=ifelse(all.cori.df$qtr<3,1,0)
eval.df$mase.lev=ifelse(eval.df$qtr<3,1,0)

all.cor.df$type="c"
all.cori.df$type="i"
all.c.df=rbind(all.cor.df,all.cori.df)

s.cor.df$type="c"
s.cori.df$type="i"
s.cori.df2$type="i"
s.c.df=rbind(s.cor.df,s.cori.df)
s.c.df2=rbind(s.cor.df,s.cori.df2)

quantile(all.cor.df$cor, c(0.025,.5,.975),na.rm = T)
quantile(all.cori.df$cor, c(0.025,.5,.975),na.rm = T)
quantile(s.cor.df$cor, c(0.025,.5,.975),na.rm = T)
quantile(s.cori.df$cor, c(0.025,.5,.975),na.rm = T)

quantile(all.cor.df$cor,na.rm = T)
quantile(all.cori.df$cor,na.rm = T)
quantile(s.cor.df$cor,na.rm = T)
quantile(s.cori.df$cor,na.rm = T)
quantile(s.cori.df2$cor,na.rm = T)

##Combine single-node and multi-node
all.c.df$sim="multi"
s.c.df$sim="single"
sim.c.df=rbind(all.c.df[,c(1,2,7,8)],s.c.df)
s.c.df2$sim="single"
sim.c.df2=rbind(all.c.df[,c(1,2,7,8)],s.c.df2)

library(RColorBrewer)

pdf(file=paste(fig.dir,"muni_cases_corr_singleimport.pdf"),width=12,height=8,useDingbats = F,pointsize = 14)
ggplot(all.c.df[all.c.df$type=="c",],aes(x=cor))+geom_histogram(binwidth = 0.05,fill="red",alpha=0.5)+theme_bw()+
  xlim(-1.05,1.05)+geom_vline(xintercept=0,linetype=2)+ylim(0,3000)
dev.off()
pdf(file=paste(fig.dir,"muni_incid_corr_singleimport.pdf"),width=12,height=8,useDingbats = F,pointsize = 14)
ggplot(all.c.df[all.c.df$type=="i",],aes(x=cor))+geom_histogram(binwidth = 0.05,fill="red",alpha=0.5)+theme_bw()+
  xlim(-1.05,1.05)+geom_vline(xintercept=0,linetype=2)+ylim(0,1500)
dev.off()
pdf(file=paste(fig.dir,"null_cases_corr.pdf"),width=12,height=8,useDingbats = F,pointsize = 14)

ggplot(s.c.df[s.c.df$type=="c",],aes(x=cor))+geom_histogram(binwidth = 0.05,fill="blue",alpha=0.5)+theme_bw()+
  xlim(-1.05,1.05)+geom_vline(xintercept=0,linetype=2)+ylim(0,3000)
dev.off()
pdf(file=paste(fig.dir,"null_incid_corr.pdf"),width=12,height=8,useDingbats = F,pointsize = 14)

ggplot(s.c.df[s.c.df$type=="i",],aes(x=cor))+geom_histogram(binwidth = 0.05,fill="blue",alpha=0.5)+theme_bw()+
  xlim(-1.05,1.05)+geom_vline(xintercept=0,linetype=2)+ylim(0,1500)
dev.off()

# The palette with black:
cbbPalette <- alpha(c("red","blue"),0.25)

# To use for fills, add
pdf(file=paste(fig.dir,"muni_incid_corr_singleimport.pdf"),width=12,height=8,useDingbats = F,pointsize = 14)
ggplot(sim.c.df[sim.c.df$type=="i",],aes(x=cor,fill=as.factor(sim)))+theme_bw()+
  geom_histogram(alpha=0.5,binwidth = 0.05,position="identity",alpha=0.5)+
  scale_fill_manual(values=alpha(c("red", "blue"), 0.1))
dev.off()

pdf(file=paste(fig.dir,"muni_cases_corr_nolegend_singleimport.pdf"),width=12,height=8,useDingbats = F,pointsize = 14)
ggplot(sim.c.df[sim.c.df$type=="c",],aes(x=cor,fill=as.factor(sim)))+theme_bw()+
  geom_histogram(alpha=0.5,binwidth = 0.05,position="identity",alpha=0.5)+xlim(-1.05,1.05)+
  scale_fill_manual(values=alpha(c("red", "blue"), 0.1),guide=FALSE)+geom_vline(xintercept=0,linetype=2)
dev.off()

##
##
## Revised figures
##
pdf(file=paste(fig.dir,"muni_incid_corr_nolegend_singleimport_sd.pdf"),width=12,height=8,useDingbats = F,pointsize = 14)
ggplot(sim.c.df2[sim.c.df2$type=="i",],aes(x=cor,fill=as.factor(sim)))+theme_bw()+
  geom_histogram(alpha=0.5,binwidth = 0.05,position="identity",alpha=0.5)+xlim(-1.05,1.05)+
  scale_fill_manual(values=alpha(c("red", "blue"), 0.1),guide=FALSE)+geom_vline(xintercept=0,linetype=2)
dev.off()

pdf(file=paste(fig.dir,"muni_incid_corr_nolegend_singleimport_multinomial.pdf"),width=12,height=8,useDingbats = F,pointsize = 14)
ggplot(sim.c.df[sim.c.df$type=="i",],aes(x=cor,fill=as.factor(sim)))+theme_bw()+
  geom_histogram(alpha=0.5,binwidth = 0.05,position="identity",alpha=0.5)+xlim(-1.05,1.05)+
  scale_fill_manual(values=alpha(c("red", "blue"), 0.1),guide=FALSE)+geom_vline(xintercept=0,linetype=2)
dev.off()



ggplot(all.c.df,aes(x=cor,fill=as.factor(type)))+
  geom_histogram(alpha=0.5,binwidth = 0.05,position="identity")+theme_bw()+
  scale_fill_discrete(name="",breaks=c(0,1),labels=c("Cases","Incidence"))

ggplot(s.c.df,aes(x=cor,fill=as.factor(type)))+
  geom_histogram(alpha=0.5,binwidth = 0.05,position="identity")+theme_bw()+
  scale_fill_discrete(name="",breaks=c(0,1),labels=c("Cases","Incidence"))



pdf(file="d:/Colombia/muni_output/figs/final_size_hist.pdf",width=10,height=8,useDingbats = F)
ggplot(all.cor.df,aes(x=cor,fill=as.factor(mase)))+
  geom_histogram(alpha=0.5,binwidth = 0.05)+theme_bw()+
  scale_fill_discrete(name="Department\nMASE",breaks=c(0,1),labels=c("Above\naverage","Below\naverage"))
dev.off()
pdf(file="d:/Colombia/muni_output/figs/final_incid_hist.pdf",width=10,height=8,useDingbats = F)
ggplot(all.cori.df,aes(x=cor,fill=as.factor(mase)))+
  geom_histogram(alpha=0.5,binwidth = 0.05)+theme_bw()+
  scale_fill_discrete(name="Department\nMASE",breaks=c(0,1),labels=c("Above\naverage","Below\naverage"))
dev.off()
## boxplot
pdf(file="d:/Colombia/muni_output/figs/final_size_boxplot.pdf",width=10,height=8,useDingbats = F)
ggplot(eval.df,aes(x=obs.cases,y=sim.cases,color=as.factor(mase.lev)))+
  geom_point()+scale_x_log10(breaks=c(1,10,100,1000,10000),labels=c(1,10,100,1000,10000),name="Observed cases")+
  scale_y_log10(breaks=c(1,10,100,1000,10000),labels=c(1,10,100,1000,10000),name="Simulated cases")+theme_bw()+
  geom_abline(intercept=0,slope=1)
dev.off()
pdf(file="d:/Colombia/muni_output/figs/final_incidence_boxplot.pdf",width=10,height=8,useDingbats = F)
ggplot(eval.df,aes(x=obs.inc,y=sim.inc,color=as.factor(mase.lev)))+
  geom_point()+scale_x_log10(breaks=c(0.1,1,10,100,1000,10000),labels=c(0.1,1,10,100,1000,10000),name="Observed reported incidence (per 1,000)")+
  scale_y_log10(breaks=c(0.1,1,10,100,1000,10000),labels=c(0.1,1,10,100,1000,10000),name="Simulated reported incidence (per 1,000)")+theme_bw()+
  geom_abline(intercept=0,slope=1)
dev.off()
# pdf(file="d:/Colombia/muni_output/final_size_hist.pdf",width=10,height=8,useDingbats = F)
# ggplot(all.cor.df,aes(x=cor,fill=as.factor(mase)))+
#   geom_histogram(alpha=0.5,binwidth = 0.05)+theme_bw()+
#   scale_fill_discrete(name="Department\nMASE",breaks=c(0,1),labels=c("Above\naverage","Below\naverage"))
# dev.off()
# ## boxplot
# pdf(file="d:/Colombia/muni_output/final_size_boxplot.pdf",width=10,height=8,useDingbats = F)
# ggplot(eval.df,aes(x=obs.cases,y=sim.cases,color=as.factor(mase.lev)))+
#   geom_point()+scale_x_log10(breaks=c(1,10,100,1000,10000),labels=c(1,10,100,1000,10000))+
#   scale_y_log10(breaks=c(1,10,100,1000,10000),labels=c(1,10,100,1000,10000))+theme_bw()+
#   geom_abline(intercept=0,slope=1)
# dev.off()
