## 
## Estimates of introduction rates for different calibs
##
library(jsonlite)
library(fields)
setwd("D:/dtk-aedes_scripts/CHIKV")

node="Colombia"
type="single"

if(type=="single"){
  fpath=paste("CHIKV_Calib_Habitat_Decay_Imports_SingleNode_",node,"/CalibManager.json",sep="")
}else{
  fpath=paste("CHIKV_Calib_Habitat_Multiple_Imports_admin2nodes_Migration_",node,"/CalibManager.json",sep="")
}

fpath="CHIKV_Calib_Habitat_Imports_Colombia/CalibManager.json"

calib.j=fromJSON(txt=fpath)

calib.params=cbind(calib.j$results$TEMPORARY_RAINFALL,calib.j$results$Temporary_Habitat_Decay_Factor,
                   eval(parse(text=paste("calib.j$results$",names(calib.j$results)[grep("Times",names(calib.j$results))],sep=""))),
                   eval(parse(text=paste("calib.j$results$",names(calib.j$results)[grep("Heights",names(calib.j$results))],sep=""))),
                   eval(parse(text=paste("calib.j$results$",names(calib.j$results)[grep("Widths",names(calib.j$results))],sep=""))))
if(type!="single"){
  calib.params=cbind(calib.params,calib.j$results$x_Regional_Migration)
}

calib.spl=Tps(calib.params,calib.j$results$total)

##get mean and peak marginal values for each parameter
params.mean=apply(calib.spl$x,2,mean)
params.max=calib.spl$x[1,]

p.mean.mat=cbind(matrix(rep(params.mean[1:4],nrow(calib.spl$x)),nrow=nrow(calib.spl$x),ncol=5,byrow=T),calib.spl$x[,6])
p.max.mat=cbind(matrix(rep(params.max[1:4],nrow(calib.spl$x)),nrow=nrow(calib.spl$x),ncol=5,byrow=T),calib.spl$x[,6])
##Mean of those better than no outbreak
noout.v=as.numeric(names(table(calib.spl$y))[which.max(table(calib.spl$y))])
nv.ind=min(which(round(calib.spl$y,1)==round(noout.v,1)))
params.mean2=apply(calib.spl$x[1:nv.ind,],2,mean)
#p.mean.mat2=cbind(matrix(rep(params.mean2[1:5],nv.ind),nrow=nv.ind,ncol=5,byrow=T),calib.spl$x[1:nv.ind,6])
p.mean.mat2=cbind(rep(params.mean2[1],nrow(calib.spl$x)),calib.spl$x[,2],rep(params.mean2[3],nrow(calib.spl$x)),rep(params.mean2[4],nrow(calib.spl$x)))
p.mean.mat4=cbind(rep(params.mean2[1],nrow(calib.spl$x)),rep(params.mean2[2],nrow(calib.spl$x)),rep(params.mean2[3],nrow(calib.spl$x)),calib.spl$x[,4])

#prof_log_lik=function(a){
#     b=(optim(1,function(z) -sum(log(dgamma(x,a,z)))))$par
#     return(-sum(log(dgamma(x,a,b))))
#}
out.p=predict(calib.spl,x=p.mean.mat2)
plot(calib.spl$x[,2],out.p,ylim=c(max(out.p)-100,max(out.p)+100),pch=19)

     ylab="Negative log likelihood",xlab="Symptomatic reporting rate")


#out.path=paste("d:/dtk-aedes_scripts/chikv/zika_output/figures/",admin.nm,"_reporting-rate_LL.pdf",sep="")
pdf(file=out.path,pointsize=12,useDingbats = F)
par(mfrow=c(1,1))
plot(calib.spl$x[,6],predict(calib.spl,x=p.mean.mat2),ylim=c(-40000,0),xlim=c(0,.25),pch=19,
     ylab="Negative log likelihood",xlab="Symptomatic reporting rate")
dev.off()

par(mfrow=c(1,3))
plot(calib.spl$x[,6],calib.spl$y,ylim=c(-10000,0))
plot(calib.spl$x[1:nv.ind,6],predict(calib.spl,x=p.mean.mat2),ylim=c(-10000,0))
plot(calib.spl$x[,6],predict(calib.spl,x=p.max.mat),ylim=c(-10000,0))  



