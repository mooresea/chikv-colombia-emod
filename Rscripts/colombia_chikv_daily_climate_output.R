

setwd("d:/dtk-aedes_scripts/CHIKV/Single_node_calibration_output_decay")
sim.dir=list.dirs(recursive=F,full.names = F)
#sim.dir=sim.dir[-which(sim.dir=="Colombia")]
sim.names=sim.dir
sim1.w=list()
for(i in 1:length(sim.dir)){
  idirs=list.dirs(path=sim.dir[i])
  idirs=idirs[-1]
  sim.infections=matrix(NA,nrow=length(idirs),ncol=obs.len)
  rownames(sim.infections)=unlist(lapply(strsplit(idirs,"_"),function(x)x[2]))
  
    wsim=fromJSON(txt=paste(idirs[1],"InsetChart.json",sep="/"))
    #sim.i=floor(as.integer(unlist(strsplit(dept.sim.files[i],"_"))[6])/2)+1
    #sample.i=as.integer(substring(unlist(strsplit(dept.sim.files[i],"_"))[7],1,1))
    itemp=wsim$Channels$`Air Temperature`$Data
    #Pop size
    irain=wsim$Channels$Rainfall$Data
    ## For now multiply cases by 10
    iweath=list(temp=itemp,rain=irain)
    sim1.w[[i]]=iweath
}

setwd("d:/dtk-aedes_scripts/CHIKV/chikv_admin2_spatial_output")
sim.dir=list.dirs(recursive=F,full.names = F)
#sim.dir=sim.dir[-which(sim.dir=="old")]
sim2.names=sim.dir
sim2.list=list()
sim2.w=list()
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
  
  wsim=fromJSON(txt=paste(idirs[1],"InsetChart.json",sep="/"))
  #sim.i=floor(as.integer(unlist(strsplit(dept.sim.files[i],"_"))[6])/2)+1
  #sample.i=as.integer(substring(unlist(strsplit(dept.sim.files[i],"_"))[7],1,1))
  itemp=wsim$Channels$`Air Temperature`$Data
  #Pop size
  irain=wsim$Channels$Rainfall$Data
  ## For now multiply cases by 10
  iweath=list(temp=itemp,rain=irain)
  sim2.w[[i]]=iweath
  print(i)
}