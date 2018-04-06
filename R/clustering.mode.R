##########################################################################
## IPCAPS Library
## Author: Kridsadakorn Chaichoompu
## Description:
##    This code is a part of Iterative Pruning to CApture Population 
##    Structure (IPCAPS) Library
##    
##Licence: GPL V3
## 
##    Copyright (C) 2016  Kridsadakorn Chaichoompu
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.

clustering.mode <- function(node,work.dir,method){
  start.time = Sys.time()
  
  cat(paste0("Node ",node,": Start clustering\n"))
  
  
  #load experiment condition
  file.name = file.path(work.dir,"RData",paste0("node",node,".RData"))
  cat(paste0("Node ",node,": Loading ",file.name,"\n"))
  load(file=file.name)
  
  if (method == "clara"){
    cl=pamk(PCs,krange=2:2,usepam=FALSE)
    cluster=cl$pamobject$clustering
  }else if (method == "pam"){
    cl=pamk(PCs,krange=2:2,usepam=TRUE)
    cluster=cl$pamobject$clustering
  }else if (method == "mixmod"){
    subPCs=as.data.frame(PCs[,1:no.significant.PC])
    mmc <- mixmodCluster(data=subPCs,nbCluster=1:3,criterion=c("BIC","ICL","NEC"),models=NULL)
    pmm = mixmodPredict(data=subPCs, classificationRule= mmc["bestResult"])
    cluster=pmm['partition']
  }else if (method == "meanshift"){
    subPCs=as.data.frame(PCs[,1:no.significant.PC])
    fit <- ms(subPCs, h=0.095)
    cluster=fit$cluster.label
  }else if (method == "apcluster"){
    subPCs=as.data.frame(PCs[,1:no.significant.PC])
    ap1 <- apcluster(negDistMat(r=5),subPCs,q=0.001)
    #ap2 <- aggExCluster(x=ap1)
    cluster=rep(2,dim(subPCs)[1])
    cluster[ap1[[2]][[1]]]=1
  }else if (method == "hclust"){
    dis=dist(PCs)
    hc=hclust(dis, method = "complete")
    cluster=cutree(hc, k=2)
  }else if (method == "rubikclust"){
    cluster=rubikClust(PCs[,1:3],min.space = 0.15)
  }else{#default Mixed clustering methods
    cluster=rubikClust(PCs[,1:3])
    if (length(unique(cluster)) == 1){ #swith to mixmod
      subPCs=as.data.frame(PCs[,1:no.significant.PC])
      mmc <- mixmodCluster(data=subPCs,nbCluster=1:3,criterion=c("BIC","ICL","NEC"),models=NULL)
      pmm = mixmodPredict(data=subPCs, classificationRule= mmc["bestResult"])
      cluster=pmm['partition']
    }
  }

  end.time = Sys.time()
  cat(paste0("Node ",node,": done for clustering\n"))
  print(end.time - start.time)
  return(cluster)
  
}



