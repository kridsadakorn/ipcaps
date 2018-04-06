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

export.groups <- function(result.dir){
  load(file.path(result.dir,"RData","leafnode.RData"))
  export.data = NULL
  for (i in 1:length(leaf.node)){
    assigned_group = which(leaf.node == leaf.node[i])
    cat(paste0("Exporting node ",leaf.node[i]," as group ",assigned_group,"\n"))
    load(file.path(result.dir,"RData",paste0("node",leaf.node[i],".RData")))
    groups = rep(assigned_group,length(index))
    nodes = rep(leaf.node[i],length(index))
    node.data = data.frame(groups,nodes,label,index)
    export.data = rbind(export.data,node.data)
  }
  colnames(export.data) = c("group","node","label","row.number")
  file.name = file.path(result.dir,"groups.txt")
  if (file.exists(file.name)){
    i = 1
    file.name = file.path(result.dir,paste0("groups",i,".txt"))
    while (file.exists(file.name)){
      i = i + 1
      file.name = file.path(result.dir,paste0("groups",i,".txt"))
    }
  }
  cat(paste0("Note: save as ",file.name,"\n"))
  write.table(export.data,file=file.name,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
  return(export.data)
}



