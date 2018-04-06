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

get.node.info <- function(cluster.obj,node){
  if (is.null(cluster.obj$output.dir)){
    cat(paste0("Incorrect parameter, please use the object returned from the function ipcaps as an input\n"))
    return(NULL)
  }else{
    file.name <- file.path(cluster.obj$output.dir,"RData",paste0("node",node,".RData"))
    if (!file.exists(file.name)){
      cat(paste0("Node ",node," doesn't exist\n"))
      return(NULL)
    }else{
      load(file.name)
      res <- list('PCs'=PCs,'eigen.fit'=eigen.fit,'index'=index,'label'=label)
      return(res)
    }
  }
}


