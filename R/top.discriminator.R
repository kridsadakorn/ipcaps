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

top.discriminator <- function(cluster.obj,group1,group2,bim.file,use.node.number=FALSE,num.top=100){
  
  if (is.null(cluster.obj)){
    cat(paste0("Incorrect parameter, please use the object returned from the function ipcaps as an input\n"))
    return(NULL)
  }
  
  raw.filename = file.path(cluster.obj$output.dir,"RData","rawdata.RData")
  if (!file.exists(raw.filename)){
    cat(paste0("Not found the rawdata file: ",raw.filename,"\n"))
    return(NULL)
  }else{
    load(raw.filename)
    if (is.null(snp.info)){
      if (!file.exists(bim.file)){
        cat(paste0("Not found the bim file: ",bim.file,"\n"))
        return(NULL)
      }else{
        snp.info <- read.table(file=bim.file, colClasses=c('factor','factor','factor','factor','factor','factor'))
      }
    }
  }
    
  
  if (num.top<0){
    cat(paste0("num.top must be more than zero\n"))
    return(NULL)
  }
  
  if (use.node.number == FALSE){
    index1 <- cluster.obj$cluster$row.number[which(cluster.obj$cluster$group == group1)]
    index2 <- cluster.obj$cluster$row.number[which(cluster.obj$cluster$group == group2)]  
  }else{
    index1 <- cluster.obj$cluster$row.number[which(cluster.obj$cluster$node == group1)]
    index2 <- cluster.obj$cluster$row.number[which(cluster.obj$cluster$node == group2)]    
  }
  
  all.fst <- fst.each.snp.hudson(raw.data,index1,index2)
  
  snp.with.fst <- cbind(snp.info,all.fst)
  colnames(snp.with.fst) <- c('chr','SNP','centimorgans','position','allele1','allele2','Fst')
  top.fst.snp <- snp.with.fst[order(-snp.with.fst$Fst),]
  ret <- head(top.fst.snp, n=num.top)
  
  return(ret)
}



