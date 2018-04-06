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

#input.data, rows = individuals, columns = markers
preprocess <- function( files, label.fname, lab.col, rdata.infile, bed.infile, cate.list, 
                        result.dir, threshold, min.fst, max.thread=NA, reanalysis=FALSE, method="mix", 
                        min.in.group=20, datatype="snp", nonlinear = FALSE ,missing.char=NA,
                        regression.file=NA,regression.col.first=NA,regression.col.last=NA, 
                        reg.method="linear", plot.as.pdf=NA,no.plot=NA){  
  
  replace.missing <- function(X,missing=NA,rep){
    if (is.na(missing)){
      idx = which(is.na(X))
      X[idx] = rep[idx]
    }else{
      idx = which(X == missing)
      X[idx] = rep[idx]
    }
    return(X)
  }
  
  do.glm = function(X,PC,method="linear"){
    if (method=="poisson"){
      ret = glm(X ~ PC ,family=poisson, na.action=na.exclude)$residuals
    } else if (method=="negative.binomial"){
      ret = glm.nb(X ~ PC, na.action=na.exclude)$residuals
    }else{
      ret = glm(X ~ PC, family= gaussian(), na.action=na.exclude)$residuals
    }
    
    return(ret)
  }
  
  if (reanalysis==FALSE){
    #New analysis
    if (file.exists(result.dir)){
      cat(paste0("Note: the directory '",result.dir,"' is existed."),"\n")
      sub.dir = "cluster_output"
      i = 0
      tmp.dir = file.path(result.dir,sub.dir)
      while (file.exists(tmp.dir)){
        i = i + 1
        tmp.dir = file.path(result.dir,paste0(sub.dir,i))
      }
      result.dir = tmp.dir
    }else{
      dir.create(file.path(result.dir), showWarnings = FALSE, recursive = TRUE)
    }
    
    
    cat("The result files will be saved to this directory:",result.dir,"\n")
    dir.create(file.path(result.dir), showWarnings = FALSE)
    img.dir = file.path(result.dir,"images")
    dir.create(file.path(img.dir), showWarnings = FALSE)
    rdata.dir = file.path(result.dir,"RData")
    dir.create(file.path(rdata.dir), showWarnings = FALSE)
    
    #create empty file
    leaf.node = c()
    file.name = file.path(result.dir,"RData","leafnode.RData")
    save(leaf.node,file=file.name)
    node.list = c()
    
  }else{
    #Old analysis
    if (file.exists(result.dir)){
      cat(paste0("Note:  directory '",result.dir,"' is existed."),"\n")
    }else{
      cat(paste0("Error:  directory '",result.dir,"' is not existed. Please refer to a result directory."),"\n")
      quit()
    }
  }
  
  #Raw data
  index = NULL
  file.name = file.path(result.dir,"RData","rawdata.RData")
  if (!is.na(files)){
    label = NULL
    raw.data = NULL
    snp.info = NULL
    ind.info = NULL
    
    if (!file.exists(file.name)){
      #import label
      label = read.table(file=label.fname, header=FALSE)
      label = label[,lab.col]
      #load genotype files
      cat(paste0("Loading the input files ... \n"))
      for (fname in files){
        #test the separator
        oneline = read.table(file=fname, header=FALSE, sep=',',nrows=1)
        #separated by comma
        if (dim(oneline)[2]>1){
          tmp_geno = read.table(file=fname, header=FALSE, sep=',')
        }else{
          #separated by white space
          tmp_geno = read.table(file=fname, header=FALSE)
        }
        
        raw.data = rbind(raw.data,tmp_geno)
      }
      # data needs to be like rows = individuals, columns = markers
      raw.data = t(raw.data) 
      
      index = seq(1,length(raw.data[,1]))
      label = unlist(label)
      label = as.factor(label)
      
      tmp = as.numeric(raw.data)
      tmp1 = matrix(tmp,nrow=length(raw.data[,1]))
      tmp2 = as.data.frame(tmp1)
      rownames(tmp2) = index
      raw.data = tmp2
      
    }else{
      cat(paste0("Loading: the file '",file.name,"' is existed."),"\n")
      #load rawdata file to prepare node 1
      load(file.name)
      index = seq(1,length(raw.data[,1]))
    }
  }else if (!is.na(rdata.infile)){
    load(rdata.infile)
    if (!exists("raw.data")){
      cat(paste0("Error: Can't find 'raw.data' in file ",rdata.infile,"\n"))
      quit()
    }
    if (!exists("label")){
      cat(paste0("Error: Can't find 'label' in file ",rdata.infile,"\n"))
      quit()
    }
    index = seq(1,length(raw.data[,1]))
    
  }else if (!is.na(bed.infile)){
    #load BED files
    prefix = gsub('.bed','',bed.infile)
    bed = paste0(prefix,".bed")
    bim = paste0(prefix,".bim")
    fam = paste0(prefix,".fam")
    bed.data = read.bed(bed,bim,fam,only.snp=FALSE)
    raw.data = bed.data$snp
    snp.info = bed.data$snp.info
    ind.info = bed.data$ind.info
    index = seq(1,length(raw.data[,1]))
    raw.data = as.data.frame(raw.data)
    rownames(raw.data) = index
    
    #import label
    label = read.table(file=label.fname, header=FALSE)
    label = label[,lab.col]
    label = unlist(label)
    label = as.factor(label)
    
  }else{
    #load CATegorical files
    raw.data = pasre.categorical.data(cate.list)
    index = seq(1,length(raw.data[,1]))
    
    #import label
    label = read.table(file=label.fname, header=FALSE)
    label = label[,lab.col]
    label = unlist(label)
    label = as.factor(label)
    
  }
  
  
  #print(dim(raw.data))
  #Resolve missing value by median
  X.median = apply(raw.data,2,median,na.rm=TRUE)
  raw.data = t(apply(raw.data,1,replace.missing,missing=missing.char,rep=X.median))
  
  #regression 
  if ((!is.na(regression.file)) && (!is.na(regression.col.first)) && !(is.na(regression.col.last))){
    PCs = read.table(file=regression.file, header=FALSE, sep='')
    PCs = data.matrix(PCs[,regression.col.first:regression.col.last])
    #print(dim(raw.data))
    cat(paste0("Correct for covariates using ",reg.method," model\n"))
    raw.data = apply(raw.data,2,do.glm,PC=PCs,method=reg.method)
    #print(dim(raw.data))
    #raw.data = lm(as.matrix(raw.data) ~ PCs, na.action=na.exclude)$residuals
  }
  
  if (exists("snp.info")){
    save(raw.data,label,snp.info,ind.info,file=file.name)
  }else{
    save(raw.data,label,file=file.name)
  }
  
  no.marker = dim(raw.data)[2]
  no.individual = dim(raw.data)[1]
  cat(paste0("Input data: ",no.individual," individuals, ",no.marker," markers\n"))
  #Save new experiment condition
  file.name = file.path(result.dir,"RData","condition.RData")
  #load some parameters to add more parameters
  save(threshold,min.fst,max.thread,no.marker,no.individual,reanalysis,result.dir,
       method,min.in.group,datatype,nonlinear,plot.as.pdf,no.plot,file=file.name)
  
  
  #Check if node 1 is existed
  file.name = file.path(result.dir,"RData","node1.RData")
  if (!file.exists(file.name)){
    save(index,file=file.name)    
  }
  
  return(result.dir)
}

#Manipulate categorical data
pasre.categorical.data <- function(files){  
  
  map.to.zero.one.list <- function(cate.data,uni.cate){
    ret = c()
    ar.uni.cate=strsplit(uni.cate,'#@')[[1]]
    for (i in 1:length(ar.uni.cate)){
      tmp = rep(0,length(cate.data))
      tmp[which(cate.data == ar.uni.cate[i])] = 1
      ret = cbind(ret,tmp)
    }

    return(ret)
  }
  
  map.to.zero.one.matrix <- function(cate.data,uni.cate){
    ret = c()
    ar.uni.cate=strsplit(uni.cate,'#@')[[1]]
    for (i in 1:length(ar.uni.cate)){
      tmp = rep(0,length(cate.data))
      tmp[which(cate.data == ar.uni.cate[i])] = 1
      ret = cbind(ret,tmp)
    }
    drop.col = which(colSums(ret) == max(colSums(ret)))
    ret = ret[,-c(drop.col)]
    return(list('test'=ret))
  }
  
  raw.data = NULL
  cat(paste0("Loading the input files ... \n"))
  
  for (fname in files){
    #test the separator
    oneline = read.table(file=fname, header=FALSE, sep=',',nrows=1)
    #separated by comma
    if (dim(oneline)[2]>1){
      tmp_geno = read.table(file=fname, header=FALSE, sep=',')
    }else{
      #separated by white space
      tmp_geno = read.table(file=fname, header=FALSE)
    }
    
    raw.data = rbind(raw.data,tmp_geno)
  }
  
  n.row = dim(raw.data)[1]
  uni.cate = apply(raw.data,2,unique)
  if (is.list(uni.cate)){
    uni.cate = mapply(paste,uni.cate,sep="#@",collapse="#@")
    raw.data = mapply(map.to.zero.one.list,cate.data=raw.data,uni.cate=uni.cate) 
  }else{
    uni.cate = apply(uni.cate,2,paste,sep="#@",collapse="#@")
    raw.data = mapply(map.to.zero.one.matrix,cate.data=raw.data,uni.cate=uni.cate)
  }
  
  raw.data = unlist(raw.data)
  raw.data = matrix(raw.data, nrow=n.row, byrow=F)
  
  raw.data = as.data.frame(raw.data)
  index = seq(1,length(raw.data[,1]))
  rownames(raw.data) = index
  return(raw.data)
}



