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


ipcaps <- function(bed=NA,rdata=NA,files=NA,label.file=NA,lab.col=1,out,plot.as.pdf=FALSE,method='mix',
                   missing=NA,covariate=NA,cov.col.first=NA,cov.col.last=NA,threshold=0.18,min.fst=0.0008,
                   min.in.group=20,no.plot = FALSE){
  
  label.fname = label.file
  label.column = lab.col
  result.dir = out
  rerun = FALSE
  rdata.infile = rdata
  bed.infile = bed
  datatype = 'snp'
  nonlinear = FALSE
  missing.char = missing
  regression.file = covariate
  regression.col.first = cov.col.first
  regression.col.last = cov.col.last
  file.list = files
  max.thread = 1

  start.time <- Sys.time()
  
  if (length(threshold)<=0){
    threshold=0.18 #work in general. The lowest value is 0.03 for the data without outlier
  }
  
  if (length(min.fst)<=0){
    min.fst=0.0008
  }
  
  usage= paste0("Usage: ?ipcaps\n")
  
  if (length(result.dir)==0){
    cat(usage)
    quit()
  }
  if ((length(label.fname)==0 || is.na(file.list)) && (length(rdata.infile)==0) && 
      (length(bed.infile)==0) && is.na(cate.list)){
    cat(usage)
    quit()
  }
  
  cat(paste0("Running ... IPCAPS \n\toutput: ",result.dir," \n"))
  
  if (length(label.fname)>0){
    cat(paste0("\tlabel file: ",label.fname,"\n"))
  }
  
  if (length(label.column)>0){
    if (is.na(as.integer(label.column))){
      #except only comma as separator, otherwise
      if (length(strsplit(label.column,',')[[1]])){
        if (anyNA(as.integer(strsplit(label.column,',')[[1]]))){
          label.column = 1
          cat(paste0("\tlabel column: ",label.column,"\n"))
        }else{
          label.column = as.integer(strsplit(label.column,',')[[1]])
        }
      }else{
        label.column = 1
        cat(paste0("\tlabel column: ",label.column,"\n"))
      }
    }else{
      label.column = as.integer(label.column)
      cat(paste0("\tlabel column: ",label.column,"\n"))
    }
  }else{
    label.column = 1
  }
  
  if (length(threshold)>0){
    cat(paste0("\tthreshold: ",threshold,"\n"))
  }
  
  if (length(min.fst)>0){
    cat(paste0("\tminimum Fst: ",min.fst,"\n"))
  }
  
  if (length(max.thread)<=0){
    max.thread=NA
  }
  if (is.na(max.thread)){
    max.thread = 1
  }else{
    if (max.thread < 1){
      max.thread = 1
    }
  }
  cat(paste0("\tmaximum thread: ",max.thread,"\n"))
  
  if (length(method)>0){
    cat(paste0("\tmethod: ",method,"\n"))
  }else{
    method="mix"
  }
  
  if (!is.na(rdata.infile)){
    if (file.exists(rdata.infile)){
      cat(paste0("\trdata: ",rdata.infile,"\n"))
    }else{
      rdata.infile = NA
    }
  }else{
    rdata.infile = NA
  }
  
  if (!is.na(bed.infile)){
    if (file.exists(bed.infile)){
      cat(paste0("\tbed: ",bed.infile,"\n"))
    }else{
      bed.infile = NA
    }
  }else{
    bed.infile = NA
  }
  
  if (!is.na(file.list)){
    cat(paste0("\tfiles: \n"))
    print(file.list)
  }
  
  if (plot.as.pdf){
    cat(paste0("\tplot.as.pdf: TRUE\n"))
    plot.as.pdf=TRUE
  }else{
    plot.as.pdf=FALSE
  }
  
  if (length(min.in.group)>0){
    if (min.in.group < 5){
      min.in.group = 5
    }
    cat(paste0("\tminimum in group: ",min.in.group,"\n"))
  }else{
    min.in.group=20
  }
  
  if (length(datatype)>0){
    cat(paste0("\tdata type: ",datatype,"\n"))
  }else{
    datatype="snp"
  }
  
  if (length(missing.char)>0){
    cat(paste0("\tmissing: ",missing.char,"\n"))
  }else{
    missing.char = NA
  }
  
  if (length(regression.file)>0 && !is.na(regression.file)){
    cat(paste0("\tcovariate file: ",regression.file,"\n"))
  }else{
    regression.file = NA
  }
  
  if (!is.na(regression.col.first) && regression.col.first>0){
    cat(paste0("\tcovariate first column: ",regression.col.first,"\n"))
  }else{
    regression.col.first = NA
  }
  
  if (!is.na(regression.col.last) && regression.col.last>0){
    cat(paste0("\tcovariate last column: ",regression.col.last,"\n"))
  }else{
    regression.col.last = NA
  }
  
  
  # read in the R functions and libraries, Don't change anything here
  # require(Matrix,quietly=TRUE)
  # require(expm,quietly=TRUE)
  # require(e1071,quietly=TRUE)
  # require(fpc)
  # require(Rmixmod,quietly=TRUE)
  # require(LPCM)
  # require(apcluster)
  # require(rARPACK)
  # require(igraph)
  #requireNamespace()
  
  
  #preprocessing step
  cat(paste0("In preprocessing step\n"))
  
  result.dir=preprocess(files=file.list, label=label.fname, lab.col=label.column, 
                        rdata.infile=rdata.infile, bed.infile=bed.infile, cate.list=cate.list, 
                        result.dir=result.dir, threshold=threshold, min.fst=min.fst, 
                        reanalysis=rerun, method=method, min.in.group=min.in.group,datatype=datatype,
                        nonlinear=nonlinear, missing.char=missing.char,
                        regression.file=regression.file, regression.col.first=regression.col.first,
                        regression.col.last=regression.col.last,plot.as.pdf=plot.as.pdf,no.plot=no.plot)
  
  #job scheduler
  cat(paste0("Start calculating\n"))
  
  #create the first task for Node 1 
  #status 0 = unprocessed, 1 = processing, 2 = done , -1 = deleted
  
  node = c(1)
  parent.node = c(0)
  status = c(0)
  
  tree = data.frame(node,parent.node,status)
  file.name = file.path(result.dir,"RData","tree.RData")
  save(tree,file=file.name)
  
  while (TRUE){
    
    file.name = file.path(result.dir,"RData","tree.RData")
    load(file.name)
    
    #check for terminate loop
    row2 = which(tree$status==2)
    row_1 = which(tree$status==-1)
    if ((length(row2)+length(row_1))==length(tree$node)){
      break
    }
    
    file.name = file.path(result.dir,"RData","condition.RData")
    load(file=file.name)
    
    #take one node to process
    which.row = which(tree$status==0)
    if (length(which.row)>0){
      exe.node.list = tree$node[which.row]
      for(idx in exe.node.list) {
        process.each.node(node=idx,work.dir=result.dir)
      }
      
    }else{
      break
    }
    
  }
  
  cat(paste0("In post process step\n"))
  cluster.tab = postprocess(result.dir=result.dir)
  
  end.time <- Sys.time()
  run.time = as.numeric(end.time-start.time, units = "secs")
  cat(paste0("Total runtime is ",run.time," sec\n"))
  file.name = file.path(result.dir,"RData","runtime.RData")
  save(run.time,file=file.name)
  ret <- list("output.dir"=result.dir,"cluster"=cluster.tab)
  return(ret)
}
# Check the result files in your output directory
# groups.txt contains the assigned groups of samples
# tree_text.html is the text result shown as binary tree
# tree_scatter_cluster.html is the scatter plots colored by clustering result
# tree_scatter_label.html is the scatter plots colored by labels
# tree_scree.html is the scree plots of eigen values			


