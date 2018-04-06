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



#The function to calculate X multipies t(X), for large matrix

XXt <- function(X,window.size= 5){
  n.col = dim(X)[2]
  window.size = as.integer(window.size)
  if (window.size < 1){
    cat(paste0("In matrices.split.by.col(), window.size must be possitive integer\n"))
    return(NULL)
  }
  if (is.null(n.col)){
    cat(paste0("In matrices.split.by.col(), n.col is null.\n"))
    return(NULL)
  }
  if (n.col < window.size){
    cat(paste0("In matrices.split.by.col(), n.col is than window.size\n"))
    return(NULL)
  }
  
  no.split = n.col %/% window.size
  
  split.matrices = lapply(1:no.split, function(col){
    X[,((col-1)*window.size+1):(col*window.size)]})
  
  if ((no.split * window.size) < n.col){
    split.matrices[[length(split.matrices)+1]] = as.matrix(X[,(no.split*window.size+1):n.col])
  }
  
  #Can be parallel
  
  partial.mul = NULL
  for (i in 1:length(split.matrices)){
    partial.mul[[i]] = split.matrices[[i]] %*% t(split.matrices[[i]])
  }
  
  mul = NULL
  for (i in 2:length(partial.mul)){
    if (is.null(mul)){
      mul = partial.mul[[1]] + partial.mul[[i]]
    }else{
      mul = mul + partial.mul[[i]]
    }
  }
  
  return(mul)
}





