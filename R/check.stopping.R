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

check.stopping <- function(eigen.value, threshold){
  
  diff.xy = function(x,y){
    ret = abs(x-y)
    return(ret)
  }
  
  cal.eigen.fit = function(eigen.value){
    I=log(eigen.value)
    X = I
    Y = I
    Y = Y[2:length(I)]
    X = X[1:(length(I)-1)]
    ret = mapply(diff.xy,X,Y)
    return(ret)
  }
  
  eigen.fit.vec = cal.eigen.fit(eigen.value) 
  eigen.fit = max(eigen.fit.vec[1:2])
  no.significant.PC = length(eigen.fit.vec[1:which(eigen.fit.vec == eigen.fit)[1]])
  if (no.significant.PC<3){
    no.significant.PC = 3
  }
  
  ret = list("status"=0,"eigen.value"=eigen.value,"eigen.fit"=eigen.fit, "threshold"=threshold, "no.significant.PC" = no.significant.PC)
  if (eigen.fit < threshold){  #case of status = 1, no more spliting, stopping criteria are met
    ret = list("status"=1,"eigen.value"=eigen.value,"eigen.fit"=eigen.fit, "threshold"=threshold)
  }
  return(ret)
}


