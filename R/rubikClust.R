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

rubikClust <- function(X, min.space=0.4,rotation=TRUE){
  
  rotate.xy <- function(xy,degree=0){
    t = degree*pi/180.0
    x = xy[1]
    y = xy[2]
    n.x = x * cos(t) - y * sin(t)
    n.y = y * cos(t) + x * sin(t)
    return(c(n.x,n.y))
  }
  
  #Initial X, Y, Z variables
  
  sort.x = sort(X[,1])
  scale.x = sort.x /(max(sort.x)-min(sort.x))
  diff.x = scale.x[2:length(scale.x)] - scale.x[1:length(scale.x)-1]
  global.border.x = which(diff.x > min.space)
  global.space.x = diff.x[global.border.x]
  
  sort.y = sort(X[,2])
  scale.y = sort.y /(max(sort.y)-min(sort.y))
  diff.y = scale.y[2:length(scale.y)] - scale.y[1:length(scale.y)-1]
  global.border.y = which(diff.y > min.space)
  global.space.y = diff.y[global.border.y]
  
  sort.z = sort(X[,3])
  scale.z = sort.z /(max(sort.z)-min(sort.z))
  diff.z = scale.z[2:length(scale.z)] - scale.z[1:length(scale.z)-1]
  global.border.z = which(diff.z > min.space)
  global.space.z = diff.z[global.border.z]
  
  theta.xy = 0
  theta.yz = 0
  theta.xz = 0
  XYZ = X
  
  if (rotation == TRUE){
    #Rotate X and Y
    
    for (th in (1:89)){
      rotated = t(apply(X[,1:2],MARGIN=1,FUN=rotate.xy,degree=th))
      
      sorted = sort(rotated[,1])
      scaled = sorted/(max(sorted)-min(sorted))
      differences = scaled[2:length(scaled)] - scaled[1:length(scaled)-1]
      border1 = which(differences > min.space)
      space1 = differences[border1]
      
      sorted = sort(rotated[,2])
      scaled = sorted/(max(sorted)-min(sorted))
      differences = scaled[2:length(scaled)] - scaled[1:length(scaled)-1]
      border2 = which(differences > min.space)
      space2 = differences[border2]
      
      #if ((length(border1)+length(border2)) > (length(global.border.x)+length(global.border.y)) ||
      #    ((length(border1)+length(border2)) == (length(global.border.x)+length(global.border.y)) &&
      #     ((sum(space1)+sum(space2)) > (sum(global.space.x)+sum(global.space.y))))){
      #if ((length(border1)+length(border2)) > (length(global.border.x)+length(global.border.y))) {
      if ((sum(space1)+sum(space2)) > (sum(global.space.x)+sum(global.space.y))){
        global.border.x = border1
        global.border.y = border2
        global.space.x = space1
        global.space.y = space2
        theta.xy = th
      }
    }
    rotated = t(apply(X[,1:2],MARGIN=1,FUN=rotate.xy,degree=theta.xy))
    XYZ = cbind(rotated,X[,3])
    
    #Here rotate Y and Z
    
    for (th in (1:89)){
      rotated = t(apply(XYZ[,2:3],MARGIN=1,FUN=rotate.xy,degree=th))
      
      sorted = sort(rotated[,1])
      scaled = sorted/(max(sorted)-min(sorted))
      differences = scaled[2:length(scaled)] - scaled[1:length(scaled)-1]
      border1 = which(differences > min.space)
      space1 = differences[border1]
      
      sorted = sort(rotated[,2])
      scaled = sorted/(max(sorted)-min(sorted))
      differences = scaled[2:length(scaled)] - scaled[1:length(scaled)-1]
      border2 = which(differences > min.space)
      space2 = differences[border2]
      
      #if ((length(border1)+length(border2)) > (length(global.border.y)+length(global.border.z)) ||
      #    ((length(border1)+length(border2)) == (length(global.border.y)+length(global.border.z)) &&
      #     ((sum(space1)+sum(space2)) > (sum(global.space.y)+sum(global.space.z))))){
      #if ((length(border1)+length(border2)) > (length(global.border.y)+length(global.border.z))){
      if ((sum(space1)+sum(space2)) > (sum(global.space.y)+sum(global.space.z))){
        global.border.y = border1
        global.border.z = border2
        global.space.y = space1
        global.space.z = space2
        theta.yz = th
      }
    }
    
    rotated = t(apply(XYZ[,2:3],MARGIN=1,FUN=rotate.xy,degree=theta.yz))
    XYZ = cbind(XYZ[,1],rotated)
    
    #Here rotate X and Z
    
    for (th in (1:89)){
      rotated = t(apply(XYZ[,c(1,3)],MARGIN=1,FUN=rotate.xy,degree=th))
      
      sorted = sort(rotated[,1])
      scaled = sorted/(max(sorted)-min(sorted))
      differences = scaled[2:length(scaled)] - scaled[1:length(scaled)-1]
      border1 = which(differences > min.space)
      space1 = differences[border1]
      
      sorted = sort(rotated[,2])
      scaled = sorted/(max(sorted)-min(sorted))
      differences = scaled[2:length(scaled)] - scaled[1:length(scaled)-1]
      border2 = which(differences > min.space)
      space2 = differences[border2]
      
      #if ((length(border1)+length(border2)) > (length(global.border.x)+length(global.border.z)) ||
      #    ((length(border1)+length(border2)) == (length(global.border.x)+length(global.border.z)) &&
      #     ((sum(space1)+sum(space2)) > (sum(global.space.x)+sum(global.space.z))))){
      #if ((length(border1)+length(border2)) > (length(global.border.x)+length(global.border.z))){
      if ((sum(space1)+sum(space2)) > (sum(global.space.x)+sum(global.space.z))){
        global.border.x = border1
        global.border.z = border2
        global.space.x = space1
        global.space.z = space2
        theta.xz = th
      }
    }
    
    rotated = t(apply(XYZ[,c(1,3)],MARGIN=1,FUN=rotate.xy,degree=theta.xz))
    XYZ = cbind(rotated[,1],XYZ[,2],rotated[,2])
    
    sort.x = sort(XYZ[,1])
    scale.x = sort.x /(max(sort.x)-min(sort.x))
    diff.x = scale.x[2:length(scale.x)] - scale.x[1:length(scale.x)-1]
    global.border.x = which(diff.x > min.space)
    
    sort.y = sort(XYZ[,2])
    scale.y = sort.y /(max(sort.y)-min(sort.y))
    diff.y = scale.y[2:length(scale.y)] - scale.y[1:length(scale.y)-1]
    global.border.y = which(diff.y > min.space)
    
    sort.z = sort(XYZ[,3])
    scale.z = sort.z /(max(sort.z)-min(sort.z))
    diff.z = scale.z[2:length(scale.z)] - scale.z[1:length(scale.z)-1]
    global.border.z = which(diff.z > min.space)
  }
  
  #Assign clusters
  groups = rep(0,length(sort.x))
  group.idx = 0
  
  if ((length(global.border.x) > 0) && (length(global.border.y) > 0) && (length(global.border.z) > 0)){
    #DO X, Y ,Z
    global.border.x = c(global.border.x,length(sort.x))
    global.border.y = c(global.border.y,length(sort.y))
    global.border.z = c(global.border.z,length(sort.z))
    start.x = min(XYZ[,1]) - 1
    for (x in 1:length(global.border.x)){
      end.x = sort.x[global.border.x[x]]
      start.y = min(XYZ[,2]) - 1
      for (y in 1:length(global.border.y)){
        end.y = sort.y[global.border.y[y]]
        start.z = min(XYZ[,3]) - 1
        for (z in 1:length(global.border.z)){
          end.z = sort.z[global.border.z[z]]
          set.x = intersect(which(XYZ[,1] > start.x),which(XYZ[,1] <= end.x))
          set.y = intersect(which(XYZ[,2] > start.y),which(XYZ[,2] <= end.y))
          set.z = intersect(which(XYZ[,3] > start.z),which(XYZ[,3] <= end.z))
          target.set = intersect(intersect(set.x,set.y),set.z)
          if (length(target.set) > 0){
            group.idx = group.idx + 1
            groups[target.set] = group.idx
          }
          start.z = end.z
        }
        start.y = end.y
      }
      start.x = end.x
    }
  }else if ((length(global.border.x) > 0) && (length(global.border.y) > 0)){
    #DO X, Y
    global.border.x = c(global.border.x,length(sort.x))
    global.border.y = c(global.border.y,length(sort.y))
    start.x = min(XYZ[,1]) - 1
    for (x in 1:length(global.border.x)){
      end.x = sort.x[global.border.x[x]]
      start.y = min(XYZ[,2]) - 1
      for (y in 1:length(global.border.y)){
        end.y = sort.y[global.border.y[y]]
        set.x = intersect(which(XYZ[,1] > start.x),which(XYZ[,1] <= end.x))
        set.y = intersect(which(XYZ[,2] > start.y),which(XYZ[,2] <= end.y))
        target.set = intersect(set.x,set.y)
        if (length(target.set) > 0){
          group.idx = group.idx + 1
          groups[target.set] = group.idx
        }
        start.y = end.y
      }
      start.x = end.x
    }
  }else if ((length(global.border.x) > 0) && (length(global.border.z) > 0)){
    #DO X, Z
    global.border.x = c(global.border.x,length(sort.x))
    global.border.z = c(global.border.z,length(sort.z))
    start.x = min(XYZ[,1]) - 1
    for (x in 1:length(global.border.x)){
      end.x = sort.x[global.border.x[x]]
      start.z = min(XYZ[,3]) - 1
      for (z in 1:length(global.border.z)){
        end.z = sort.z[global.border.z[z]]
        set.x = intersect(which(XYZ[,1] > start.x),which(XYZ[,1] <= end.x))
        set.z = intersect(which(XYZ[,3] > start.z),which(XYZ[,3] <= end.z))
        target.set = intersect(set.x,set.z)
        if (length(target.set) > 0){
          group.idx = group.idx + 1
          groups[target.set] = group.idx
        }
        start.z = end.z
      }
      start.x = end.x
    }
  }else if ((length(global.border.y) > 0) && (length(global.border.z) > 0)){
    #Do Y, Z
    global.border.y = c(global.border.y,length(sort.y))
    global.border.z = c(global.border.z,length(sort.z))
    start.y = min(XYZ[,2]) - 1
    for (y in 1:length(global.border.y)){
      end.y = sort.y[global.border.y[y]]
      start.z = min(XYZ[,3]) - 1
      for (z in 1:length(global.border.z)){
        end.z = sort.z[global.border.z[z]]
        set.y = intersect(which(XYZ[,2] > start.y),which(XYZ[,2] <= end.y))
        set.z = intersect(which(XYZ[,3] > start.z),which(XYZ[,3] <= end.z))
        target.set = intersect(set.y,set.z)
        if (length(target.set) > 0){
          group.idx = group.idx + 1
          groups[target.set] = group.idx
        }
        start.z = end.z
      }
      start.y = end.y
    }
  }else if ((length(global.border.x) > 0)){
    #Do X
    global.border.x = c(global.border.x,length(sort.x))
    start.x = min(XYZ[,1]) - 1
    for (x in 1:length(global.border.x)){
      end.x = sort.x[global.border.x[x]]
      target.set = intersect(which(XYZ[,1] > start.x),which(XYZ[,1] <= end.x))
      if (length(target.set) > 0){
        group.idx = group.idx + 1
        groups[target.set] = group.idx
      }
      start.x = end.x
    }
  }else if ((length(global.border.y) > 0)){
    #Do y
    global.border.y = c(global.border.y,length(sort.y))
    start.y = min(XYZ[,2]) - 1
    for (y in 1:length(global.border.y)){
      end.y = sort.y[global.border.y[y]]
      target.set = intersect(which(XYZ[,2] > start.y),which(XYZ[,2] <= end.y))
      if (length(target.set) > 0){
        group.idx = group.idx + 1
        groups[target.set] = group.idx
      }
      start.y = end.y
    }
  }else if ((length(global.border.z) > 0)){
    #Do Z
    global.border.z = c(global.border.z,length(sort.z))
    start.z = min(XYZ[,3]) - 1
    for (z in 1:length(global.border.z)){
      end.z = sort.z[global.border.z[z]]
      target.set = intersect(which(XYZ[,3] > start.z),which(XYZ[,3] <= end.z))
      if (length(target.set) > 0){
        group.idx = group.idx + 1
        groups[target.set] = group.idx
      }
      start.z = end.z
    }
  }else{
    #Only 1 group
    groups = rep(1,length(sort.z))
  }
  
  return(groups)
  
  #For debugging
  #ret = list("XYZ"=XYZ,"borders"=c(global.border.x,global.border.y,global.border.z),"groups"=groups,
  #           "theta.xy"=theta.xy,"theta.yz"=theta.yz,"theta.xz"=theta.xz,
  #           "space.x"=global.space.x,"space.y"=global.space.y,"space.z"=global.space.z)
  #return(ret)
}



