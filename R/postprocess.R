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

postprocess <- function( result.dir, reanalysis=FALSE){

  file.name = file.path(result.dir,"RData","leafnode.RData")
  load(file=file.name)
  file.name = file.path(result.dir,"RData","condition.RData")
  load(file=file.name)

  if (length(leaf.node) == 0){
    leaf.node = c(1)
    save(leaf.node,file=file.name, compress = 'bzip2')
  }

  #Generate HTML output file
  cluster.tab = export.groups(result.dir)
  save.html(result.dir)
  if (no.plot == FALSE){
    save.plots(result.dir)
  }
  cat("The result files were saved at: ",result.dir,"\n")
  return (cluster.tab)

}


