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

save.plots.label.html <- function(output.dir){
  load(file.path(output.dir,"RData","leafnode.RData"))
  load(file.path(output.dir,"RData","tree.RData"))
  load(file.path(output.dir,"RData","rawdata.RData"))
  load(file.path(output.dir,"RData","condition.RData"))
  global.label = label
  node.list = sort(tree$node)
  
  txt_leafnode = ""
  txt_data = ""
  
  test.dir=file.path(output.dir,"images.new")
  img.dir="images"
  if (file.exists(test.dir)){
    img.dir="images.new"
  }
  
  for (i in node.list){
    parent_node = ""
    if (i > 1){
      parent_node = tree$parent.node[which(tree$node == i)]
    }
    PCs.file = file.path(output.dir,"RData",paste0("node",i,".RData"))
    load(PCs.file)
    label = global.label[index]
    list.sum = c()
    u.label = sort(unique(label))
    for (j in 1:length(u.label)){
      co = length(label[label == u.label[j]])
      strout = paste0(u.label[j],"&nbsp;(",co,")")
      list.sum = c(list.sum,strout)
    }
    content = ""
    for (s in list.sum){
      content = paste0(content,s,"<br>")
    }
    if (file.exists(file.path(output.dir,"images",paste0("scatter_by_label_preview",i,".jpg")))){
      txt_data = paste0(txt_data,"[{v:'",i,"', f:'<div class=\"box_class\" align=\"center\"><p class=\"head_class\">Node ",i,"</p><a href=\"",img.dir,"/scatter_by_label",i,".png\" target=\"_blank\"><img src=\"",img.dir,"/scatter_by_label_preview",i,".jpg\" width=200 height=200><br />view</a></div>'}, '",parent_node,"', '']")
    }else if (file.exists(file.path(output.dir,"images",paste0("scatter_by_label_preview",i,".pdf")))){
      txt_data = paste0(txt_data,"[{v:'",i,"', f:'<div class=\"box_class\" align=\"center\"><p class=\"head_class\">Node ",i,"</p><a href=\"",img.dir,"/scatter_by_label",i,".pdf\" target=\"_blank\"><embed src=\"",img.dir,"/scatter_by_label_preview",i,".pdf\" width=200 height=200><br />view</a></div>'}, '",parent_node,"', '']")
    }else{
      txt_data = paste0(txt_data,"[{v:'",i,"', f:'<div class=\"box_class\"><p class=\"head_class\">Node ",i,"</p><p class=\"subhead_class\">under cutoff (<",min.in.group,")</p><br>",content,"</div>'}, '",parent_node,"', '']")
    }
    if (!(i == node.list[length(node.list)])){
      txt_data = paste0(txt_data,",\n")
    }else{
      txt_data = paste0(txt_data,"\n")
    }
    
    no_idx = which(i == node.list) - 1
    if (i %in% leaf.node){
      txt_leafnode = paste0(txt_leafnode,"data.setRowProperty(",no_idx,", 'style', 'border: 3px solid #DB6E6E; background-color:#FFE1E1');\n")
    }
  }
  
  txt_title = "Scatterplots colored by predefined labels"
  txt_body = "Scatterplots colored by predefined labels"
 
  txt_html = output.template$template
  txt_html[output.template$lno_data] = txt_data
  txt_html[output.template$lno_leafnode] = txt_leafnode
  txt_html[output.template$lno_body] = txt_body
  txt_html[output.template$lno_title] = txt_title

  fo = file(file.path(output.dir,"tree_scatter_label.html"),"w")
  for (i in txt_html){ write(i,fo)}
  close(fo)
  
}



