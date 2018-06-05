
#' Detecting top discriminators between two groups
#'
#' @description This function detects top discriminators that contribute to
#' group separation based on the fixation index (Fst).
#'
#' @param cluster.obj The object which is returned from \code{\link{ipcaps}}.
#' @param group1 To specify the first group number to be compared. (also see
#' \code{use.node.number})
#' @param group2 To specify the second group number to be compared. (also see
#' \code{use.node.number})
#' @param bim.file Option: In case that SNP information is not provided to
#' \code{\link{ipcaps}}, an absolute path of SNP information file is required.
#' It needs to be in PLINK format (bim). See more details at:
#' \url{http://zzz.bwh.harvard.edu/plink/data.shtml}.
#' @param use.node.number To specify whether a group number or a node number
#' is be used. If TRUE, a node nubmer is used instead. Default = FALSE.
#' @param num.top A number of top Fst SNPs to be returned. Default = 100.
#'
#' @return The returned value is a data.frame of SNP information sorting by Fst
#' in descending order, which contains 7 columns, chr, SNP, centimorgans,
#' position, allele1, allele2, and Fst. The column 1-6 are SNP information from
#' the bim file. The column Fst contains estimated Fst between group1 and group2.
#'
#' @export
#'
#' @examples
#'
#' \donttest{
#' # Importantly, bed file, bim file, and fam file are required
#' # Use the example files embedded in the package
#' BED.file <- system.file("extdata","simSNP.bed",package="IPCAPS")
#' LABEL.file <- system.file("extdata","simSNP_individuals.txt",package="IPCAPS")
#' my.cluster <- ipcaps(bed=BED.file,label.file=LABEL.file,lab.col=2,out=tempdir())
#' table(my.cluster$cluster$label,my.cluster$cluster$group)
#' # 1 2 3 4 5 6
#' # outlier4 5 4 1 0 0 0
#' # pop1 0 0 0 0 250 0
#' # pop2 0 0 0 0 0 250
#' # pop3 0 0 0 250 0 0
#'
#' #Identify top discriminators between groups, for example, group 4 and group 5
#' top.snp <-top.discriminator(my.cluster,4,5)
#' #or, specify the bim file
#' #top.snp <-top.discriminator(my.cluster,4,5,bim.file="simSNP.bim")
#' head(top.snp)
#' # chr SNP centimorgans position allele1 allele2 Fst
#' #V5452 1 marker5452 0 54520000 A T 0.11337260
#' #V2348 1 marker2348 0 23480000 A T 0.11194490
#' #V8244 1 marker8244 0 82440000 A T 0.09556580
#' #V5972 1 marker5972 0 59720000 A T 0.08747794
#' #V3561 1 marker3561 0 35610000 A T 0.08725860
#' #V8419 1 marker8419 0 84190000 A T 0.08293494
#'
#' #Alternatively, it is possible to refer to node numbers instead of group numbers
#' table(my.cluster$cluster$label,my.cluster$cluster$node)
#' # 2 4 6 7 8 9
#' # outlier4 5 4 1 0 0 0
#' # pop1 0 0 0 0 250 0
#' # pop2 0 0 0 0 0 250
#' # pop3 0 0 0 250 0 0
#'
#' #Identify top discriminators between groups, for example, node 7 and node 8
#' top.snp2 <-top.discriminator(my.cluster,7,8,use.node.number=TRUE)
#' head(top.snp2)
#' # chr SNP centimorgans position allele1 allele2 Fst
#' #V5452 1 marker5452 0 54520000 A T 0.11337260
#' #V2348 1 marker2348 0 23480000 A T 0.11194490
#' }

top.discriminator <- function(cluster.obj,group1,group2,bim.file,use.node.number=FALSE,num.top=100){

  raw.data <- NULL

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



