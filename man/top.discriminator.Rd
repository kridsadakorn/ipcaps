\name{top.discriminator}
\alias{top.discriminator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
A function to detect top discriminators between two groups.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function detects top discriminators that contribute to group separation based on the fixation index (Fst).
}
\usage{
top.discriminator(cluster.obj,group1,group2,bim.file,use.node.number=FALSE,num.top=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cluster.obj}{
%%     ~~Describe \code{x} here~~
The object which is returned from the function \emph{ipcaps}.
}
  \item{group1}{
%%     ~~Describe \code{x} here~~
To specify the first group number to be compared. (also see \emph{use.node.number})
}
  \item{group2}{
%%     ~~Describe \code{x} here~~
To specify the second group number to be compared. (also see \emph{use.node.number})
}
  \item{bim.file}{
%%     ~~Describe \code{x} here~~
Option: In case that SNP information is not provided to the \emph{ipcaps} function, an absolute path of SNP information file is required. It needs to be in PLINK format (.\emph{bim}). See more details at: 

http://zzz.bwh.harvard.edu/plink/data.shtml
}
  \item{use.node.number}{
%%     ~~Describe \code{x} here~~
To specify whether a group number or a node number will be used. If \emph{TRUE}, a node nubmer is used instead. Default = \emph{FALSE}.
}
  \item{num.top}{
%%     ~~Describe \code{x} here~~
A number of top Fst SNPs to be returned. Default = 100.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~

}
\value{

  The returned value is a \emph{data.frame} of SNP information sorting by Fst in descending  order, which contains 7 columns, \emph{chr}, \emph{SNP}, \emph{centimorgans}, \emph{position}, \emph{allele1}, \emph{allele2}, and \emph{Fst}. The column 1-6 are SNP information from the \emph{bim} file. The column \emph{Fst} contains estimated Fst between \emph{group1} and \emph{group2}.

}
\references{

}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# Importantly, bed file, bim file, and fam file are required
# Use the example files embedded in the package

BED.file <- system.file("extdata","simSNP.bed",package="IPCAPS")
LABEL.file <- system.file("extdata","simSNP_individuals.txt",package="IPCAPS")

my.cluster <- ipcaps(bed=BED.file,label.file=LABEL.file,lab.col=2,out=getwd())
table(my.cluster$cluster$label,my.cluster$cluster$group)
#             1   2   3   4   5   6
#  outlier4   5   4   1   0   0   0
#  pop1       0   0   0   0 250   0
#  pop2       0   0   0   0   0 250
#  pop3       0   0   0 250   0   0

#Identify top discriminators between groups, for example, group 4 and group 5
top.snp <-top.discriminator(my.cluster,4,5)
#or, specify the bim file
#top.snp <-top.discriminator(my.cluster,4,5,bim.file="simSNP.bim")
head(top.snp)
#      chr        SNP centimorgans position allele1 allele2        Fst
#V5452   1 marker5452            0 54520000       A       T 0.11337260
#V2348   1 marker2348            0 23480000       A       T 0.11194490
#V8244   1 marker8244            0 82440000       A       T 0.09556580
#V5972   1 marker5972            0 59720000       A       T 0.08747794
#V3561   1 marker3561            0 35610000       A       T 0.08725860
#V8419   1 marker8419            0 84190000       A       T 0.08293494

#Alternatively, it is possible to refer to node numbers instead of group numbers
table(my.cluster$cluster$label,my.cluster$cluster$node)
#             2   4   6   7   8   9
#  outlier4   5   4   1   0   0   0
#  pop1       0   0   0   0 250   0
#  pop2       0   0   0   0   0 250
#  pop3       0   0   0 250   0   0

#Identify top discriminators between groups, for example, node 7 and node 8
top.snp2 <-top.discriminator(my.cluster,7,8,use.node.number=TRUE)
head(top.snp2)
#      chr        SNP centimorgans position allele1 allele2        Fst
#V5452   1 marker5452            0 54520000       A       T 0.11337260
#V2348   1 marker2348            0 23480000       A       T 0.11194490
#V8244   1 marker8244            0 82440000       A       T 0.09556580
#V5972   1 marker5972            0 59720000       A       T 0.08747794
#V3561   1 marker3561            0 35610000       A       T 0.08725860
#V8419   1 marker8419            0 84190000       A       T 0.08293494


}

