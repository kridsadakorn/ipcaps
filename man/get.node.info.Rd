\name{get.node.info}
\alias{get.node.info}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
A function to obtain the information for specified node.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used to obtain the information for specified node from the output object of ipcaps().
}
\usage{
get.node.info(cluster.obj,node)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cluster.obj}{
%%     ~~Describe \code{x} here~~
The object which is returned from the function \emph{ipcaps}.
}
  \item{node}{
%%     ~~Describe \code{x} here~~
To specify the node number as shown in the HTML output files.}

}
\details{
%%  ~~ If necessary, more details than the description above ~~

}
\value{

  The returned object is a \emph{list} of node information including by \emph{PCs}, \emph{eigen.fit}, \emph{index}, \emph{label}. \emph{PCs} represent the pricipal components of this node. \emph{eigen.fit} is the EigenFit value of this node. \emph{index} is a vector of row number (individuals) of \emph{raw.data} (input data). \emph{label} is the vector of labels of all individuals that belongs to this node.

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

#Here, to obtain the information of specified node, for example, node 3
node.info <- get.node.info(my.cluster,3)
ls(node.info)
#[1] "PCs"       "eigen.fit" "index"     "label"

}

