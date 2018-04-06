\name{simSNP}
\alias{simSNP}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Synthetic dataset containing single nucleotide polymorphisms (SNP) and population labels.
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
The \emph{simSNP} is the simulated data set which consists of 10,000 independent SNPs and 760 individuals belonging to one of three populations (250 individuals each) and 10 outlying individuals. The pairwise genetic distance between populations is set to Fst=0.005 [1].  
}
\usage{data("simSNP")}
\format{

}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
  The \emph{simSNP} data set consists of:
  \describe{
    \item{\emph{raw.data}}{A numeric matrix of 760 rows and 10,000 columns representing individuals and SNPs respectively. The matrix \emph{raw.data} contains the number 0, 1, and 2 representing SNP in additive coding.}
    \item{\emph{label}}{A character vector of 760 elements containing labels or populations of 760 individuals which they belong. Three populations and outliers were labeled as "pop1", "pop2", "pop3", and "outlier4".}
    \item{\emph{PC}}{A numeric matrix of 760 rows and 10 columns containing the first ten principal components.}
    

}
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
[1] Balding DJ, Nichols RA. A method for quantifying differentiation between populations at multi-allelic loci and its implications for investigating identity and paternity. Genetica. 1995 Jun;96(1-2):3-12. 
}
\examples{

data(simSNP)
dim(raw.data)
#[1]   760 10000

raw.data[1:5,1:5]
#     [,1] [,2] [,3] [,4] [,5]
#[1,]    1    2    0    1    0
#[2,]    1    1    0    0    2
#[3,]    0    1    0    1    1
#[4,]    1    2    0    1    1
#[5,]    1    1    0    0    1

length(label)
#[1] 760

unique(label)
#[1] 760  10

dim(PC)
#[1] 760  10
}
\keyword{data sets}
