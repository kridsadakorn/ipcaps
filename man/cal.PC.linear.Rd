\name{cal.PC.linear}
\alias{cal.PC.linear}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
A function for linear principal component analysis (PCA)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function is able to handle 2 types of data; linear and SNP data set in additive coding (0, 1, and 2).
}
\usage{
cal.PC.linear(X, PCscore=TRUE, no.pc=NA, data.type="linear", XXT=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{x} here~~
A data matrix which rows represent samples and columns represent features.
}
  \item{PCscore}{
%%     ~~Describe \code{x} here~~
To specify whether scaled PCs will be calculated or not. If \emph{FALSE}, eigenvectors are returned instead. Default = \emph{TRUE}.
}
  \item{no.pc  }{
%%     ~~Describe \code{x} here~~
A number of PCs to be calculated. If \emph{no.pc} is set, PCs are patially calculated. Otherwise all PCs are obtained after calculation. Default = \emph{NA}.
}
  \item{data.type}{
%%     ~~Describe \code{x} here~~
To specify a type of data matrix \emph{X}. It can be set to "linear" and "snp". Default = "linear".
}
  \item{XXT}{
%%     ~~Describe \code{x} here~~
To specify how pricipal components (PCs) are calculated. If \emph{TRUE}, PCs are calculated from \emph{X.t(X)}, otherwise \emph{X} is used directly. Default = \emph{TRUE}.
}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
\emph{XXT} is useful option especially an input matrix \emph{X} contains many columns. Enabling this option, it helps to reduce computation complexity. Regardless the option \emph{XXT} is enable or not, optained PCs are the same.
}
\value{

  The returned value is a list with components:
  \describe{
    \item{\emph{PC}}{a PC matrix which rows represent samples and columns represent PCs.}
    \item{\emph{evalue}}{a vector of eigen values.}
    }

}
\references{
%% ~put references to the literature/web site here ~
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

data(simSNP) #to load simulated data set

PCs <- cal.PC.linear(raw.data)
summary(PCs)
#       Length Class  Mode   
#PC     577600 -none- numeric
#evalue    760 -none- numeric

print(PCs$PC[1:5,1:3])
#           [,1]       [,2]        [,3]
#[1,] -1.0326907  0.4161087 -1.43709887
#[2,] -0.5855313 -0.1809882 -0.53587760
#[3,]  1.2100050 -1.5454890  3.22493101
#[4,] -2.5581070  0.7142676  0.04434291
#[5,]  0.2211739 -1.8262511 -1.31319919

print(PCs$evalue[1:3])
#[1] 205508.24  37222.45  35198.47

plot.3views(PCs$PC[,1:3],label)

PCs <- cal.PC.linear(raw.data,PCscore = FALSE)
summary(PCs)
#       Length Class  Mode   
#PC     577600 -none- numeric
#evalue    760 -none- numeric

print(PCs$PC[1:5,1:3])
#              [,1]          [,2]          [,3]
#[1,] -0.0022780102  0.0021567724 -0.0076599296
#[2,] -0.0012916223 -0.0009380973 -0.0028562995
#[3,]  0.0026691474 -0.0080105706  0.0171893146
#[4,] -0.0056429226  0.0037021882  0.0002363536
#[5,]  0.0004878869 -0.0094658149 -0.0069995278

print(PCs$evalue[1:3])
#[1] 205508.24  37222.45  35198.47

plot.3views(PCs$PC[,1:3],label)

PCs <- cal.PC.linear(raw.data,no.pc=3)
summary(PCs)
#       Length Class  Mode   
#PC     2280   -none- numeric
#evalue    3   -none- numeric

print(PCs$PC[1:5,])
#           [,1]       [,2]        [,3]
#[1,]  1.0326907  0.4161087 -1.43709887
#[2,]  0.5855313 -0.1809882 -0.53587760
#[3,] -1.2100050 -1.5454890  3.22493101
#[4,]  2.5581070  0.7142676  0.04434291
#[5,] -0.2211739 -1.8262511 -1.31319919

print(PCs$evalue)
#[1] 205508.24  37222.45  35198.47

plot.3views(PCs$PC[,1:3],label)
}

