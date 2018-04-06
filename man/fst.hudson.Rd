\name{fst.hudson}
\alias{fst.hudson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
A function for average fixation index (Fst) calculation between two groups.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Fixation index (Fst) calculation was implemented using Hudson method [1,2].

}
\usage{
fst.hudson(raw.data, idx.p1, idx.p2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raw.data}{
%%     ~~Describe \code{x} here~~
A matrix contains the number 0, 1, and 2 representing SNP in additive coding. Rows represent individuals and columns represent SNP.
}
  \item{idx.p1}{
%%     ~~Describe \code{x} here~~
An integer vector contains the row indices of first population in \emph{raw.data}.
}
  \item{idx.p2}{
%%     ~~Describe \code{x} here~~
An integer vector contains the row indices of second population in \emph{raw.data}.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~

}
\value{
The function returns an average Fst value between 2 specified groups.
}
\references{
%% ~put references to the literature/web site here ~
[1] Bhatia,G. et al. (2013) Estimating and interpreting Fst: The impact of rare variants. Genome Res., 23, 1514:1521.

[2] Hudson,R.R. et al. (1992) Estimation of levels of gene flow from DNA sequence data. Genetics, 132, 583:589.

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
idx1 <- which(label == 'pop1')
idx2 <- which(label == 'pop2')
fst <- fst.hudson(raw.data,idx1,idx2)

print(fst)
#[1] 0.004843083
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

